#!/usr/bin/env python3
from pathlib import Path

p=Path("external/SWIFT/src/cosmology.c")
s=p.read_text()

needle='#include <gsl/gsl_interp.h>\n#endif\n'
insert='''#include <gsl/gsl_interp.h>\n#include <gsl/gsl_errno.h>\n#endif\n\n#ifdef HAVE_LIBGSL\nstatic int k01d2_qag_diag(const int line, gsl_function *F, double a, double b,\n                          double epsabs, double epsrel, size_t limit, int key,\n                          gsl_integration_workspace *workspace,\n                          double *result, double *abserr) {\n  static int handler_disabled = 0;\n  if (!handler_disabled) {\n    gsl_set_error_handler_off();\n    handler_disabled = 1;\n  }\n  const int status = gsl_integration_qag(F, a, b, epsabs, epsrel, limit, key,\n                                         workspace, result, abserr);\n  if (status != GSL_SUCCESS) {\n    fprintf(stderr,\n            "K01D2_QAG_FAIL line=%d status=%d bounds=[%.17g,%.17g] epsabs=%.17g epsrel=%.17g result=%.17g abserr=%.17g\\n",\n            line, status, a, b, epsabs, epsrel,\n            result ? *result : NAN, abserr ? *abserr : NAN);\n    fflush(stderr);\n    exit(86);\n  }\n  return status;\n}\n#define gsl_integration_qag(F,A,B,EA,ER,L,K,W,R,E) \\\n  k01d2_qag_diag(__LINE__,F,A,B,EA,ER,L,K,W,R,E)\n#endif\n'''

if needle not in s:
    raise SystemExit("include anchor not found")
s=s.replace(needle,insert,1)
p.write_text(s)
print("K01D2_D1_SWIFT_QAG_DIAGNOSTIC_PATCH_APPLIED")
