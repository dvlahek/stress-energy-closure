#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

static double f1(double x, void *param) {
  double y = *(double *)param;
  double numerator = x * x * hypot(x, y);
  if (x > 20 + log(numerator)) return numerator * exp(-x);
  return numerator / (1.0 + exp(x));
}

static double f2(double w, void *param) {
  if (w == 0.0) return 0.0;
  return f1(1.0 / w, param) / (w * w);
}

int main(void) {
  gsl_set_error_handler_off();

  const double kb_ev_k = 8.617333262e-5;
  const double Tnu_K = 1.951757805;
  const double Tnu_eV = Tnu_K * kb_ev_k;
  const double m = 0.06;
  const double a_safe = 0.5e-7 / m;
  const double y_min = a_safe * m / Tnu_eV;
  const double y_max = m / Tnu_eV;
  const int N = 512;
  const double tolerances[] = {1e-13,1e-12,1e-11,1e-10};
  const int NT = 4;

  printf("Y_MIN %.17g\n", y_min);
  printf("Y_MAX %.17g\n", y_max);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);
  if (!w) return 2;

  for (int it=0; it<NT; ++it) {
    const double tol = tolerances[it];
    int fail1=0, fail2=0;
    int first_status1=0, first_status2=0;
    double first_y1=NAN, first_y2=NAN;
    double max_rel_abserr1=0.0, max_rel_abserr2=0.0;

    for (int i=0; i<N; ++i) {
      double t = (double)i/(double)(N-1);
      double y = exp(log(y_min) + t*(log(y_max)-log(y_min)));
      double r=0.0,e=0.0;
      gsl_function F1 = {&f1,&y};
      gsl_function F2 = {&f2,&y};

      int s1 = gsl_integration_qag(&F1,0.0,1.0,0.0,tol,100000,GSL_INTEG_GAUSS61,w,&r,&e);
      if (isfinite(r) && r!=0.0) {
        double re=fabs(e/r); if (re>max_rel_abserr1) max_rel_abserr1=re;
      }
      if (s1 != GSL_SUCCESS) {
        fail1++;
        if (!first_status1) { first_status1=s1; first_y1=y; }
      }

      r=0.0; e=0.0;
      int s2 = gsl_integration_qag(&F2,0.0,1.0,0.0,tol,100000,GSL_INTEG_GAUSS61,w,&r,&e);
      if (isfinite(r) && r!=0.0) {
        double re=fabs(e/r); if (re>max_rel_abserr2) max_rel_abserr2=re;
      }
      if (s2 != GSL_SUCCESS) {
        fail2++;
        if (!first_status2) { first_status2=s2; first_y2=y; }
      }
    }

    printf("TOL %.1e FAIL_F1 %d FAIL_F2 %d FIRST_STATUS_F1 %d FIRST_Y_F1 %.17g FIRST_STATUS_F2 %d FIRST_Y_F2 %.17g MAX_REL_ABSERR_F1 %.17g MAX_REL_ABSERR_F2 %.17g\n",
      tol,fail1,fail2,first_status1,first_y1,first_status2,first_y2,max_rel_abserr1,max_rel_abserr2);
  }

  gsl_integration_workspace_free(w);
  return 0;
}
