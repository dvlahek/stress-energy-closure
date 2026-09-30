#!/usr/bin/env python3
"""E35: positive spherical profile low-k theorem and toy counterexample, offline."""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUT=HERE/'e35_offline_lowk_profile_bound_result.json'
checks=0

def check(ok,label):
    global checks
    if not ok:
        raise AssertionError(label)
    checks+=1

def sinc(x):
    return math.sin(x)/x if x else 1.

def u(k,shells):
    return sum(float(w)*sinc(k*r) for w,r in shells)

def r2(shells):
    return sum(float(w)*r*r for w,r in shells)

def simpson(f,n=2048):
    if n%2: raise ValueError('even n required')
    h=1./n
    acc=f(0)+f(1)
    for i in range(1,n):
        acc+=(4 if i%2 else 2)*f(i*h)
    return acc*h/3

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

source=[(Fraction(1,5),0.3),(Fraction(1,2),0.8),(Fraction(3,10),1.0)]
test=[(Fraction(7,10),0.4),(Fraction(3,10),0.9)]
check(sum(w for w,r in source)==1,'positive source normalized')
check(sum(w for w,r in test)==1,'positive test normalized')
check(all(w>0 and r>=0 for w,r in source+test),'positive nonnegative radii')
r2s,r2t=r2(source),r2(test)
rows=[]
for k in (0.,0.05,0.1,0.173,0.5,1.,4.,8.):
    us,ut=u(k,source),u(k,test)
    point_error=abs(1-us*ut)
    bound=k*k*(r2s+r2t)/6
    clipped=min(2.,bound)
    check(abs(us)<=1+1e-14,'source normalized characteristic function')
    check(abs(ut)<=1+1e-14,'test normalized characteristic function')
    check(-1e-14<=1-us<=k*k*r2s/6+1e-14,'single-source second-moment bound')
    check(-1e-14<=1-ut<=k*k*r2t/6+1e-14,'single-test second-moment bound')
    check(point_error<=bound+1e-14,'source test product bound')
    check(point_error<=clipped+1e-14,'clipped product bound')
    rows.append({'k_dimensionless':k,'u_source':us,'u_test':ut,
                 'actual_abs_product_minus_point':point_error,
                 'second_moment_bound':bound,'clipped_bound':clipped})

# E34 profiles are mass-conserving, unit total mass; calculate actual integral
# and its rigorously integrated point-profile replacement bound.
R0,EPS=1.,0.1
K=lambda t:math.exp(t-1.)
g=lambda t:math.sin(math.pi*t)**2
histories={}
for sgn,name in ((+1,'outward'),(-1,'inward')):
    each=[]
    for k in (0.,0.1,0.5,1.,2.,4.,8.):
        rt=lambda t:R0+sgn*EPS*g(t)
        J=sinc(k*R0)*simpson(lambda t:K(t)*sinc(k*rt(t)))
        Jpoint=simpson(K)
        bd=(k*k/6)*simpson(lambda t:K(t)*(R0*R0+rt(t)**2))
        error=abs(J-Jpoint)
        check(error<=bd+1e-13,'retarded point-mass error bound')
        check(abs(sinc(k*R0))<=1+1e-14,'test form factor <=1')
        each.append({'k_dimensionless':k,'retarded_product':J,
                     'point_mass_retarded':Jpoint,
                     'actual_absolute_replacement_error':error,
                     'proven_integral_bound':bd})
    histories[name]=each

# Independent counterexample: same endpoint radii and velocities, but the
# interior shell expands strongly. Endpoint r2 cannot bound intermediate r2.
kbad=0.15
r_endpoint=1.
r_mid=10.
false_endpoint_bound=kbad*kbad*(r_endpoint*r_endpoint+r_endpoint*r_endpoint)/6
actual_mid_error=abs(1-sinc(kbad*r_endpoint)*sinc(kbad*r_mid))
check(actual_mid_error>false_endpoint_bound,'endpoint-only moment bound must FAIL')
check(abs(1+9*math.sin(0)**2-r_endpoint)<1e-15,'transient first endpoint')
check(abs(1+9*math.sin(math.pi)**2-r_endpoint)<1e-14,'transient second endpoint')
check(abs(1+9*math.sin(math.pi/2)**2-r_mid)<1e-14,'transient midpoint')
# The correct midpoint moment bound remains valid.
correct_mid_bound=kbad*kbad*(r_endpoint*r_endpoint+r_mid*r_mid)/6
check(actual_mid_error<=correct_mid_bound+1e-14,'correct time-local moment bound')
# Negative control: signed mass weights would not be a positive distribution.
bad=[(Fraction(3,2),0.3),(Fraction(-1,2),1.0)]
check(any(w<0 for w,r in bad),'signed-weight negative control rejected')
parents={name:sha(HERE/name) for name in (
    'E34_E35_OFFLINE_PROTOCOL.md',
    'e34_offline_causal_history_result.json',
    'E32_OFFLINE_RADIAL_PROFILE_ENCLOSURE.md',
    'E33_OFFLINE_M200C_PSEUDOEVOLUTION.md',
    'e29_20_halo_two_epoch_pair_catalog.json')}
report={
    'classification':'E35_OFFLINE_POSITIVE_SPHERICAL_LOW_K_SECOND_MOMENT_BOUND',
    'source_files_sha256':parents,
    'theorem':{
        'single_profile':'0<=1-u(k)<=k^2 * <r^2>/6 and |u(k)|<=1 for positive spherical normalized mass',
        'source_test':'|1-u_src(k,t)*u_test(k,T)|<=k^2*(<r^2>src(t)+<r^2>test(T))/6, also <=2',
        'retarded':'|A_profile(k)-A_point(k)| <= k^2/6 * integral dt |K(k,t)| m(t) [<r^2>src(t)+<r^2>test(T)]',
        'essential_condition':'A valid bound on the source second moment at every contributing time, not endpoints alone; positive m(t).',
        'no_claim':'Not the E28 actual force, F+/F- contrast, UV convergence, anisotropic full-profile bound, or higher-order Born bound.'
    },
    'toy_source_shells':[{'weight':float(w),'radius':r} for w,r in source],
    'toy_test_shells':[{'weight':float(w),'radius':r} for w,r in test],
    'source_second_moment':r2s,'test_second_moment':r2t,
    'dimensionless_k_grid':rows,
    'retarded_history_toy':histories,
    'endpoint_only_counterexample':{
        'endpoints_r':[r_endpoint,r_endpoint],
        'transient_midpoint_r':r_mid,
        'k_dimensionless':kbad,
        'invalid_bound_from_endpoint_moments':false_endpoint_bound,
        'actual_midpoint_abs_product_minus_point':actual_mid_error,
        'valid_bound_from_midpoint_moment':correct_mid_bound,
        'interpretation':'Endpoint mass profiles alone do not supply uniform temporal second moments.'
    },
    'qa_checks_passed':checks,
    'physical_total_halo_drag_certified':False,
    'observed_odd_read':False,
    'a03':'PHYSICAL_UNCERTIFIED','a04':'BLOCKED'
}
enc=json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+'\n'
if OUT.exists():
    if OUT.read_text(encoding='utf-8')!=enc:
        raise RuntimeError('IMMUTABLE_E35_OUTPUT_EXISTS_WITH_DIFFERENT_CONTENT')
else:
    OUT.write_text(enc,encoding='utf-8')
print('E35_OFFLINE_LOW_K_PROFILE_QA_PASS checks='+str(checks))
print('R2_SOURCE',f'{r2s:.9g}','R2_TEST',f'{r2t:.9g}')
print('ILLUSTRATIVE_K1_BOUND',f'{rows[5]["clipped_bound"]:.12g}',
      'ACTUAL',f'{rows[5]["actual_abs_product_minus_point"]:.12g}')
print('ENDPOINT_ONLY_COUNTEREXAMPLE actual=',f'{actual_mid_error:.12g}',
      'false_bound=',f'{false_endpoint_bound:.12g}')
print('OUTPUT',OUT)
