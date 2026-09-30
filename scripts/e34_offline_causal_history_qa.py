#!/usr/bin/env python3
"""E34: fully offline mathematical witness. No user data arrays or external deps."""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / 'e34_offline_causal_history_result.json'
R0, EPS, T, KNUM = 1.0, 0.1, 1.0, 1.0
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def sinc(x):
    return math.sin(x)/x if x else 1.0


def g(t):
    return math.sin(math.pi*t/T)**2


def radius(t, sgn):
    return R0 + sgn*EPS*g(t)


def vr(t, sgn):
    return sgn*EPS*math.pi/T*math.sin(2*math.pi*t/T)


def kern(t):
    return math.exp(-(T-t)/T)


def simpson(func, panels):
    if panels <= 0 or panels % 2:
        raise ValueError('Simpson panels must be positive and even')
    step=T/panels
    total=func(0.) + func(T)
    for i in range(1, panels):
        total += (4 if i % 2 else 2)*func(i*step)
    return total*step/3


def J(sgn, panels):
    return sinc(KNUM*R0)*simpson(lambda t: kern(t)*sinc(KNUM*radius(t,sgn)), panels)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for sign in (+1,-1):
    for t in (0., T):
        check(abs(radius(t,sign)-R0)<1.e-14, 'same endpoint radius')
        check(abs(vr(t,sign))<1.e-14, 'zero endpoint radial velocity')
    check(radius(T/2,sign)>0, 'positive shell radius')
    check(abs(sinc(KNUM*radius(0.,sign))-sinc(KNUM*R0))<1.e-14, 'same endpoint profile')
    # Weak radial continuity: d/dt E[phi(r)] == E[phi\'(r) vr].
    for t in (0.2,0.38,0.7):
        for power in (1,2,3):
            h=1.e-5
            lhs=(radius(t+h,sign)**power-radius(t-h,sign)**power)/(2*h)
            rhs=power*radius(t,sign)**(power-1)*vr(t,sign)
            check(abs(lhs-rhs)<1.e-8, 'shell weak continuity for polynomial test')

values = {str(n):{'J_plus':J(+1,n),'J_minus':J(-1,n),'signed_gap_minus_plus':J(-1,n)-J(+1,n)}
          for n in (64,128,256,512,1024)}
for n in (64,128,256,512,1024):
    v=values[str(n)]
    check(v['J_minus']>v['J_plus'], 'strict causal gap')
    check(v['signed_gap_minus_plus']>0, 'gap positive')
check(abs(values['512']['signed_gap_minus_plus']-values['1024']['signed_gap_minus_plus'])<1.e-10,
      'Simpson gap convergence 512 to 1024')
check(abs(values['256']['J_plus']-values['1024']['J_plus'])<1.e-10,
      'Simpson plus convergence 256 to 1024')
check(abs(values['256']['J_minus']-values['1024']['J_minus'])<1.e-10,
      'Simpson minus convergence 256 to 1024')
for t in (0.05,0.25,0.50,0.75,0.95):
    check(sinc(KNUM*radius(t,+1)) < sinc(KNUM*radius(t,-1)), 'sinc sign over shell range')
# Lipschitz bound for each u history: |sinc'(x)|<=1/2, |dr/dt|<=eps*pi/T.
L=KNUM*EPS*math.pi/(2*T)
interp_integral=simpson(lambda t:kern(t)*2*L*t*(T-t)/T, 1024)
u0=sinc(KNUM*R0)
endpoint_response=sinc(KNUM*R0)*u0*simpson(kern, 1024)
for s in (+1,-1):
    check(abs(J(s,1024)-endpoint_response) <= abs(sinc(KNUM*R0))*interp_integral+1.e-14,
          'temporal Lipschitz integrated error bound')
# Explicit negative control: incorrect direction of strict ordering must be rejected.
check(not (J(+1,1024)>J(-1,1024)), 'reverse-signed-gap negative control')
parents={name:sha(HERE/name) for name in (
    'E30_OFFLINE_WIND_SIGN_AND_TWO_TRACER_ESTIMAND.md',
    'E32_OFFLINE_RADIAL_PROFILE_ENCLOSURE.md',
    'E33_OFFLINE_M200C_PSEUDOEVOLUTION.md',
    'e29_20_halo_two_epoch_pair_catalog.json',
    'E34_E35_OFFLINE_PROTOCOL.md')}
report={
    'classification':'E34_OFFLINE_TOY_CAUSAL_PROFILE_HISTORY_NONIDENTIFIABILITY',
    'source_files_sha256':parents,
    'construction':{
        'time_domain':'dimensionless [0,1]', 'r0':R0, 'eps':EPS,
        'k_dimensionless':KNUM, 'r_plus':'1+0.1*sin(pi*t)^2',
        'r_minus':'1-0.1*sin(pi*t)^2',
        'source_and_test_mass':'same positive constant M',
        'endpoint_shell_radii':[R0,R0], 'endpoint_radial_velocities':[0,0],
        'source_u':'sinc(k*r(t))', 'common_test_u':u0,
        'illustrative_kernel':'exp(-(1-t)); not E28 physical kernel',
        'continuity':'distributional shell transport, identical mass and fixed center'
    },
    'simpson':values,
    'endpoint_only_prediction':endpoint_response,
    'conditional_lipschitz_L':L,
    'lipschitz_integral_bound_before_test_factor':interp_integral,
    'theorem':{
        'claim':'Endpoint positions, masses, radial velocities and profiles do not identify retarded source-profile integral.',
        'rigorous_sign':'On [0.9,1.1], sinc(k*r) decreases in radius for k=1; K>0 and g>0 in (0,1), so J_minus>J_plus.',
        'conditional_stability':'Given independently justified |du/dt|<=L, |u(t)-u_lin(t)| <= 2L*t*(T-t)/T and integrated error <=2L/T*integral |K|*t*(T-t)dt.',
        'signed_kernel_caveat':'This fixed toy K is positive; no numerical nonzero E28 actual kernel gap is claimed.'
    },
    'qa_checks_passed':checks,
    'physical_total_halo_drag_certified':False,
    'observed_odd_read':False,
    'a03':'PHYSICAL_UNCERTIFIED', 'a04':'BLOCKED'
}
encoded=json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+'\n'
if OUT.exists():
    if OUT.read_text(encoding='utf-8') != encoded:
        raise RuntimeError('IMMUTABLE_E34_OUTPUT_EXISTS_WITH_DIFFERENT_CONTENT')
else:
    OUT.write_text(encoded, encoding='utf-8')
print('E34_OFFLINE_CAUSAL_HISTORY_QA_PASS checks='+str(checks))
print('J_PLUS',f'{values["1024"]["J_plus"]:.15g}')
print('J_MINUS',f'{values["1024"]["J_minus"]:.15g}')
print('SIGNED_GAP_MINUS_PLUS',f'{values["1024"]["signed_gap_minus_plus"]:.15g}')
print('SIMSPON_512_1024_GAP_DIFF',f'{abs(values["512"]["signed_gap_minus_plus"]-values["1024"]["signed_gap_minus_plus"]):.3g}')
print('OUTPUT',OUT)
