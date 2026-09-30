#!/usr/bin/env python3
"""E36: positive kinetic fields with equal 0th-2nd moments; free-streaming witness."""
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'e36_offline_hidden_kinetic_initial_data_result.json'
EPS=0.1
A=math.exp(1.5)/4
C=-3*math.exp(-0.5)/4
N=0

def check(ok,label):
    global N
    if not ok:raise AssertionError(label)
    N+=1

def F0(v):return math.exp(-v*v/2)/math.sqrt(2*math.pi)
def h(v):return math.cos(v)-A*math.cos(2*v)+C
def f(sign,x,v,t):return F0(v)*(1+sign*EPS*math.cos(x-v*t)*h(v))
def gaussian_cos(q):return math.exp(-q*q/2)
def gaussian_v2_cos(q):return (1-q*q)*gaussian_cos(q)
def G(s):return -0.5*math.exp(-(s*s+1)/2)*(math.cosh(s)-1)**2

def g_direct(s):
    return .5*(gaussian_cos(1-s)+gaussian_cos(1+s)) - A/2*(gaussian_cos(2-s)+gaussian_cos(2+s)) + C*gaussian_cos(s)

def integrate(fun, panels=8000, vmax=10):
    a,b=-vmax,vmax
    step=(b-a)/panels
    total=fun(a)+fun(b)
    for i in range(1,panels):total+=(4 if i%2 else 2)*fun(a+i*step)
    return total*step/3

check(EPS*(1+A+abs(C))<1,'global positivity bound')
check(abs(gaussian_cos(1)-A*gaussian_cos(2)+C)<1.e-15,'analytic zeroth mode cancellation')
check(abs(gaussian_v2_cos(1)-A*gaussian_v2_cos(2)+C)<1.e-15,'analytic second mode cancellation')
check(abs(G(0))<1.e-15,'initial hidden mode zero density')
for s in (0.1,0.5,1.0,2.0):
    check(abs(G(s)-g_direct(s))<1.e-14,'closed-form Gaussian free-streaming transform')
    check(G(s)<0,'strict transient density mode sign')
    check(abs(G(-s)-G(s))<1.e-14,'even initial velocity modulation')
check(abs((G(0.001)/(0.001**4))-(-math.exp(-0.5)/8))<2.e-7,'fourth-order temporal onset')
for t in (0.,.2,.5,1.):
    for sign in (+1,-1):
        for x in (0.,math.pi/3,math.pi):
            for v in (-7.,-1.,0.,1.,7.):
                check(f(sign,x,v,t)>0,'pointwise positivity sample')
for sign in (+1,-1):
    for x in (0.,math.pi/3,math.pi):
        for moment,expected in ((0,1.),(1,0.),(2,1.)):
            q=integrate(lambda v:f(sign,x,v,0)*v**moment)
            check(abs(q-expected)<2.e-11,'equal initial 0,1,2 moments')
    for t in (0.,.5,1.):
        for x in (0.,math.pi/3,math.pi):
            rho=integrate(lambda v:f(sign,x,v,t))
            prediction=1+sign*EPS*math.cos(x)*G(t)
            check(abs(rho-prediction)<2.e-11,'free-streamed density numerical-analytic agreement')
    check(abs(integrate(lambda v:f(sign,0.,v,0)*v**2)-1)<2.e-11,
          'initial pressure same at x zero')

positive_lower_bound=1-EPS*(1+A+abs(C))
check(positive_lower_bound>0,'explicit global positive density')
check(abs((1+EPS*G(1))+(1-EPS*G(1))-2)<1e-14,'same baseline density and reversal')
check(abs((1+EPS*G(1))-(1-EPS*G(1)))>1e-3,'later densities distinguishable')
check(not (G(1)>0),'wrong evolution sign rejected')
parents={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in (
    'E36_OFFLINE_PROTOCOL.md',
    'E34_OFFLINE_TWO_SNAPSHOT_CAUSAL_HISTORY_NONIDENTIFIABILITY.md',
    'E35_OFFLINE_LOWK_PROFILE_ERROR_BOUND.md',
    'e34_offline_causal_history_result.json',
    'e35_offline_lowk_profile_bound_result.json')}
report={
 'classification':'E36_OFFLINE_TOY_HIDDEN_INITIAL_KINETIC_MODE_FREE_STREAMING',
 'parent_sha256':parents,
 'kinetic_setup':{'domain':'1D x periodic [0,2pi), v in R','background':'standard normal F0(v)',
   'amplitude_eps':EPS,'h':'cos(v)-exp(3/2)/4*cos(2v)-3*exp(-1/2)/4',
   'global_lower_density_factor_bound':positive_lower_bound,
   'initial_equal_moments':{'rho':1,'momentum':0,'second_velocity_moment':1},
   'initial_full_velocity_distribution_equal':False,
   'equation':'free collisionless streaming only: df/dt+v*dxf=0'},
 'mode_formula':'G(s)=-(1/2)*exp(-(s^2+1)/2)*(cosh(s)-1)^2',
 'density_at_t1_x0':{'plus':1+EPS*G(1),'minus':1-EPS*G(1),'signed_plus_minus':2*EPS*G(1)},
 'mode_at_t1':G(1),
 'fourth_order_coefficient_limit':-math.exp(-.5)/8,
 'qa_checks_passed':N,
 'scope_limit':'Not self-gravitating or Einstein-Vlasov constraint-matched; not halo forcing or eBOSS galaxy signal.',
 'observed_odd_read':False,
 'a03':'PHYSICAL_UNCERTIFIED','a04':'BLOCKED'}
enc=json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+'\n'
if OUT.exists():
 if OUT.read_text(encoding='utf-8')!=enc:raise RuntimeError('IMMUTABLE_E36_OUTPUT_DIFFERS')
else:OUT.write_text(enc,encoding='utf-8')
print('E36_OFFLINE_HIDDEN_KINETIC_MODE_PASS checks='+str(N))
print('GLOBAL_POSITIVITY_FACTOR_LOWER_BOUND',f'{positive_lower_bound:.12g}')
print('G(1)',f'{G(1):.15g}')
print('PLUS_MINUS_DENSITY_DIFFERENCE_t1_x0',f'{2*EPS*G(1):.15g}')
print('FOURTH_ORDER_COEFFICIENT',f'{-math.exp(-.5)/8:.15g}')
print('OUTPUT',OUT)
