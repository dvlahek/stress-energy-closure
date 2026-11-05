import math
T=1.0; f0=1.0; f1=1.2; d=(f1-f0)/T; L=.8
def ell(t): return f0+d*t
def up(t): return min(f0+L*t, f1+L*(T-t))
def lo(t): return max(f0-L*t, f1-L*(T-t))
def ep(t): return min((L-d)*t,(L+d)*(T-t))
def em(t): return min((L+d)*t,(L-d)*(T-t))
Bmax=(L*L-d*d)*T/(2*L)
n=0
for i in range(101):
    t=i/100
    assert up(t)>=ell(t)-1e-14>=lo(t)-1e-14; n+=1
    assert up(t)-ell(t)<=ep(t)+1e-12; n+=1
    assert ell(t)-lo(t)<=em(t)+1e-12; n+=1
    assert max(ep(t),em(t))<=Bmax+1e-12; n+=1
print("E37_PASS", n, "BMAX", Bmax)
