k=.6
R20=.20; R21=.32; R2t=.18
d0=min(2.,k*k*R20/6); d1=min(2.,k*k*R21/6); dt=min(2.,k*k*R2t/6)
assert 0<=d0<=2 and 0<=d1<=2 and 0<=dt<=2
for i in range(101):
    s=i/100
    dl=(1-s)*d0+s*d1
    assert min(d0,d1)-1e-15 <= dl <= max(d0,d1)+1e-15
print("E38_PASS", d0, d1, dt)
