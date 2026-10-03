import math
Cp=[4328422.243894628, 11606037.777232643, 19862015.365340475, 33329438.921102364]
Cm=[4297706.35960601, 11550362.014643177, 19824645.207705196, 33268689.17565131]

def dot(a,b): return sum(x*y for x,y in zip(a,b))
def norm(a): return math.sqrt(dot(a,a))
c=dot(Cp,Cm)/(norm(Cp)*norm(Cm))
theta=math.acos(max(-1,min(1,c)))
alpha=dot(Cm,Cp)/dot(Cm,Cm)
res=[Cp[i]-alpha*Cm[i] for i in range(4)]
rf=norm(res)/norm(Cp)
assert abs(rf-math.sin(theta))<1e-10
target=[1.,2.,3.,4.]
for C in (Cp,Cm):
    chi=[target[i]/C[i] for i in range(4)]
    assert all(abs(chi[i]*C[i]-target[i])<1e-12 for i in range(4))
print("E40_PASS", theta, rf)
