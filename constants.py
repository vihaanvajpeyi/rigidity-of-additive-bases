import math
from fractions import Fraction as Fr
# Lambda table
for s in range(3,9):
    c=s*(s-3); L=min(Fr(i)+Fr(c,i) for i in range(1,50))
    l=s+Fr(1,2)-L/2; u=s+Fr(1,2)+L/2; sq=math.sqrt(s*s-3*s+0.25)
    print(s,L,float(l),round(s+.5-sq,4),float(u),round(s+.5+sq,4), float(L)<=2*sq+1e-12)
# Lambda_s < 2 sqrt for s in 4..2000
bad=[s for s in range(4,3000) if not min(i+s*(s-3)/i for i in range(1,s+1)) < 2*math.sqrt(s*s-3*s+.25)-1e-12]
print("strict fails",bad)
# c_k
ck=lambda k:2*k-1-math.sqrt((2*k-1)**2-1)
for k in (2,4,6,8,16,24): print("c",k,ck(k),(1-ck(k))**2-4*(k-1)*ck(k))
# perfect square claim: floor improves iff K square
diff=[]
for K in range(2,20001):
    Ks=2*(K//2); c=ck(Ks)
    a=(math.sqrt(K)-1)**2; b=(math.sqrt(K)-math.sqrt(1+c))**2
    fa=math.floor(a+1e-12); fb=math.floor(b+1e-12)
    sq=int(math.isqrt(K))**2==K
    if (fb<fa)!=sq: diff.append((K,a,b))
print("perfect-square claim exceptions up to 20000:",diff[:10],len(diff))
for K in (4,9,16,25):
    c=ck(2*(K//2)); print(K,(math.sqrt(K)-math.sqrt(1+c))**2,(math.sqrt(K)-1)**2)
# Cor 7
c=ck(6)
disc=4*(1+c)**2-64*c
u1=(-2*(1+c)-math.sqrt(disc))/8; u2=(-2*(1+c)+math.sqrt(disc))/8
print("u range",u1,u2,"M",4+2*u1,4+2*u2,"P2max",(1+c)**2/16-c)
import numpy as np
us=np.linspace(u1,u2,200001)
P1=0.5*((1+c)*(1+us/2)+us**2-us); P3=0.5*((1+c)*(1+us/2)+us**2+us)
print("P1min",P1.min(),"P3min",P3.min())
# threshold
print("7-4sqrt3",7-4*math.sqrt(3), "root of (1+c)^2=16c", (14-math.sqrt(192))/2, "kappa3 root", ck(4))
# Thm 1.5 P bounds
us=np.linspace(-0.5,0,100001); print("P2max thm",(-us/2-us**2).max(), "P3min",(0.5*(1+1.5*us+us**2)).min())
