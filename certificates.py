import math, itertools, random
from fractions import Fraction as Fr
def R(A,m):
    r=[0]*m
    for a in A:
        for b in A: r[(a+b)%m]+=1
    return r
def D(A,m):
    d=[0]*m
    for a in A:
        for b in A: d[(a-b)%m]+=1
    return d
def stats(A,m):
    r=R(A,m); d=D(A,m); k=len(A)
    N=[r.count(i) for i in range(max(r)+1)]
    t=sum(1 for x in r if x%2)
    E=sum((d[g]-3)*(d[g]-4) for g in range(1,m))
    return r,d,k,N,t,E
cert={101:[0,2,9,11,17,31,42,52,56,61,68,69,72,74,79,80,97,98,100],
103:[1,8,9,11,19,34,39,42,55,61,63,68,72,79,80,94,95,98,100],
104:[2,4,5,8,9,23,26,28,32,34,39,45,49,54,61,65,90,99,100],
105:[4,23,27,28,35,40,44,46,60,61,73,75,79,82,83,89,98,101,103],
107:[0,4,6,7,10,11,18,19,21,28,41,43,56,61,70,81,82,90,95],
108:[11,13,17,18,27,29,33,35,38,46,56,57,70,80,83,89,90,93,95]}
for m,A in cert.items():
    assert len(set(A))==len(A)==19 and all(0<=a<m for a in A)
    r,d,k,N,t,E=stats(A,m)
    N=N+[0]*(8-len(N))
    lhs=E+4*N[4]; rhs=8*(N[1]+N[7])+7*k-3*t-12
    print(m,"min/max R",min(r),max(r),"N",N[1:7],"t",t,"E",E,"id",lhs,rhs, "N4/m",round(N[4]/m,3))
# Z43 N4=0 set
A=[1,4,7,8,10,16,17,21,22,23,25,30,33,35,40]; m=43
r,d,k,N,t,E=stats(A,m); print("Z43",min(r),max(r),N, "E",E)
# m=37 example
A=[0,4,14,21,22,25,26,32,34,35]; r=R(A,37); print("Z37 ex",min(r),max(r))
