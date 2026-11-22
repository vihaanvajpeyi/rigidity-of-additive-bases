# Sidon sets in cyclic groups: Singer (M=q^2+q+1), Bose (M=q^2-1), Ruzsa (M=p(p-1)), q prime.
import itertools, math
def factor(n):
    f=[];d=2
    while d*d<=n:
        if n%d==0:
            f.append(d)
            while n%d==0:n//=d
        d+=1
    if n>1:f.append(n)
    return f
def prim_root(p):
    fs=factor(p-1)
    for g in range(2,p):
        if all(pow(g,(p-1)//f,p)!=1 for f in fs): return g
class Ext:
    # F_q[x]/(poly), poly monic of degree n given by low coefficients c[0..n-1]: x^n = -sum c_i x^i
    def __init__(s,q,c): s.q=q; s.c=c; s.n=len(c)
    def mul(s,a,b):
        q,n=s.q,s.n; r=[0]*(2*n-1)
        for i,ai in enumerate(a):
            if ai:
                for j,bj in enumerate(b): r[i+j]=(r[i+j]+ai*bj)%q
        for k in range(2*n-2,n-1,-1):
            t=r[k]
            if t:
                r[k]=0
                for i in range(n): r[k-n+i]=(r[k-n+i]-t*s.c[i])%q
        return tuple(r[:n])
    def pw(s,a,e):
        r=tuple([1]+[0]*(s.n-1))
        while e:
            if e&1:r=s.mul(r,a)
            a=s.mul(a,a);e>>=1
        return r
def prim_poly(q,n):
    N=q**n-1; fs=factor(N); one=tuple([1]+[0]*(n-1)); x=tuple([0,1]+[0]*(n-2))
    for c in itertools.product(range(q),repeat=n):
        if c[0]==0: continue
        F=Ext(q,list(c))
        if F.pw(x,N)!=one: continue
        if all(F.pw(x,N//f)!=one for f in fs): return F
def singer(q):
    F=prim_poly(q,3); M=q*q+q+1; x=(0,1,0); cur=(1,0,0); L=[]
    for i in range(M):
        if cur[2]==0: L.append(i)
        cur=F.mul(cur,x)
    assert len(L)==q+1; return M,sorted(L)
def bose(q):
    F=prim_poly(q,2); M=q*q-1; x=(0,1); cur=(1,0); L=[]
    for i in range(M):
        if cur[1]==1: L.append(i)
        cur=F.mul(cur,x)
    assert len(L)==q; return M,sorted(L)
def ruzsa(p):
    g=prim_root(p); M=p*(p-1); L=[]
    for i in range(p-1):
        y=pow(g,i,p)
        # CRT: n = i mod (p-1), n = y mod p
        n=(i*p*pow(p,-1,p-1) + y*(p-1)*pow(p-1,-1,p))%M
        L.append(n)
    return M,sorted(L)
def is_sidon(M,L):
    s=set()
    for i,a in enumerate(L):
        for b in L[i:]:
            t=(a+b)%M
            if t in s: return False
            s.add(t)
    return True
if __name__=="__main__":
    for q in (5,7,11,13):
        for nm,fn in (("singer",singer),("bose",bose),("ruzsa",ruzsa)):
            M,L=fn(q); print(nm,q,M,len(L),is_sidon(M,L))
