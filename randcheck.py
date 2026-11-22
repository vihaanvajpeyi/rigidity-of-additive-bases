# Random verification of the identities/inequalities in general finite abelian groups.
import itertools, random, sys
from collections import Counter
random.seed(12345)
def group(mods):
    els=list(itertools.product(*[range(n) for n in mods]))
    idx={e:i for i,e in enumerate(els)}
    add=[[idx[tuple((a+b)%n for a,b,n in zip(x,y,mods))] for y in els] for x in els]
    neg=[idx[tuple((-a)%n for a,n in zip(x,mods))] for x in els]
    return len(els), add, neg
def stats(m, add, neg, A):
    R=[0]*m; D=[0]*m
    for a in A:
        ra=add[a]; na=add[a]
        for b in A:
            R[ra[b]]+=1; D[add[a][neg[b]]]+=1
    return R,D
def run(mods, smax, trials, want):
    m,add,neg=group(mods); zero=0
    found=0; viol=0
    for _ in range(trials):
        k=random.randint(max(2,int((1.6*m)**.5)), int(((2*smax+1)*m)**.5))
        A=random.sample(range(m),k)
        R,D=stats(m,add,neg,A)
        if min(R)<1 or max(R)>2*smax+1: continue
        found+=1
        N=Counter(R); t=sum(v%2 for v in R)
        # general inequality for all j
        for j in range(0,3*smax+3):
            if 2*(smax-j)*k*k+2*(smax+j)*k < (4*smax-j*(j+1))*m + j*(j+1): viol+=1
        if smax==3:
            E=sum((D[d]-3)*(D[d]-4) for d in range(m) if d!=zero)
            if E+4*N[4]!=8*(N[1]+N[7])+7*k-3*t-12: viol+=1
            if not (k-3 <= N[4]+E/4 <= 3*k-3): viol+=1
            if not (3*m+3-5*k <= k*k <= 4*m-10+7*k): viol+=1
        if found>=want: break
    return m,found,viol
tot=0
for mods in ([29],[31],[33],[36],[37],[40],[45],[6,6],[2,2,3,3],[4,8],[3,3,5],[5,7],[2,2,2,5],[7,7]):
    m,f,v=run(mods,3,400000,1500); tot+=f
    print("G=Z_"+"xZ_".join(map(str,mods)),"m=%d sets with 1<=R<=7: %d  violations: %d"%(m,f,v),flush=True)
print("total s=3:",tot)
for s in (4,5):
    for mods in ([40],[6,6],[60],[4,4,4]):
        m,f,v=run(mods,s,100000,1000)
        print("s=%d G=Z_"%s+"xZ_".join(map(str,mods)),"m=%d sets: %d violations: %d"%(m,f,v),flush=True)
