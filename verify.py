import sys, ast, re, math
from collections import Counter
def load():
    rows=[]
    for line in open('rm_U6.txt'):
        p=line.split(' ',4)
        if len(p)>=5 and p[2]=='True': rows.append((int(p[0]), ast.literal_eval(p[4])))
    for fn in ('found_103.txt','s3a.txt','s3b.txt','found_105.txt','found_s5.txt'):
        try:
            for line in open(fn):
                if line.startswith('FOUND'):
                    rows.append((int(re.search(r'm=(\d+)',line).group(1)), ast.literal_eval(re.search(r'A=(\[.*\])',line).group(1))))
        except FileNotFoundError: pass
    return rows
ok=True; worst=0
for m,A in load():
    k=len(A)
    R=[0]*m; D=[0]*m
    for a in A:
        for b in A: R[(a+b)%m]+=1; D[(a-b)%m]+=1
    assert min(R)>=1 and max(R)<=7
    N=Counter(R); t=sum(v%2 for v in R); E=sum((D[d]-3)*(D[d]-4) for d in range(1,m))
    # Theorem 1
    assert E+4*N[4]==8*(N[1]+N[7])+7*k-3*t-12
    assert k-3 <= N[4]+E/4 <= 3*k-3
    # Corollary: sizes
    assert math.sqrt(3*m)-2.5 < k < 2*math.sqrt(m)+4
    assert N[4]+t <= 4*k-3 and sum(1 for d in range(1,m) if D[d] not in (3,4)) <= 6*k-6
    assert N[6] >= (m-16*k+9)/4 and N[6] <= (m+4*k-5)/2 and N[2] >= (m-9*k+8)/2 and N[2] <= (3*m+16*k-9)/4
    # Autocorrelation proposition
    S=[1 if R[g]>=5 else 0 for g in range(m)]; T=[1 if (g!=0 and D[g]>=4) else 0 for g in range(m)]
    for x in range(1,m):
        cS=sum(S[n]*S[(n+x)%m] for n in range(m)); rT=sum(T[d]*T[(x-d)%m] for d in range(m))
        err=abs(16*cS-(2*k*k-5*m)-rT)
        worst=max(worst, err/k**1.5)
        assert err <= 151*k**1.5
    # exact convolution identity
    for x in (1,2,5):
        assert sum(R[n]*R[(n+x)%m] for n in range(m))==sum(D[d]*D[(x-d)%m] for d in range(m))
print("all checks passed on",len(load()),"certificates; worst autocorrelation error / k^1.5 =",round(worst,3))
