import sys, ast, re, cmath, math
sys.path.insert(0,'cyc')
from collections import Counter
from sidon import singer
def prof(m, A):
    k=len(A)
    R=Counter((a+b)%m for a in A for b in A); D=Counter((a-b)%m for a in A for b in A)
    Rv=[R.get(n,0) for n in range(m)]
    N=Counter(Rv); t=sum(v%2 for v in Rv)
    E=sum((D.get(d,0)-3)*(D.get(d,0)-4) for d in range(1,m))
    fmax=max(abs(sum(cmath.exp(2j*math.pi*a*x/m) for a in A)) for x in range(1,m))
    return dict(m=m,k=k,min=min(Rv),max=max(Rv),N=[N[i] for i in range(1,8)],t=t,E=E,
                lhs=E+4*N[4],rhs=8*(N[1]+N[7])+7*k-3*t-12,fmax=fmax,fb=math.sqrt(m/2+7*k-10))
print("== Singer one-exception sets A = L u (beta - L) ==")
for q in (7,11,13,17,19,23):
    M,L=singer(q); Ls=set(L); LL=set((a+b)%M for a in L for b in L)
    beta=next(b for b in range(M) if b not in LL)
    A=sorted(Ls|set((beta-l)%M for l in L)); assert len(A)==2*(q+1)
    R=Counter((a+b)%M for a in A for b in A)
    Rv=[R.get(n,0) for n in range(M)]
    off=[Rv[n] for n in range(M) if n!=beta]
    c=Counter(off)
    print(f"q={q} m={M} |A|={len(A)}: R(beta)={Rv[beta]}; off beta: min {min(off)} max {max(off)}; #R=4: {c[4]} ({c[4]/M:.3f} m), #R=2: {c[2]}, #R=6: {c[6]}; 3|A|-3={3*len(A)-3}")
print("== certificates with max R = 6 ==")
rows=[]
for fn in ('rm_U6.txt',):
    for line in open(fn):
        p=line.split(' ',4)
        if len(p)<5 or p[2]!='True': continue
        rows.append(prof(int(p[0]), ast.literal_eval(p[4])))
for fn in ('found_103.txt','s3a.txt','s3b.txt','found_105.txt','found_s5.txt'):
    try:
        for line in open(fn):
            if line.startswith('FOUND'):
                mm=int(re.search(r'm=(\d+)',line).group(1)); A=ast.literal_eval(re.search(r'A=(\[.*\])',line).group(1))
                rows.append(prof(mm,A))
    except FileNotFoundError: pass
bad=[r for r in rows if r['lhs']!=r['rhs'] or r['min']<1 or r['max']>6 or r['fmax']>r['fb']]
print("certificates:",len(rows),"violations of identity/Fourier bound:",len(bad))
print("N4/m range: %.3f .. %.3f"%(min(r['N'][3]/r['m'] for r in rows),max(r['N'][3]/r['m'] for r in rows)))
print("(N4+E/4) vs [k-3, 3k-3]:", all(r['k']-3 <= r['N'][3]+r['E']/4 <= 3*r['k']-3 for r in rows))
print("max |A^(chi)|/sqrt(m/2+7k-10): %.3f"%max(r['fmax']/r['fb'] for r in rows), " max |A^|/sqrt(m/2): %.3f"%max(r['fmax']/math.sqrt(r['m']/2) for r in rows))
for r in rows[-3:]: print(r)
