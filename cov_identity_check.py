import numpy as np
rng=np.random.default_rng(7)
Nmax=6000
A=[a for a in range(Nmax) if rng.random()<min(1,1.1/np.sqrt(a+1))]
f=np.zeros(Nmax); f[A]=1.0
delta=0.02
# phi coefficients
n=np.arange(Nmax); phi=-(delta/(1+delta))*(1/(1+delta))**n; phi[0]+=1
def poly_mul(a,b,L):
    s=len(a)+len(b)-1; S=1<<int(np.ceil(np.log2(s)))
    return np.real(np.fft.ifft(np.fft.fft(a,S)*np.fft.fft(b,S)))[:L]
for name,coef in (("F",f),("phiF",poly_mul(phi,f,Nmax))):
  for N in (150,400):
    x=1-1/N; rho=np.sqrt(x)
    c=coef.copy(); deg=len(c)
    psi=poly_mul(c,c,2*deg-1)
    nn=np.arange(deg); n2=np.arange(2*deg-1)
    w=x**nn; w2=x**n2
    m=[np.sum(nn**j*c**2*w) for j in range(3)]
    E=[np.sum(n2**j*psi**2*w2) for j in range(3)]
    lam=-m[1]/m[0]
    # LHS via grid
    Mg=1<<int(np.ceil(np.log2(4*deg+8)))
    cz=c*rho**nn
    Phi=np.fft.ifft(cz,Mg)*Mg
    G=np.fft.ifft(cz*(nn+lam),Mg)*Mg
    P2=np.abs(Phi)**2; G2=np.abs(G)**2
    lhs=np.mean((P2-m[0])*G2)
    # RHS lemma with fixed lambda: (theta/2+lam)^2 (E - m^2) - 1/2 (m theta^2 m - (theta m)^2)
    # theta^j of m^2: m^2 -> 2 m m1, 2 m1^2+2 m m2
    f0=E[0]-m[0]**2; f1=E[1]-2*m[0]*m[1]; f2=E[2]-2*m[1]**2-2*m[0]*m[2]
    rhs=0.25*f2+lam*f1+lam**2*f0-0.5*(m[0]*m[2]-m[1]**2)
    # corollary
    sig2=(m[0]*m[2]-m[1]**2)/m[0]**2
    y=E[0]/m[0]**2-1
    y2=E[2]/m[0]**2-4*E[1]*m[1]/m[0]**3-2*E[0]*m[2]/m[0]**3+6*E[0]*m[1]**2/m[0]**4
    lhs2=np.mean((P2-m[0])*(G2-m[0]*sig2))
    rhs2=-0.5*(1-y)*m[0]**2*sig2+0.25*m[0]**2*y2
    kappa=np.mean((G2-np.mean(G2))**2)/(m[0]*sig2)**2
    print(name,N,"lemma rel err %.2e"%(abs(lhs-rhs)/abs(rhs)),"cor rel err %.2e"%(abs(lhs2-rhs2)/abs(rhs2)),
          "intG2/(m sig2)-1 %.1e"%(np.mean(G2)/(m[0]*sig2)-1),"y %.3f kappa %.3f"%(y,kappa),
          "pointwise holds",(1-y)*sig2-0.5*y2 <= 2*sig2*np.sqrt(kappa*y)+1e-9)
