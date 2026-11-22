import numpy as np, sys
def run(A, Np, delta, K=None, verbose=True):
    A=np.asarray(A); cut=A.max()+1
    x=1-1/Np; rho=np.sqrt(x)
    L=1<<int(np.ceil(np.log2(4*cut+10)))
    ind=np.zeros(L); ind[A]=1
    # coefficients of Phi = phi*F, phi=1-(delta/(1+delta)) sum (z/(1+delta))^k
    q=1/(1+delta); S=np.zeros(L); s=0.0
    for nn in range(L):
        s=ind[nn]+q*s; S[nn]=s
    Phi=ind-(delta*q)*S
    n=np.arange(L); xn=x**n
    def mom(c2): return (c2*xn).sum(),(n*c2*xn).sum(),(n*n*c2*xn).sum()
    m0,m1,m2=mom(Phi**2)
    PP=np.real(np.fft.ifft(np.fft.fft(Phi)**2))   # coefficients of Phi^2 (L big enough to avoid wrap for relevant n)
    E0,E1,E2=mom(PP**2)
    lam=-m1/m0
    c=Phi*rho**n; cz=n*Phi*rho**n
    F=np.fft.ifft(c)*L; ZF=np.fft.ifft(cz)*L; G=ZF+lam*F
    u=np.abs(F)**2; W=np.abs(G)**2
    D0=E0-m0**2; D1=E1-2*m0*m1; D2=E2-(2*m1**2+2*m0*m2)
    deficit=0.5*(m0*m2-m1**2)
    rhs=0.25*D2+lam*D1+lam**2*D0-deficit
    cov=((u-u.mean())*W).mean()
    y=D0/m0**2
    varw=W.var(); Ew=W.mean()
    # Abel statistics of R for comparison
    R=np.rint(np.real(np.fft.ifft(np.fft.fft(ind)**2))).astype(np.int64)
    Mx=(1-x)*(R*xn).sum(); Qx=(1-x)*(R*R*xn).sum(); V=Qx-Mx*Mx
    if verbose:
        print("Np=%g delta=%.3g  m_Phi/F(x)=%.4f  cov=%.6g rhs=%.6g relerr=%.1e"%(Np,delta,m0/(x**A).sum(),cov,rhs,abs(cov-rhs)/abs(rhs)))
        print("  y_Phi=%.4f   V/M-1=%.4f   ratio (theta/2+lam)^2Delta / deficit = %.4f"%(y,V/Mx-1,(rhs+deficit)/deficit))
        print("  Var(w)/(Ew)^2=%.4f   E w^2/(Ew)^2=%.4f   CS check: deficit-term=%.4g <= sqrt(D0*Varw)=%.4g"%(varw/Ew**2,(W**2).mean()/Ew**2,-rhs,np.sqrt(D0*varw)))
    return dict(y=y,V=V,M=Mx,varw=varw/Ew**2)
if __name__=="__main__":
    rng=np.random.default_rng(int(sys.argv[1])); Np=float(sys.argv[2]); kappa=2.2
    cut=int(40*Np); a=np.arange(1,cut)
    A=np.concatenate([[0],a[rng.random(cut-1)<np.minimum(1,kappa/(2*np.sqrt(a)))]])
    for d in [Np**-0.5*3, Np**-0.6, Np**-0.75]:
        run(A,Np,d)
