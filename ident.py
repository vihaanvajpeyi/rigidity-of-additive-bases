import numpy as np, sys
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 1)
Np=float(sys.argv[2]) if len(sys.argv)>2 else 2000.0
kappa=2.2
cut=int(40*Np)
# random 'basis-like' set: P(a in A)=kappa/(2 sqrt a)
a=np.arange(1,cut); A=np.concatenate([[0],a[rng.random(cut-1)<np.minimum(1,kappa/(2*np.sqrt(a)))]])
x=1-1/Np; rho=np.sqrt(x)
L=1<<int(np.ceil(np.log2(2*cut+10)))
c=np.zeros(L); c[A]=rho**A
F=np.fft.ifft(c)*L                 # F(rho e(j/L)) = sum c_a e(a j/L)
cz=np.zeros(L); cz[A]=A*rho**A
ZF=np.fft.ifft(cz)*L
# moments
w=x**A.astype(float); m0=w.sum(); m1=(A*w).sum(); m2=(A*A*w).sum()
lam=-m1/m0
G=ZF+lam*F
u=np.abs(F)**2; W=np.abs(G)**2
m=u.mean()
lhs=((u-m)*W).mean()
# R and E moments via convolution of indicator
ind=np.zeros(L); ind[A]=1
R=np.rint(np.real(np.fft.ifft(np.fft.fft(ind)**2))).astype(np.int64)[:2*cut]
n=np.arange(len(R)); xn=x**n
E0=(R**2*xn).sum(); E1=(n*R**2*xn).sum(); E2=(n*n*R**2*xn).sum()
D0=E0-m0**2; D1=E1-2*m0*m1; D2=E2-(2*m1**2+2*m0*m2)
rhs=0.25*D2+lam*D1+lam**2*D0 - 0.5*(m0*m2-m1**2)
print("m(circle)=%.6g  m(sum)=%.6g"%(m,m0))
print("lhs=%.10g rhs=%.10g relerr=%.2e"%(lhs,rhs,abs(lhs-rhs)/abs(rhs)))
print("E|G|^2=%.6g  m*sigma2=%.6g"%(W.mean(), m0*(m2/m0-(m1/m0)**2)))
print("E|G|^4/(E|G|^2)^2 = %.4f"%((W**2).mean()/W.mean()**2))
M=(1-x)*m0**2; Q=(1-x)*E0; print("M=%.4f Q=%.4f V=Q-M^2=%.4f  V/M=%.4f"%(M,Q,Q-M*M,(Q-M*M)/M))
print("D0(deficit)=%.6g   (theta/2+lam)^2 Delta=%.6g"%(0.5*(m0*m2-m1**2), 0.25*D2+lam*D1+lam**2*D0))
