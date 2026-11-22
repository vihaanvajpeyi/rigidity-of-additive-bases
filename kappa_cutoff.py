import numpy as np
for seed in (1,2,3):
  rng=np.random.default_rng(seed); L=80000; N=2000; x=1-1/N; delta=N**-0.6
  a=np.arange(L); f=(rng.random(L)<np.minimum(1,1.1/np.sqrt(np.maximum(a,1)))).astype(float)
  S=1<<19
  n=np.arange(L); phi=-(delta/(1+delta))*(1/(1+delta))**n; phi[0]+=1
  c=np.real(np.fft.ifft(np.fft.fft(f,S)*np.fft.fft(phi,S)))[:L]
  w=x**n; m0=np.sum(c**2*w); m1=np.sum(n*c**2*w); m2=np.sum(n**2*c**2*w)
  lam=-m1/m0; sig2=(m0*m2-m1**2)/m0**2
  cz=c*np.sqrt(x)**n; G=np.fft.ifft(cz*(n+lam),S)*S; G2=np.abs(G)**2
  P=np.fft.ifft(cz,S)*S; P2=np.abs(P)**2
  kap=np.mean((G2-G2.mean())**2)/(m0*sig2)**2; y=np.mean((P2-m0)**2)/m0**2
  print(seed,"kappa(phiF) %.3f  y %.3f"%(kap,y))
