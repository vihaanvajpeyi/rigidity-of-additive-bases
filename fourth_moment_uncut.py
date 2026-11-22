import numpy as np
def run(N,seed):
  rng=np.random.default_rng(seed); L=40*N; x=1-1/N
  a=np.arange(L); f=(rng.random(L)<np.minimum(1,1.1/np.sqrt(np.maximum(a,1)))).astype(float)
  S=1<<int(np.ceil(np.log2(4*L)))
  n=a; w=x**n; m0=np.sum(f*w); m1=np.sum(n*f*w); m2=np.sum(n**2*f*w)
  lam=-m1/m0; sig2=(m0*m2-m1**2)/m0**2
  cz=f*np.sqrt(x)**n; G=np.fft.ifft(cz*(n+lam),S)*S; G2=np.abs(G)**2
  return np.mean(G2**2)/np.mean(G2)**2
for N in (250,1000,4000,16000):
  r=[run(N,s) for s in range(4)]
  print(N,"int|G_F|^4/(int|G_F|^2)^2:",np.round(r,3))
