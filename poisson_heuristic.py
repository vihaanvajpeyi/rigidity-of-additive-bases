# Crude Poisson model of Section 8 ("Why small moduli are not representative").
# Unordered representation counts are independent Poisson(k^2/2m), required to lie in
# {1,2,3} off 2A and in {0,1,2} on 2A. Prints log of the expected number of k-subsets
# with 1<=R_A<=6 (maximised over k) and log(m*phi(m)), the number of affine maps.
import math
def lg(m,k):
    l=k*k/(2*m); P=lambda r: math.exp(-l)*l**r/math.factorial(r)
    pe=P(1)+P(2)+P(3); po=P(0)+P(1)+P(2)
    return math.lgamma(m+1)-math.lgamma(k+1)-math.lgamma(m-k+1)+(m-k)*math.log(pe)+k*math.log(po)
phi=lambda n: sum(1 for i in range(1,n+1) if math.gcd(i,n)==1)
for m in range(100,451,50):
    b=max(range(2,m),key=lambda k:lg(m,k))
    print(m,b,round(lg(m,b),2),round(math.log(m*phi(m)),2))
