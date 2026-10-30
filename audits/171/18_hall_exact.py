# EXACT rational evaluation (Fractions) of the Hall union bound of sec 18, lines 1253-1299,
# at modest n. e, e^{.01} replaced by rational upper bounds. Generic constant K taken = 10,
# e^{O(r)} slack taken = n^20, M = 2^n, L = floor(K M / C_n), t0 = floor(0.3 n).
from fractions import Fraction as F
from math import factorial
E01 = F(1010051,1000000)  # > e^{0.01}
def small_term(n,C,K,p,poly):
    M=2**n; L=(K*M)//C
    return poly*F(L**p,factorial(p))*4**(p-1)*p**(p-1)*(p+1)**4*2*F(K**(p-2),M**(p+1))*E01**(2*n)*K**p
def large_term(n,C,K,p,poly):
    M=2**n; L=(K*M)//C
    return poly*F(L**p,factorial(p))*4**(p-1)*p**(p-1)*2*M*F(K**p,M**p)*K**p
def lg2(x): return x.numerator.bit_length()-x.denominator.bit_length()
K=10
for n in (200,300):
    poly=n**20; t0=int(0.3*n)
    for C in (10**3,10**4,10**5,10**6):
        s=sum(small_term(n,C,K,p,poly) for p in range(3,t0))
        l=large_term(n,C,K,t0,poly)
        # exact ratio check for geometric tail of small-set series
        rmax=max(small_term(n,C,K,p+1,1)/small_term(n,C,K,p,1) for p in range(3,t0-1))
        print(f"n={n} C_n={C}: small-sum<1? {s<1} (log2~{lg2(s)}), large<1? {l<1} (log2~{lg2(l)}), max ratio term(p+1)/term(p) = {float(rmax):.4f}")
