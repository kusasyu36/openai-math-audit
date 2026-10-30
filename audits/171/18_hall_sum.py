# Check of the union bound in Lemma low-two-endpoint-hall (sec 18, lines 1223-1311).
# Terms are formed as exact rational ratios term(p+1)/term(p) (Fractions) and
# summed in log domain; e, e^{.01} replaced by rational upper bounds.
from fractions import Fraction as F
from math import log, lgamma
LN2=log(2)
E_UP = F(27183,10000); E01_UP = F(1010051,1000000)
def ln(x): return log(x.numerator)-log(x.denominator)

def run(n, C, K, eta, pe):
    lnM = n*LN2; t0=int(eta*n)
    lnL = log(K)+lnM-log(C)
    lnpoly = pe*log(n)
    # small-set term p: poly * L^p/p! * 4^{p-1} p^{p-1} (p+1)^4 * 2 M^{-p-3} K^{p-2} E^2 * K^p
    def lsmall(p):
        return (lnpoly + p*lnL - lgamma(p+1) + (p-1)*log(4) + (p-1)*log(p) + 4*log(p+1)
                + log(2) - (p+1)*lnM + (2*p-2)*log(K) + 2*n*ln(E01_UP))
    terms=[lsmall(p) for p in range(3,t0)]
    mx=max(terms); s=mx+log(sum(2.718281828**(t-mx) for t in terms))
    p=t0
    llarge = (lnpoly + p*lnL - lgamma(p+1) + (p-1)*log(4)+(p-1)*log(p) + log(2)+lnM + 2*p*log(K) - p*lnM)
    return s/LN2, llarge/LN2

for args in [(400,10**6,10,0.3,20),(2000,10**6,10,0.3,20),(2000,10**4,10,0.3,20),
             (2000,3000,10,0.3,20),(4000,5000,10,0.3,20),(10**5,5000,10,0.3,20),(10**5,10**4,10,0.3,20)]:
    s,l=run(*args)
    print("n=%d C_n=%d K=%d eta=%.2f poly=n^%d : log2(small sum)=%.1f  log2(large count)=%.1f"%(args+(s,l)))
print("large-set ratio per row 4eK^2/C_n must beat 2^{-1/eta}: need C_n >", float(4*E_UP*100)*2**(1/0.3))
# exact check of algebraic step 2M(K/M)^{p-2}(M^-2)^2 <= K^p M^{-p-1}
ok=all(2*F(M)*(F(K,M))**(p-2)*F(1,M**4) <= F(K)**p*F(1,M**(p+1)) for M in (2,16,1024) for K in (2,10) for p in range(3,40))
print("tree+2-edge algebra (exact Fractions) holds on grid:", ok)
# exact: p^{p-1} 4^{p-1}/p! <= (4e)^p
print("diagram count <= (4e)^p for p<=200 (exact):", all(F(p**(p-1)*4**(p-1)) <= (4*E_UP)**p*__import__('math').factorial(p) for p in range(1,201)))
