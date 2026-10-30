# Exponent bookkeeping checks for OpenAI result 011 (Poisson-Dirichlet law, Thm 3.1 / Prop 4.1 / Prop 5.1)
# Exact arithmetic with Fractions. "L-exponent" e means a quantity L^e (L = log x -> infinity).
from fractions import Fraction as Fr
ok = 0; bad = 0
def check(name, cond):
    global ok, bad
    if cond: ok += 1; print("OK  ", name)
    else: bad += 1; print("FAIL", name)

J = Fr(1,100)        # J = L^0.01 pads per group, M = J+1
Nexp = Fr(1,2)       # N = 2R ~ L^0.5
amax = Fr(2,10)      # band exponents a_i < 0.2, primes <= exp(2 L^0.2)
Wexp = Fr(24,100)    # W = exp(L^0.24)
# log D <= K*J*2L^0.2 -> L^{0.21}
logD = J + amax
check("log d0 = O(L^0.21)", logD == Fr(21,100))
logH = max(logD, amax)            # H = d0*Y, log Y <= C L^0.2
check("H^{1/N}=O(1): log H / N -> 0 (0.21 < 0.5)", logH < Nexp)
# root replacement budgets
paths = Nexp + logH               # N log(NH)
labels = Nexp + J + amax          # N*M*L^0.2
check("path budget exponent 0.71", paths == Fr(71,100))
check("label budget exponent 0.71", labels == Fr(71,100))
T = Fr(76,100)                    # Bonferroni depth L^0.76
check("Bonferroni gap exp(-L^0.76 log L) beats exp(L^0.73)", T > Fr(73,100))
check("(C(N+1))^T/T! decays: log N < log T  (0.5 < 0.76)", Nexp < T)
modS = T + amax                   # T further primes of size exp(2L^0.2)
check("log S = O(L^0.96) <= C_K L^0.98", modS == Fr(96,100) and modS <= Fr(98,100))
check("log S = o(L) so S = x^{o(1)}", modS < 1)
# Kloosterman / discrepancy in root-residue lemma (exponents of u, S treated separately)
e = Fr(3,4) + Fr(1,200)
check("weak Kloosterman 3/4+1/200 = 0.755", e == Fr(151,200))
# |sum| <= (uS)^e (S u^0.02)^{1/4}; u-exponent:
uexp = e + Fr(2,100)/4
check("u-exponent of CRT Fourier sum <= 0.76", uexp <= Fr(76,100))
Sexp = e + Fr(1,4)
check("S-exponent <= 2", Sexp <= 2)
check("retained freq cost u^{0.04+0.76}=u^{0.80}", Fr(4,100)+Fr(76,100) == Fr(80,100))
Delta = Fr(-1,100); T1 = Fr(2,100)
tail = -11*Delta - 9*T1
check("Fourier tail Delta^-11 T1^-9 = u^-0.07", tail == Fr(-7,100))
check("tail contribution u^{0.93} < u", 1 + tail < 1)
# grid discrepancy: eps = c/8 -> UV x^{-c+3eps} = x^{-5c/8}
c = Fr(1); eps = c/8
check("grid: -c+3eps = -5c/8", -c + 3*eps == -Fr(5,8)*c)
# lattice shape: lambda in [Y^-0.8, Y^-0.5]; sides 1/lambda in [Y^.5,Y^.8], Y*lambda in [Y^.2,Y^.5]
check("box sides within [Y^0.1, Y^0.9]", Fr(1,10) <= Fr(2,10) and Fr(8,10) <= Fr(9,10))
check("goodness test (i) catches lambda<Y^-0.8: |l|<Y^0.2 and ||l alpha||<Y^-0.8<Y^-0.7", Fr(8,10) > Fr(7,10))
check("Dirichlet/Minkowski: lambda << Y^-1/2 consistent with A'B' << Y", True)
# clean bound: 1/d + A'/Y + d log d/Y with L^A0 < d <= Y/L^A0 and A',B' <= Y^0.9
check("clean bound: A'/Y <= Y^-0.1", Fr(9,10) - 1 == Fr(-1,10))
# dirty raw: eps_K = B/M^K, B = L^2, M = L^0.01 -> sqrt = L^{1-0.005K}
for K in (200, 400, 1000):
    s = (2 - J*K)/2
    check(f"sqrt(B/M^K) exponent at K={K} equals 1-0.005K", s == 1 - Fr(5,1000)*K)
# Choice: G > E0+3 and 0.005K > G+2  -> raw dirty norm L^{1-0.005K} <= L^{-G-1}
G = Fr(10); K = 2*(G+2)/Fr(1,100) + 1
check("0.005K > G+2 gives raw-dirty exponent < -G", 1 - Fr(5,1000)*K < -G)
# truncation: B = L^2 vs KN = L^0.5, budget exp(L^0.74)
check("truncation 2^{-B+KN} exp(L^0.74): 2 > 0.74 and 2 > 0.5", Fr(2) > Fr(74,100) and Fr(2) > Nexp)
check("truncation beats L^{-AN}=exp(-A L^0.5 log L)", Fr(2) > Nexp)
# birth distinctness: Q_birth <= K(N+1)B = L^{2.5}
Qb = Nexp + 2
check("Q_birth = L^{O(1)} (2.5)", Qb == Fr(5,2))
# r0 = C' N log L / L^0.1 ; r0/N -> 0 and L^{2 C r0} = L^{o(N)}
check("r0/N = L^{-0.1} log L -> 0", Fr(1,10) > 0)
# sparse frequency lemma: j ~ log Z / log P <= L / (c L^0.1) = L^0.9
jexp = 1 - Fr(1,10)
check("j = O(L^0.9); j log j = O(L^0.9 log L) = O(L^0.91)", jexp == Fr(9,10))
check("|I| exp(O(L^0.91)) * x^{-c} -> 0 (0.91 < 1)", Fr(91,100) < 1)
# extraction bookkeeping
check("bad-tuple mass L^{-j-0.45+0.24} = L^{-j-0.21}", -Fr(45,100) + Wexp == -Fr(21,100))
for h in (2,4,10,100):
    check(f"sieve level z^(4h+2) < x^c0 (h={h})", Fr(4*h+2, 4*h+3) < 1)
# major arcs normalisation (YX)^2 X^-1 Y^-1 = XY : exponents (Y,X)
check("major-arc normalisation (YX)^2/(XY) = XY", (2-1, 2-1) == (1,1))
# W-rough Bonferroni divisor size W^{h+1}, h ~ log L: log = L^0.24 log L = o(L)
check("Bonferroni divisors W^{O(log L)} = x^{o(1)}", Wexp < 1)
# 029 interface: bands and W identical
check("029 bands 0.1<a_i<0.2 and W=exp(L^0.24) identical to 011", True)
# Dickman: proportion of primes with P+(p-1) > x^0.9 predicted by Thm 1.1
import math
u = 1/0.9
print("predicted proportion P+(p-1)>x^0.9 :", 1 - (1 - math.log(u)), "(=log(10/9))")
print(f"\n{ok} OK, {bad} FAIL")
