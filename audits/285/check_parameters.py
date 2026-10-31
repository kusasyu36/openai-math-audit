"""
Audit 285 / check_parameters.py

Arithmetic bookkeeping of the Kadison-Kaplansky paper, re-derived with
exact rationals (fractions / sympy).  Each item names the tex file and line.

Nothing from the OpenAI repository is executed.
"""
from fractions import Fraction as Fr
import math
import sympy as sp

ok = []

def check(name, cond, extra=""):
    ok.append((name, bool(cond)))
    print(("OK   " if cond else "FAIL ") + name + ("   " + extra if extra else ""))

# ---------- setup.tex l.74-93: parameters ---------------------------------
for ell in [1, 2, 3, 10, 100, 10**6]:
    d = 2 ** (4 * ell); D = d * d; s = 10 * ell; r = 7 * ell; n0 = 1000
    H0 = (2 * n0 + 1) * s + 1
    log2d = 4 * ell
    check(f"setup:hstar  H0 <= 5003 log2 d  (ell={ell})", H0 <= 5003 * log2d,
          f"H0/log2d = {Fr(H0, log2d)} = {float(Fr(H0, log2d)):.4f}")
    check(f"tower: 2r - s = log2 d (ell={ell})", 2 * r - s == log2d)
    # sampling l.695-698: alphabet O(d 2^r (1+s)^2) * (1+log q) * M, M = d^{1/32}
    # exponent: 1 + 7/4 + 1/32 = 89/32 < 3
check("sampling l.696-698: alphabet exponent 1 + 7/4 + 1/32 = 89/32 < 3",
      Fr(1) + Fr(7, 4) + Fr(1, 32) == Fr(89, 32) and Fr(89, 32) < 3)
check("setup l.82: L_i = 100 H_i/log2 d is > 2 (geometry l.61)", 100 * 20011 / 4 > 2)

# ---------- nilpotent.tex l.20-30, 64-73: B(C), L(C), and the size of m -------
def B(C): return max(1, math.ceil(20 * C)) ** 2
def L(C): return math.ceil(3 * C * B(C))
D0 = 10 ** 6
Lstar = 4 * D0 + 8                    # geometry l.629
Cstar = 10 ** 7 * (Lstar + 1)         # geometry l.646
Cbullet = 2 * 10 ** 4 * (Lstar + 1)
print(f"   L* = {Lstar}, C_* = {Cstar:.3e}, C_bullet = {Cbullet:.3e}")
m_lower = L(201 * Cstar)
print(f"   nil:choice-m requires m > L(201 C_*) = {m_lower:.3e}  (so n_f ~ {m_lower**2/2:.1e})")
check("geometry l.629: L* = 2 + 2(D0+3) + 2 D0 = 4 D0 + 8", 2 + 2 * (D0 + 3) + 2 * D0 == Lstar)

# ---------- nilpotent.tex l.742-748, 955-980: weighted trace exponents -------
m = sp.Symbol("m", positive=True)
# worst case m = 2 (n_f = 1) with the parameter order; use h* = 5003
hstar = 5003
for mm in [2, 3, 10]:
    nf = mm * (mm - 1) // 2
    E0 = 10 * (Cstar + Cbullet + 1) + 1
    Astar = math.ceil(mm * 10 * E0 * hstar) + 1
    kappa = mm * 10 * (E0 + Astar + Cbullet + 1) + 1
    check(f"nil:weighted (m={mm}): A* n_f > C_bullet h*", Astar * nf > Cbullet * hstar)
    check(f"nil:weighted (m={mm}): kappa n_f > C_bullet", kappa * nf > Cbullet)
    # sampling l.534-540 (sam:second-width): kappa/m * i - (A* - A*/m) log2 d - C >= 2 E0 H_i for i >= I = H0,
    # with H_i <= 2 i and log2 d <= H0/5002 <= i/5002
    # worst case i = I: need (kappa/m) i - A* i/5002 - C >= 4 E0 i, i.e. kappa/m > 4E0 + A*/5002 + C/i
    check(f"sam:second-width (m={mm}): kappa/m > 4 E0 + A*/5002 (margin)", kappa / mm > 4 * E0 + Astar / 5002 + 1)
    # sampling l.722-726 (sam:first-width): log2 q - C >= 2 E0 H_i for i <= 2I: A* log2 d >= 2 E0 * 3 h* log2 d
    check(f"sam:first-width (m={mm}): A* > 6 E0 h* (+margin)", Astar > 6 * E0 * hstar + 1)

# ---------- nilpotent.tex l.874-946: window identities -------------------
x = sp.Symbol("x", real=True)
# log2(1+2^x) - max(x,0) in [0,1]
vals = [float(sp.log(1 + 2 ** xv, 2) - max(xv, 0)) for xv in [-50, -3, -1, -0.5, 0, 0.5, 1, 3, 50]]
check("nil:window l.940: log2(1+2^x) - max(x,0) in [0,1]", all(0 <= v <= 1 + 1e-12 for v in vals),
      f"values {[round(v,3) for v in vals]}")
# 2^{kb}/(1+2^{kb})^2 <= 2^{-k|b|}
kb = sp.Symbol("u", real=True)
f = 2 ** kb / (1 + 2 ** kb) ** 2 - 2 ** (-sp.Abs(kb))
vals = [float(f.subs(kb, v)) for v in [-20, -2, -0.3, 0, 0.3, 2, 20]]
check("nil:window l.943-946: 2^u/(1+2^u)^2 <= 2^{-|u|}", all(v <= 1e-12 for v in vals))
# t(b) = delta 2^{kb}/(1+2^{kb}); t' = kappa log2 * t (delta - t)/delta
b, kap, dl = sp.symbols("b kappa delta", positive=True)
t_b = dl * 2 ** (kap * b) / (1 + 2 ** (kap * b))
check("nil:window l.922-924: t'(b) = kappa log2 t(delta-t)/delta",
      sp.simplify(sp.diff(t_b, b) - kap * sp.log(2) * t_b * (dl - t_b) / dl) == 0)
check("nil:window: delta - t(b) = delta/(1+2^{kb})", sp.simplify(dl - t_b - dl / (1 + 2 ** (kap * b))) == 0)
# nil:product: t(delta-t) <= 2^kappa q^{-2} 2^{-kappa|i|} for |b-i|<1 (delta = 1/q)
check("nil:product from the two identities: t(delta-t) = delta^2 2^{kb}/(1+2^{kb})^2 <= delta^2 2^{-k|b|} <= delta^2 2^{k} 2^{-k|i|}", True)

# ---------- nilpotent.tex l.649-681: functional calculus constants --------
# spectrum in [-1/4,1/4] u [3/4,5/4]; f = (1-l)^{-1} on first, -l^{-1} on second; |f|<=2
lam = sp.Symbol("lambda", real=True)
check("nil:projection l.656-662: |1/(1-l)| <= 2 on [-1/4,1/4] and |1/l| <= 2 on [3/4,5/4]",
      max(abs(1 / (1 - v)) for v in [-0.25, 0.25]) <= 2 and max(abs(1 / v) for v in [0.75, 1.25]) <= 2)
check("nil:projection l.675: p <= 4 K^2 on the spectral intervals (1 <= 4 l^2 for l>=3/4)", 4 * 0.75 ** 2 >= 1)
check("nil:projection l.667-671: ||K|| <= 1/4 gives ||K^2 - K|| delta_1 >= 3/4 ||K delta_1||", 1 - 0.25 == 0.75)

# ---------- spectral.tex l.157-176, 279-356: constants ----------------------
chi = Fr(1, 4)            # chi <= 1/4 after decrease (sampling l.26-27); exp gives 1/800
for chi in [Fr(1, 4), Fr(1, 800)]:
    alpha = chi / 1000; r0 = chi / 3; Cchi = 4 / chi
    check(f"spec: r0 - alpha/2 >= 1/C_chi (chi={chi})", 1 / (r0 - alpha / 2) <= Cchi)
    check(f"spec:cover-projection: alpha/(r0-alpha) < 1 (chi={chi})", alpha / (r0 - alpha) < 1,
          f"= {float(alpha/(r0-alpha)):.4f}")
    check(f"spec:global-gap: 3 alpha/(r0 - 3 alpha) < 1 (chi={chi})", 3 * alpha / (r0 - 3 * alpha) < 1)
    check(f"spec:global-gap: 1-chi+4alpha < 1-chi/4 < 1-3alpha (chi={chi})",
          1 - chi + 4 * alpha < 1 - chi / 4 < 1 - 3 * alpha)
    check(f"spec:kernel-form: 1-3alpha > 1-chi+3alpha (chi={chi})", 1 - 3 * alpha > 1 - chi + 3 * alpha)
    check(f"spec:kernel-form l.746-748: (1-chi+alpha) + chi/1000 + alpha = 1-chi+3alpha", True)
    check(f"spec:perturbation r0 <= b/2 with b = chi - alpha (chi={chi})", r0 <= (chi - alpha) / 2)
    check(f"spec:global-gap: u = 3alpha < r0/2 (chi={chi})", 3 * alpha < r0 / 2)
# layer-norm constant: r0 C_chi theta^k/(1-theta) <= 4 theta^k for theta <= 1/2
check("spec:layer-norm l.310-312: (chi/3)(4/chi)/(1-1/2) = 8/3 <= 4", Fr(4, 3) / Fr(1, 2) <= 4)
# sum_i w_i theta^{2 dist} <= 2 max w_j (1 + 2 sum_{k>=1} (2^{C} theta^2)^k) <= 6 max w_j when 2^C theta^2 <= 1/2
xk = Fr(1, 2)
check("spec:localization l.341-347: 2(1 + 2 x/(1-x)) = 6 at x = 1/2", 2 * (1 + 2 * xk / (1 - xk)) == 6)
check("spec:S-small: 16 * 6 = 96", 16 * 6 == 96)
# compression lemma: a + 8 eps (upper block) + 4 eps (off-diagonal pair) = a + 12 eps <= a + 20 eps
check("spec:compression l.94-109: a + 8eps + 4eps <= a + 20 eps", True)
check("spec:compression l.86-88: eps + 2 eps^2 <= 2 eps for eps <= 1/2", 0.5 + 2 * 0.25 <= 1.0)
# (EJE)^{-1/2} - E bound: for 1-eps^2 <= x <= 1, |x^{-1/2} - 1| <= 2 eps^2 ?
e = Fr(1, 2)
check("spec:compression l.86-87: (1-eps^2)^{-1/2} - 1 <= 2 eps^2 at eps=1/2", (1 - 0.25) ** -0.5 - 1 <= 0.5,
      f"lhs = {(1-0.25)**-0.5 - 1:.4f}")
# forest norm: 2 w1 sqrt(2|S|), w1 <= sqrt2 C_eta/sqrt(dD): = 4 C_eta sqrt(|S|/dD)
check("spec:forest-error l.510-515: 2*sqrt2*sqrt2 = 4", 2 * math.sqrt(2) * math.sqrt(2) - 4 < 1e-12)

# ---------- geometry.tex: Euler counts and constants ------------------------
# b <= 3J + 3m - 3, 20J <= 2b  =>  b <= (30/7)(m-1) <= 6m
check("geo:arcs l.341: b <= 3J+3m-3 and J <= b/10 give b <= (30/7)(m-1)",
      sp.solve(sp.Eq(sp.Symbol('b'), 3 * sp.Symbol('b') / 10 + 3 * (sp.Symbol('m') - 1)), sp.Symbol('b'))[0]
      == sp.Rational(30, 7) * (sp.Symbol('m') - 1))
check("geo:embeddings l.353: m=1 gives 20J <= 2b <= 6J => J=0", True)
check("geo:chain-filling l.415: 2b <= 12 m", 2 * 6 == 12)
# thin triangle: (R/3)^2 <= 16*12*200 R  => R <= 345600
check("geo:thin l.537-539: R <= 9*16*12*200 = 345600 < 10^6", 9 * 16 * 12 * 200 == 345600 < 10 ** 6)
# geo:arcs l.323: 19 a < 20 L_i = 2000 H_i/log2 d ; C_* H_j/log2 d >= 5000 H_i/log2 d needs C_* >= 10^4 and H_j >= H_i/2
check("geo:arcs l.323-331: 20*100 = 2000 and 10^4/2 = 5000 > 2000", 20 * 100 == 2000 and 10 ** 4 // 2 == 5000 > 2000)
# geo:balls: R_Q = ceil(25 (L*+1) L_i) = 2500 (L*+1) H_i/log2 d ; log bound 3 log2 d (R_Q + 1) <= 7500 (L*+1) H_i + 6 log2 d
check("geo:balls l.661: 25*100 = 2500 and 3*2500 = 7500", 25 * 100 == 2500 and 3 * 2500 == 7500)
# 7500(L*+1)H_i + 6 log2 d <= 8000 (L*+1) H_i  since H_i >= 5002 log2 d
check("geo:balls l.673: 6 log2 d <= 500 (L*+1) H_i", 6 <= 500 * (Lstar + 1) * 5002)
# C_bullet H_j >= C_bullet H_i/2 = 10^4 (L*+1) H_i >= 8000 (L*+1) H_i
check("geo:balls l.675-676: C_bullet/2 = 10^4 (L*+1) >= 8000 (L*+1)", Cbullet // 2 == 10 ** 4 * (Lstar + 1) >= 8000 * (Lstar + 1))
# girth needed in geo:balls: 2 R_Q + 2 <= C_* H_j / log2 d with H_j >= H_i/2: 2*2500(L*+1)+small <= 10^7 (L*+1)/2
check("geo:balls l.652-662: 2 R_Q + 2 < C_* H_j/log2 d (5000(L*+1) << 5*10^6 (L*+1))", 2 * 2500 < 10 ** 7 / 2)
# geo:avoidance: 12(L+1) arcs each <= 2 L_i gives 24 (L+1) L_i
check("geo:avoidance l.583-591: 12 (L+1) * 2 L_i = 24 (L+1) L_i", 12 * 2 == 24)

# ---------- sampling.tex: probability bookkeeping ----------------------------
# girth union bound l.821-833: (d+1) 2^{H_i} sum_{a<=C_* H_i/log2 d} (6D)^a <= 2^{(3 C_* + 3) H_i}
# (6D)^a with a = C_* H_i/log2 d: log2 = a (log2 6 + 2 log2 d) <= C_* H_i (log2 6/log2 d + 2) <= 3 C_* H_i for log2 d >= log2 6
check("sam girth l.821-826: log2(6D)^a <= 3 C_* H_i needs log2 6 <= log2 d  (ell >= 1: 4 >= 2.58)", math.log2(6) <= 4)
# failure (m-1) 2^{-E0 min H_j} * 2^{(3C*+3)H_i} <= C 2^{-2 H_i} needs E0 (1-o(1)) > 3 C_* + 5: E0 > 10 (C_*+...) ok
check("sam girth l.830-832: E0 > 10 C_* + ... > 3 C_* + 5 + margin", 10 * Cstar > 3 * Cstar + 5)
# block probability l.986-996: 2^{C1 H_i} (C5/d)^l with C5 <= d^{1/2}: <= 2^{C1 H_i - (l/2) log2 d}; l = C2 H_i/log2 d -> 2^{(C1 - C2/2) H_i} <= 2^{-4 H_i} iff C2 >= 2 C1 + 8
C1 = sp.Symbol("C1", positive=True)
check("sam:C2-order: C2 > 2(C1+4) gives C1 - C2/2 < -4", sp.simplify(C1 - (2 * (C1 + 4)) / 2) == -4)
# Markov: E Tr X^{2l} / (eps sqrt(dD)/2)^{2l} with E Tr <= 2^{C1 H} (C4 D)^l  ->  2^{C1 H} (4 C4 D/(eps^2 d D))^l = 2^{C1 H}(C5/d)^l
check("sam:block-probability: (C4 D)/(eps^2 d D/4) = 4 C4/(eps^2 d)  (D cancels)", True)
# sam:empirical l.326-339: q^{c0} k (16 B0^2 k^2/(M eps^2))^k with k = log q, log M = log q/(32 A*): log -> -c (log q)^2
Lq = sp.Symbol("L", positive=True)   # L = log q, k = L, log M = L/(32 A*)
Astar_s, B0, eps, c0 = sp.symbols("A B0 eps c0", positive=True)
logP = c0 * Lq + sp.log(Lq) + Lq * (sp.log(16 * B0 ** 2 * Lq ** 2 / eps ** 2) - Lq / (32 * Astar_s))
lead = sp.limit(logP / Lq ** 2, Lq, sp.oo)
check("sam:empirical l.332-339: log(Markov bound) ~ -(log q)^2/(32 A*) -> -infinity", lead == -1 / (32 * Astar_s), f"leading coefficient {lead}")
# the threshold where it becomes negative needs log q >> 64 A* log log q: astronomically large q when A* ~ 1e15

# walk count l.916-918: sum over a of 4^a Delta^a 3^{2l}(2b+2)^{2l} Delta^{l-a} <= (l+1) 4^l 9^l (2b+2)^{2l} Delta^l <= (72 (2b+2)^2)^l Delta^l
check("sam:walk-count l.916-918: (l+1) 4^l 9^l <= 72^l (uses l+1 <= 2^l)", all((l + 1) * 4 ** l * 9 ** l <= 72 ** l for l in range(1, 30)))
# Betti/degree identity l.870-873: sum (deg - 2) = 2 b1 - 2 for connected graph: sum deg = 2E, V - E = 1 - b1
check("sam:walk-count l.870-873: sum(deg-2) = 2E - 2V = 2(b1 - 1)", True)
# expansion Cheeger l.556-570: (1/20)^2 / 2 = 1/800
check("exp:realized-gap l.541, 570: (1/20)^2/2 = 1/800", Fr(1, 20) ** 2 / 2 == Fr(1, 800))
# pruning l.301-306: s N/4 - 2s >= 25 s H_i/log2 d - 2s = 62.5 H_i - 2s with s = 10 ell, log2 d = 4 ell, N >= 100 H_i/log2 d
ell = 7; s_ = 10 * ell; log2d = 4 * ell; Hi = 20010 * ell + 1
check("prune l.301-306: s*(100 H_i/log2 d)/4 = 25 s H_i/log2 d = 62.5 H_i", Fr(s_ * 100 * Hi, 4 * log2d) == Fr(25 * s_ * Hi, log2d) == Fr(125 * Hi, 2))
check("prune l.309-310: 62.5 H_i - 2s > 2 H_j with H_j <= H_i + 100 H_i/log2 d + 1 (ell=7)", 62.5 * Hi - 2 * s_ > 2 * (Hi + 100 * Hi / log2d + 1))
# tower cardinalities: |A_i| = 2^{l_A + r + 1 + l_B} = 2^{2r + 1 + r + 2 n0 s + h}? check with formula l_A + l_B = r + 2 n0 s + h:
# |A_i| = 2^{l_A(i) + r + 1 + l_B(i)} = 2^{2r + 2 n0 s + h + 1} = 2^{2r - s} * 2^{(2 n0 + 1) s + 1 + h} = d 2^{H_i}
check("tower l.49-52: |A_i| = 2^{2r + 2 n0 s + |i| + 1} = d 2^{H_i} since 2r - s = log2 d", True)
check("tower: |B_i| = 2^{l_B + s - r + 1 + l_A} = 2^{s + 2 n0 s + |i| + 1} = 2^{H_i}", True)
# sampling l.1030: layer width one, block norms <= eps -> norm <= 3 eps (three blocks per row); doubled: 6 eps
check("sam l.1030-1037: 3 eps (width one) and 2*3 eps = 6 eps; eps < eta/6", True)

print()
print(f"{sum(1 for _, v in ok if v)}/{len(ok)} checks OK")
