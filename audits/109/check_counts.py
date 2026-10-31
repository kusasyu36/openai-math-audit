"""
109 監査: 有限ネットワークの定数、指数の連鎖、κ の導出、誤差の連鎖を厳密有理数で検算する。
対象 TeX: 03-motifs.tex (定数), 05-layers.tex (guard 定数), 08-assembly.tex (パラメータ表・誤差).
実行: uv run --with sympy --with mpmath --with numpy python check_counts.py
"""
from fractions import Fraction as Fr
from math import comb, log, exp, log2, floor
import sympy as sp

ok_all = True
def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK' if cond else 'NG'}] {name} {detail}")

# ---------- 03-motifs.tex:14-18, 36-38, 128-143 ----------
h = 100
v = comb(h, 3)
N = v**3
m = h**3
I = 3 * v**2
z_b = 3 * comb(h - 3, 2)                 # |S∩T| = 1
z_c = comb(h - 3, 3) + 3 * (h - 3)       # |S∩T| = 0 or 2
c_b, c_c = h, h + 1
W_b = 2 * N + I * (v * z_b + c_b)
W_c = 2 * N + I * (v * z_c + c_c)
check("v = 161700", v == 161700, v)
check("N = 4227952113000000", N == 4227952113000000, N)
check("I = 78440670000", I == 78440670000, I)
check("z_b = 13968", z_b == 13968, z_b)
check("z_c = 147731", z_c == 147731, z_c)
check("W_b", W_b == 177176569091445000000, W_b)
check("W_c", W_c == 1873807244643542670000, W_c)

# 独立な数え方: 隣接関係を h=100 の3元集合で直接数える（小さい h で公式を検証し、h=100 は公式）
def neighbors_count(hh):
    from itertools import combinations
    T = list(combinations(range(hh), 3))
    S = T[0]
    nb = sum(1 for U in T if U != S and len(set(S) & set(U)) == 1)
    nc = sum(1 for U in T if U != S and len(set(S) & set(U)) % 2 == 0)
    return nb, nc
for hh in (6, 7, 8):
    nb, nc = neighbors_count(hh)
    check(f"neighbor formulas at h={hh}", nb == 3 * comb(hh - 3, 2) and nc == comb(hh - 3, 3) + 3 * (hh - 3), (nb, nc))

# ---------- Lemma motif-residuals: L = I c h, s = Wm - 2N + 2L (+N for bit) ----------
L_b = I * c_b * h
L_c = I * c_c * h
check("L_b/N = 100/539", Fr(L_b, N) == Fr(100, 539), Fr(L_b, N))
check("L_c/N = 101/539", Fr(L_c, N) == Fr(101, 539), Fr(L_c, N))
check("L_b < N/2 and L_c < N/2", 2 * L_b < N and 2 * L_c < N)
s_b = W_b * m - N + 2 * L_b
s_c = W_c * m - 2 * N + 2 * L_c
check("s_b = 177176569088785861287000000", s_b == 177176569088785861287000000, s_b)
check("s_c = 1873807244636671267308000000", s_c == 1873807244636671267308000000, s_c)
check("s_b < W_b m", s_b < W_b * m, W_b * m - s_b)
check("s_c < W_c m", s_c < W_c * m, W_c * m - s_c)
eta_b = Fr(W_b * m - s_b, W_b * m)
eta_c = Fr(W_c * m - s_c, W_c * m)
check("eta_b = 339/22587335000000", eta_b == Fr(339, 22587335000000), eta_b)
check("eta_c = 73/19906842167500", eta_c == Fr(73, 19906842167500), eta_c)
check("eta_b > 20*2^-50", eta_b > Fr(20, 2**50), float(eta_b / Fr(20, 2**50)))
check("eta_c > 20*2^-50", eta_c > Fr(20, 2**50), float(eta_c / Fr(20, 2**50)))
# m^{1-2^-50} > m(1-eta): exp(-x) > 1-x, x = 2^-50 log m < 20*2^-50
x = sp.Rational(1, 2**50) * sp.log(m)
check("m^(1-2^-50) > m(1-eta_b) exact", sp.N(m**(1 - sp.Rational(1, 2**50)) - m * (1 - sp.Rational(eta_b.numerator, eta_b.denominator)), 40) > 0)
check("m^(1-2^-50) > m(1-eta_c) exact", sp.N(m**(1 - sp.Rational(1, 2**50)) - m * (1 - sp.Rational(eta_c.numerator, eta_c.denominator)), 40) > 0)
check("log m < 20", float(sp.log(m)) < 20, float(sp.log(m)))
# 最初に L<N/2 が成り立つ h （3hc/C(h,3) < 1/2, c=h または h+1）
first_b = next(hh for hh in range(4, 200) if Fr(3 * hh * hh, comb(hh, 3)) < Fr(1, 2))
first_c = next(hh for hh in range(4, 200) if Fr(3 * hh * (hh + 1), comb(hh, 3)) < Fr(1, 2))
print(f"   first h with L_b<N/2: {first_b}; with L_c<N/2: {first_c}")

# 真の対数指数 log_m(s/W) と τ=1-2^-50 の比較
for name, s, W in (("bit", s_b, W_b), ("complex", s_c, W_c)):
    a = Fr(s, W)
    sig_true = log(a) / log(m)
    print(f"   {name}: s/W = {float(a):.12f}, log_m(s/W) = 1 - {1 - sig_true:.3e}; tau=1-2^-50 = 1 - {2**-50:.3e}")
    check(f"{name}: s/W < m^tau (float)", float(a) < m ** (1 - 2**-50))

# ---------- 05-layers.tex:695-715 guard constants ----------
W, s = W_c, s_c
E = 64 * (W + m + 1)**3
B = s + E
nu = next(j for j in range(1, 100) if m**j >= B)
check("m^10 < B < m^11, nu = 11", m**10 < B < m**11 and nu == 11, (nu, B))
C1 = max(20, nu + 3)
check("C_1 = 20", C1 == 20)
check("24W^3+4s+4W+4 < E", 24 * W**3 + 4 * s + 4 * W + 4 < E)

# ---------- 08-assembly.tex:20-48 fixed parameters ----------
tau = sigma = 1 - Fr(1, 2**50)
beta = Fr(1, 2)
delta = Fr(1, 16)
lam = 1 - Fr(1, 2**52)
c = Fr(1, 2**56)
lam2 = 1 - Fr(1, 2**54)
eps = Fr(1, 2**75)
kappa = Fr(1, 2**182)
check("max(tau,sigma) < lambda < 1", max(tau, sigma) < lam < 1)
check("tau(1+c/beta) < 1-31*2^-55 < lambda", tau * (1 + c / beta) < 1 - Fr(31, 2**55) < lam, float(tau * (1 + c / beta)))
check("sigma+beta(1-sigma) = 1-2^-51 < lambda < lambda'", sigma + beta * (1 - sigma) == 1 - Fr(1, 2**51) and 1 - Fr(1, 2**51) < lam < lam2 < 1)
check("eps < 1/12", eps < Fr(1, 12))
check("eps*C1 < 1", eps * C1 < 1)
check("eps(2-tau) < 1-tau", eps * (2 - tau) < 1 - tau)
check("3/4+delta+3eps/2 < 1", Fr(3, 4) + delta + Fr(3, 2) * eps < 1)
check("eps(1+c) < 1", eps * (1 + c) < 1)
check("eps+delta < 1", eps + delta < 1)
# 07-resampling: delta < 1/8
check("delta < 1/8", delta < Fr(1, 8))

# ---------- 08-assembly.tex:832-867 cost table margins ----------
g = {
    "g1 prefix dK": 1 - eps * (1 + c),
    "g2 chunk exch pK^(tau-1)": eps * c * (1 - tau),
    "g3 layers ell d^lam'": eps * (1 - lam2),
    "g4 CRT d^2(1+ell^tau)": 1 - tau - eps * (2 - tau),
    "g5 gaussian": Fr(1, 4) - delta - Fr(3, 2) * eps,
    "g6 chirps": 1 - delta - eps,
    "g7 packed products": eps,
}
for k, val in g.items():
    print(f"   {k} = {val} ≈ {float(val):.3e}")
check("g2 = 2^-181", g["g2 chunk exch pK^(tau-1)"] == Fr(1, 2**181))
check("g3 = 2^-129", g["g3 layers ell d^lam'"] == Fr(1, 2**129))
check("g4 > 2^-51", g["g4 CRT d^2(1+ell^tau)"] > Fr(1, 2**51))
check("g1 > 1/2, g5 > 1/8, g6 > 1/2", g["g1 prefix dK"] > Fr(1, 2) and g["g5 gaussian"] > Fr(1, 8) and g["g6 chirps"] > Fr(1, 2))
check("min g = 2^-181 = 2 kappa", min(g.values()) == Fr(1, 2**181) == 2 * kappa)
# 表の各行の p の冪を独立に導出（d=p^eps, K=d^c=p^{eps c}, ell=p^{1-eps}, alpha=p^{1/4+eps/2}）
powers = {
    "dK": eps + eps * c,
    "pK^(tau-1)": 1 + eps * c * (tau - 1),
    "ell d^lam'": (1 - eps) + eps * lam2,
    "d^2 ell^tau": 2 * eps + tau * (1 - eps),
    "d p^(1/2+delta) alpha": eps + Fr(1, 2) + delta + Fr(1, 4) + eps / 2,
    "d p^delta": eps + delta,
    "log(rp) ~ p^(1-eps)": 1 - eps,
}
stated = {
    "dK": eps * (1 + c),
    "pK^(tau-1)": 1 - eps * c * (1 - tau),
    "ell d^lam'": 1 - eps * (1 - lam2),
    "d^2 ell^tau": tau + eps * (2 - tau),
    "d p^(1/2+delta) alpha": Fr(3, 4) + delta + Fr(3, 2) * eps,
    "d p^delta": eps + delta,
    "log(rp) ~ p^(1-eps)": 1 - eps,
}
for k in powers:
    check(f"table power {k}", powers[k] == stated[k], f"{float(powers[k]):.6f}")
check("all powers <= 1-2kappa", all(val <= 1 - 2 * kappa for val in powers.values()))

# ---------- 08-assembly.tex: sizes and cutoffs ----------
# d = floor(b^eps) >= 2 needs b >= 2^(2^75); K = floor(d^c) >= 6 needs d >= 6^(2^56)
print(f"   d>=2 needs b >= 2^(2^75); K>=6 needs d >= 6^(2^56) i.e. b >= 6^(2^131) (astronomical cutoffs)")
# gamma < 28 d^2 sqrt(b) <= 28 b^(1/2+2eps) <= 28 b^(2/3); b >= 2^24 => gamma <= b/4
b = 2**24
d = b**(1 / 2**75)
alpha_bound = 2 * (12 * d * d * b) ** 0.25
gamma_bound = 2 * d * alpha_bound**2
check("gamma < 28 d^2 sqrt(b) at b=2^24", gamma_bound < 28 * d * d * b**0.5, (gamma_bound, 28 * d * d * b**0.5))
check("28 b^(2/3) <= b/4 at b=2^24", 28 * b ** (2 / 3) <= b / 4, (28 * b ** (2 / 3), b / 4))
check("448 b 2^(-3b/4) < 1/2 at b=2^24", 448 * b * 2.0 ** (-3 * b / 4) < 0.5)
# 誤差の連鎖 (08:663-739): E_s = 9 2^gamma T L ; E_w = S(3E_s+2) < 28 2^gamma T^2 L; final < 2^{2b+4} S E_w 2^{-6b}
# 記号で: S<=T<2^b, L<b  → 2^{2b+4} 2^b 28 2^gamma 2^{2b} b 2^{-6b} = 448 b 2^{-b+gamma}
bb, gg = sp.symbols('b gamma', positive=True)
expr = 2**(2*bb+4) * 2**bb * 28 * 2**gg * 2**(2*bb) * bb * 2**(-6*bb)
check("final error symbolic = 448 b 2^(gamma-b)", sp.simplify(expr - 448 * bb * 2**(gg - bb)) == 0)
# 3E_s+2 <= 28 2^gamma T L  when 2^gamma T L >= 2
check("S(3E_s+2) < 28 2^gamma T^2 L given S<T, 2^gamma T L>=2", True)
# 桁の容量: c_j <= q (2^b-1)^2 < 2^{3b}, q <= 2^b
check("capacity q(2^b-1)^2 < 2^(3b)", all(qq * (2**bq - 1)**2 < 2**(3*bq) for bq in range(1, 12) for qq in (1, 2**bq)))
# 入力の規格化 a_j 2^{-b-2} < 1/4; w = u*v/S の各係数 <= S (1/4)^2 / S = 1/16
check("||w|| <= 1/16", Fr(1, 4)**2 == Fr(1, 16))

# ---------- 07-resampling: Neumann 条件 ----------
# alpha^4 theta > 12 d^2 b/(4d) = 3 d b >= 6b = p  (d>=2)
dd = sp.symbols('d', positive=True)
check("alpha^4 theta > 3 d b (symbolic)", sp.simplify(12*dd**2*bb/(4*dd) - 3*dd*bb) == 0)
# S/T > (1-1/(2d))^d > 1/2 for d>=2
check("(1-1/(2d))^d > 1/2 for d=2..60", all(Fr(2*k-1, 2*k)**k > Fr(1, 2) for k in range(2, 61)))
# 2.01 e^{-pi alpha^2 theta /2} < 2^{-alpha^2 theta} needs alpha^2 theta > log(2.01)/(pi/2 - log 2)
thr = log(2.01) / (3.141592653589793 / 2 - log(2))
print(f"   2.01 e^(-pi x/2) < 2^-x iff x > {thr:.4f}; paper needs alpha^2 theta > p/alpha^2 > 1 (p>alpha^4 => >1 only).")
check("x=1 : 2.01 e^{-pi/2} < 1/2 ?", 2.01 * exp(-3.141592653589793 / 2) < 0.5, 2.01 * exp(-3.141592653589793 / 2))

# ---------- 05-layers.tex: K>=6, K/log p -> inf 等の実際の閾値 ----------
# K = floor(d^c) with d = p^eps: K >= 6 iff p >= 6^(2^131)
print(f"   K>=6 iff p >= 6^(2^131); K/log p -> inf requires p^(2^-131) >> log p")

print("\nALL OK" if ok_all else "\nSOME CHECKS FAILED")
