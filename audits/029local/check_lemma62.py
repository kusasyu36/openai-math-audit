# -*- coding: utf-8 -*-
"""
029 Lemma 6.2 (Local effect of reflection), 04-reflection.tex:339-381,
at the prime p = (zeta12 - 2) of F = Q(zeta12), residue field F_13.

SCOPE (written before running).
  Tested:   the finite algebra that leads from the paper's stated local rules
            to the claims of the lemma, with every sum over F_13 enumerated
            and every equality decided exactly in Z[zeta_156][1/13].
  Assumed, NOT tested (inputs of the lemma):
            - the lift identity  w n(H) = n(1/H) a(H^-2) k_v   (03-theta.tex:114-123)
            - the torus law      C(p^e1 u, p^e2 v) = chi(u)^(2 e2) chi(v)^(-2 e1)
                                                     (04-reflection.tex:394-399)
            - the spherical values W(a(p^e u))       (03-theta.tex:147-157)
            - psi_p unramified, W(n(x) g) = psi_p(x) W(g)
  So this is a consistency check of the lemma given its local inputs.  It says
  nothing about the metaplectic cover, Prop 5.1, Lemma 6.3 or the theorem.

RELATION TO THE EARLIER CHECK.
  audits/029core/local_checks.py Part D already ran the same lemma on Z/13 in
  floating point with one character and one additive character.  New here:
  (i) the character and the prime are built from Q(zeta12) (residue symbol,
  generator, injectivity of mu_12); (ii) exact arithmetic; (iii) all 12
  nontrivial additive characters of F_13 and both orientations of mu_6 in C
  (we cannot pin down the paper's "allowed generator", so we run every choice);
  (iv) each claim is a separately named check; (v) negative controls.
  NOTE added after the first run: the 24 conventions are Galois conjugates of
  one another (check G), so for the exact equalities they carry the information
  of a single run.  The vanishing for e < 0 is a restatement of the rule
  W(a(p^e u)) = 0 and is not counted as a check.

CRITERIA (fixed before running).
  Equalities: exact (zero in Q(zeta_156)).  Inequalities: float, must hold with
  |lhs| <= rhs * (1 + 1e-12).  Negative controls N1, N2 must FAIL for every
  convention.  N3 is exploratory (reported, not a pass/fail criterion).
"""
import cmath, math, sys
from fractions import Fraction as Fr

P = 13
N = 156            # lcm(12, 13): zeta6 = z^26, e(1/13) = z^12
E_MAX = 9

# ---------------------------------------------------------------- polynomials
def polymul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for k, y in enumerate(b):
                r[i + k] += x * y
    return r
def polydiv_exact(a, b):          # b monic, low-to-high coefficients
    a = a[:]; q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1]
        q[i] = c
        for k, y in enumerate(b):
            a[i + k] -= c * y
    assert not any(a), "inexact division"
    return q
def mobius(n):
    m, r, p = n, 1, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
def cyclotomic(n):
    num, den = [1], [1]
    for d in range(1, n + 1):
        if n % d == 0:
            f = [-1] + [0] * (d - 1) + [1]
            mu = mobius(n // d)
            if mu == 1: num = polymul(num, f)
            elif mu == -1: den = polymul(den, f)
    return polydiv_exact(num, den)
PHI = cyclotomic(N)
DEG = len(PHI) - 1
assert DEG == 48
PHI12 = cyclotomic(12)
assert PHI12 == [1, 0, -1, 0, 1]

# ---------------------------------------------------------------- Q(zeta_156)
class El:
    """element of Q[Z/156] read in Q(zeta_156): sum v[k] z^k, v[k] in Q"""
    __slots__ = ("v",)
    def __init__(self, v): self.v = v
    @staticmethod
    def mono(k, c=1):
        v = [0] * N; v[k % N] = c; return El(v)
    @staticmethod
    def const(c): return El.mono(0, c)
    def __add__(s, o): return El([a + b for a, b in zip(s.v, o.v)])
    def __sub__(s, o): return El([a - b for a, b in zip(s.v, o.v)])
    def scale(s, c): return El([a * c for a in s.v])
    def __mul__(s, o):
        r = [0] * N
        for i, x in enumerate(s.v):
            if x:
                for k, y in enumerate(o.v):
                    if y: r[(i + k) % N] += x * y
        return El(r)
    def shift(s, k):               # multiply by z^k
        k %= N
        return El(s.v[-k:] + s.v[:-k]) if k else El(s.v[:])
    def conj(s): return El([s.v[(-k) % N] for k in range(N)])
    def is_zero(s):
        a = s.v[:]
        for i in range(N - 1, DEG - 1, -1):
            c = a[i]
            if c:
                for k, y in enumerate(PHI):
                    if y: a[i - DEG + k] -= c * y
        return not any(a)
    def __eq__(s, o): return (s - o).is_zero()
    def val(s): return sum(complex(float(c)) * cmath.exp(2j * cmath.pi * k / N) for k, c in enumerate(s.v) if c)
ZERO, ONE = El.const(0), El.const(1)

# ---------------------------------------------------------------- field data
results = []
def check(name, ok, extra=""):
    results.append(bool(ok))
    print(("OK   " if ok else "FAIL ") + name + (("   [" + extra + "]") if extra else ""))
    sys.stdout.flush()

print("== K0: the prime (zeta12 - 2) of Q(zeta12) ==")
check("Phi_12(2) = 13, so N(zeta12 - 2) = prod over embeddings (2 - zeta12^k) = 13 and (zeta12 - 2) is a prime of degree 1 containing 13",
      sum(c * 2 ** i for i, c in enumerate(PHI12)) == 13)
# numerical cross-check of the norm from the four embeddings
nrm = 1
for k in (1, 5, 7, 11): nrm *= (2 - cmath.exp(2j * cmath.pi * k / 12))
check("norm recomputed from the 4 complex embeddings = 13", abs(nrm - 13) < 1e-12, "%.12f" % nrm.real)
ords = sorted({pow(2, k, P) for k in range(12)})
check("zeta12 -> 2 mod p; 2 has order 12 mod 13, so mu_12 reduces injectively", len(ords) == 12 and pow(2, 12, P) == 1)
check("q = 13 = 1 mod 12 and chi_p(-1) = (-1)^((q-1)/6) = 1  (02-arithmetic.tex:406-408)", P % 12 == 1 and pow(P - 1, (P - 1) // 6, P) == 1)
# residue symbol: chi_p(x) is the element of mu_6 = <zeta12^2> congruent to x^((q-1)/6) = x^2; zeta12^2 -> 4
MU6 = {pow(4, k, P): k for k in range(6)}
check("mu_6 = <zeta12^2> reduces to {4^k}: 6 distinct residues, and x^2 lies in it for every unit x",
      len(MU6) == 6 and all(pow(x, 2, P) in MU6 for x in range(1, P)))
def chi_log(x):                    # chi_p(x) = (zeta12^2)^chi_log(x), from the residue symbol
    return MU6[pow(x, (P - 1) // 6, P)]
LOG2 = {pow(2, k, P): k for k in range(12)}
check("the residue symbol equals the character of 029core Part D (2^k -> zeta6^k): the character itself is not new",
      all(chi_log(x) == LOG2[x] % 6 for x in range(1, P)))
dphi = sum(i * c * 2 ** (i - 1) for i, c in enumerate(PHI12) if i) % P
c_gen = (-dphi) % P     # 1/(zeta-2) = -g(2)/13 with g(2) = Phi12'(2) mod 13  ->  v/p = -Phi12'(2) v / 13 mod Z_13
print("info: for the generator p = zeta12 - 2,  v/p = %d*v/13 mod Z_13 (Phi12'(2) = %d mod 13); with psi_p(x) = e(-+{x}) this is c = %d or %d below."
      % (c_gen, dphi, (-c_gen) % P, c_gen))
print("info: the paper's allowed generator may differ by a unit, so every c in F_13^x is run.")

# ---------------------------------------------------------------- the lemma, one convention
def run(o, c, verbose=False):
    """o = +-1: which primitive 6th root of unity in C the generator of mu_6 goes to.
       c in F_13^x: psi_p(v/p) = e(c v / 13)."""
    fails = {}
    def rec(key, ok, info=None):
        if not ok: fails.setdefault(key, []).append(info)
    def chi_idx(x, j):             # exponent of z for chi(x)^j, x a unit
        return (26 * o * j * chi_log(x % P)) % N
    def psi_idx(a): return (12 * c * (a % P)) % N
    inv = lambda x: pow(x, -1, P)
    # Gauss sums
    tau = {j: sum((El.mono(chi_idx(v, -j) + psi_idx(v)) for v in range(1, P)), ZERO) for j in range(1, 6)}
    for j in range(1, 6):
        rec("K1 tau*conj(tau) = 13", tau[j] * tau[j].conj() == El.const(13), j)
    tauinv = {j: tau[j].conj().scale(Fr(1, 13)) for j in range(1, 6)}   # valid once K1 holds
    G2 = sum((El.mono(chi_idx(x, 2) + psi_idx(x)) for x in range(1, P)), ZERO)   # sqrt(q) * g_2
    rec("K1 |g_2| = 1", G2 * G2.conj() == El.const(13), "g2")
    def alpha_pref(j): return tauinv[j] if j else El.const(Fr(-1, 13))
    def W(e):                      # W_e, exact
        if e < 0: return ZERO
        k, r = divmod(e, 3)
        if r == 0: return El.const(Fr(1, 13 ** k))
        if r == 1: return G2.scale(Fr(1, 13 ** (k + 1)))
        return ZERO
    # K2 mask identity, all t in F_13 (eq:reflection-mask-identity)
    for j in range(6):
        for t in range(P):
            if j:
                lhs = tauinv[j] * sum((El.mono(chi_idx(v, -j) + psi_idx(t * v)) for v in range(1, P)), ZERO)
                want = El.mono(chi_idx(t, j)) if t else ZERO
            else:
                lhs = El.const(Fr(P - 1, P)) + sum((El.mono(psi_idx(t * v), Fr(-1, P)) for v in range(1, P)), ZERO)
                want = ONE if t else ZERO
            rec("K2 mask identity", lhs == want, (j, t))
    # the digit sum from the stated rules
    def M_direct(j, e, t, coc_o=1, vpow=None):
        """sum over unit digits v of  chi(v)^-j [j != 0] * psi(p^(e-1) t / v) * C(p^(e-2) t, p^2 v^-2) * chi(t v^-2)^(2e)"""
        acc = [0] * N
        for v in range(1, P):
            vi = inv(v); u2 = vi * vi % P
            idx = chi_idx(v, -j) if j else 0
            if e == 0: idx += psi_idx(t * vi)                  # e >= 1: argument integral, psi = 1
            idx += coc_o * (chi_idx(t, 2 * 2) + chi_idx(u2, -2 * (e - 2)))   # torus law
            idx += chi_idx(t * u2, 2 * e)                      # unit rule of W
            acc[idx % N] += 1
        return El(acc)
    def R_direct(j, e, t, coc_o=1):
        if e < 0: return ZERO      # W(a(p^e ...)) = 0
        return alpha_pref(j) * W(e) * M_direct(j, e, t, coc_o)
    def R_formula(j, e, t, vexp=4):                            # eq:reflection-local-formula
        acc = [0] * N
        for v in range(1, P):
            idx = (chi_idx(v, -j) if j else 0) + chi_idx(v, vexp)
            if e == 0: idx += psi_idx(t * inv(v))
            acc[idx % N] += 1
        return (alpha_pref(j) * W(e) * El(acc)).shift(chi_idx(t, 4 + 2 * e))
    Q3 = sum((El.mono(chi_idx(x, 3) + psi_idx(x)) for x in range(1, P)), ZERO)
    obs = {}
    n1_fail = n2_fail = n3_fail = False
    for j in range(6):
        for e in range(0, E_MAX + 1):
            Rt = {t: R_direct(j, e, t) for t in range(1, P)}
            s = Rt[1]
            b = (8 - j) % 6 if e == 0 else (4 + 2 * e) % 6
            for t in range(1, P):
                rec("K3 definition via rules = eq:reflection-local-formula", Rt[t] == R_formula(j, e, t), (j, e, t))
                rec("K4 monomial R(t) = s chi(t)^b, b_{j,0} = 8-j, b_{j,e} = 4+2e", Rt[t] == s.shift(chi_idx(t, b)), (j, e, t))
                if not (Rt[t] == R_formula(j, e, t, vexp=2)): n1_fail = True
                wrong_b = (-j) % 6 if e == 0 else (2 * e) % 6
                if not (Rt[t] == s.shift(chi_idx(t, wrong_b))): n2_fail = True
                if not (R_direct(j, e, t, coc_o=-1) == R_formula(j, e, t)): n3_fail = True
            sz = s.is_zero()
            if e > 0 and j != 4: rec("K5 e>0: s = 0 unless j = 4", sz, (j, e))
            if e > 0 and e % 3 == 2: rec("K5 e = 2 mod 3: s = 0", sz, (j, e))
            if e > 0 and j == 1: rec("K6 j = 1: vanishes for e > 0", sz, (j, e))
            if j == 4 and e > 0 and e % 3 != 2: obs[("j4 nonzero", e)] = not sz
            a = abs(s.val())
            if e == 0:
                rec("K7 |s_{j,0}| <= 1", a <= 1 + 1e-12, (j, a))
                want2 = Fr(1) if j in (1, 2, 3, 5) else Fr(1, 13)
                rec("K7x |s_{j,0}|^2 = 1 (j=1,2,3,5), 1/13 (j=0,4), exact", s * s.conj() == El.const(want2), j)
            else:
                we = abs(W(e).val())
                k, r = divmod(e, 3)
                beta = 13.0 ** (-k) if r == 1 else (13.0 ** (0.5 - k) if r == 0 else 0.0)
                rec("K7 e>0: |s| <= q^(1/2) |W_e| = beta_p(e)", a <= math.sqrt(13) * we * (1 + 1e-12) and abs(math.sqrt(13) * we - beta) <= 1e-12 * max(beta, 1e-300) + 0.0, (j, e, a, we, beta))
    # K6: omega
    omega = R_direct(1, 0, 1)
    rec("K6 R_{p,1,0}(t) = omega chi(t)", all(R_direct(1, 0, t) == omega.shift(chi_idx(t, 1)) for t in range(1, P)))
    rec("K6 omega = (sum chi^3(x) psi(x/p)) / tau_1", omega * tau[1] == Q3)
    rec("K6 |omega| = 1 exactly", omega * omega.conj() == ONE)
    return fails, dict(n1=n1_fail, n2=n2_fail, n3=n3_fail), obs, omega

print("== K1..K7 over all conventions: o in {+1,-1}, c in 1..12; j = 0..5, e = 0..%d, all 12 units t ==" % E_MAX)
KEYS = ["K1 tau*conj(tau) = 13", "K1 |g_2| = 1", "K2 mask identity",
        "K3 definition via rules = eq:reflection-local-formula",
        "K4 monomial R(t) = s chi(t)^b, b_{j,0} = 8-j, b_{j,e} = 4+2e",
        "K5 e>0: s = 0 unless j = 4", "K5 e = 2 mod 3: s = 0",
        "K6 j = 1: vanishes for e > 0", "K6 R_{p,1,0}(t) = omega chi(t)",
        "K6 omega = (sum chi^3(x) psi(x/p)) / tau_1", "K6 |omega| = 1 exactly",
        "K7 |s_{j,0}| <= 1", "K7x |s_{j,0}|^2 = 1 (j=1,2,3,5), 1/13 (j=0,4), exact",
        "K7 e>0: |s| <= q^(1/2) |W_e| = beta_p(e)"]
allf = {k: [] for k in KEYS}
neg = {"n1": [], "n2": [], "n3": []}
obs_all, omegas, omega_exact = [], {}, {}
for o in (1, -1):
    for c in range(1, P):
        f, n, ob, omE = run(o, c); om = omE.val(); omega_exact[(o, c)] = omE
        for k, v in f.items(): allf.setdefault(k, []).append(((o, c), v[:2]))
        for k in neg: neg[k].append(n[k])
        obs_all.append(all(ob.values())); omegas[(o, c)] = om
        print("  convention o=%+d c=%2d done: %d failing keys, omega = %+.6f%+.6fi" % (o, c, len(f), om.real, om.imag)); sys.stdout.flush()
for k in allf:
    check(k + "  (24 conventions)", not allf[k], str(allf[k][:2]) if allf[k] else "")
print("== G: the 24 conventions are Galois conjugates of one ==")
def galois(el, m):
    v = [0] * N
    for k, x in enumerate(el.v): v[(k * m) % N] = x
    return El(v)
gal_ok = True
for (o, c), omE in omega_exact.items():
    m = next(m for m in range(1, N) if m % 6 == o % 6 and m % 13 == c and math.gcd(m, N) == 1)
    if not (galois(omega_exact[(1, 1)], m) == omE): gal_ok = False
check("G: omega(o,c) = sigma_m(omega(1,1)) for the automorphism z -> z^m, m = o mod 6, m = c mod 13 (so the exact equalities in the 24 runs are ONE computation and its conjugates, not 24 independent confirmations)", gal_ok)
print("== negative controls ==")
check("N1: chi(v)^4 -> chi(v)^2 in eq:reflection-local-formula is detected (K3 fails) in every convention", all(neg["n1"]))
check("N2: exponents b -> b - 2 (i.e. -j at e=0, 2e at e>0) are detected (K4 fails) in every convention", all(neg["n2"]))
print("info N3 (exploratory): torus law with the opposite orientation relative to the W unit rule -> K3 fails in %d of 24 conventions" % sum(neg["n3"]))
print("info: for j = 4 and e > 0 with e != 2 mod 3 the value s is nonzero in every convention: %s (the lemma only claims 'zero unless')" % all(obs_all))
print("info: omega takes %d distinct values over the 24 conventions (it depends on the convention; only |omega| = 1 and the formula are claimed)"
      % len({(round(v.real, 9), round(v.imag, 9)) for v in omegas.values()}))
print()
print("SUMMARY: %d checks, %d OK, %d FAIL" % (len(results), sum(results), len(results) - sum(results)))
