# -*- coding: utf-8 -*-
"""Independent numerical checks for the audit of 029 Theorem 1.2 (core).
Self-written; does not import or execute any code from the paper's repository.
Each check prints OK/FAIL with the tex location it refers to.
"""
import cmath, math, itertools, sys
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp

OK = 0
FAIL = 0
def check(name, cond, extra=""):
    global OK, FAIL
    if cond:
        OK += 1
        print("OK   " + name + (("  " + extra) if extra else ""))
    else:
        FAIL += 1
        print("FAIL " + name + (("  " + extra) if extra else ""))

def close(a, b, tol=1e-8):
    return abs(a - b) <= tol * (1 + abs(a) + abs(b))

Z6 = cmath.exp(2j * cmath.pi / 6)

# ---------------------------------------------------------------------------
# Part A.  Finite fields F_q, q = p or p^2, q = 1 mod 12.
#   02-arithmetic.tex:438-594  (Gauss/Jacobi identities, Davenport-Hasse lift,
#   sextic duplication, sign of gamma_3)
# ---------------------------------------------------------------------------
class GF:
    """F_q for q = p^f, f in {1,2}.  Elements are ints 0..q-1 (a + p*b for a+b*alpha)."""
    def __init__(self, p, f):
        self.p, self.f, self.q = p, f, p ** f
        if f == 2:
            # non-residue n mod p, alpha^2 = n
            n = next(n for n in range(2, p) if pow(n, (p - 1) // 2, p) == p - 1)
            self.n = n
        # find generator and logs
        q = self.q
        for g in range(2, q):
            x, seen = 1, 1
            ok = True
            # order test via full enumeration (q small)
            logs = {1: 0}
            x = 1
            for k in range(1, q - 1):
                x = self.mul(x, g)
                if x in logs:
                    ok = False
                    break
                logs[x] = k
            if ok and len(logs) == q - 1:
                self.g, self.log = g, logs
                break
        self.exp = [0] * (q - 1)
        for x, k in self.log.items():
            self.exp[k] = x
    def mul(self, x, y):
        p = self.p
        if self.f == 1:
            return x * y % p
        a, b = x % p, x // p
        c, d = y % p, y // p
        return (a * c + self.n * b * d) % p + p * ((a * d + b * c) % p)
    def add(self, x, y):
        p = self.p
        if self.f == 1:
            return (x + y) % p
        return (x % p + y % p) % p + p * ((x // p + y // p) % p)
    def neg(self, x):
        p = self.p
        if self.f == 1:
            return (-x) % p
        return (-(x % p)) % p + p * ((-(x // p)) % p)
    def tr(self, x):
        p = self.p
        return x % p if self.f == 1 else (2 * (x % p)) % p
    def psi(self, x):
        return cmath.exp(2j * cmath.pi * self.tr(x) / self.p)
    def chi(self, x, j=1):
        """sextic character chi^j, zero-extended (also for j = 0 mod 6)."""
        if x == 0:
            return 0
        return Z6 ** ((self.log[x] * j) % 6)
    def norm_to_prime(self, x):
        # x^(p+1) for f=2 lies in F_p
        if x == 0:
            return 0
        k = (self.log[x] * (self.p + 1)) % (self.q - 1)
        y = self.exp[k]
        assert y < self.p
        return y
    def gauss(self, j):
        return sum(self.chi(x, j) * self.psi(x) for x in range(1, self.q))
    def jacobi(self, a, b):
        s = 0
        one = 1
        for x in range(self.q):
            y = self.add(one, self.neg(x))
            s += self.chi(x, a) * self.chi(y, b)
        return s

def eisenstein_coords(z):
    """write complex z = a + b*omega with a,b (nearly) integers"""
    w = cmath.exp(2j * cmath.pi / 3)
    b = z.imag / w.imag
    a = z.real - b * w.real
    return a, b

print("== Part A: finite-field identities (02-arithmetic.tex:438-594) ==")
fields = [(13, 1), (37, 1), (61, 1), (5, 2), (7, 2), (11, 2), (13, 2)]
store = {}
for p, f in fields:
    K = GF(p, f)
    Q = K.q
    assert Q % 12 == 1
    G = {j: K.gauss(j) for j in range(1, 6)}
    J = K.jacobi(2, 2)           # rho = chi^2 cubic
    four = K.add(K.add(1, 1), K.add(1, 1))
    c4 = K.chi(four)
    chim1 = K.chi(K.neg(1))
    name = "F_%d" % Q
    check(name + ": chi(-1)=1", close(chim1, 1))
    check(name + ": G2*G4=Q   (eq:cubic-jacobi)", close(G[2] * G[4], Q))
    check(name + ": G2^2=J(rho,rho)*G4", close(G[2] ** 2, J * G[4]))
    check(name + ": (G2/sqrtQ)^3 = J/sqrtQ", close((G[2] / math.sqrt(Q)) ** 3, J / math.sqrt(Q)))
    a, b = eisenstein_coords(J)
    ai, bi = round(a), round(b)
    check(name + ": J in Z[omega], |J|^2=Q, J = -1 mod 3",
          close(a, ai) and close(b, bi) and ai * ai - ai * bi + bi * bi == Q
          and (ai + 1) % 3 == 0 and bi % 3 == 0, "J=%d%+d*omega" % (ai, bi))
    g3 = G[3] / math.sqrt(Q)
    check(name + ": gamma_3 = +-1", close(g3, 1) or close(g3, -1), "gamma_3=%+.0f" % g3.real)
    check(name + ": G1*G4 = conj(chi(4))*G2*G3  (sextic duplication)", close(G[1] * G[4], c4.conjugate() * G[2] * G[3]))
    gam = {j: G[j] / math.sqrt(Q) for j in G}
    check(name + ": gamma1*gamma2 = conj(chi(4))*gamma3*gamma2^3", close(gam[1] * gam[2], c4.conjugate() * gam[3] * gam[2] ** 3))
    check(name + ": gamma_j*gamma_{-j}=1 (j=1,2)", close(gam[1] * gam[5], 1) and close(gam[2] * gam[4], 1))
    store[(p, f)] = (K, J)
    if f == 2 and p % 3 == 2:
        # inert in Q(omega): pi_+ = -p, so the paper predicts J = -pi_+ = p and gamma_2^3 = +1
        check(name + ": inert case J = p = -pi_+  (eq:eisenstein-jacobi with pi_+=-p)", close(J, p))
        check(name + ": inert case gamma_2^3 = +1 = -Lambda", close(gam[2] ** 3, 1))
    if f == 2 and p % 3 == 1:
        # split in Q(omega): residue degree 2 of F-prime over the Q(omega)-prime.
        # Davenport-Hasse: J_{p^2}(rho o N, rho o N) = -J_p(rho,rho)^2.
        # Build rho on F_p as x -> (chi^2 restricted through norm): rho0(N x) = chi^2(x).
        rho0 = {}
        consistent = True
        for x in range(1, Q):
            nx = K.norm_to_prime(x)
            v = K.chi(x, 2)
            if nx in rho0 and not close(rho0[nx], v):
                consistent = False
            rho0[nx] = v
        check(name + ": chi^2 factors through the residue norm (lift of a cubic character of F_p)", consistent and len(rho0) == p - 1)
        J0 = sum(rho0[x] * rho0[(1 - x) % p] for x in range(2, p))
        check(name + ": Davenport-Hasse  J_{p^2} = (-1)^{f-1} J_p^f = -J_p^2", close(J, -J0 ** 2), "J_p=%s" % (complex(round(J0.real, 6), round(J0.imag, 6)),))
        check(name + ": hence gamma_2^3 = -(pi_+/|pi_+|)^2 with pi_+=-J_p", close(gam[2] ** 3, -((-J0) / abs(J0)) ** 2))

# ---------------------------------------------------------------------------
# Part B.  The local Poisson sums C_{t,k}(J)  (05-poisson.tex:498-595) on Z/13^n,
#          computed from the raw definition (eq:C-local-definition) and also the
#          un-factored finite Fourier transform (eq:finite-Poisson-transform).
# ---------------------------------------------------------------------------
print("== Part B: C_{t,k}(J) from the definition vs eq:C-local-values (05-poisson.tex:498-595) ==")
P = 13
# sextic character on (Z/13)^*: 2 is a primitive root mod 13
LOG13 = {pow(2, k, P): k for k in range(12)}
def chi13(x, j):
    x %= P
    if x == 0:
        return 0.0
    return Z6 ** ((LOG13[x] * j) % 6)
CHI = {j: np.array([chi13(x, j) for x in range(P)]) for j in range(-12, 13)}
def chivec(arr, j):
    return CHI[((j + 12) % 6)][arr % P] if (j % 6) != 0 else (arr % P != 0).astype(complex)
def e_(num, den):
    return np.exp(2j * np.pi * (num % den) / den)

def C_def(t, k, J):
    """raw definition eq:C-local-definition with psi^o(x)=exp(2 pi i x) on p^{-n}Z/Z"""
    pk, pt, pJ = P ** k, P ** t, P ** J
    tot = 0.0
    d = np.arange(pk)
    v = np.arange(pt)
    cond = ((pJ - pt * d) % pk == 0)
    dd = d[cond]
    if k > 0:
        cd = chivec(dd, k)
    else:
        cd = np.ones(len(dd), dtype=complex)
    if t > 0:
        cv = chivec(v, t)
    else:
        cv = np.ones(1, dtype=complex)
    den = P ** (k + t)
    for di, cdi in zip(dd, cd):
        if cdi == 0:
            continue
        num = (pJ - pt * int(di)) * v
        tot += cdi * np.sum(cv * e_(num, den))
    return tot

def T_(vexp, y):
    """T_v(y) = sum_{x mod p} chi(x)^v psi^o(y x / p), unit mask kept for v=0 mod 6"""
    x = np.arange(P)
    return np.sum(chivec(x, vexp) * e_((int(y) % P) * x, P))
TAU = {v: T_(v, 1) for v in range(0, 12)}
def C_formula(t, k, J):
    Q = P
    if t == 0:
        if k == 0:
            return 1.0
        return 1.0 if J == 0 else 0.0
    if k == 0:
        return (Q ** (t - 1)) * T_(t, P ** (J - t + 1)) if J >= t - 1 else 0.0
    if k == 1:
        return TAU[1] * (Q ** (t - 1)) * T_(t - 1, P ** (J - t)) if J >= t else 0.0
    return 0.0
bad = []
count = 0
for t in range(0, 4):
    for k in range(0, 3):
        if t + k > 4:
            continue
        for J in range(0, 7):
            a = C_def(t, k, J)
            b = C_formula(t, k, J)
            count += 1
            if not close(a, b, 1e-7):
                bad.append((t, k, J, a, b))
check("C_{t,k}(J): raw definition = eq:C-local-values for t<=3,k<=2,t+k<=4,J<=6 (%d cases, p=13)" % count, not bad, str(bad[:3]))
check("tau_0=-1, |tau_v|=sqrt(Q), tau_1*tau_5=Q", close(TAU[0], -1) and all(close(abs(TAU[v]), math.sqrt(P)) for v in range(1, 6)) and close(TAU[1] * TAU[5], P))

# the un-factored finite transform (eq:finite-Poisson-transform) for s=p^k, A=p^t (unit parts 1), h=p^J:
def finite_transform(t, k, J):
    pk, pt = P ** k, P ** t
    mod = P ** (k + t)
    m = np.arange(mod)
    # g_s(-m) = sum_{d mod s} chi_s(d) e(-m d / s)
    if k > 0:
        d = np.arange(pk)
        cd = chivec(d, k)
        gs = np.array([np.sum(cd * e_(-int(mm) * d, pk)) for mm in m])
    else:
        gs = np.ones(mod, dtype=complex)
    cA = chivec(m, t) if t > 0 else np.ones(mod, dtype=complex)
    return np.sum(gs * cA * e_((P ** J) * m, mod))
bad = []
for t, k in [(1, 0), (1, 1), (2, 1), (0, 1), (3, 1), (1, 2), (2, 0), (3, 0)]:
    for J in range(0, 5):
        lhs = finite_transform(t, k, J)           # = q_s * F(s,A,h) before the 1/sqrt(q_s) factors
        rhs = (P ** k) * C_formula(t, k, J)
        if not close(lhs, rhs, 1e-6):
            bad.append((t, k, J, lhs, rhs))
check("sum_{m mod sA} g_s(-m) chi_A(m) e(hm/(sA)) = q_s*C_{t,k}(J)  (eq:finite-Poisson-transform, one-prime model)", not bad, str(bad[:3]))

# closed Euler factor eq:Euler-good-exact
print("== Part B2: closed form of the good Euler factor (05-poisson.tex:597-654) ==")
rng = np.random.default_rng(20261111)
g1 = TAU[1] / math.sqrt(P)          # gamma_1 in the psi^o sign; chi(-1)=1 so sign-free
g2 = T_(2, 1) / math.sqrt(P)
def euler_series(xi, w, z, vartheta, theta_p, L=14, M=14):
    Q = P
    tot = 0.0
    for e0 in (0, 1):
        for l in range(0, L):
            t = e0 + 3 * l
            for k in (range(0, 80) if t == 0 else (0, 1, 2)):   # t=0 admits every k>=1 (C_{0,k}=1_{J=0})
                for mp_ in range(0, M):
                    J = 6 * mp_
                    c = C_formula(t, k, J)
                    if c == 0:
                        continue
                    tot += (g2 ** e0) * (vartheta ** (e0 + 3 * l)) * (theta_p ** (k - e0 - 3 * l)) * complex(c) \
                        * cmath.exp(math.log(Q) * (-e0 * (0.5 + xi) - l * (1 + 3 * xi) - k * w - 6 * mp_ * z))
    return tot
bad = []
for trial in range(6):
    xi = 1.3 + rng.random() + 1j * rng.normal()
    w = 0.8 + rng.random() + 1j * rng.normal()
    z = 0.5 + rng.random() + 1j * rng.normal()
    vartheta = cmath.exp(2j * cmath.pi * rng.random())
    theta_p = Z6 ** int(rng.integers(0, 6))
    eta = -g1 * g2 * vartheta            # eq:gauss-corrections: gamma1*gamma2*vartheta = -eta
    Q = P
    V = Q ** (-6 * z); W = theta_p * Q ** (-w)
    D1 = eta * theta_p.conjugate() * Q ** (-xi)
    E1 = vartheta ** 6 * Q ** (4 - 6 * xi) * V
    closed = 1 / (1 - V) + W / (1 - W) + E1 / (1 - E1) * ((1 - 1 / Q) / (1 - V) + W) \
        - D1 / (1 - E1) * (1 + (Q - 1) * W * V / (1 - V))
    ser = euler_series(xi, w, z, vartheta, theta_p)
    if not close(ser, closed, 1e-7):
        bad.append((trial, ser, closed))
check("F_{u,p} series (eq:exact-Euler-factor, j=0) = closed form eq:Euler-good-exact, 6 random points", not bad, str(bad[:2]))
check("|gamma_1*gamma_2| = 1 so |eta|=1 in the test", close(abs(g1 * g2), 1))

# ---------------------------------------------------------------------------
# Part C.  Appendix: shells of the Jacquet integral (10-whittaker-calculation.tex:295-457)
# ---------------------------------------------------------------------------
print("== Part C: Jacquet shells and the exceptional exponent (10-whittaker-calculation.tex:331-457) ==")
def I_shell_bruteforce(e, l):
    """int_{O^x} chi(v)^{-8l} conj(psi(p^{e-l} v)) dv, psi of conductor O, vol(O)=1"""
    n = l - e
    if n <= 0:
        # additive character trivial: (1/p) * sum over units mod p
        v = np.arange(P)
        return np.sum(chivec(v, -8 * l)) / P
    mod = P ** n
    v = np.arange(mod)
    return np.sum(chivec(v, -8 * l) * np.conj(e_(v, mod))) / mod
def g_norm(j):
    x = np.arange(P)
    return np.sum(chivec(x, j) * e_(x, P)) / math.sqrt(P)
def I_shell_formula(e, l):
    Q = P
    if 1 <= l <= e:
        return (1 - 1 / Q) if l % 3 == 0 else 0.0
    if l == e + 1:
        return (Q ** -0.5) * g_norm(-2 * l) if l % 3 else -1 / Q
    return 0.0
bad = []
for e in range(0, 4):
    for l in range(1, e + 4):
        if l - e > 3:
            continue
        a = I_shell_bruteforce(e, l); b = I_shell_formula(e, l)
        if not close(a, b, 1e-9):
            bad.append((e, l, a, b))
check("I_{e,l} brute force = eq:app-shell-integrals (e<=3, l<=e+3)", not bad, str(bad[:3]))

def J_vector(e, s, lmax=None):
    """Vector Jacquet value in the basis (e_0,e_1,e_2) for inducing exponent -s on a(y)
    (the paper is s=1/6).  Derived independently:
      common factor Q^{-e+e/2-e*s};  unit ball -> e_{-e};
      shell l>=1 -> Q^{2 l s} * I_{e,l} * e_{2l-e}   (measure Q^l, modular Q^{-l}, torus Q^{2ls})."""
    Q = P
    vec = np.zeros(3, dtype=complex)
    vec[(-e) % 3] += 1.0
    for l in range(1, e + 2):
        vec[(2 * l - e) % 3] += Q ** (2 * l * s) * I_shell_formula(e, l)
    return Q ** (-e + e / 2 - e * s) * vec
s0 = 1 / 6
J0 = J_vector(0, s0)
g2n, g4n = g_norm(2), g_norm(4)
exp_vals = {}
okall = True
for e in range(0, 10):
    Je = J_vector(e, s0)
    k, r = divmod(e, 3)
    pred = {0: P ** (-k), 1: P ** (-k - 0.5) * g2n, 2: 0.0}[r]
    if not all(close(Je[i], pred * J0[i], 1e-10) for i in range(3)):
        okall = False
check("s=1/6: J_e = W_e * J_0 with W_{3k}=Q^-k, W_{3k+1}=Q^{-k-1/2} g_2, W_{3k+2}=0 (e<=9)  (eq:app-J-values, eq:app-pure-values)", okall)
check("J_0 = e_0 + Q^{-1/6} g_4 e_2", close(J0[0], 1) and close(J0[1], 0) and close(J0[2], P ** (-1 / 6) * g4n))
# away from the exceptional exponent the vectors J_e are NOT all proportional to J_0
def aligned(s):
    J0 = J_vector(0, s)
    for e in (3, 4, 6):
        Je = J_vector(e, s)
        M = np.array([J0, Je])
        if np.linalg.matrix_rank(M, tol=1e-9) > 1:
            return False
    return True
check("alignment of J_e with J_0 holds at s=1/6 and fails at s=0.1, 0.2, 0.3 (the exceptional point is detected)",
      aligned(1 / 6) and not aligned(0.1) and not aligned(0.2) and not aligned(0.3))
# orientation consistency: Gauss character chi^2 <-> unit rule chi(u)^{2e}
u = 5
x = np.arange(P)
lhs = np.sum(chivec(x, 2) * e_(x * pow(u, -1, P), P)) / P
rhs = chi13(u, 2) * np.sum(chivec(x, 2) * e_(x, P)) / P
check("sum chi(x)^2 psi(x/(pi u)) = chi(u)^2 sum chi(x)^2 psi(x/pi)  (10-whittaker-calculation.tex:483-495)", close(lhs, rhs))

# ---------------------------------------------------------------------------
# Part D.  Local Weyl lemma (04-reflection.tex:339-431) re-derived from the stated rules
# ---------------------------------------------------------------------------
print("== Part D: local reflection lemma (04-reflection.tex:78-97, 339-431) ==")
def tau_pj(j):
    v = np.arange(1, P)
    return np.sum(chivec(v, -j) * e_(v, P))
def alpha(j, v):
    v %= P
    if j % 6 != 0:
        return 0.0 if v == 0 else chi13(v, -j) / tau_pj(j)
    return (1 - 1 / P) if v == 0 else -1 / P
bad = []
for j in range(6):
    for t in range(P):
        s = sum(alpha(j, v) * cmath.exp(2j * cmath.pi * t * v / P) for v in range(P))
        want = chi13(t, j) if j else (0.0 if t == 0 else 1.0)
        if not close(s, want, 1e-10):
            bad.append((j, t))
check("mask identity sum_v alpha_{p,j}(v) psi(tv/p) = chi(t)^j incl. j=0 mask (eq:reflection-mask-identity)", not bad)

def W_pure(e):
    if e < 0: return 0.0
    k, r = divmod(e, 3)
    return {0: P ** (-k), 1: P ** (-k - 0.5) * g2n, 2: 0.0}[r]
def C_torus(e1, u1, e2, u2):
    # C(p^{e1} u1, p^{e2} u2) = chi(u1)^{2 e2} chi(u2)^{-2 e1}
    return chi13(u1, 2 * e2) * chi13(u2, -2 * e1)
def R_direct(j, e, t):
    """sum over unit digits v of alpha(v) W(a(p^{e-2} t) w n(v/p)) using
       w n(H) = n(1/H) a(-H^-2) k,  H=v/p :  psi(p^{e-1} t / v) * C(p^{e-2}t, p^2 v^-2) * W(a(p^e t v^-2))."""
    tot = 0.0
    for v in range(1, P):
        vinv = pow(v, -1, P)
        add = cmath.exp(2j * cmath.pi * (t * vinv % P) / P) if e == 0 else (1.0 if e >= 1 else None)
        if e < 0:
            continue
        cocyc = C_torus(e - 2, t, 2, vinv * vinv % P)
        unit = t * vinv * vinv % P
        Wval = chi13(unit, 2 * e) * W_pure(e)
        tot += alpha(j, v) * add * cocyc * Wval
    return tot
bad = []
omega_p = None
for j in range(6):
    for e in range(0, 7):
        s_je = R_direct(j, e, 1)
        b = (8 - j) % 6 if e == 0 else (4 + 2 * e) % 6
        for t in range(1, P):
            val = R_direct(j, e, t)
            if not close(val, s_je * chi13(t, b) if b else s_je, 1e-9):
                bad.append(("monomial", j, e, t))
        if e > 0 and j != 4 and abs(s_je) > 1e-9:
            bad.append(("vanish", j, e))
        if e > 0 and e % 3 == 2 and abs(s_je) > 1e-9:
            bad.append(("vanish e=2 mod 3", j, e))
        beta = 1.0 if e == 0 else (P ** (-(e - 1) // 3) if e % 3 == 1 else (P ** (0.5 - e // 3) if e % 3 == 0 else 0.0))
        if abs(s_je) > beta + 1e-9:
            bad.append(("beta bound", j, e, abs(s_je), beta))
        if j == 1 and e == 0:
            omega_p = s_je
check("R_{p,j,e}(t) = s_{p,j,e} chi(t)^{b_{j,e}}, b_{j,0}=8-j, b_{j,e}=4+2e; vanishing for e>0,j!=4 and e=2 mod 3; |s|<=beta_p(e)  (eq:reflection-local-monomial, eq:reflection-frozen-weight)", not bad, str(bad[:3]))
check("|omega_p| = 1 and row exponent at j=1,e=0 is 1  (eq:reflection-valuation-one)", close(abs(omega_p), 1))
q3 = np.sum(chivec(np.arange(P), 3) * e_(np.arange(P), P))
check("omega_p = (quadratic Gauss sum)/tau_{p,1}", close(omega_p, q3 / tau_pj(1)))
check("twist involution j -> -(j+2) mod 6 fixes the pairing 1<->3 (sextic <-> quadratic), 0<->4, 2<->2, 5<->5",
      [(-(j + 2)) % 6 for j in range(6)] == [4, 3, 2, 1, 0, 5])

# ---------------------------------------------------------------------------
# Part E.  Exponent bookkeeping (02a:69-76,131-166; 06-zero-free.tex:371-683)
# ---------------------------------------------------------------------------
print("== Part E: exponents ==")
a, b = Fr(9, 10), Fr(133, 1000)
M0, h0 = a + b, 1 - a
refl = (2 * M0 - 1 - b) / 2
Cstar1 = a / 2 + h0 / 6
check("reflection exponent (2M0-1-b)/2 = 933/2000", refl == Fr(933, 2000))
check("C_*(1) = a/2 + h0/6 = 7/15; margin = 1/6000", Cstar1 == Fr(7, 15) and Cstar1 - refl == Fr(1, 6000))
check("margin identity: C_*(1) - refl = (4 - 4a - 3b)/6 (so reflection needs 4a+3b<4)", Cstar1 - refl == (4 - 4 * a - 3 * b) / 6, "4a+3b=%s" % (4 * a + 3 * b))
def Estar(sx, sw, sz):
    return a / 2 + sx - 1 + h0 * sz + b * (sw - 1)
d1 = Fr(1, 10 ** 7)
s1 = 1 + d1
rel = lambda x: x - Cstar1
Rw = rel(Estar(s1, Fr(95, 100), Fr(1, 6) + d1)); Rz = rel(Estar(s1, 1, Fr(16, 100)))
check("R_w, R_z relative exponents < -0.0005", Rw < Fr(-5, 10 ** 4) and Rz < Fr(-5, 10 ** 4), "%.9f %.9f" % (Rw, Rz))
base = rel(Estar(s1, Fr(3, 1000), Fr(4, 10)))
small = base + Fr(9, 100) * Fr(11, 10)
exc = base + Fr(10001, 100000) * Fr(6, 10)
good = rel(Estar(Fr(998, 1000), Fr(3, 1000), Fr(4, 10))) + Fr(10001, 100000) * Fr(11, 10)
check("rows U<=Z^.09: %.9f; exceptional: %.9f; good central: %.9f; all < -0.0004" % (small, exc, good),
      small < Fr(-4, 10 ** 4) and exc < Fr(-4, 10 ** 4) and good < Fr(-4, 10 ** 4))
# the same nonprincipal bound WITHOUT the zero-density shift (xi stays at 1+): is it below the principal scale?
no_shift = base + Fr(10001, 100000) * Fr(11, 10)
check("WITHOUT the xi-shift to .998 the good-row bound is ABOVE C_*(1) (%.9f > 0): the zero-density input is indispensable" % no_shift, no_shift > 0)
# general window: need  1.5*h0 - 0.997*b - 0.002 (+1e-7, +1.1e-5 slack) < h0/6 - 0.0004   and   4a+3b<4
lhs_np = Fr(11, 10) * Fr(10001, 100000) + h0 * (Fr(4, 10) - Fr(1, 6)) - b * Fr(997, 1000) - Fr(2, 1000)
check("good-row margin recomputed = -377/300000 (06-zero-free.tex:577-580, without the 1e-7)", lhs_np == Fr(-377, 300000))
R0 = 350000
far = Fr(3849, 1000) + h0 * R0 + Fr(10001, 100000) * (1 - R0)
check("far rows exponent = .44901 < C_*(1) - .01", far == Fr(44901, 100000) and far < Cstar1 - Fr(1, 100))
check("1/6000 - 1e-5 > 1/20000 > 1e-6", Fr(1, 6000) - Fr(1, 10 ** 5) > Fr(1, 20000) > Fr(1, 10 ** 6))
check("zero detector: 10(1/2-.996)+1/2+1 = -3.46; 11*(2-2*.996)=.088; .088+.1 < 1/2",
      10 * (Fr(1, 2) - Fr(996, 1000)) + Fr(3, 2) == Fr(-346, 100) and 11 * (2 - 2 * Fr(996, 1000)) == Fr(88, 1000))
th = math.log(1.002 / 0.999) / math.log(1.003 / 0.999)
check("three circles: theta=%.4f<1; 2-1.0035=.9965>.996; 2-1.002=.998" % th, th < 1)
# Euler-region exponents (05-poisson.tex:673-689)
r1 = (Fr(998, 1000), Fr(3, 1000), Fr(4, 10)); r2 = (Fr(998, 1000), Fr(95, 100), Fr(16, 100))
def worst(f):
    return max(f(*r1), f(*r2))
check("Euler regions: max exponents of VW, D1V, D1W, Q*D1WV, E1 are all < -1",
      worst(lambda x, w, z: -6 * z - w) < -1 and worst(lambda x, w, z: -x - 6 * z) < -1 and worst(lambda x, w, z: -x - w) <= Fr(-1001, 1000)
      and worst(lambda x, w, z: 1 - x - w - 6 * z) <= Fr(-1908, 1000) and worst(lambda x, w, z: 4 - 6 * x - 6 * z) <= Fr(-2948, 1000),
      "D1W worst = %s" % worst(lambda x, w, z: -x - w))
# sensitivity: how much room is there in (a,b)?
def admissible(aa, bb):
    # Only TWO of the paper's constraints, with the paper's contour constants frozen
    # (w=.003, z=.4, xi=.998, row threshold U<=Z^{h0+1e-5}); other constraints are not tested here.
    hh = 1 - aa
    refl_ok = (4 - 4 * aa - 3 * bb) / 6 > Fr(1, 20000)            # reflection beats the principal scale with the stated loss
    good_ok = Fr(11, 10) * (hh + Fr(1, 100000)) + hh * (Fr(4, 10) - Fr(1, 6)) - bb * Fr(997, 1000) - Fr(2, 1000) < Fr(-4, 10 ** 4)
    return refl_ok and good_ok and bb < aa
cnt = 0; tot = 0
for ai in range(850, 1000):
    for bi in range(1, 300):
        tot += 1
        if admissible(Fr(ai, 1000), Fr(bi, 1000)):
            cnt += 1
print("     (a,b) on the 1/1000 grid, .85<=a<1, 0<b<.3: %d of %d satisfy these two margins (paper's contour constants frozen; other constraints not tested)" % (cnt, tot))
check("(a,b)=(.9,.133) satisfies these two margins (computed)", admissible(a, b))

# ---------------------------------------------------------------------------
# Part F.  Archimedean formulas (03-theta.tex:263-290, 680-725; 02a:31-36)
# ---------------------------------------------------------------------------
print("== Part F: archimedean formulas ==")
mp.mp.dps = 30
for s in (mp.mpf("0.3"), mp.mpf("1.7"), mp.mpc("0.4", "1.1")):
    c = mp.mpf("2.3")
    lhs = mp.quad(lambda t: t ** s * mp.besselk(mp.mpf(1) / 3, c * mp.sqrt(t)), [0, 1, 10, mp.inf])
    rhs = 2 ** (2 * s + 1) * c ** (-2 * s - 2) * mp.gamma(s + mp.mpf(5) / 6) * mp.gamma(s + mp.mpf(7) / 6)
    check("Mellin of t*K_{1/3}(c sqrt t) at s=%s  (eq:theta-bessel-mellin)" % mp.nstr(s, 3), abs(lhs - rhs) < mp.mpf(10) ** -15 * (1 + abs(rhs)))
for xi in (mp.mpf("0.7"), mp.mpc("-0.3", "0.9")):
    lhs = mp.quad(lambda t: mp.exp(-mp.log(t) ** 2 / 4) / (2 * mp.sqrt(mp.pi)) * t ** (xi - 1), [0, 1, mp.inf])
    check("Mellin of P_G at xi=%s equals exp(xi^2)  (eq:reflection-gaussian)" % mp.nstr(xi, 3), abs(lhs - mp.exp(xi ** 2)) < mp.mpf(10) ** -15)
nu0 = mp.mpf(2) / 3; c = mp.mpf("1.9"); rho = mp.mpf("0.8")
lap = mp.pi * rho ** (2 * nu0) / mp.gamma(2 * nu0) * mp.quad(lambda t: t ** (2 * nu0 - 2) * mp.exp(-t * rho ** 2 - c ** 2 / (4 * t)), [0, 1, mp.inf])
bes = mp.pi / mp.gamma(2 * nu0) * 2 * (c / 2) ** (2 * nu0 - 1) * rho * mp.besselk(2 * nu0 - 1, c * rho)
check("Laplace form of the complex Jacquet integral = const * rho * K_{1/3}(c rho)  (03-theta.tex:277-286)", abs(lap - bes) < mp.mpf(10) ** -15)
# direct 2-D integral (radial reduction: int_{R^2} e^{-i c x1} (rho^2+|x|^2)^{-2nu0} dx = 2 pi int_0^inf J_0(c r) r (rho^2+r^2)^{-2 nu0} dr)
direct = rho ** (2 * nu0) * 2 * mp.pi * mp.quadosc(lambda r: mp.besselj(0, c * r) * r * (rho ** 2 + r ** 2) ** (-2 * nu0), [0, mp.inf], omega=c)
check("direct oscillatory integral agrees with the Bessel form", abs(direct - bes) < mp.mpf(10) ** -8, "rel.err=%s" % mp.nstr(abs(direct - bes) / abs(bes), 3))
# decay rate of the smoothing kernel's Mellin transform: |1/Gamma|^{2r} ~ exp(pi r |t|)
r = 2
t = mp.mpf(40)
val = abs(1 / (mp.gamma(mp.mpc(0.5, t) + mp.mpf(5) / 6) * mp.gamma(mp.mpc(0.5, t) + mp.mpf(7) / 6))) ** r
check("|Gamma(s+5/6)Gamma(s+7/6)|^{-r} <= C(1+|t|)^C exp(pi r |t|) sanity at t=40, r=2 (03-theta.tex:737-742)",
      mp.log(val) < mp.pi * r * t + 5 * mp.log(1 + t))

print("TOTAL OK=%d FAIL=%d" % (OK, FAIL))
