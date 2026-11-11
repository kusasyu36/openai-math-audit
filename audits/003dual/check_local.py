# -*- coding: utf-8 -*-
"""
003 (The Quasi-Riemann Hypothesis), paper2.tex -- stage 1 of audits/003dual.

Exhaustive check, at small primary primes p of O = Z[omega], of the finite
identities in the paragraph "The local character transformation":

  eq:ray-fourier          paper2.tex:2992   C_{p,j}(h) in closed form (4 cases)
  eq:theta-local-fourier  paper2.tex:2948   finite Fourier inversion of chi_p^j
  eq:ray-local-transform  paper2.tex:3232   sum_{h != 0} C_{p,j}(h) chi_p(a)^-2 e(eps h^-1 x / p)
                                            = chi_p(sigma)^-2 omega_{p,j} B_{p,j}(x)
  omega_{p,j}             paper2.tex:3239   three cases, |omega_{p,j}| = 1
  B_{p,j}                 paper2.tex:1735   three cases; B_{p,1} = chi_p^3 quadratic

SCOPE (written before running).
  Tested:  the identities above as statements about the finite field O/(p), with
           a := sigma*h, for EVERY sigma, eps in (O/p)^x, EVERY x in O/p (x = 0
           included, where the characters are extended by zero) and j = 0..5.
           Definitions are taken from the tex: chi_p = sextic residue symbol
           (footnote at :256), e(z) = exp(4 pi i Im z / sqrt 3) = exp(2 pi i Tr(z/lambda))
           (:306), gamma_j(p) = N^-1/2 sum chi_p^j(x) e(x/p) (:310), primary = 1 mod 3.
  NOT tested here (inputs or other steps of the proof of Prop. lem:reflection):
           - that the actual a, sigma_p, eps_p of the proof are of this form
             (a = sigma_p h_p mod p, eq:ray-additive-crt): see check_global_sample.py,
             which is a SAMPLE and is not part of the stage-1 gate;
           - eq:ray-multiplier (Kubota character), the cusp expansions d_sigma,
             the archimedean transform, the constant C, the proposition itself.
  So a pass means: the finite algebra from (ray-fourier, a = sigma h, the CRT
  phase) to (ray-local-transform) is right at the primes tried.  Nothing more.

HOW MUCH INFORMATION THERE IS (so the count is not inflated).
  sigma enters both sides only through chi_p(sigma)^-2, and (eps, x) enter the
  left side only through the residue eps*x.  So among the 6 (N-1)^2 N triples per
  prime, the informative cases are about 6 N (pairs (j, t = eps x)).  The two
  primary primes above a split rational prime are complex conjugates; one is run.
  All equalities are identities in Q(omega, zeta_q), so the choice of embedding
  into C carries no information either.

ARITHMETIC.
  Every equality is decided exactly.  A value is an integer array c[a][b]
  standing for sum c[a][b] zeta6^a zeta_q^b (zeta6 = 1 + omega in O, zeta_q =
  exp(2 pi i / q), q the rational prime under p).  It is zero iff the q
  coefficients in Z[omega] are all equal (Phi_q is irreducible over Q(omega), q != 3).
  sqrt N is avoided: with g_j = sqrt N gamma_j, every identity is multiplied by N.
  After that scaling the identity only sees the TOTAL power of N; the split of
  N^-1/2 between omega_{p,j} and B_{p,j} is pinned by G1 (|g_j|^2 = N), which is
  therefore part of the gate.

CRITERIA (fixed before running).  Gate of stage 1 = all of the following, at
every prime in PRIMES:
  K0  field data: norm, primary, mu_6 injective mod p, chi_p multiplicative,
      chi_p(zeta6^k) consistent, chi_p^3 = quadratic symbol, e(./p) a well defined
      nontrivial additive character of O/p, chi_p(-1) = (-1)^((N-1)/6)
  G1  g_j * conj(g_j) = N                      (j = 1..5)   [=> |gamma_j| = 1, |omega| = 1]
  G2  g_1 * g_5 = chi_p(-1) N                  (eq:g, paper2.tex:844)
  F1  N C_{p,j}(h) = closed form eq:ray-fourier, all j, all h (4 cases)
  F2  sum_y chi_p^j(y) e(-y/p) = chi_p(-1)^j g_j            (display before eq:ray-fourier)
  F3  N chi_p^j(x) = sum_h N C_{p,j}(h) e(h x / p), all j, all x (x = 0 included)
  J4  sum_{y != 0} e(x y / p) = -1 + N [p | x], all x
  T1  eq:ray-local-transform, j in {1,2,3,5}, all sigma, eps, x
  T2  same, j = 4
  T3  same, j = 0 (active)
  B1  B_{p,1}(x) = chi_p(x)^-3 = chi_p(x)^3 = quadratic symbol, all x
  Negative controls (each must behave exactly as stated, else the gate fails):
  N1  chi_p(a)^+2 in place of chi_p(a)^-2 on the left (forgetting to conjugate
      kappa): must FAIL for every prime and every j
  N2  drop chi_p(-1)^j from omega_{p,j} (first case): must FAIL exactly when
      chi_p(-1) = -1 and j is odd (j = 1,3,5); must still hold otherwise
  N3  omega_{p,0} with the sign + instead of -: must FAIL for every prime
  Reported, part of the gate only if it runs (needs scipy for the import):
  A   chi_p and gamma_j agree with audits/029reflect/kubota.py (Prime.sext exactly,
      gauss_gamma to 1e-10).  This anchors the conventions OUTSIDE the identity:
      both sides of T1-T3 are built from the same tables, so a convention slip
      would cancel there.
  D0  implementation self-check (not an independent test): for N = 7, 13 the
      left side is recomputed by plain Python loops straight from the
      definitions, without memoisation on eps*x, and compared with the tables.

Run (about 1-2 minutes):
  UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy \
      python -u audits/003dual/check_local.py > audits/003dual/check_local.out
"""
import cmath, math, os, sys, time
import numpy as np

UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]     # zeta6^k, zeta6 = 1 + omega
U0 = np.array([u[0] for u in UNITS], dtype=np.int64)
U1 = np.array([u[1] for u in UNITS], dtype=np.int64)


def mul(x, y):
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)
def conj(x): return (x[0] - x[1], -x[1])
def norm(x): return x[0] * x[0] - x[0] * x[1] + x[1] * x[1]
def omcoef(x): return x[1]                    # Tr(x / lambda) for x = a + b omega


def split_prime(q):
    """the primary prime a + b omega of norm q with b > 0 (the other one is its conjugate)"""
    for b in range(3, 4 * q, 3):
        for a in range(-4 * q, 4 * q):
            if a % 3 == 1 and a * a - a * b + b * b == q:
                return (a, b)
    raise ValueError(q)


# split primes above 7, 13, 19, 31, 37, 43 and the inert primes -5, -11
PRIMES = [split_prime(q) for q in (7, 13, 19, 31, 37, 43)] + [(-5, 0), (-11, 0)]


# ------------------------------------------------------------ exact values
def iszero(D):
    """D[..., a, b] ~ sum D zeta6^a zeta_q^b ; True where the value is 0"""
    P0 = np.tensordot(D, U0, axes=([-2], [0]))       # [..., b]
    P1 = np.tensordot(D, U1, axes=([-2], [0]))
    return (P0 == P0[..., :1]).all(axis=-1) & (P1 == P1[..., :1]).all(axis=-1)

def fmul(A, B):
    q = A.shape[1]; R = np.zeros((6, q), dtype=np.int64)
    for a in range(6):
        for b in range(q):
            c = int(A[a, b])
            if c: R += c * np.roll(np.roll(B, a, axis=0), b, axis=1)
    return R

def fconj(A):
    q = A.shape[1]
    return A[(-np.arange(6)) % 6][:, (-np.arange(q)) % q]

def const(v, q):
    R = np.zeros((6, q), dtype=np.int64); R[0, 0] = v; return R

def embed(A):
    q = A.shape[1]; tot = 0
    for a in range(6):
        for b in range(q):
            if A[a, b]:
                tot += int(A[a, b]) * cmath.exp(1j * cmath.pi * a / 3) * cmath.exp(2j * cmath.pi * b / q)
    return tot


# ------------------------------------------------------------ the finite field O/(p)
class Local:
    def __init__(self, elt):
        self.elt = elt; N = norm(elt); self.N = N
        if elt[1] == 0:
            self.inert = True; q = abs(elt[0]); self.q = q
            self.reps = [(a, b) for a in range(q) for b in range(q)]
        else:
            self.inert = False; q = N; self.q = q
            self.r = (-elt[0] * pow(elt[1], -1, q)) % q          # omega = r mod p
            self.reps = [(k, 0) for k in range(q)]
        N_, q_ = N, q
        idx = self.idx
        self.MUL = np.array([[idx(mul(x, y)) for y in self.reps] for x in self.reps], dtype=np.int64)
        self.ADD = np.array([[idx((x[0] + y[0], x[1] + y[1])) for y in self.reps] for x in self.reps], dtype=np.int64)
        self.NEG = np.array([idx((-x[0], -x[1])) for x in self.reps], dtype=np.int64)
        one = idx((1, 0))
        INV = np.zeros(N, dtype=np.int64)
        for i in range(1, N):
            w = np.nonzero(self.MUL[i] == one)[0]
            assert len(w) == 1
            INV[i] = w[0]
        self.INV = INV
        # e(z/p) = exp(2 pi i omcoef(z conj p) / N) = zeta_q^E(z)
        cp = conj(elt); E = np.zeros(N, dtype=np.int64)
        for i, x in enumerate(self.reps):
            E[i] = self.Eraw(x)
        self.E = E
        # chi_p(x) = zeta6^K(x), the sixth root of unity congruent to x^((N-1)/6)
        z6 = [idx(u) for u in UNITS]; self.z6 = z6
        K = np.full(N, -1, dtype=np.int64)
        for i in range(1, N):
            t = self.powi(i, (N - 1) // 6)
            K[i] = z6.index(t)
        self.K = K
        self.nz = np.arange(1, N, dtype=np.int64)
        self.km1 = int(K[idx((-1, 0))])                           # chi_p(-1) = zeta6^km1
        # g_j = sum_y chi_p^j(y) e(y/p)
        self.g = []
        for j in range(6):
            flat = ((j * K[self.nz]) % 6) * q + E[self.nz]
            self.g.append(np.bincount(flat, minlength=6 * q).reshape(6, q))

    def idx(self, x):
        if self.inert: return (x[0] % self.q) * self.q + x[1] % self.q
        return (x[0] + x[1] * self.r) % self.q

    def Eraw(self, x):
        v = omcoef(mul(x, conj(self.elt))); d = self.N // self.q
        assert v % d == 0
        return (v // d) % self.q

    def powi(self, i, e):
        res = self.idx((1, 0)); base = i
        while e:
            if e & 1: res = int(self.MUL[res, base])
            base = int(self.MUL[base, base]); e >>= 1
        return res

    def NC(self, j):
        """NC[h] = N * C_{p,j}(h) = sum_y chi_p^j(y) e(-h y / p), from the definition"""
        N, q, K, E = self.N, self.q, self.K, self.E; nz = self.nz
        out = np.zeros((N, 6, q), dtype=np.int64)
        a6 = (j * K[nz]) % 6
        for h in range(N):
            flat = a6 * q + E[self.NEG[self.MUL[h, nz]]]
            out[h] = np.bincount(flat, minlength=6 * q).reshape(6, q)
        return out

    def lhs_table(self, j, sign):
        """T[si, t] = N * sum_{h != 0} C_{p,j}(h) chi_p(sigma h)^sign e(t h^-1 / p),
        as the double sum over (h, y) of monomials, from the definitions.
        sign = -2 is the paper; t stands for eps * x."""
        N, q, K, E = self.N, self.q, self.K, self.E; nz = self.nz; n = N - 1
        A6 = (j * K[nz])[None, :]                                  # chi^j(y)
        EQ = E[self.NEG[self.MUL[nz][:, nz]]]                      # e(-h y / p), [h, y]
        INVh = self.INV[nz]
        out = np.zeros((n, N, 6, q), dtype=np.int64)
        for si in range(n):
            s6 = (sign * K[self.MUL[nz[si], nz]]) % 6              # chi(a)^sign, a = sigma h
            base = ((A6 + s6[:, None]) % 6) * q
            for t in range(N):
                u = E[self.MUL[t, INVh]]                           # e(t h^-1 / p)
                flat = base + (EQ + u[:, None]) % q
                out[si, t] = np.bincount(flat.ravel(), minlength=6 * q).reshape(6, q)
        return out

    def rhs_coef(self, j, variant=None):
        """N * omega_{p,j} * (N-power part of B), without the monomial in sigma, eps, x"""
        q = self.q
        if j in (1, 2, 3, 5):
            c = fmul(self.g[j], self.g[(j + 2) % 6])
            if variant != "N2":
                c = np.roll(c, (j * self.km1) % 6, axis=0)         # chi_p(-1)^j
            return c
        if j == 4: return self.g[4]
        if j == 0:
            return self.g[2] if variant == "N3" else -self.g[2]

    def compare(self, j, T, variant=None):
        """number of (sigma, eps, x) with  T[sigma, eps x] != N chi(sigma)^-2 omega B(x)"""
        N, q, K = self.N, self.q, self.K; nz = self.nz; n = N - 1
        coef = self.rhs_coef(j, variant)
        shifts = np.stack([np.roll(coef, s, axis=0) for s in range(6)])     # times zeta6^s
        xs = np.arange(N); Kx = K.copy(); Kx[0] = 0
        bad = 0; tot = 0
        for si in range(n):
            ks = int(K[nz[si]])
            for ei in range(n):
                eps = int(nz[ei]); ke = int(K[eps])
                Lrows = T[si, self.MUL[eps, xs]]                             # [x, 6, q]
                if j == 4:
                    R = np.repeat(shifts[(-2 * ks) % 6][None], N, axis=0) * (-1)
                    R[0] = shifts[(-2 * ks) % 6] * (N - 1)                   # -1 + N [p | x]
                else:
                    m = j + 2 if j != 0 else 2
                    ex = (-2 * ks - m * ke - m * Kx) % 6
                    R = shifts[ex]
                    R[0] = 0                                                 # chi extended by zero
                ok = iszero(Lrows - R)
                bad += int((~ok).sum()); tot += N
        return bad, tot


# ------------------------------------------------------------ checks
RESULTS = []
def report(name, ok, detail=""):
    RESULTS.append((name, ok))
    print("  %-4s %s  %s" % (name, "PASS" if ok else "FAIL", detail)); sys.stdout.flush()


def check_K0(L):
    N, q, K, E = L.N, L.q, L.K, L.E; elt = L.elt
    ok = True
    ok &= (elt[0] % 3 == 1 and elt[1] % 3 == 0)                    # primary
    ok &= (L.idx(elt) == 0) and (N % 6 == 1)
    ok &= (N == q) != L.inert and (N == q * q) == L.inert
    ok &= len(set(L.z6)) == 6                                      # mu_6 injective mod p
    ok &= all(int(K[L.z6[k]]) == (k * (N - 1) // 6) % 6 for k in range(6))
    nz = L.nz
    Km = K[L.MUL[nz][:, nz]]
    ok &= bool(((K[nz][:, None] + K[nz][None, :]) % 6 == Km).all())     # multiplicative
    ok &= bool(((E[:, None] + E[None, :]) % q == E[L.ADD]).all())      # additive
    ok &= bool((E != 0).any()) and int(E[0]) == 0
    # e(z/p) depends only on z mod p: raw definition on a box of O
    for a in range(-6, 7):
        for b in range(-6, 7):
            ok &= (L.Eraw((a, b)) == int(E[L.idx((a, b))]))
    # chi^3 = quadratic symbol
    one = L.idx((1, 0))
    for i in range(1, N):
        ok &= ((3 * int(K[i])) % 6 == 0) == (L.powi(i, (N - 1) // 2) == one)
    ok &= L.km1 == (3 * ((N - 1) // 6)) % 6
    report("K0", bool(ok), "N=%d q=%d %s chi_p(-1)=%+d" % (N, q, "inert" if L.inert else "split", 1 if L.km1 == 0 else -1))


def check_G(L):
    N, q = L.N, L.q
    ok = all(bool(iszero(fmul(L.g[j], fconj(L.g[j])) - const(N, q))) for j in range(1, 6))
    report("G1", ok, "g_j conj(g_j) = N, j=1..5")
    rhs = np.roll(const(N, q), L.km1, axis=0)
    report("G2", bool(iszero(fmul(L.g[1], L.g[5]) - rhs)), "g_1 g_5 = chi_p(-1) N")


def check_F(L, NCs):
    N, q, K, E = L.N, L.q, L.K, L.E
    ok1 = True; ok2 = True; ok3 = True
    for j in range(6):
        NC = NCs[j]
        cj = np.roll(L.g[j], (j * L.km1) % 6, axis=0)              # chi(-1)^j g_j
        for h in range(N):
            if j != 0 and h != 0: R = np.roll(cj, (-j * int(K[h])) % 6, axis=0)
            elif j != 0: R = const(0, q)
            elif h != 0: R = const(-1, q)
            else: R = const(N - 1, q)
            ok1 &= bool(iszero(NC[h] - R))
        if j != 0:
            ok2 &= bool(iszero(NC[L.idx((1, 0))] - cj))
        # inversion: N chi^j(x) = sum_h NC[h] e(h x / p)
        for x in range(N):
            S = np.zeros((6, q), dtype=np.int64)
            for h in range(N):
                S += np.roll(NC[h], int(E[L.MUL[h, x]]), axis=1)
            R = const(0, q) if x == 0 else np.roll(const(N, q), (j * int(K[x])) % 6, axis=0)
            ok3 &= bool(iszero(S - R))
    report("F1", ok1, "eq:ray-fourier, %d pairs (j,h)" % (6 * N))
    report("F2", ok2, "sum chi^j(y) e(-y/p) = chi(-1)^j g_j, j=1..5")
    report("F3", ok3, "Fourier inversion, %d pairs (j,x)" % (6 * N))
    ok4 = True
    for x in range(N):
        flat = E[L.MUL[x, L.nz]]
        S = np.zeros((6, q), dtype=np.int64); S[0] = np.bincount(flat, minlength=q)
        ok4 &= bool(iszero(S - const(-1 + (N if x == 0 else 0), q)))
    report("J4", ok4, "sum_{y!=0} e(xy/p) = -1 + N[p|x]")


def check_T(L):
    N = L.N; n = N - 1
    res = {}
    for j in range(6):
        T = L.lhs_table(j, -2)
        res[j] = L.compare(j, T)
        res[("N2", j)] = L.compare(j, T, "N2") if j in (1, 2, 3, 5) else None
        res[("N3", j)] = L.compare(j, T, "N3") if j == 0 else None
        if N in (7, 13): res[("D0", j)] = direct_selfcheck(L, j, T)
        del T
        T1 = L.lhs_table(j, +2)
        res[("N1", j)] = L.compare(j, T1)
        del T1
    b = sum(res[j][0] for j in (1, 2, 3, 5)); t = sum(res[j][1] for j in (1, 2, 3, 5))
    report("T1", b == 0, "j=1,2,3,5: %d mismatches / %d triples (sigma,eps,x)" % (b, t))
    report("T2", res[4][0] == 0, "j=4: %d / %d" % res[4])
    report("T3", res[0][0] == 0, "j=0: %d / %d" % res[0])
    # B1: B_{p,1}(x) = chi(x)^-3 = chi(x)^3 = quadratic symbol (the latter is in K0)
    K = L.K
    report("B1", all((-3 * int(K[i])) % 6 == (3 * int(K[i])) % 6 for i in range(1, N)), "chi^-3 = chi^3")
    # negative controls
    n1 = [res[("N1", j)][0] for j in range(6)]
    report("N1", all(v > 0 for v in n1), "chi(a)^+2: mismatches per j=0..5: %s (must all be > 0)" % n1)
    n2 = {j: res[("N2", j)][0] for j in (1, 2, 3, 5)}
    want = {j: (L.km1 == 3 and j % 2 == 1) for j in (1, 2, 3, 5)}
    report("N2", all((n2[j] > 0) == want[j] for j in n2),
           "no chi(-1)^j: mismatches %s ; expected to fail for j in %s" % (n2, [j for j in want if want[j]]))
    report("N3", res[("N3", 0)][0] > 0, "omega_{p,0} with +: %d mismatches (must be > 0)" % res[("N3", 0)][0])
    if N in (7, 13):
        report("D0", all(res[("D0", j)] for j in range(6)), "plain-loop recomputation of the left side agrees with the table")
    return 6 * n * n * N, 6 * N


def direct_selfcheck(L, j, T):
    """left side by plain Python loops from the definitions, (eps h^-1) x order, no memoisation"""
    N, q, K, E, MUL, NEG, INV = L.N, L.q, L.K, L.E, L.MUL, L.NEG, L.INV
    for sigma in range(1, N):
        for eps in range(1, N):
            for x in range(N):
                cnt = np.zeros((6, q), dtype=np.int64)
                for h in range(1, N):
                    a = int(MUL[sigma, h])
                    u = int(E[MUL[MUL[eps, INV[h]], x]])
                    for y in range(1, N):
                        cnt[(j * int(K[y]) - 2 * int(K[a])) % 6, (int(E[NEG[MUL[h, y]]]) + u) % q] += 1
                if not (cnt == T[sigma - 1, int(MUL[eps, x])]).all(): return False
    return True


def check_A(L):
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, os.path.join(here, "..", "029reflect"))
        import kubota
    except Exception as ex:                                        # noqa
        print("  A    SKIP  (cannot import kubota.py: %s)" % type(ex).__name__); return
    pr = kubota.get_prime(L.elt)
    ok = all(pr.sext(L.reps[i]) == int(L.K[i]) for i in range(1, L.N)) and pr.sext(L.reps[0]) is None
    dmax = 0.0
    for j in range(1, 6):
        dmax = max(dmax, abs(kubota.gauss_gamma(pr, j) - embed(L.g[j]) / math.sqrt(L.N)))
    report("A", ok and dmax < 1e-10, "chi_p = kubota.Prime.sext (all residues); max |gamma_j - kubota.gauss_gamma| = %.1e" % dmax)


def main():
    t0 = time.time()
    print("003 stage 1: local character transformation, exhaustive over O/(p)")
    print("criteria: see the header of this file (fixed before the first run)\n")
    tot_tr = 0; tot_inf = 0
    for elt in PRIMES:
        L = Local(elt)
        print("p = %d%+d*omega   N(p) = %d   [%.0fs]" % (elt[0], elt[1], L.N, time.time() - t0)); sys.stdout.flush()
        check_K0(L)
        check_A(L)
        check_G(L)
        check_F(L, [L.NC(j) for j in range(6)])
        a, b = check_T(L)
        tot_tr += a; tot_inf += b
        print()
    names = sorted(set(n for n, _ in RESULTS))
    print("summary by check (over %d primes):" % len(PRIMES))
    for nme in names:
        r = [ok for n, ok in RESULTS if n == nme]
        print("  %-4s %d / %d" % (nme, sum(r), len(r)))
    allok = all(ok for _, ok in RESULTS)
    print("\ntriples (j, sigma, eps, x) compared in T1-T3: %d ; informative cases (j, t = eps x): about %d" % (tot_tr, tot_inf))
    print("checks: %d passed / %d   =>  STAGE 1 GATE: %s" % (sum(ok for _, ok in RESULTS), len(RESULTS), "PASS" if allok else "FAIL"))
    print("elapsed %.0fs" % (time.time() - t0))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
