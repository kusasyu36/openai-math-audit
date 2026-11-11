# -*- coding: utf-8 -*-
"""Kubota-Patterson cubic theta function over Q(omega): self-written tool kit.

Source of the formulas (transcribed by me from Dunn-Radziwill, arXiv:2109.07463v3):
  (1.7)  g~(c) = |c|^{-1} sum_{x mod c} (x/c)_3 e^(x/c),   e^(z) = exp(2 pi i (z + conj z))
  (5.1)  action of SL2(C) on H^3
  (5.4)  chi(gamma) = (c/a)_3 for gamma in Gamma_1(3), c != 0;  chi = 1 on SL2(Z)
  (5.6)-(5.8)  theta(w) = sigma v^{2/3} + sum_{mu in lambda^{-3}O} tau(mu) v K_{1/3}(4 pi |mu| v) e^(mu z)
  (5.13),(5.14) tau_1, tau_2;  Appendix A: d_10(mu) = omega tau_2(omega mu) e^(mu),
                                           d_19(mu) = omega^2 tau_1(omega^2 mu) e^(mu)
  gamma_10 = [[1,0],[omega,1]],  gamma_19 = [[1,0],[omega^2,1]].
Whether these transcriptions (and the formulas themselves) are right is NOT assumed:
run_gate.py tests them numerically.

Elements of Z[omega] are pairs (a, b) = a + b*omega of Python ints.
"""
import cmath, math, os, pickle
import numpy as np
from scipy.special import kv

W3 = cmath.exp(2j * cmath.pi / 3)
Z6 = cmath.exp(1j * cmath.pi / 3)          # = -omega^2 = 1 + omega
E9 = cmath.exp(2j * cmath.pi / 9)
SIGMA = 3 ** 2.5 / 2
LAM = (1, 2)                               # lambda = 1 + 2 omega = sqrt(-3)
LAMC = 1j * math.sqrt(3)
ONE, OMEGA, OMEGA2 = (1, 0), (0, 1), (-1, -1)
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]   # zeta6^k, k=0..5


def norm(x): return x[0] * x[0] - x[0] * x[1] + x[1] * x[1]
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def neg(x): return (-x[0], -x[1])
def mul(x, y):
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)
def conj(x): return (x[0] - x[1], -x[1])
def cx(x): return x[0] + x[1] * W3
def powe(x, n):
    r = (1, 0)
    for _ in range(n): r = mul(r, x)
    return r
def exact_div(y, x):
    """y / x if x | y else None"""
    n = norm(x); a, b = mul(y, conj(x))
    if a % n == 0 and b % n == 0: return (a // n, b // n)
    return None
def divmod_e(y, x):
    """Euclidean division: y = q x + r with N(r) < N(x)"""
    n = norm(x); a, b = mul(y, conj(x))
    q = ((2 * a + n) // (2 * n), (2 * b + n) // (2 * n))
    r = sub(y, mul(q, x))
    best = (norm(r), q, r)
    for da in (-1, 0, 1):
        for db in (-1, 0, 1):
            q2 = (q[0] + da, q[1] + db); r2 = sub(y, mul(q2, x))
            if norm(r2) < best[0]: best = (norm(r2), q2, r2)
    return best[1], best[2]
def xgcd(a, c):
    """return (g, s, t) with s a + t c = g  (g a gcd)"""
    r0, r1 = a, c; s0, s1 = (1, 0), (0, 0); t0, t1 = (0, 0), (1, 0)
    while r1 != (0, 0):
        q, r = divmod_e(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, sub(s0, mul(q, s1))
        t0, t1 = t1, sub(t0, mul(q, t1))
    return r0, s0, t0
def is_primary(x): return x[0] % 3 == 1 and x[1] % 3 == 0
def primary_assoc(x):
    for u in UNITS:
        y = mul(x, u)
        if is_primary(y): return y
    raise ValueError(x)
def tr(x): return 2 * x[0] - x[1]          # x + conj(x)
def omcoef(x): return x[1]                  # Tr(x / lambda) = b for x = a + b omega


# ---------------------------------------------------------------- primes
_PR = {}
class Prime:
    """primary prime pi of Z[omega], prime to 3"""
    def __init__(self, elt):
        assert is_primary(elt)
        self.elt = elt; self.N = norm(elt)
        if elt[1] == 0:
            self.inert = True; self.q = abs(elt[0]); self.rat = self.q
            assert self.N == self.q ** 2
            self.z6 = [(1, 0)]
            for _ in range(5):
                t = mul(self.z6[-1], (1, 1)); self.z6.append((t[0] % self.q, t[1] % self.q))
        else:
            self.inert = False; p = self.N; self.p = p; self.rat = p
            self.r = (-elt[0] * pow(elt[1], -1, p)) % p            # omega = r mod pi
            assert (self.r * self.r + self.r + 1) % p == 0
            self.z6 = [pow((1 + self.r) % p, k, p) for k in range(6)]
    def _pow(self, x, e):
        if self.inert:
            q = self.q; a, b = x[0] % q, x[1] % q
            if a == 0 and b == 0: return None
            res = (1, 0); base = (a, b)
            while e:
                if e & 1:
                    res = mul(res, base); res = (res[0] % q, res[1] % q)
                base = mul(base, base); base = (base[0] % q, base[1] % q)
                e >>= 1
            return res
        v = (x[0] + x[1] * self.r) % self.p
        if v == 0: return None
        return pow(v, e, self.p)
    def cub(self, x):
        """exponent k (value omega^k) of the cubic residue symbol (x/pi)_3; None if pi | x"""
        t = self._pow(x, (self.N - 1) // 3)
        if t is None: return None
        return {self.z6[0]: 0, self.z6[2]: 1, self.z6[4]: 2}[t]
    def sext(self, x):
        """exponent k (value zeta6^k, zeta6 = -omega^2) of the sextic symbol; None if pi | x"""
        assert (self.N - 1) % 6 == 0
        t = self._pow(x, (self.N - 1) // 6)
        if t is None: return None
        return self.z6.index(t)

def get_prime(elt):
    if elt not in _PR: _PR[elt] = Prime(elt)
    return _PR[elt]

_SPLIT = {}
def split_primes(p):
    """the two primary primes above a rational prime p = 1 mod 3"""
    if p not in _SPLIT:
        found = None
        for b in range(1, int(math.sqrt(4 * p / 3)) + 2):
            disc = b * b - 4 * (b * b - p)
            if disc < 0: continue
            s = math.isqrt(disc)
            if s * s == disc and (b + s) % 2 == 0:
                found = ((b + s) // 2, b); break
        assert found and norm(found) == p
        pi = primary_assoc(found)
        _SPLIT[p] = (get_prime(pi), get_prime(primary_assoc(conj(pi))))
    return _SPLIT[p]

def rat_factor(n):
    out = []; d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0: n //= d; e += 1
            out.append((d, e))
        d += 1 if d == 2 else 2
    if n > 1: out.append((n, 1))
    return out

def factor(x, ratfac=None):
    """x = unit * lambda^k * prod pi^e (pi primary primes).  returns (unit, k, [(Prime, e)])"""
    assert x != (0, 0)
    x0 = x; k = 0
    while norm(x) % 3 == 0:
        x = exact_div(x, LAM); k += 1
    fac = []
    for P, e in (ratfac if ratfac is not None else rat_factor(norm(x))):
        if P == 3: continue
        if P % 3 == 2:
            pr = get_prime(primary_assoc((P, 0))); v = 0
            while True:
                y = exact_div(x, pr.elt)
                if y is None: break
                x = y; v += 1
            if v: fac.append((pr, v))
        else:
            for pr in split_primes(P):
                v = 0
                while True:
                    y = exact_div(x, pr.elt)
                    if y is None: break
                    x = y; v += 1
                if v: fac.append((pr, v))
    assert norm(x) == 1, (x0, x)
    return x, k, fac

def cubic_symbol(x, a):
    """exponent of (x/a)_3 for a = 1 mod 3 (any, via factorisation of a); None if not coprime"""
    assert is_primary(a), a
    u, k, fac = factor(a)
    assert u == (1, 0) and k == 0
    tot = 0
    for pr, e in fac:
        c = pr.cub(x)
        if c is None: return None
        tot += c * e
    return tot % 3


# ---------------------------------------------------------------- Gauss sums
def gauss_tilde(pr):
    """g~(pi) = N^{-1/2} sum_x (x/pi)_3 e^(x/pi), DR (1.7)"""
    if pr.inert:
        q = pr.q; s = pr.elt[0]            # pi = s = +-q
        A, B = np.meshgrid(np.arange(q, dtype=np.int64), np.arange(q, dtype=np.int64), indexing="ij")
        ra, rb = np.ones_like(A), np.zeros_like(A); ba, bb = A.copy(), B.copy()
        e = (q * q - 1) // 3
        while e:
            if e & 1:
                ra, rb = (ra * ba - rb * bb) % q, (ra * bb + rb * ba - rb * bb) % q
            ba, bb = (ba * ba - bb * bb) % q, (2 * ba * bb - bb * bb) % q
            e >>= 1
        k = np.full(A.shape, -1)
        for kk, z in ((0, pr.z6[0]), (1, pr.z6[2]), (2, pr.z6[4])):
            k[(ra == z[0]) & (rb == z[1])] = kk
        nz = (A != 0) | (B != 0)
        assert np.all(k[nz] >= 0)
        ph = np.exp(2j * np.pi * (k / 3.0 + (2 * A - B) / float(s)))
        return np.sum(ph[nz]) / q
    p = pr.p; a, b = pr.elt
    xs = np.arange(1, p, dtype=np.int64)
    v = xs                                    # residues are rational integers
    res = np.ones(p - 1, dtype=np.int64); base = v.copy(); e = (p - 1) // 3
    while e:
        if e & 1: res = res * base % p
        base = base * base % p
        e >>= 1
    k = np.where(res == 1, 0, np.where(res == pr.z6[2], 1, 2))
    assert np.all((res == 1) | (res == pr.z6[2]) | (res == pr.z6[4]))
    t = (2 * a - b) % p                       # Tr(x/pi) = x Tr(conj pi)/p = x (2a-b)/p
    return np.sum(np.exp(2j * np.pi * (k / 3.0 + (xs * t % p) / p))) / math.sqrt(p)

def gauss_brute(c_elt, fac):
    """direct g~(c) for squarefree primary c with prime list fac (test only)"""
    N = norm(c_elt); g = math.gcd(c_elt[0], c_elt[1]); cc = cx(c_elt); tot = 0
    for i in range(N // g):
        for j in range(g):
            x = (i, j); val = 1
            for pr in fac:
                k = pr.cub(x)
                if k is None: val = 0; break
                val *= W3 ** k
            if val == 0: continue
            z = cx(x) / cc
            tot += val * cmath.exp(4j * cmath.pi * z.real)
    return tot / math.sqrt(N)


# ---------------------------------------------------------------- enumeration of primary A = c d^3
class ThetaData:
    """all primary A = c d^3 (c squarefree) with N(A) <= Nmax, with
       G(A) = conj(g~(c)) |d|,  kl = exponent of (lambda/c)_3,  kw = exponent of (omega/c)_3"""
    def __init__(self, Nmax, cache_dir=None, verbose=True):
        self.Nmax = Nmax
        fn = os.path.join(cache_dir, "thetadata_%d.pkl" % Nmax) if cache_dir else None
        if fn and os.path.exists(fn):
            with open(fn, "rb") as f: d = pickle.load(f)
            self.A, self.G, self.kl, self.kw, self.checks = d
            return
        spf = np.zeros(Nmax + 1, dtype=np.int64)
        for i in range(2, Nmax + 1):
            if spf[i] == 0:
                spf[i::i] = np.where(spf[i::i] == 0, i, spf[i::i])
        gt = {}
        def g_of(pr):
            if pr.elt not in gt:
                c = conj(pr.elt)
                if (not pr.inert) and c in gt: gt[pr.elt] = gt[c].conjugate()
                else: gt[pr.elt] = gauss_tilde(pr)
            return gt[pr.elt]
        A_list, G_list, kl_list, kw_list = [], [], [], []
        R = int(math.sqrt(4 * Nmax / 3)) + 3
        checks = {"mult": [], "cube": 0.0, "abs": 0.0}
        for a in range(-R, R + 1):
            if a % 3 != 1: continue
            for b in range(-R - (-R % 3), R + 1, 3):
                N = norm((a, b))
                if N > Nmax or N == 0: continue
                rf = []; m = N
                while m > 1:
                    P = int(spf[m]); e = 0
                    while m % P == 0: m //= P; e += 1
                    rf.append((P, e))
                u, k, fac = factor((a, b), rf)
                assert u == (1, 0) and k == 0
                if any(e % 3 == 2 for _, e in fac): continue
                sq = [pr for pr, e in fac if e % 3 == 1]
                dabs = 1.0
                for pr, e in fac: dabs *= pr.N ** ((e // 3) / 2)
                g = 1.0 + 0j
                for pr in sq: g *= g_of(pr)
                for i in range(len(sq)):
                    for j in range(i + 1, len(sq)):
                        g *= W3 ** (sq[i].cub(sq[j].elt) + sq[j].cub(sq[i].elt))
                if len(sq) == 2 and len(checks["mult"]) < 6:
                    c_elt = mul(sq[0].elt, sq[1].elt)
                    if norm(c_elt) < 700:
                        checks["mult"].append((c_elt, abs(gauss_brute(c_elt, sq) - g)))
                A_list.append((a, b)); G_list.append(g.conjugate() * dabs)
                kl_list.append(sum(pr.cub(LAM) for pr in sq) % 3)
                kw_list.append(sum(pr.cub(OMEGA) for pr in sq) % 3)
        # g~(pi)^3 = -pi/|pi|  (DR (1.8)), |g~| = 1
        for elt, val in gt.items():
            z = cx(elt)
            checks["cube"] = max(checks["cube"], abs(val ** 3 + z / abs(z)))
            checks["abs"] = max(checks["abs"], abs(abs(val) - 1))
        checks["nprimes"] = len(gt)
        self.A = np.array(A_list, dtype=np.int64); self.G = np.array(G_list)
        self.kl = np.array(kl_list); self.kw = np.array(kw_list); self.checks = checks
        if fn:
            with open(fn, "wb") as f: pickle.dump((self.A, self.G, self.kl, self.kw, self.checks), f)


class Series:
    """F(z, v) = const * v^{2/3} + sum_k coef_k v K_{1/3}(4 pi |mu_k| v) e^(mu_k z).
       num_k = (integer pair) lam^den * mu_k,  den = 3 or 4."""
    def __init__(self, num, coef, den, const=0.0):
        num = np.asarray(num, dtype=np.int64); coef = np.asarray(coef, dtype=complex)
        mu = (num[:, 0] + num[:, 1] * W3) / LAMC ** den
        o = np.argsort(np.abs(mu), kind="stable")
        self.num, self.coef, self.mu, self.den, self.const = num[o], coef[o], mu[o], den, const
        self.amu = np.abs(self.mu)
    def twisted(self, factor_array, const=None):
        s = Series.__new__(Series)
        s.num, s.mu, s.amu, s.den = self.num, self.mu, self.amu, self.den
        s.coef = self.coef * factor_array
        s.const = self.const if const is None else const
        return s
    def conj(self):
        """series of w -> conj(F(w)): coefficient at mu is conj(coef(-mu))"""
        return Series(-self.num, np.conj(self.coef), self.den, np.conj(self.const))
    def abs_sum(self, v, T=36.0, deriv=None):
        """sum of the absolute values of the terms (the scale against which cancellation is measured)"""
        n = np.searchsorted(self.amu, T / (4 * math.pi * v), side="right")
        t = np.abs(self.coef[:n]) * v * kv(1.0 / 3, 4 * math.pi * self.amu[:n] * v)
        if deriv: t = t * 2 * math.pi * self.amu[:n]
        return abs(self.const) * v ** (2.0 / 3) * (0 if deriv else 1) + np.sum(t)
    def __call__(self, z, v, T=36.0, deriv=None):
        """deriv = 'zbar' -> d/d(zbar) ;  'z' -> d/dz  (termwise)"""
        n = np.searchsorted(self.amu, T / (4 * math.pi * v), side="right")
        if n == len(self.amu) and n > 0 and 4 * math.pi * self.amu[-1] * v < T - 1e-9:
            raise RuntimeError("series table too short: need |mu| up to %.3f, have %.3f"
                               % (T / (4 * math.pi * v), self.amu[-1]))
        mu = self.mu[:n]
        terms = self.coef[:n] * v * kv(1.0 / 3, 4 * math.pi * self.amu[:n] * v) * np.exp(4j * math.pi * (mu * z).real)
        if deriv == "zbar": return np.sum(terms * 2j * math.pi * np.conj(mu))
        if deriv == "z": return np.sum(terms * 2j * math.pi * mu)
        return self.const * v ** (2.0 / 3) + np.sum(terms)


def build_series(td, Mmax, swap9=False):
    """theta (DR (5.7)), F_10* and F_19* (DR (5.13),(5.14), App. A) for |mu| <= Mmax"""
    A = td.A; NA = A[:, 0] ** 2 - A[:, 0] * A[:, 1] + A[:, 1] ** 2
    G, kl, kw = td.G, td.kl, td.kw
    def times(elt):
        e0, e1 = elt
        return np.stack([A[:, 0] * e0 - A[:, 1] * e1, A[:, 0] * e1 + A[:, 1] * e0 - A[:, 1] * e1], axis=1)
    num, coef = [], []
    lim3 = 27 * Mmax ** 2 * (1 + 1e-12)        # N(lam^3 mu) <= 27 Mmax^2
    # case 4:  mu = +-lam^{3n-3} A, n >= 0
    n = 0
    while 27 ** n <= lim3:
        sel = NA * 27 ** n <= lim3
        if n == 0: assert lim3 <= td.Nmax * (1 + 1e-9), "ThetaData too short for Mmax"
        L = powe(LAM, 3 * n)
        for sgn in (1, -1):
            num.append(times((sgn * L[0], sgn * L[1]))[sel]); coef.append((G * 3 ** (n / 2 + 2.5))[sel])
        n += 1
    # cases 1-3:  mu = +- omega^k lam^{3n-4} A, n >= 1
    n = 1
    while 3 ** (3 * n - 1) <= lim3:
        sel = NA * 3 ** (3 * n - 1) <= lim3
        L = powe(LAM, 3 * n - 1)
        for k, pref in ((0, 1.0), (1, E9 if swap9 else E9.conjugate()), (2, E9.conjugate() if swap9 else E9)):
            Lk = mul(L, [ONE, OMEGA, OMEGA2][k])
            c = pref * W3 ** ((k * kw + 2 * kl) % 3) * G * 3 ** (n / 2 + 2)
            for sgn in (1, -1):
                num.append(times((sgn * Lk[0], sgn * Lk[1]))[sel]); coef.append(c[sel])
        n += 1
    theta = Series(np.concatenate(num), np.concatenate(coef), 3, SIGMA)
    # cusp series: mu in lam^{-4} O, numerators m = lam^4 mu = 9 mu
    lim4 = 81 * Mmax ** 2 * (1 + 1e-12)
    sel = NA <= lim4          # complete only up to |mu| <= sqrt(Nmax)/9 ; Series.__call__ guards the range
    def ecek(nm):                                # e^(mu) = exp(2 pi i Tr(m)/9)
        return np.exp(2j * np.pi * (2 * nm[:, 0] - nm[:, 1]) / 9.0)
    base = [G * W3 ** ((2 * kl) % 3), E9.conjugate() * G * W3 ** ((kw + 2 * kl) % 3), E9 * G * W3 ** ((2 * kw + 2 * kl) % 3)]
    # d_10(mu) = omega tau_2(omega mu) e^(mu);  tau_2 lives on -lam^{-4} A {1, omega, omega^2}
    n10, c10 = [], []
    for unit_mu, c in ((neg(OMEGA2), 9 * W3 * W3 ** 2 * base[0]),      # omega mu = -lam^-4 A
                       (neg(ONE), 9 * W3 * base[1]),                    # omega mu = -omega lam^-4 A
                       (neg(OMEGA), 9 * W3 * W3 * base[2])):            # omega mu = -omega^2 lam^-4 A
        nm = times(unit_mu)[sel]; n10.append(nm); c10.append(c[sel] * ecek(nm))
    F10 = Series(np.concatenate(n10), np.concatenate(c10), 4, 0.0)
    # d_19(mu) = omega^2 tau_1(omega^2 mu) e^(mu);  tau_1 lives on lam^{-4} A {1, omega, omega^2}
    n19, c19 = [], []
    for unit_mu, c in ((OMEGA, 9 * W3 ** 2 * W3 * base[0]),            # omega^2 mu = lam^-4 A
                       (OMEGA2, 9 * W3 ** 2 * W3 ** 2 * base[1]),      # omega^2 mu = omega lam^-4 A
                       (ONE, 9 * W3 ** 2 * base[2])):                   # omega^2 mu = omega^2 lam^-4 A
        nm = times(unit_mu)[sel]; n19.append(nm); c19.append(c[sel] * ecek(nm))
    F19 = Series(np.concatenate(n19), np.concatenate(c19), 4, 0.0)
    return theta, F10, F19


# ---------------------------------------------------------------- SL2 action and coset machinery
def mobius(g, z, v):
    """DR (5.1); g = ((a,b),(c,d)) complex numbers"""
    (a, b), (c, d) = g
    den = abs(c * z + d) ** 2 + abs(c) ** 2 * v * v
    return ((a * z + b) * (c * z + d).conjugate() + a * c.conjugate() * v * v) / den, v / den

def mat_mul(X, Y):
    return ((add(mul(X[0][0], Y[0][0]), mul(X[0][1], Y[1][0])), add(mul(X[0][0], Y[0][1]), mul(X[0][1], Y[1][1]))),
            (add(mul(X[1][0], Y[0][0]), mul(X[1][1], Y[1][0])), add(mul(X[1][0], Y[0][1]), mul(X[1][1], Y[1][1]))))
def mat_inv(X): return ((X[1][1], neg(X[0][1])), (neg(X[1][0]), X[0][0]))
def mat_det(X): return sub(mul(X[0][0], X[1][1]), mul(X[0][1], X[1][0]))
def mat_cx(X): return ((cx(X[0][0]), cx(X[0][1])), (cx(X[1][0]), cx(X[1][1])))
def mat_mod3(X): return tuple(tuple((e[0] % 3, e[1] % 3) for e in row) for row in X)
I2 = ((ONE, (0, 0)), ((0, 0), ONE))
def Tmat(t): return ((ONE, t), ((0, 0), ONE))
EMAT = (((0, 0), (-1, 0)), (ONE, (0, 0)))
PMAT = ((OMEGA, (0, 0)), ((0, 0), OMEGA2))
G10 = ((ONE, (0, 0)), (OMEGA, ONE))
G19 = ((ONE, (0, 0)), (OMEGA2, ONE))
COSET_REPS = {"I": I2, "P": PMAT, "P2": mat_mul(PMAT, PMAT), "g10": G10, "g19": G19}

_SL2Z3 = None
def sl2z_lifts():
    """integer lifts of the 24 elements of SL2(Z/3), by BFS on T and E"""
    global _SL2Z3
    if _SL2Z3 is None:
        gens = [Tmat(ONE), EMAT, mat_inv(Tmat(ONE)), mat_inv(EMAT)]
        seen = {mat_mod3(I2): I2}; frontier = [I2]
        while frontier:
            nxt = []
            for X in frontier:
                for g in gens:
                    Y = mat_mul(X, g); k = mat_mod3(Y)
                    if k not in seen: seen[k] = Y; nxt.append(Y)
            frontier = nxt
        assert len(seen) == 24
        _SL2Z3 = seen
    return _SL2Z3

def kubota(g1):
    """exponent of chi(g1) = (c/a)_3 for g1 = I mod 3  (DR (5.4))"""
    assert mat_mod3(g1) == mat_mod3(I2), g1
    a, c = g1[0][0], g1[1][0]
    if c == (0, 0): return 0
    k = cubic_symbol(c, a)
    assert k is not None
    return k

def decompose(g):
    """g in SL2(Z[omega]).  Find  g = g1 * sigma * rep * T_t  with g1 = I mod 3, sigma in SL2(Z),
       rep in COSET_REPS, t in Z[omega].  Then theta(g w) = chi(g1) * theta(rep (w + t)).
       returns (kexp, repname, t)"""
    assert mat_det(g) == ONE
    lifts = sl2z_lifts()
    for name, rep in COSET_REPS.items():
        ri = mat_inv(rep)
        for t0 in range(3):
            for t1 in range(3):
                t = (t0, t1)
                M = mat_mul(mat_mul(g, Tmat(neg(t))), ri)
                k = mat_mod3(M)
                if all(e[1] == 0 for row in k for e in row):
                    sigma = lifts[k]
                    g1 = mat_mul(M, mat_inv(sigma))
                    return kubota(g1), name, t
    raise RuntimeError("no coset found")


class Kubota:
    """theta and its translates through cusp expansions"""
    def __init__(self, Nmax, Mmax, cache_dir=None):
        self.td = ThetaData(Nmax, cache_dir)
        self.theta, self.F10, self.F19 = build_series(self.td, Mmax)
        self._dec = {}
        self.conj_kappa = False          # negative control switch
    def F(self, name, z, v, **kw):
        if name == "I": return self.theta(z, v, **kw)
        if name == "P": return self.theta(W3 ** 2 * z, v, **kw)       # P w = (omega^2 z, v)
        if name == "P2": return self.theta(W3 * z, v, **kw)
        if name == "g10": return self.F10(z, v, **kw)
        if name == "g19": return self.F19(z, v, **kw)
    def cusp_data(self, a, c):
        """for a reduced fraction a/c: matrix g with g(inf) = a/c and its decomposition"""
        key = (a, c)
        if key not in self._dec:
            g0, s, t = xgcd(a, c)                 # s a + t c = g0 (unit)
            u = exact_div(ONE, g0)
            d, b = mul(s, u), neg(mul(t, u))      # a d - b c = 1
            g = ((a, b), (c, d))
            assert mat_det(g) == ONE
            self._dec[key] = (g,) + decompose(g)
        return self._dec[key]
    def translate(self, a, c, z, v):
        """theta(z + a/c, v) evaluated via  theta(g w') = chi(g1) F_rep(w' + t),  w' = g^{-1}(z + a/c, v)"""
        g, kexp, name, t = self.cusp_data(a, c)
        zz, vv = mobius(mat_cx(mat_inv(g)), z + cx(a) / cx(c), v)
        if self.conj_kappa: kexp = -kexp
        return W3 ** kexp * self.F(name, zz + cx(t), vv)
    def translate_parts(self, a, c):
        g, kexp, name, t = self.cusp_data(a, c)
        return g, kexp, name, t


def reduce_frac(a, c):
    g0, _, _ = xgcd(a, c)
    return exact_div(a, g0), exact_div(c, g0)


def sextic_table(pr, j):
    """tab[a % P, b % P] = chi_p(a + b omega)^j (0 on multiples of p; j = 0 gives the unit indicator)"""
    P = pr.rat; tab = np.zeros((P, P), dtype=complex)
    for a in range(P):
        for b in range(P):
            k = pr.sext((a, b))
            if k is not None: tab[a, b] = Z6 ** ((k * j) % 6)
    return tab

def gauss_gamma(pr, j):
    """gamma_j(p) = N^{-1/2} sum_y chi_p(y)^j e(y/p),  e(x) = e^(x/lambda)   (003 eq:intro-gauss-sum)"""
    P = pr.rat; tot = 0; pc = conj(pr.elt)
    seen = set()
    for a in range(P):
        for b in range(P):
            if pr.inert: key = (a, b)
            else:
                key = (a + b * pr.r) % P
                if key in seen: continue
                seen.add(key)
            k = pr.sext((a, b))
            if k is None: continue
            tot += Z6 ** ((k * j) % 6) * cmath.exp(2j * cmath.pi * omcoef(mul((a, b), pc)) / pr.N)
    return tot / math.sqrt(pr.N)
