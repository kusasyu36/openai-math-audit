# -*- coding: utf-8 -*-
"""
Common tool kit for stages 2-4 of audits/003dual (003, paper2.tex, Proposition lem:reflection over Q(omega)).
No check is run by importing this file.

Imports (read only, nothing in audits/029reflect is modified; run python with -B):
  check_global_sample (cgs): exact arithmetic of Z[omega] on pairs (a, b) = a + b*omega
  check_local (cl)         : class Local = chi_p exponent table K, exact Gauss sums g_j, embed()
  kubota (kb)              : theta series tau (DR (5.7)), cusp series F10 = theta(gamma_10 w), F19 = theta(gamma_19 w)

Conventions (tex lines of paper2.tex):
  e(z) = exp(2 pi i Tr(z/lambda))  (:306),   e~(z) = exp(2 pi i Tr z) = e(lambda z)  (:2900)
  primary = 1 mod 3;  S = primes above 2 and 3 (:665);  chi_p = sextic symbol, zeta6 = 1 + omega (:256)
  phi(n) = chi_n(lambda)^2 Psi_0(n) for n = 1 mod 3, (n, 2) = 1, else 0   (:2909)
  gamma_+ = [[1,0],[omega,1]] = kubota G10,  gamma_- = [[1,0],[omega^2,1]] = kubota G19   (:1773)
  d_sigma(l) = conj t_sigma(-l) = coefficients of conj theta(gamma_sigma w)   (:3118)
  UNIT CONVENTION (the lead auditor's reading, tested in stage3a): c is normalised as in :3047
  (a = 1 mod 3 if lambda | c, else c = 1 mod 3) and c_0 := c / r, so that c = c_0 r exactly.
"""
import sys, os, math, cmath
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "029reflect"))
import check_global_sample as cgs
import check_local as cl
import kubota as kb

CACHE = os.path.join(HERE, "cache")
mul, add, sub, neg, conj, norm = cgs.mul, cgs.add, cgs.sub, cgs.neg, cgs.conj, cgs.norm
powe, prod, exact_div, divides, xgcd, gcd = cgs.powe, cgs.prod, cgs.exact_div, cgs.divides, cgs.xgcd, cgs.gcd
inv_mod, val, is_primary, mod3, crt = cgs.inv_mod, cgs.val, cgs.is_primary, cgs.mod3, cgs.crt
UNITS, ZERO, ONE, LAM, TWO = cgs.UNITS, cgs.ZERO, cgs.ONE, cgs.LAM, cgs.TWO
LAM2, LAM3, LAM4 = powe(LAM, 2), powe(LAM, 3), powe(LAM, 4)
W3 = cmath.exp(2j * cmath.pi / 3)
Z6 = cmath.exp(1j * cmath.pi / 3)
SIGMA0 = 3 ** 2.5 / 2


def emod(x, m):
    return cgs.divmod_e(x, m)[1]

def cxf(x):
    """complex value of a + b omega (python ints of any size -> float)"""
    return complex(x[0] - x[1] / 2.0, x[1] * math.sqrt(3) / 2)

def esum(xs):
    s = ZERO
    for x in xs: s = add(s, x)
    return s


# ---------------------------------------------------------------- cubic residue symbol for large arguments
# (al/be)_3 for primary be, by cubic reciprocity and the supplementary laws.  The supplementary laws are NOT
# taken from the literature: the tables SUPW, SUPL (exponents of (omega/be)_3, (lambda/be)_3 as functions of
# be mod 9) are calibrated by brute force against kubota.cubic_symbol (factorisation + Euler criterion), and
# jacobi3 is validated against kubota.cubic_symbol in stage3a (check J0).
SUPW, SUPL = {}, {}
def _calibrate():
    bad = 0
    for a in range(-40, 41):
        for b in range(-39, 40, 3):
            n = (a, b)
            if a % 3 != 1 or norm(n) > 1500 or norm(n) == 1: continue
            key = (a % 9, b % 9)
            kw, kl = kb.cubic_symbol(kb.OMEGA, n), kb.cubic_symbol(kb.LAM, n)
            if SUPW.setdefault(key, kw) != kw or SUPL.setdefault(key, kl) != kl: bad += 1
    assert bad == 0 and len(SUPW) == 9 and len(SUPL) == 9, (bad, len(SUPW))
_calibrate()

def unit_lam_primary(x):
    """x != 0  ->  (t, k, y) with x = zeta6^t lambda^k y, y primary"""
    k = 0
    while True:
        y = exact_div(x, LAM)
        if y is None: break
        x = y; k += 1
    for t, u in enumerate(UNITS):
        y = mul(x, conj(u))
        if is_primary(y): return t, k, y
    raise ValueError(x)

def jacobi3(al, be):
    """exponent k of (al/be)_3 = omega^k for primary be; None if not coprime; (x/1)_3 = 1 for every x"""
    assert is_primary(be), be
    e = 0
    while True:
        if norm(be) == 1: return e % 3
        al = emod(al, be)
        if al == ZERO: return None
        t, k, y = unit_lam_primary(al)
        key = (be[0] % 9, be[1] % 9)
        e += 2 * t * SUPW[key] + k * SUPL[key]         # zeta6 = -omega^2, (-1/be)_3 = 1
        al, be = be, y                                 # (y/be)_3 = (be/y)_3, both primary


# ---------------------------------------------------------------- local data at p
_LOC = {}
class Loc:
    """chi_p, gamma_j(p), C_{p,j}, B_{p,j}, omega_{p,j} as complex numbers; built on check_local.Local"""
    def __init__(self, elt):
        self.elt = elt; self.L = cl.Local(elt); self.N = self.L.N
        self.gam = [cl.embed(self.L.g[j]) / math.sqrt(self.N) for j in range(6)]     # gamma_j(p), j = 1..5 (index 0 unused)
        self.chim1 = Z6 ** self.L.km1
        self._C = {}
    def k(self, x):
        """sextic exponent of chi_p(x) (zeta6^k), None if p | x"""
        v = int(self.L.K[self.L.idx(x)])
        return None if v < 0 else v
    def chi(self, x, j=1):
        k = self.k(x)
        return 0.0 if k is None else Z6 ** ((k * j) % 6)
    def C(self, j):
        """C_{p,j}(h) for h in L.reps order, FROM THE DEFINITION :2945 (not from eq:ray-fourier)"""
        if j not in self._C:
            NC = self.L.NC(j)
            self._C[j] = np.array([cl.embed(NC[h]) for h in range(self.N)]) / self.N
        return self._C[j]
    def omega(self, j, eps):
        """omega_{p,j} of :3239-3244 (eps = eps_p as an element of O)"""
        ke = self.k(eps); g = self.gam
        if j == 4: return g[4]
        if j == 0: return -g[2] * Z6 ** ((-2 * ke) % 6)
        return self.chim1 ** j * g[j] * g[(j + 2) % 6] * Z6 ** ((-(j + 2) * ke) % 6)
    def idx_arr(self, m):
        L = self.L
        if L.inert: return (m[:, 0] % L.q) * L.q + (m[:, 1] % L.q)
        return (m[:, 0] + m[:, 1] * L.r) % L.q
    def Btab(self, j):
        """B_{p,j} of :1735-1742 as a table over L.reps order"""
        N = self.N; K = self.L.K
        out = np.zeros(N, dtype=complex)
        if j == 4:
            out[:] = -1.0 / math.sqrt(N); out[0] = (N - 1) / math.sqrt(N)
        elif j == 0:
            out[1:] = Z6 ** ((-2 * K[1:]) % 6) / math.sqrt(N)
        else:
            out[1:] = Z6 ** ((-(j + 2) * K[1:]) % 6)
        return out
    def chitab(self, j):
        K = self.L.K; out = np.zeros(self.N, dtype=complex)
        out[1:] = Z6 ** ((j * K[1:]) % 6)
        return out

def loc(elt):
    if elt not in _LOC: _LOC[elt] = Loc(elt)
    return _LOC[elt]


# ---------------------------------------------------------------- exact additive phases
def T_frac(num, den):
    """Tr(num/(den*lambda)) mod 1 as a Fraction:  e(num/den) = exp(2 pi i T_frac)"""
    return cgs.T(num, den)

class AddChar:
    """m -> e(num * m / den) for integer pairs m (vectorised, exact reduction before the float exp)"""
    def __init__(self, num, den):
        self.Nd = norm(den)
        w = mul(num, conj(den))
        self.w0, self.w1 = w[0] % self.Nd, w[1] % self.Nd
        assert self.Nd < 2 ** 31
    def __call__(self, m):
        m = np.asarray(m, dtype=np.int64)
        m0, m1 = m[:, 0] % self.Nd, m[:, 1] % self.Nd
        t = ((self.w0 * m1) % self.Nd + (self.w1 * ((m0 - m1) % self.Nd)) % self.Nd) % self.Nd
        return np.exp(2j * np.pi * t / float(self.Nd))


# ---------------------------------------------------------------- Psi_0 and phi
def psi0_exp(name, n):
    """exponent k (value omega^k) of Psi_0(n) for primary n prime to 2.  Candidates (multiplicativity and
    periodicity are CHECKED in stage2, not assumed):
      A: trivial                       B: conj (lambda/n)_3
      C: conj (lambda/n)_3 * omega^t,  t = ((n-1)/3 mod lambda) in F_3
      D: rho_2(n) = the cube root of unity congruent to n mod 2
      E: C * rho_2
      F: conj (lambda/n)_3 * omega^(-t)   (added after table P2 of stage2 showed that C reaches only
         sigma = '+' in the case v_lambda(c) = 1; F is used in stage 3 only)"""
    def rho2(n): return {(1, 0): 0, (0, 1): 1, (1, 1): 2}[(n[0] % 2, n[1] % 2)]
    def lam(n): return jacobi3(LAM, n)
    def t3(n):
        x = exact_div(sub(n, ONE), (3, 0)); return (x[0] + x[1]) % 3
    if name == "A": return 0
    if name == "B": return (-lam(n)) % 3
    if name == "C": return (-lam(n) + t3(n)) % 3
    if name == "D": return rho2(n)
    if name == "E": return (-lam(n) + t3(n) + rho2(n)) % 3
    if name == "F": return (-lam(n) - t3(n)) % 3
    raise ValueError(name)

def phi_exp(name, n, with2=True):
    """exponent of phi(n) (:2909), or None where phi(n) = 0.  with2=False drops the condition at 2
    (the defect of the earlier phi_W, used only as a negative control)."""
    if not is_primary(n): return None
    if with2 and (n[0] % 2 == 0 and n[1] % 2 == 0): return None
    if not with2:
        assert name in ("A", "B", "C", "F")
    return (jacobi3(LAM, n) + psi0_exp(name, n)) % 3

def phi_table(name, Lint, with2=True):
    """phi on O/(Lint), Lint a rational integer: array [x0, x1] complex.  Representatives are shifted to
    be nonzero primary elements; periodicity mod Lint is checked in stage2 (check P1)."""
    tab = np.zeros((Lint, Lint), dtype=complex)
    for x0 in range(Lint):
        for x1 in range(Lint):
            k = phi_exp(name, (x0, x1), with2)
            if k is not None: tab[x0, x1] = W3 ** k
    return tab

def phi_hat(tab, Lint):
    """phi-hat(h_0) = N(L)^-1 sum_{x mod L} phi(x) e(-h_0 x / L)  (:2936), L = Lint; returns dict h0 -> value.
    e(-h x / L) = exp(-2 pi i (h x)_omega-coefficient / L)."""
    xs0, xs1 = np.meshgrid(np.arange(Lint), np.arange(Lint), indexing="ij")
    out = {}
    for h0 in range(Lint):
        for h1 in range(Lint):
            # (h x)[1] = h0 x1 + h1 x0 - h1 x1
            t = (h0 * xs1 + h1 * xs0 - h1 * xs1) % Lint
            out[(h0, h1)] = complex(np.sum(tab * np.exp(-2j * np.pi * t / Lint)) / Lint ** 2)
    return out


# ---------------------------------------------------------------- the construction of the proof (:3008-3230)
class Setup:
    def __init__(self, Lelt):
        self.L = Lelt
        n = norm(Lelt)
        while n % 2 == 0: n //= 2
        while n % 3 == 0: n //= 3
        assert n == 1, "L must be supported on S = {2, lambda}"
        assert divides(TWO, Lelt) and divides(LAM, Lelt)
        self.M = mul(powe(LAM, 12), powe(Lelt, 4)); self.M2 = mul(self.M, self.M)
        self.VM = {LAM: val(self.M, LAM), TWO: val(self.M, TWO)}

    def build(self, h0, primes, ks, raw=False, u0shift=ZERO, c0_unit=ONE):
        """h_p = M^2 k_p (raw=True: h_p = k_p, negative control).  Returns a dict with everything.
        c0_unit != 1 implements the ALTERNATIVE unit reading c_0' = c_0 * unit (control only)."""
        L, M2 = self.L, self.M2
        r = prod(primes)
        hp = list(ks) if raw else [mul(M2, k) for k in ks]
        num = mul(LAM2, add(mul(h0, r), mul(L, esum([mul(h, exact_div(r, p)) for h, p in zip(hp, primes)]))))
        den = mul(L, r)
        g = gcd(num, den)
        a, c = exact_div(num, g), exact_div(den, g)
        if divides(LAM, c): u = [u for u in UNITS if is_primary(mul(a, u))]
        else: u = [u for u in UNITS if is_primary(mul(c, u))]
        assert len(u) == 1
        a, c = mul(a, u[0]), mul(c, u[0])
        c0 = exact_div(c, r)
        assert c0 is not None, "r does not divide c"
        c0 = mul(c0, c0_unit)
        # delta' (:3056-3060)
        pairs = []
        for q in (LAM, TWO):
            if divides(q, c0):
                m = powe(q, self.VM[q] + val(c0, q)); pairs.append((inv_mod(a, m), m))
            else:
                pairs.append((ZERO, powe(q, self.VM[q])))
        if norm(r) > 1: pairs.append((inv_mod(a, r), r))
        dl = crt(pairs)
        bg = exact_div(sub(mul(a, dl), ONE), c)
        assert bg is not None
        vl = val(c, LAM)
        if vl >= 2:
            case = "3|c"; u0 = None; g1 = (a, bg, c, dl); sigma = "0"
        elif vl == 1:
            case = "v=1"
            u0 = [u for u in (LAM, neg(LAM)) if mod3(sub(u, c)) == (0, 0)]
            assert len(u0) == 1; u0 = u0[0]
            g1 = (sub(a, mul(bg, u0)), bg, sub(c, mul(dl, u0)), dl)
            sigma = "-" if u0 == LAM else "+"
        else:
            case = "coprime"
            u0 = add(mod3(a), mul((3, 0), u0shift))
            g1 = (neg(bg), add(a, mul(bg, u0)), neg(dl), add(c, mul(dl, u0)))
            sigma = {0: "0", 1: "-", 2: "+"}[u0[1] % 3]
        D0 = mul(LAM3, c0)
        rinv = inv_mod(r, D0) if norm(D0) > 1 else ONE
        sig, eps = [], []
        for p in primes:
            cp = exact_div(c, p)
            s = mul(LAM2, cp); sig.append(s)
            eps.append(neg(inv_mod(mul(mul(LAM3, cp), s), p)))
        return dict(h0=h0, primes=list(primes), r=r, hp=hp, a=a, c=c, c0=c0, dl=dl, bg=bg, case=case, u0=u0,
                    g1=g1, sigma=sigma, D0=D0, rinv=rinv, sig=sig, eps=eps)


def kappa_literal(cfg):
    """exponent of kappa(g_1) = (c_1/a_1)_3 (:3152-3154), from the matrix g_1 itself"""
    a1, b1, c1, d1 = cfg["g1"]
    assert [mod3(e) for e in cfg["g1"]] == [(1, 0), (0, 0), (0, 0), (1, 0)], "g_1 != I mod 3"
    if c1 == ZERO: return 0
    k = jacobi3(c1, a1)
    assert k is not None
    return k

def kappa_formula(cfg):
    """(exponent of kappa(g_1), exponent of kappa_0) from the right side of eq:ray-multiplier (:3175-3183)"""
    a, c0, r, bg, u0 = cfg["a"], cfg["c0"], cfg["r"], cfg["bg"], cfg["u0"]
    ar = jacobi3(a, r)
    if cfg["case"] == "3|c":
        k0 = jacobi3(c0, a)
    elif cfg["case"] == "v=1":
        k0 = jacobi3(neg(u0), sub(a, mul(u0, bg))) + jacobi3(exact_div(c0, u0), a)
    else:
        # (a/c_0)_3 depends only on the ideal (c_0); with c_0 = c/r it is primary already.  The primary
        # associate is taken only so that the unit-reading controls (c_0 * unit) are defined in this case.
        k0 = jacobi3(a, unit_lam_primary(c0)[2])
    return (k0 + ar) % 3, k0 % 3


def Cprime(cfg, phihat, jdict, inactive=()):
    """the scalar of :3420-3424 WITHOUT the factor -(i/81) conj(alpha(c))^2:
       phi-hat(h_0) conj(kappa_0) prod_{inactive}(1 - 1/N) prod_{active} chi_p(sigma_p)^-2 omega_{p,j_p}"""
    _, k0 = kappa_formula(cfg)
    val_ = phihat * W3 ** ((-k0) % 3)
    for p in inactive: val_ *= 1 - 1.0 / norm(p)
    for p, s, e in zip(cfg["primes"], cfg["sig"], cfg["eps"]):
        lp = loc(p)
        val_ *= Z6 ** ((-2 * lp.k(s)) % 6) * lp.omega(jdict[p], e)
    return val_

def alpha(c):
    z = cxf(c); return z / abs(z)

def Cfull(cfg, phihat, jdict, inactive=()):
    """C of :3420-3424"""
    return -1j / 81 * (alpha(cfg["c"]).conjugate()) ** 2 * Cprime(cfg, phihat, jdict, inactive)

def psi_char(cfg):
    """psi(x) = e(-delta'_0 r^-1 x / D_0)  (:3225)"""
    return AddChar(neg(mul(cfg["dl"], cfg["rinv"])), cfg["D0"])


# ---------------------------------------------------------------- theta series (floating point)
class Ser:
    """sum_m coef(m) v K_{1/3}(4 pi |l| v) e~(l z),  l = m / lambda^4 (lambda^4 = 9), plus const v^{2/3}"""
    def __init__(self, m, coef, const=0.0):
        m = np.asarray(m, dtype=np.int64); coef = np.asarray(coef, dtype=complex)
        ell = (m[:, 0] + m[:, 1] * W3) / 9.0
        o = np.argsort(np.abs(ell), kind="stable")
        self.m, self.coef, self.ell, self.const = m[o], coef[o], ell[o], const
        self.aell = np.abs(self.ell)
    def cut(self, v, T):
        return int(np.searchsorted(self.aell, T / (4 * math.pi * v), side="right"))
    def with_coef(self, f, const=None):
        s = Ser.__new__(Ser)
        s.m, s.ell, s.aell = self.m, self.ell, self.aell
        s.coef = self.coef * f
        s.const = self.const if const is None else const
        return s
    def __call__(self, z, v, T=40.0, complete_to=None, phase=None, deriv=None, absum=False):
        """phase: optional AddChar applied to m (an exact rational shift).  complete_to: the |l| up to which
        the table is complete; an error is raised if the truncation T needs more."""
        from scipy.special import kv
        n = self.cut(v, T)
        if complete_to is not None and T / (4 * math.pi * v) > complete_to * (1 + 1e-12):
            raise RuntimeError("table too short: need |l| <= %.3f, complete to %.3f" % (T / (4 * math.pi * v), complete_to))
        ell = self.ell[:n]
        terms = self.coef[:n] * v * kv(1.0 / 3, 4 * math.pi * self.aell[:n] * v)
        if absum:
            t = np.abs(terms)
            if deriv: t = t * 2 * math.pi * self.aell[:n]
            return (0.0 if deriv else abs(self.const) * v ** (2.0 / 3)) + float(np.sum(t))
        terms = terms * np.exp(4j * math.pi * (ell * z).real)
        if phase is not None: terms = terms * phase(self.m[:n])
        if deriv == "zbar": return complex(np.sum(terms * 2j * math.pi * np.conj(ell)))
        if deriv == "z": return complex(np.sum(terms * 2j * math.pi * ell))
        return self.const * v ** (2.0 / 3) + complex(np.sum(terms))

def _to_ser(s):
    """kubota.Series -> Ser (index m = lambda^4 l)"""
    num = s.num
    if s.den == 3:
        m = np.stack([num[:, 0] * LAM[0] - num[:, 1] * LAM[1], num[:, 0] * LAM[1] + num[:, 1] * LAM[0] - num[:, 1] * LAM[1]], axis=1)
    else:
        m = num
    return Ser(m, s.coef, s.const)

class Theta:
    """theta and the three sequences d_sigma (sigma in '0', '+', '-') as Ser objects.
       complete['theta'], complete['cusp']: the |l| up to which the tables are complete."""
    def __init__(self, Nmax):
        Mmax = math.sqrt(Nmax / 27.0)
        self.K = kb.Kubota(Nmax, Mmax, CACHE)
        self.Nmax = Nmax
        self.theta = _to_ser(self.K.theta)
        self.d = {"0": _to_ser(self.K.theta.conj()), "+": _to_ser(self.K.F10.conj()), "-": _to_ser(self.K.F19.conj())}
        self.complete = {"0": Mmax, "+": math.sqrt(Nmax) / 9.0, "-": math.sqrt(Nmax) / 9.0}
    def thetabar_at(self, z, v, shift=None, T=40.0):
        """conj theta(z + A/C, v) straight from the series of conj theta; shift = (A, C) exact (reduced or not)"""
        ph = None
        if shift is not None:
            A, C = shift
            ph = AddChar(A, mul(LAM3, C))             # e~(l A/C) = e(m A/(lambda^3 C)),  m = lambda^4 l
        return self.d["0"](z, v, T=T, complete_to=self.complete["0"], phase=ph)


def cusp_point(c, z, v):
    """the part of g^-1(z + a/c, v) that does not depend on delta' (:3166-3170):
       z'' = -conj(z)/(c^2 (v^2+|z|^2)),  v' = v/(N(c)(v^2+|z|^2))"""
    cc = cxf(c); q = v * v + abs(z) ** 2
    return -z.conjugate() / (cc * cc * q), v / (abs(cc) ** 2 * q)


def dual_value(th, cfg, jd, z2, v2, T=36.0, variant=None, deriv=None):
    """sum over l of d_sigma(l) psi(lambda^4 l) prod_{p active} B_{p,j_p}(lambda^4 l) v2 K_{1/3}(4 pi |l| v2) e~(l z2)
    plus, for sigma = 0, the constant mode (3^{5/2}/2) v2^{2/3} prod_p B_{p,j_p}(0)   [the l = 0 term of the same
    local identity; it is nonzero only if every active prime has j = 4].
    variant (negative controls): 'swap' = d_+ <-> d_- ; 'psi_no_rinv' = psi without r^-1 ; 'Bsign' = chi^{+(j+2)}
    in B for j != 0, 4.   deriv = 'z': termwise d/dz (the constant mode drops)."""
    from scipy.special import kv
    sg = cfg["sigma"]
    if variant == "swap": sg = {"0": "0", "+": "-", "-": "+"}[sg]
    s = th.d[sg]
    need = T / (4 * math.pi * v2)
    if need > th.complete[sg] * (1 + 1e-12):
        raise RuntimeError("table too short: need |l| <= %.3f, complete to %.3f (sigma %s)" % (need, th.complete[sg], sg))
    n = s.cut(v2, T); m = s.m[:n]
    if variant == "psi_no_rinv": ps = AddChar(neg(cfg["dl"]), cfg["D0"])
    else: ps = psi_char(cfg)
    coef = s.coef[:n] * ps(m); c0 = SIGMA0 if sg == "0" else 0.0
    for p in cfg["primes"]:
        lp = loc(p); j = jd[p]; bt = lp.Btab(j)
        if variant == "Bsign" and j not in (0, 4): bt = np.conj(bt)
        coef = coef * bt[lp.idx_arr(m)]; c0 = c0 * bt[0]
    terms = coef * v2 * kv(1.0 / 3, 4 * math.pi * s.aell[:n] * v2) * np.exp(4j * math.pi * (s.ell[:n] * z2).real)
    if deriv == "z": return complex(np.sum(terms * 2j * math.pi * s.ell[:n]))
    return c0 * v2 ** (2.0 / 3) + complex(np.sum(terms))
