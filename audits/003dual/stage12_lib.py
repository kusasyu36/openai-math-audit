# -*- coding: utf-8 -*-
"""
Stage 12, common code: a dual-tail bound that uses the STRUCTURE of the coefficients d_sigma instead of the pointwise
worst case 27 N(l)^{1/6}.  Importing runs nothing.  Written by the stage-12 auditor (Claude Opus 5.5) on 2026-11-19, AFTER
the stage-5 results (stage5_a_compact.out, stage5_a2_tight.out) had been seen.

PREMISE (not proved here, stated in the report): the coefficients d_sigma(l) for ALL l (also beyond the table) are given
by the formulas implemented in kubota.build_series (DR (5.7), (5.13), (5.14), App. A; paper2.tex :3104-3128, :1858-1867):
  sigma = +,- :  d(l) != 0 only if  lambda^4 l = u A,  u in a fixed set of THREE units (+-{1, omega, omega^2}, one sign
                 per sigma), A primary;   |d(l)| = 9 |G(A)|.
  sigma = 0   :  l = +- lambda^{3k-3} A (k >= 0),   |d(l)| = 3^{k/2 + 5/2} |G(A)|      [two units]
                 l = u lambda^{3k-4} A (k >= 1, u any of the 6 units),  |d(l)| = 3^{k/2 + 2} |G(A)|
  where A is primary (so coprime to lambda) and |G(A)| = N(b)^{1/2} if A = n b^3 with n squarefree (n, b primary),
  G(A) = 0 if some prime divides A to an exponent = 2 mod 3   (kubota.ThetaData: dabs = prod N(pr)^{(e//3)/2}).
  [paper :3123-3128 states |t_0| <= 3^{k/2+2}|b|, 3^{k/2+5/2}|b| and |t_+-| <= 9|b| at m = -4, with l = u lambda^m n b^3.]
  The lambda-adic valuation of l is fixed by (family, k), so the map (family, k, unit, A) -> l is injective.

DERIVATION (all steps are inequalities between sums of non-negative terms).
  D_sigma(s, R) := sum_{N(l) > R} |d_sigma(l)| N(l)^-s .  With N(l) = N(A)/81 for sigma = +-:
     D_+-(s, R) = 3 * 9 * 81^s * P(s, 81 R)
     D_0(s, R)  = sum_{k>=0} 2 * 3^{k/2+5/2} * 3^{-(3k-3)s} P(s, R/3^{3k-3}) + sum_{k>=1} 6 * 3^{k/2+2} * 3^{-(3k-4)s} P(s, R/3^{3k-4})
  P(s, M) := sum_{A primary, N(A) > M} |G(A)| N(A)^-s
          <= sum_{b primary} N(b)^{1/2 - 3s} Z(s, M / N(b)^3)          [A -> (n, b) injective; squarefreeness of n dropped]
  Z(s, M') := sum_{n primary, N(n) > M'} N(n)^-s
          <= sum_{n primary, M' < N(n) <= T0} N(n)^-s  (exact enumeration)  +  S(s, max(M', T0); d = 3)
  S(s, R; d) = bound (S)/(T3) of stage4_bounds (Abel summation with the Voronoi count of a lattice translate of minimal
  distance d; 1 + 3O is a translate of 3O, d = 3).  The b-sum: exact for N(b) <= B0, and for N(b) > B0
  Z(s, .) <= Zfull(s) := Z(s, 1/2) and sum_{N(b) > B0} N(b)^{1/2-3s} <= S(3s - 1/2, B0; 3).  The k-sums of D_0 are
  cut at k = KMAX; the rest is bounded by Zfull-type terms (every P <= Pfull(s) := P(s, 1/2)) times a geometric series.
  Group tail (beyond the table, N(l) > R):
     |sum_{N(l) > R} C d(l) alpha psi prod B V^sharp(N(l) X/N(c)^2)| <= |C| Bmax min_A (c X/N(c)^2)^-A Mt(A) D_sigma(1/2 + A, R)
  using |alpha psi| = 1, |prod B| <= Bmax, |V^sharp(x)| <= (c x)^-A Mt(A)  [(V') of stage5_bump, A > 0].
NON-RIGOROUS STEPS (evaluation, not derivation): Mt(A) and the envelope are numerical quadratures (stage5_bump.py); see
  Mt_upper below for a one-sided Riemann sum over the tau grid; the tau-integral beyond the grid end and the 800-node
  Gauss-Legendre integral inside the envelope are NOT bounded.  Floating-point sums are not bounded.
"""
import math
import numpy as np
import stage4_bounds as sb

T0 = 2.0e4      # exact enumeration of primary n up to this norm
B0 = 2.0e4      # exact enumeration of primary b up to this norm
KMAX = 40

def primaries(T):
    """norms of all primary elements (a = 1 mod 3, b = 0 mod 3) with norm <= T, sorted"""
    R = int(math.sqrt(4 * T / 3)) + 3
    a = np.arange(-R, R + 1); a = a[a % 3 == 1]
    b = np.arange(-R - (-R % 3), R + 1, 3)
    A, B = np.meshgrid(a, b, indexing="ij")
    N = (A * A - A * B + B * B).ravel()
    return np.sort(N[(N <= T) & (N > 0)]).astype(float)

_NP = primaries(max(T0, B0))

_SUF = {}
def _suffix(s):
    """suffix sums of N(n)^-s over primary n with N(n) <= T0 (sorted ascending): suf[i] = sum_{j >= i}"""
    key = round(s, 12)
    if key not in _SUF:
        Np = _NP[_NP <= T0]
        v = Np ** (-s)
        _SUF[key] = np.concatenate([np.cumsum(v[::-1])[::-1], [0.0]])
    return _SUF[key]

def Z_vec(s, Mp):
    """upper bound for sum_{n primary, N(n) > Mp} N(n)^-s, vectorised over Mp"""
    Mp = np.atleast_1d(np.asarray(Mp, dtype=float))
    Np = _NP[_NP <= T0]; suf = _suffix(s)
    i = np.searchsorted(Np, Mp, side="right")          # first index with N > Mp
    ex = suf[i]
    Rr = np.maximum(Mp, T0)
    rho = 3.0 / math.sqrt(3); cov = math.sqrt(3) * 9.0 / 2
    S = math.pi / cov * s * (Rr ** (1 - s) / (s - 1) + 2 * rho * Rr ** (0.5 - s) / (s - 0.5) + rho * rho * Rr ** (-s) / s)
    return ex + S

def Z(s, Mp):
    return float(Z_vec(s, Mp)[0])

_NB = None
def P(s, M):
    """upper bound for sum_{A primary, N(A) > M} |G(A)| N(A)^-s"""
    global _NB
    if _NB is None: _NB = np.unique(_NP[_NP <= B0], return_counts=True)
    un, cnt = _NB
    tot = float(np.sum(cnt * un ** (0.5 - 3 * s) * Z_vec(s, M / un ** 3)))
    tot += sb.S_pow(3 * s - 0.5, B0, 3.0) * Z(s, 0.5)
    return tot

def D_pm(s, R):
    return 3 * 9 * 81.0 ** s * P(s, 81 * R)

def D_0(s, R):
    tot = 0.0; Pf = P(s, 0.5)
    for k in range(0, KMAX + 1):
        tot += 2 * 3 ** (k / 2 + 2.5) * 3.0 ** (-(3 * k - 3) * s) * P(s, R / 3.0 ** (3 * k - 3))
    for k in range(1, KMAX + 1):
        tot += 6 * 3 ** (k / 2 + 2) * 3.0 ** (-(3 * k - 4) * s) * P(s, R / 3.0 ** (3 * k - 4))
    # k > KMAX: P <= Pf; ratio of consecutive terms is 3^{1/2 - 3s} < 1
    q = 3.0 ** (0.5 - 3 * s)
    k = KMAX + 1
    tot += Pf * (2 * 3 ** (k / 2 + 2.5) * 3.0 ** (-(3 * k - 3) * s) + 6 * 3 ** (k / 2 + 2) * 3.0 ** (-(3 * k - 4) * s)) / (1 - q)
    return tot

def D_sigma(sg, s, R):
    return D_0(s, R) if sg == "0" else D_pm(s, R)

_DC = {}
def D_cached(sg, s, R):
    key = (sg, round(s, 9), round(R, 6))
    if key not in _DC: _DC[key] = D_sigma(sg, s, R)
    return _DC[key]

def struct_tail(mtf, X, sg, Nc, R, Cabs, Bmax, A_grid=None):
    """structured tail bound of one group; mtf(A) = the Mt bound to use; returns (bound, A*)"""
    best = (float("inf"), None)
    for A in (A_grid or bump_A_grid()):
        Mv = mtf(A)
        if Mv == float("inf"): continue
        b = math.exp(-A * math.log(sb.CK * X / Nc ** 2)) * Mv * D_cached(sg, 0.5 + A, R)
        if b < best[0]: best = (b, A)
    return (Cabs * Bmax * best[0], best[1])

def bump_A_grid():
    return [0.25 * k for k in range(4, 49)]

class MtUp:
    """cached Mt_upper for one Bump"""
    def __init__(self, bump): self.b = bump; self.c = {}
    def __call__(self, A):
        if A not in self.c: self.c[A] = Mt_upper(self.b, A)
        return self.c[A]

def Mt_upper(bump, A):
    """one-sided Riemann sum of (a/pi) int Env(A,t) |Q(A+it)| dt over the tau grid of the Bump object:
    on [t_i, t_{i+1}] uses Env(t_i) (the computed Env is non-increasing in tau: every Gauss-Legendre term is) and
    max(|Q(t_i)|, |Q(t_{i+1})|) (|Q| checked monotone separately).  First cell [0, t_1] uses Env(0) max|Q|.
    Beyond the last grid point: NOT bounded (the integrand there is checked to be < e^-69 of its maximum)."""
    le = bump.log_env(A); t = bump.taus
    lq = np.array([sb.logQ(complex(A, x)).real for x in t])
    lqm = np.maximum(lq[:-1], lq[1:])
    f = le[:-1] + lqm + math.log(bump.a / math.pi)
    if le[-1] + lq[-1] > (le + lq).max() - 69.0: return float("inf")
    mx = f.max()
    if mx > 700: return float("inf")
    return float(np.sum(np.exp(f - mx) * np.diff(t))) * math.exp(mx)


# ---------------------------------------------------------------- added for stage12_main2 (after stage12_main.out was seen)
# (1) boundary: E.groups keeps table terms with |l| <= sqrt(R) by a floating-point comparison, so a term with N(l) = R
#     exactly may be in neither the table sum nor D_sigma(s, R) (which counts N(l) > R).  struct_tail2 evaluates D at
#     R (1 - 1e-9): every l with N(l) >= R is then counted by the bound (double counting only makes the bound larger).
# (2) D_paper: the same bound using ONLY the printed statement of the paper (:1858-1867): l = u lambda^m n b^3, m >= -4,
#     u any of the 6 units, |d(l)| <= 27 * 3^{m/6} |b|, N(l) = 3^m N(n b^3).  Reference column only.
def D_paper(s, R, mmax=120):
    tot = 0.0
    for m in range(-4, mmax + 1):
        tot += 6 * 27 * 3.0 ** (m / 6.0) * 3.0 ** (-m * s) * P(s, R / 3.0 ** m)
    q = 3.0 ** (1 / 6.0 - s); m = mmax + 1
    tot += 6 * 27 * 3.0 ** (m / 6.0) * 3.0 ** (-m * s) * P(s, 0.5) / (1 - q)
    return tot

_DP = {}
def struct_tail2(mtf, X, sg, Nc, R, Cabs, Bmax, paper=False, A_grid=None):
    Re = R * (1 - 1e-9)
    best = (float("inf"), None)
    for A in (A_grid or bump_A_grid()):
        Mv = mtf(A)
        if Mv == float("inf"): continue
        if paper:
            key = (round(0.5 + A, 9), round(Re, 6))
            if key not in _DP: _DP[key] = D_paper(0.5 + A, Re)
            Dv = _DP[key]
        else:
            Dv = D_cached(sg, 0.5 + A, Re)
        b = math.exp(-A * math.log(sb.CK * X / Nc ** 2)) * Mv * Dv
        if b < best[0]: best = (b, A)
    return (Cabs * Bmax * best[0], best[1])
