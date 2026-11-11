# -*- coding: utf-8 -*-
"""
Stage 4: error bounds used by stages 2-3, written as UPPER BOUNDS (not estimates), and the weight transform
V_*^sharp of eq:theta-weight (paper2.tex:1752-1757).  Importing this file runs nothing; running it prints
the tables used to choose (X, w, cut-offs) for stage3_d BEFORE any comparison is made.

WEIGHT.  W(y) = exp(-(log y)^2/(4w)),  V_*(y) = sqrt(y) W(y)  (:1293).  DEVIATION FROM THE TEXT: the paper
assumes W smooth with COMPACT support in (0, inf) (:670, :1869); the log-Gaussian is not compactly supported.
The proof (:3284-3436) only uses that V_*-hat decays fast on vertical lines, which holds here, but this is
my reading, and the statement tested in stage3_d is the identity for this W, not for a compactly supported one.
  V_*-hat(s) = int V_*(x) x^s dx/x = sqrt(4 pi w) exp(w (s + 1/2)^2)                         [exact]
  V_*^sharp(x) = (1/2 pi i) int_(0) V_*-hat(-t) Q(t) (c x)^-t dt,   c = (2 pi)^4/27,
  Q(t) = Gamma(7/6+t) Gamma(5/6+t) / (Gamma(7/6-t) Gamma(5/6-t)),   |Q(i tau)| = 1.

BOUNDS (all are absolute errors).
 (K)  K_{1/3}(x) <= K_{1/2}(x) = sqrt(pi/(2x)) e^-x  (x > 0; K_nu increases with nu >= 0).
 (D)  |d_sigma(l)| <= 27 N(l)^{1/6} for all three sequences (derived from the formulas implemented in
      kubota.build_series: 3^{n/2+5/2}|d|, 3^{n/2+2}|d|, 9|d| with |g~| = 1; cf. eq:theta-coefficient-bound :1867).
      CHECKED on the tables below (max of |d(l)|/N(l)^{1/6}); it is an input, proved by hand, verified on the table.
 (L)  lattice count: a lattice with hexagonal cell of minimal distance d has at most
      U(t) = pi (sqrt t + d/sqrt 3)^2 / (sqrt 3 d^2/2) points l with N(l) = |l|^2 <= t.
      (each point owns its Voronoi hexagon, of circumradius d/sqrt 3, inside the disc of radius sqrt t + d/sqrt 3)
      d = 3 for primary elements 1 + 3O,  d = 1/9 for lambda^-4 O,  d = 1/sqrt 27 for lambda^-3 O.
 (S)  Abel summation: for f >= 0 decreasing on [R, inf):  sum_{N(l) > R} f(N(l)) <= int_R^inf U(t) (-f'(t)) dt
      = U(R) f(R) + int_R^inf U'(t) f(t) dt.
 (T1) tail of a theta-type series at height v beyond 4 pi |l| v > T (T >= 2):
      sum |d(l)| v K_{1/3}(4 pi |l| v) <= 27 v * sum_{N(l) > R} N^{1/6} sqrt(pi/(8 pi v sqrt N)) exp(-4 pi v sqrt N),
      R = (T/(4 pi v))^2, bounded by (S); the summand is a constant times N^{-1/12} exp(-4 pi v sqrt N), decreasing in N.
      The integral in (S) is cut at 400 R, i.e. 4 pi |l| v = 20 T, where the integrand is below e^{-700}.
 (T2) direct sum of stage3_d: sum_{N(A) > Nmax} |a(A)| W(N(A)/X) X^-1/2, |a(A)| = |d| <= N(A)^{1/6}; by (S) with
      f(t) = t^{1/6} exp(-(log(t/X))^2/(4w)), decreasing for t > X e^{w/3}.
 (V)  for every A > 0 (no pole of the integrand in 0 <= Re t <= A; the poles of Q are at t = -5/6 - k, -7/6 - k):
      |V_*^sharp(x)| <= (c x)^-A M(A),  M(A) = (sqrt(4 pi w)/(2 pi)) int exp(w((1/2 - A)^2 - tau^2)) |Q(A + i tau)| d tau.
      M(A) is computed by numerical quadrature (scipy quad, relative accuracy requested 1e-10): this number is
      an ESTIMATE of a quantity that enters an upper bound; it is the only non-rigorous step of (V).
 (T3) dual sum of a group beyond N(l) > R:  |C| Bmax 27 sum_{N(l) > R} N^{1/6 - 1/2} |V^sharp(N X / N(c)^2)|
      <= |C| Bmax 27 min_A { (c X/N(c)^2)^-A M(A) * S(1/3 + A, R) },  S(s, R) = bound (S) for f(t) = t^-s:
      S(s, R) = (pi/cov) s [ R^{1-s}/(s-1) + 2 rho R^{1/2-s}/(s-1/2) + rho^2 R^-s / s ],  s > 1, rho = d/sqrt 3.
 (Q)  quadrature of V^sharp by the trapezoidal rule with step h on [-Tm, Tm]:
      truncation <= (sqrt(4 pi w)/(2 pi)) e^{w/4} * 2 int_Tm^inf e^{-w tau^2} d tau     (|Q| = 1 on the line)
      discretisation <= [ (c x)^{a} M(-a) + (c x)^{-a} M(a) ] / (e^{2 pi a/h} - 1)  for 0 < a < 5/6
      (the integrand is analytic in |Im tau| < 5/6; Trefethen-Weideman, SIAM Rev. 56 (2014), Thm 5.1).
 (R)  floating-point rounding is NOT bounded; it is ESTIMATED as 1e-13 times the sum of the absolute values of
      the terms (labelled 'rounding allowance (estimate)' wherever it is used).
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage4_bounds.py > audits/003dual/stage4_bounds.out 2>/dev/null
"""
import math, sys
import numpy as np
from scipy.special import loggamma
from scipy.integrate import quad

CK = (2 * math.pi) ** 4 / 27

def logQ(t):
    return loggamma(7.0 / 6 + t) + loggamma(5.0 / 6 + t) - loggamma(7.0 / 6 - t) - loggamma(5.0 / 6 - t)

_M = {}
def M_of(A, w):
    """M(A) of (V); also valid for -5/6 < A < 0 (used in (Q))"""
    key = (round(A, 9), w)
    if key not in _M:
        pre = math.sqrt(4 * math.pi * w) / (2 * math.pi)
        # the integrand is normalised by its value at tau = 0 to avoid overflow; M = inf if exp overflows
        l0 = w * (0.5 - A) ** 2 + logQ(complex(A, 0.0)).real
        f = lambda tau: math.exp(min(700.0, -w * tau * tau + logQ(complex(A, tau)).real - logQ(complex(A, 0.0)).real))
        Tm = math.sqrt(800.0 / w)
        val, err = quad(f, 0, Tm, epsabs=0, epsrel=1e-10, limit=2000)
        _M[key] = 2 * pre * val * math.exp(l0) if l0 < 700 else float("inf")
    return _M[key]

A_GRID = [0.25 * k for k in range(1, 161)]
def vsharp_bound(x, w, extra_pow=0.0):
    """min over A of (c x)^-A M(A)   [bound (V)] ; returns (bound, A)"""
    best = (float("inf"), None)
    for A in A_GRID:
        Mv = M_of(A, w)
        if Mv == float("inf"): continue
        b = math.exp(-A * math.log(CK * x)) * Mv
        if b < best[0]: best = (b, A)
    return best

def vsharp(x, w, h=0.02):
    """V_*^sharp(x) for an array x, by the trapezoidal rule on Re t = 0 (bound (Q) for its error)"""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    Tm = math.sqrt(40.0 / w) + 1.0
    n = int(Tm / h) + 1; tau = np.arange(-n, n + 1) * h
    t = 1j * tau
    g = np.exp(w * (0.5 - t) ** 2 + logQ(t)) * (math.sqrt(4 * math.pi * w) / (2 * math.pi)) * h
    ly = np.log(CK * x); out = np.empty(len(x), dtype=complex)
    for i0 in range(0, len(x), 4000):
        out[i0:i0 + 4000] = np.exp(-1j * np.outer(ly[i0:i0 + 4000], tau)) @ g
    return out

def quad_error_bound(x, w, h=0.02, a=0.8):
    """bound (Q) at one x: (truncation, discretisation)"""
    Tm = math.sqrt(40.0 / w) + 1.0
    pre = math.sqrt(4 * math.pi * w) / (2 * math.pi)
    trunc = pre * math.exp(w / 4) * math.sqrt(math.pi / w) * math.erfc(math.sqrt(w) * Tm)
    y = CK * x
    disc = (y ** a * M_of(-a, w) + y ** (-a) * M_of(a, w)) / math.expm1(2 * math.pi * a / h)
    return trunc, disc

def U_count(t, d):
    return math.pi * (math.sqrt(t) + d / math.sqrt(3)) ** 2 / (math.sqrt(3) * d * d / 2)
def dU(t, d):
    return math.pi * (1 + d / math.sqrt(3) / math.sqrt(t)) / (math.sqrt(3) * d * d / 2)

def abel_tail(f, R, d, upper=None):
    """bound (S): U(R) f(R) + int_R^inf U'(t) f(t) dt  (f decreasing on [R, inf) is the caller's duty)"""
    up = upper if upper is not None else R * 1e8
    val, err = quad(lambda u: dU(math.exp(u), d) * f(math.exp(u)) * math.exp(u), math.log(R), math.log(up), epsabs=0, epsrel=1e-9, limit=2000)
    return U_count(R, d) * f(R) + val

def direct_tail_bound(X, w, Nmax):
    """bound (T2)"""
    assert Nmax > X * math.exp(w / 3)
    f = lambda t: t ** (1.0 / 6) * math.exp(-(math.log(t / X)) ** 2 / (4 * w))
    return abel_tail(f, Nmax, 3.0) / math.sqrt(X)

def S_pow(s, R, d):
    rho = d / math.sqrt(3); cov = math.sqrt(3) * d * d / 2
    assert s > 1
    return math.pi / cov * s * (R ** (1 - s) / (s - 1) + 2 * rho * R ** (0.5 - s) / (s - 0.5) + rho * rho * R ** (-s) / s)

def dual_tail_bound(X, w, Nc, R, Cabs, Bmax):
    """bound (T3) for one group; returns (bound, A)"""
    best = (float("inf"), None)
    for A in A_GRID:
        if A <= 0.75: continue
        Mv = M_of(A, w)
        if Mv == float("inf"): continue
        b = math.exp(-A * math.log(CK * X / Nc ** 2)) * Mv * S_pow(1.0 / 3 + A, R, 1.0 / 9)
        if b < best[0]: best = (b, A)
    return Cabs * Bmax * 27 * best[0], best[1]

def theta_tail_bound(v, T, d=1.0 / 9):
    """bound (T1): tail of sum |d(l)| v K_{1/3}(4 pi |l| v) over 4 pi |l| v > T, any of the three sequences
    (lambda^-4 O contains lambda^-3 O, so d = 1/9 covers theta too)"""
    R = (T / (4 * math.pi * v)) ** 2
    f = lambda N: 27 * v * N ** (1.0 / 6) * math.sqrt(math.pi / (8 * math.pi * v * math.sqrt(N))) * math.exp(-4 * math.pi * v * math.sqrt(N))
    return abel_tail(f, R, d, upper=R * 400)


if __name__ == "__main__":
    print(__doc__)
    print("(V) M(A) and the bound min_A (c x)^-A M(A), against |V^sharp(x)| computed by quadrature")
    print("    (the quadrature value has a rounding floor of about 3e-16 ABSOLUTE: it is a sum of ~2000 terms of size <= 4e-4;")
    print("     the check is  |V^sharp| <= bound(V) + (Q) + 2e-15,  the 2e-15 being a rounding allowance (estimate))")
    for w in (0.05, 0.1, 0.2):
        print("  w = %.2f:  M(A): " % w + "  ".join("A=%g: %.2e" % (A, M_of(A, w)) for A in (0.5, 1, 2, 4, 8, 12, 16, 20)))
        xs = np.array([1e-2, 1e-1, 1, 10, 1e2, 1e3, 3e3, 1e4, 3e4, 1e5])
        vs = np.abs(vsharp(xs, w))
        for x, vv in zip(xs, vs):
            b, A = vsharp_bound(x, w); tr, di = quad_error_bound(x, w)
            print("     x = %-8g  |V^sharp| = %.3e   bound (V) = %.3e (A = %g)   quadrature error bound (Q): truncation %.1e, discretisation %.1e  %s"
                  % (x, vv, b, A, tr, di, "" if vv <= b * (1 + 1e-9) + tr + di + 2e-15 else "  <-- VIOLATION"))
    print("\n(T2) direct tail bound, Nmax = 1e6")
    for w in (0.05, 0.1, 0.2):
        print("  w = %.2f: " % w + "  ".join("X=%g: %.1e" % (X, direct_tail_bound(X, w, 1e6)) for X in (4000, 8000, 16000, 24000, 32000, 48000, 64000)))
    print("\n(T3) dual tail bound per unit |C| Bmax, cut at R = N(l) limit of the tables (192.4^2 for sigma = 0, 111.1^2 for the cusps)")
    for w in (0.05, 0.1, 0.2):
        for Nc in (36, 63, 252, 28, 84, 52, 100, 156, 468):
            print("  w = %.2f  N(c) = %4d : " % (w, Nc) + "  ".join("X=%g: %.1e/%.1e" % (X, dual_tail_bound(X, w, Nc, 192.4 ** 2, 1, 1)[0], dual_tail_bound(X, w, Nc, 111.1 ** 2, 1, 1)[0]) for X in (4000, 8000, 16000, 32000, 64000)))
    print("\n(T1) theta tail bound at T = 36 and T = 40 (to be compared with the size of the functions, O(1)/sqrt v):")
    print("  (stage2 uses T = 40 with 0.05 <= v <= 0.35; layers B, D, E0 use T = 40 with v >= 0.97/|c| >= 0.0169 and v' >= 1/(1.051 |c|) >= 0.0166)")
    for v in (0.0166, 0.02, 0.05, 0.1, 0.3, 1.0):
        print("  v = %-5g  T=36: %.2e   T=40: %.2e" % (v, theta_tail_bound(v, 36.0), theta_tail_bound(v, 40.0)))
