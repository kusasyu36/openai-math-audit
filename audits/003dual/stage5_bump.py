# -*- coding: utf-8 -*-
"""
Stage 5 (a), tools: a COMPACTLY SUPPORTED smooth weight, its transform V_*^sharp (eq:theta-weight, paper2.tex:1752-1757)
and error bounds in the manner of stage4_bounds.py.  Importing runs nothing; running prints the planning tables
(bounds only, no comparison) from which the parameter a and the runs of stage5_a_compact.py are chosen.

WEIGHT (inside the hypothesis W in C_c^infty((0, inf)), :670, :1293):
  W(y) = b(log(y)/a),  b(u) = exp(-(1-u^2)^-p) for |u| < 1, 0 otherwise (p = 1, 2, 3, ...);  support [e^-a, e^a].
  V_*(y) = sqrt(y) W(y).   p = 1 is the standard bump.  [p > 1 was introduced after the first planning table (p = 1 only)
  showed dual tail bounds >= 5e-2 per unit |C| for every candidate run: the transform of the p = 1 bump decays like
  exp(-const |tau|^(1/2)) and V^sharp like exp(-const x^(1/8)).  For general p the rates are |tau|^(p/(p+1)), x^(p/(4(p+1))).
  All these functions are in C_c^infty; the choice of (p, a) is made from the bounds alone, before any comparison.]
  V_*-hat(s) = int V_*(x) x^s dx/x = a int_{-1}^{1} b(u) exp((s + 1/2) a u) du                                [exact formula]
  V_*^sharp(x) = (1/2 pi) int V_*-hat(-i tau) Q(i tau) (c x)^(-i tau) d tau,  c = (2 pi)^4/27,  Q as in stage4_bounds.

NUMERICS.  V_*-hat(-i tau): trapezoidal rule in u with NU nodes (b is smooth and flat at +-1; the aliasing error is the
  transform itself at frequency >= pi NU - a|tau|, far below rounding).  V^sharp: trapezoidal rule in tau, step H, |tau| <= TM.

BOUNDS (absolute errors).
 (ENV) rigorous envelope of the transform.  b extends holomorphically off the real axis; by Cauchy's theorem the u-integral
      may be taken along z = u - i eta (1 - u^2), sign(tau) eta > 0, for any eta (the integrand tends to 0 at +-1 inside the
      sector |arg(1 -+ z)| < pi/(2p); on the path arg(1 - z^2) -> arctan(2 eta) at the ends, so eta < tan(pi/(2p))/2 is
      required and the eta-grid is cut at 0.9 of that).  Hence for t = A + i tau
        |V_*-hat(-t)| <= a * Env(A, |tau|),
        Env(A, w) = min_eta int_{-1}^{1} exp( -Re (1 - z^2)^-p + (1/2 - A) a u - a w eta (1 - u^2) ) sqrt(1 + 4 eta^2 u^2) du.
      The remaining integral (smooth, positive) is done by Gauss-Legendre with 800 nodes; eta runs over a grid (any eta
      gives a valid bound, the grid only affects sharpness).  Non-rigorous step: that 1-D quadrature.
 (V') |V_*^sharp(x)| <= (c x)^-A Mt(A),  Mt(A) = (a/pi) int_0^inf Env(A, tau) |Q(A + i tau)| d tau   (A > 0; no pole crossed).
      Mt is integrated numerically on a log grid in tau up to TAUMAX (trapezoid in log tau), the integrand being below
      1e-30 of its maximum at TAUMAX (checked; otherwise Mt = inf).  Non-rigorous step: this quadrature.
 (T3') dual tail of a group beyond N(l) > R:  |C| Bmax 27 min_A { (c X/N(c)^2)^-A Mt(A) S(1/3 + A, R) }  (S as in stage4_bounds).
 (Q')  quadrature of V^sharp: truncation <= (a/pi) int_TM^inf Env(0, tau) d tau ;
      discretisation <= [ (c x)^a' Mt(-a') + (c x)^-a' Mt(a') ] / (e^{2 pi a'/H} - 1), a' = 0.8  (as (Q) of stage4_bounds).
 The direct sum has NO truncation error when X e^a <= Nmax: W vanishes beyond.
 Rounding is not bounded, only estimated (as in stage 4).
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage5_bump.py > audits/003dual/stage5_bump.out 2>/dev/null
"""
import math, sys
import numpy as np
import stage4_bounds as sb

CK = sb.CK
_GL = np.polynomial.legendre.leggauss(800)
A_GRID = [0.25 * k for k in range(4, 49)]          # 1.0 .. 12.0

class Bump:
    def __init__(self, a, p=1, NU=8001, H=0.04):
        self.a = a; self.p = p; self.NU = NU; self.H = H
        emax = 3.0 if p == 1 else 0.9 * math.tan(math.pi / (2 * p)) / 2
        self.ETAS = ETAS = np.concatenate([[0.0], np.exp(np.linspace(math.log(1e-3), math.log(emax), 60))])
        u, wu = _GL
        self.u = u; self.wu = wu
        one = 1 - u * u
        # for every eta: Re (1 - z^2)^-p, z = u - i eta (1-u^2):  1 - z^2 = one * (1 + eta^2 one + 2 i u eta)
        self.re = np.array([((one * (1 + e * e * one + 2j * u * e)) ** (-p)).real for e in ETAS])     # [eta, u]
        self.jac = np.array([np.sqrt(1 + 4 * e * e * u * u) for e in ETAS])
        self._env = {}; self._Mt = {}
        self.taus = np.concatenate([[0.0], np.exp(np.linspace(math.log(1e-2), math.log(3.0e5 / a), 700))])
        self.TM = None

    def W(self, y):
        uu = np.log(np.asarray(y, dtype=float)) / self.a
        out = np.zeros_like(uu); m = np.abs(uu) < 1
        out[m] = np.exp(-(1 - uu[m] ** 2) ** (-float(self.p)))
        return out

    def log_env(self, A):
        """log Env(A, tau) on self.taus"""
        if A not in self._env:
            a = self.a; u = self.u; one = 1 - u * u
            base = -self.re + (0.5 - A) * a * u[None, :] + np.log(self.jac * self.wu[None, :])        # [eta, u]
            out = np.empty(len(self.taus))
            for i, w in enumerate(self.taus):
                ex = base - a * w * self.ETAS[:, None] * one[None, :]
                mx = ex.max(axis=1)
                val = mx + np.log(np.exp(ex - mx[:, None]).sum(axis=1))
                out[i] = val.min()
            self._env[A] = out
        return self._env[A]

    def Mt(self, A):
        if A not in self._Mt:
            le = self.log_env(A); t = self.taus
            lq = np.array([sb.logQ(complex(A, x)).real for x in t])
            f = le + lq + math.log(self.a / math.pi)
            mx = f.max()
            if f[-1] > mx - 69.0: self._Mt[A] = float("inf")
            else:
                g = np.exp(f - mx)
                val = g[0] * t[1] + np.sum(0.5 * (g[1:-1] * t[1:-1] + g[2:] * t[2:]) * np.diff(np.log(t[1:])))     # [0, t1] + log-trapezoid
                self._Mt[A] = float(val) * math.exp(mx) if mx < 700 else float("inf")
        return self._Mt[A]

    def vsharp_bound(self, x):
        best = (float("inf"), None)
        for A in A_GRID:
            b = math.exp(-A * math.log(CK * x)) * self.Mt(A)
            if b < best[0]: best = (b, A)
        return best

    def dual_tail_bound(self, X, Nc, R, Cabs, Bmax):
        best = float("inf")
        for A in A_GRID:
            Mv = self.Mt(A)
            if Mv == float("inf"): continue
            best = min(best, math.exp(-A * math.log(CK * X / Nc ** 2)) * Mv * sb.S_pow(1.0 / 3 + A, R, 1.0 / 9))
        return Cabs * Bmax * 27 * best

    def setup_quadrature(self):
        """choose TM: the smallest grid tau with (a/pi) int_tau^inf Env(0, .) <= 1e-17, and tabulate G(tau) = V-hat(-i tau)"""
        le = self.log_env(0.0); t = self.taus
        g = np.exp(le) * self.a / math.pi
        tail = np.concatenate([np.cumsum((0.5 * (g[1:-1] * t[1:-1] + g[2:] * t[2:]) * np.diff(np.log(t[1:])))[::-1])[::-1], [0.0]])    # int_{t_i}^{tmax}, i >= 1
        i = int(np.argmax(tail <= 1e-17)); self.TM = float(t[1 + i]); self.trunc = float(tail[i])
        n = int(self.TM / self.H) + 1; self.tau = np.arange(-n, n + 1) * self.H
        uu = np.linspace(-1, 1, self.NU)[1:-1]; hu = 2.0 / (self.NU - 1)
        bb = np.exp(-(1 - uu * uu) ** (-float(self.p)) + 0.5 * self.a * uu) * hu * self.a
        G = np.empty(len(self.tau), dtype=complex)
        for i0 in range(0, len(self.tau), 500):
            G[i0:i0 + 500] = np.exp(-1j * self.a * np.outer(self.tau[i0:i0 + 500], uu)) @ bb
        self.G = G
        env = self.a * np.exp(np.interp(np.abs(self.tau), t, le)); ok = env > 1e-9 * env.max()     # below that G is rounding noise
        self.env_check = float(np.max(np.abs(G[ok]) / env[ok]))    # must be <= 1 (up to interpolation of the envelope in tau)

    def vsharp(self, x, invert=False):
        if self.TM is None: self.setup_quadrature()
        x = np.atleast_1d(np.asarray(x, dtype=float))
        lq = np.array([sb.logQ(1j * tt) for tt in self.tau])
        g = self.G * np.exp(-lq if invert else lq) * (self.H / (2 * math.pi))
        ly = np.log(CK * x); out = np.empty(len(x), dtype=complex)
        for i0 in range(0, len(x), 1000):
            out[i0:i0 + 1000] = np.exp(-1j * np.outer(ly[i0:i0 + 1000], self.tau)) @ g
        return out

    def quad_error_bound(self, x, ap=0.8):
        if self.TM is None: self.setup_quadrature()
        y = CK * x
        disc = (y ** ap * self.Mt(-ap) + y ** (-ap) * self.Mt(ap)) / math.expm1(2 * math.pi * ap / self.H)
        return self.trunc, disc


if __name__ == "__main__":
    print(__doc__)
    NMAX = 1.0e6
    # (table limit R, N(c) values, sum over groups of |C| Bmax (rough: 36 groups * (1/36..1/12)/81)) for candidate runs
    CAND = [("B {}", 111.11 ** 2, 4), ("B {p7:1}", 111.11 ** 2, 28), ("F {} / C {}", 111.11 ** 2, 12), ("F {p7:1}", 111.11 ** 2, 84),
            ("A {}", 192.45 ** 2, 36), ("A {p7:1}", 192.45 ** 2, 252), ("B {p13:1}", 111.11 ** 2, 52)]
    for p, a in ((1, 3.0), (2, 2.0), (2, 3.0), (2, 4.0), (3, 2.0), (3, 3.0), (3, 4.0), (4, 3.0), (6, 3.0)):
        bp = Bump(a, p); bp.setup_quadrature()
        X = math.floor(NMAX * math.exp(-a) * 0.999)
        print("p = %d, a = %g : largest X with X e^a <= 10^6: %d ; TM = %.1f (truncation bound %.1e), %d tau nodes ; max |V-hat| / envelope where envelope > 1e-9 max = %.3f (must be <= 1)"
              % (p, a, X, bp.TM, bp.trunc, len(bp.tau), bp.env_check))
        print("   Mt(A): " + "  ".join("A=%g: %.2e" % (A, bp.Mt(A)) for A in (1, 2, 3, 4, 6, 8, 10, 12)))
        xs = np.array([1.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8])
        vs = np.abs(bp.vsharp(xs))
        print("   x: |V^sharp| / bound (V') : " + "  ".join("%g: %.1e / %.1e" % (x, v, bp.vsharp_bound(x)[0]) for x, v in zip(xs, vs)))
        for name, R, Nc in CAND:
            print("   %-12s N(c) = %3d: dual tail bound per unit |C| Bmax at X = %d : %.1e ; at X/4 : %.1e" % (name, Nc, X, bp.dual_tail_bound(X, Nc, R, 1, 1), bp.dual_tail_bound(X / 4, Nc, R, 1, 1)))
        sys.stdout.flush()
