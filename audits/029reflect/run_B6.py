# -*- coding: utf-8 -*-
"""Task B, part 6: the compactly supported bump weight of the old script (audits/029core/theta_sum_check.py).
   U(x) = exp(-1/(1-t^2)), t = (x-1.25)/0.75 on [0.5, 2].   T(X) = X^{-1/2} sum_A a_W(A) conj(alpha(A)) U(N(A)/X).
   Dual formula (validated in run_B.py for log-Gaussian weights) with the Mellin transform of U computed numerically.
   Purpose: does the dual side reproduce the 'noise' values 0.01-0.08 that the old script printed for X <= 32000?"""
import os, math, cmath
import numpy as np
from scipy.special import loggamma
from kubota import *

CACHE = os.environ.get("TMPDIR", "/tmp"); NMAX = 200000
K = Kubota(NMAX, math.sqrt(NMAX / 27.0), CACHE)
lam2 = mul(LAM, LAM)
vals = {}
for a in range(-40, 41):
    for b in range(-40, 41):
        if is_primary((a, b)) and norm((a, b)) <= 1500: vals.setdefault((a % 9, b % 9), cubic_symbol(LAM, (a, b)))
def phiW(n):
    key = (n[0] % 9, n[1] % 9); return W3 ** vals[key] if key in vals else 0.0
rows = []
for h0 in range(9):
    for h1 in range(9):
        s = 0
        for n0 in range(9):
            for n1 in range(9):
                f = phiW((n0, n1))
                if f != 0: s += f * cmath.exp(-2j * cmath.pi * omcoef(mul((n0, n1), (h0, h1))) / 9)
        s /= 81
        if abs(s) > 1e-13:
            a, cc = reduce_frac(mul(lam2, (h0, h1)), (9, 0)); rows.append((s, a, cc) + K.translate_parts(a, cc))
agg = {}
for c, a, cc, g, kexp, rep, t in rows:
    assert rep == "I" and norm(cc) == 9
    S = K.theta; cc_c = cx(cc); d = cx(g[1][1]); shift = -d / cc_c + cx(t)
    Fz = S.coef * 2j * math.pi * np.conj(S.mu) * np.exp(4j * math.pi * (S.mu * shift).real)
    ek = c * (-1.0 / cc_c ** 2) * np.conj(W3 ** kexp * Fz) / 9
    beta = 4 * math.pi * S.amu / 9
    for n0, n1, val in zip(S.num[:, 0].tolist(), S.num[:, 1].tolist(), (ek / beta ** 3).tolist()):
        agg[(n0, n1)] = agg.get((n0, n1), 0.0) + val
keys = [k for k, v in agg.items() if abs(v) > 1e-13]
E = np.array([agg[k] for k in keys]); NL = np.array([(k[0] ** 2 - k[0] * k[1] + k[1] ** 2) / 27.0 for k in keys])
B2 = (4 * math.pi) ** 2 * NL / 81
NLMAX = K.theta.amu[-1] ** 2

# Mellin transform of U on the line Re s = sigma, by the trapezoid rule in y = log t (U is C^infty with compact support)
SIG = -0.5
ny = 8192; y = np.linspace(math.log(0.5), math.log(2.0), ny); dy = y[1] - y[0]
tt = (np.exp(y) - 1.25) / 0.75
Uy = np.zeros(ny); m = np.abs(tt) < 1; Uy[m] = np.exp(-1.0 / (1 - tt[m] ** 2))
TM = 1500.0; dtau = 0.01; tau = np.arange(-TM, TM + dtau / 2, dtau)
Uh = np.empty(len(tau), dtype=complex)
for i0 in range(0, len(tau), 4000):
    tb = tau[i0:i0 + 4000]
    Uh[i0:i0 + 4000] = (Uy[None, :] * np.exp(np.outer(SIG + 1j * tb, y))).sum(axis=1) * dy
s = SIG + 1j * tau
logR = (1 - 2 * s) * math.log(2) + loggamma(5.0 / 3 - s) + loggamma(4.0 / 3 - s) + math.log(4 * 3 ** -2.5) \
       + s * math.log((2 * math.pi) ** 2 / 27) - loggamma(s + 1.0 / 3) - loggamma(s + 2.0 / 3)
base = Uh * np.exp(logR)
print("Mellin transform of the bump times the gamma quotient, |U^(s) R(s)| on Re s = %.1f:" % SIG)
for T0 in (0, 10, 30, 100, 300, 1000, 1500):
    i = np.argmin(np.abs(tau - T0)); print("   tau=%6.0f  |U^|=%.3e  |U^ R|=%.3e" % (T0, abs(Uh[i]), abs(base[i])))
def H(u):
    out = np.empty(len(u), dtype=complex); lu = np.log(u)
    for i0 in range(0, len(u), 500):
        out[i0:i0 + 500] = (np.exp(np.outer(lu[i0:i0 + 500], s)) * base[None, :]).sum(axis=1) * dtau / (2 * math.pi)
    return out
xs = np.array([1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9])
print("dual kernel of the bump, |H_U|(x), x = (2 pi)^4 N(ell) X/(27 N(c)^2):")
print("   " + "  ".join("x=%g: %.1e" % (x, h) for x, h in zip(xs, np.abs(H(xs * 27 / math.pi ** 2)))))

A = K.td.A; NA = (A[:, 0] ** 2 - A[:, 0] * A[:, 1] + A[:, 1] ** 2).astype(float)
Ac = A[:, 0] + A[:, 1] * W3; alb = np.conj(Ac / np.abs(Ac)); aW = np.conj(K.td.G) * W3 ** K.td.kl
def U(x):
    t = (x - 1.25) / 0.75; out = np.zeros_like(x); mm = np.abs(t) < 1; out[mm] = np.exp(-1.0 / (1 - t[mm] ** 2)); return out
old = {125: 5.413e-02, 250: 4.486e-03, 500: 8.200e-02, 1000: 5.619e-02, 2000: 2.800e-02, 4000: 2.196e-02, 8000: 1.170e-02, 16000: 4.108e-02, 32000: 1.136e-02}
print("\n      X     T direct (re)      T dual (re)        |diff|     old script |T_ang|   dual freq used   largest x used   |H_U| there")
for X in (1000, 2000, 4000, 8000, 16000, 32000, 64000, 100000):
    td_ = np.sum(aW * alb * U(NA / X)) / math.sqrt(X)
    u = X * B2; x = u * math.pi ** 2 / 27
    uu, inv = np.unique(np.round(u, 9), return_inverse=True)
    Hh = H(uu)
    tdu = np.sum(E * Hh[inv]) / math.sqrt(X)
    print("  %6d   % .8e   % .8e   %.1e   %s            %6d        %.2e      %.1e" % (X, td_.real, tdu.real, abs(td_ - tdu),
          ("%.3e" % old[X]) if X in old else "    -    ", len(E), x.max(), abs(Hh[-1])))
print("\n(The dual table stops at N(ell) = %.0f; |diff| is the truncation of the dual series at 'largest x used'." % NLMAX)
print(" It shrinks as X grows because the same table then reaches larger x.)")
