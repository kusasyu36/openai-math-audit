# -*- coding: utf-8 -*-
"""Task B, part 7: turn the dual-formula PREDICTION at X = 32000 into a direct check.
   T(X) = X^{-1/2} sum_A a(A) conj(alpha(A)) exp(-(log(N(A)/X))^2/(4w)),  a_P = g~(c)|d|,  a_W = gamma_2(c)|d|.
   run_B.py is not modified; the dual-side code below is copied from it.

   CRITERION, FIXED BEFORE THE COMPUTATION:
   - direct sum over all primary A with N(A) <= NMAX_DIRECT = 1,000,000.
   - truncation error estimate of the direct sum:  E_trunc = X^{-1/2} * rho * m * int_{NMAX}^{inf} W(t/X) dt,
       rho = 0.403 (density of primary elements per unit norm, 2 pi/(9 sqrt 3)),
       m   = mean of |a(A)| over NMAX/2 < N(A) <= NMAX (measured; includes the zeros of non-theta-supported A in rho).
     This is an estimate without cancellation (an upper-bound heuristic, not a rigorous bound).
   - dual side: terms with x <= xcut kept (kernel below 1e-15 beyond, see run_B.out), quadrature error ~1e-13.
   - PASS iff |direct - dual| <= TOL := max(1e-9, 10 * E_trunc).   Prediction printed in run_B.out: w=0.1: |T_P| = 7.171e-06,
     |T_W| = 3.078e-04 ; w=0.05: |T_P| = 5.088e-04, |T_W| = 1.604e-03.
"""
import os, math, cmath, time
import numpy as np
from scipy.special import loggamma
from scipy.integrate import quad
from kubota import *

CACHE = os.environ.get("TMPDIR", "/tmp")
X = 32000.0
NMAX_DIRECT = 1000000
print(__doc__)
t0 = time.time()
K = Kubota(200000, math.sqrt(200000 / 27.0), CACHE)        # dual side: only small frequencies are needed
print("dual-side tool kit ready (%.0f s)" % (time.time() - t0)); t0 = time.time()

lam2 = mul(LAM, LAM); vals = {}
for a in range(-40, 41):
    for b in range(-40, 41):
        if is_primary((a, b)) and norm((a, b)) <= 1500: vals.setdefault((a % 9, b % 9), cubic_symbol(LAM, (a, b)))
def phiW(n):
    key = (n[0] % 9, n[1] % 9); return W3 ** vals[key] if key in vals else 0.0
def phiP(n): return 1.0 if (n[0] % 3 == 1 and n[1] % 3 == 0) else 0.0
def rows_of(phi, q):
    rows = []
    for h0 in range(q):
        for h1 in range(q):
            s = 0
            for n0 in range(q):
                for n1 in range(q):
                    f = phi((n0, n1))
                    if f != 0: s += f * cmath.exp(-2j * cmath.pi * omcoef(mul((n0, n1), (h0, h1))) / q)
            s /= q * q
            if abs(s) > 1e-13:
                a, cc = reduce_frac(mul(lam2, (h0, h1)), (q, 0)); rows.append((s, a, cc) + K.translate_parts(a, cc))
    return rows
def dual_terms(rows, nl_max):
    agg = {}
    for c, a, cc, g, kexp, rep, t in rows:
        cc_c = cx(cc); d = cx(g[1][1]); Nc = norm(cc); shift = -d / cc_c + cx(t)
        if rep in ("I", "P", "P2"):
            S = K.theta; rot = {"I": ONE, "P": OMEGA2, "P2": OMEGA}[rep]
        else:
            S = K.F10 if rep == "g10" else K.F19; rot = ONE
        n = np.searchsorted(S.amu, math.sqrt(nl_max) * 1.0001, side="right")
        num = np.stack([S.num[:n, 0] * rot[0] - S.num[:n, 1] * rot[1], S.num[:n, 0] * rot[1] + S.num[:n, 1] * rot[0] - S.num[:n, 1] * rot[1]], axis=1)
        ell = S.mu[:n] * cx(rot)
        Fz = S.coef[:n] * 2j * math.pi * np.conj(ell) * np.exp(4j * math.pi * (ell * shift).real)
        ek = c * (-1.0 / (cc_c ** 2)) * np.conj(W3 ** kexp * Fz) / Nc
        beta = 4 * math.pi * np.abs(ell) / Nc
        d_ = agg.setdefault((Nc, S.den), {})
        for n0, n1, val in zip(num[:, 0].tolist(), num[:, 1].tolist(), (ek / beta ** 3).tolist()):
            d_[(n0, n1)] = d_.get((n0, n1), 0.0) + val
    E, B2 = [], []
    for (Nc, den), d_ in agg.items():
        for (n0, n1), val in d_.items():
            if abs(val) < 1e-13: continue
            nl = (n0 * n0 - n0 * n1 + n1 * n1) / 3.0 ** den
            E.append(val); B2.append((4 * math.pi) ** 2 * nl / Nc ** 2)
    return np.array(E), np.array(B2)
def H_kernel(u, w, sigma=-0.5):
    Tm = math.sqrt(75.0 / w); n = int(2 * Tm / 0.004) + 1
    tau = np.linspace(-Tm, Tm, n); s = sigma + 1j * tau; dt = tau[1] - tau[0]
    logR = (1 - 2 * s) * math.log(2) + loggamma(5.0 / 3 - s) + loggamma(4.0 / 3 - s) + math.log(4 * 3 ** -2.5) \
           + s * math.log((2 * math.pi) ** 2 / 27) - loggamma(s + 1.0 / 3) - loggamma(s + 2.0 / 3)
    base = np.exp(0.5 * math.log(4 * math.pi * w) + w * s * s + logR)
    out = np.empty(len(u), dtype=complex); lu = np.log(u)
    for i0 in range(0, len(u), 2000):
        out[i0:i0 + 2000] = (np.exp(np.outer(lu[i0:i0 + 2000], s)) * base[None, :]).sum(axis=1) * dt / (2 * math.pi)
    return out
XCUT = {0.1: 3.0e5, 0.05: 1.2e6}
ROWS = {"P": (rows_of(phiP, 3), 1), "W": (rows_of(phiW, 9), 81)}
pred = {}
for w in (0.1, 0.05):
    for lab in ("P", "W"):
        rows, nc2 = ROWS[lab]
        nl_max = XCUT[w] * 27 * nc2 / ((2 * math.pi) ** 4 * X)       # N(ell) needed for x <= xcut
        E, B2 = dual_terms(rows, nl_max * 1.01)
        u = X * B2; x = u * math.pi ** 2 / 27; sel = x <= XCUT[w]
        uu, inv = np.unique(np.round(u[sel], 9), return_inverse=True)
        pred[(w, lab)] = np.sum(E[sel] * H_kernel(uu, w)[inv]) / math.sqrt(X)
        print("dual prediction  w=%.2f  %s :  % .10e %+.3ei   (|T| = %.4e, %d frequencies, x in [%.3g, %.3g])"
              % (w, lab, pred[(w, lab)].real, pred[(w, lab)].imag, abs(pred[(w, lab)]), int(sel.sum()), x[sel].min(), x[sel].max()))
print("dual side done (%.0f s).  Now the direct side: enumerating A and Gauss sums up to N = %d ..." % (time.time() - t0, NMAX_DIRECT)); t0 = time.time()

td = ThetaData(NMAX_DIRECT, CACHE)
print("direct-side data ready (%.0f s): %d theta-supported primary A, %d primes ; max| |g~|-1 | = %.1e, max|g~^3 + pi/|pi|| = %.1e"
      % (time.time() - t0, len(td.A), td.checks["nprimes"], td.checks["abs"], td.checks["cube"]))
A = td.A; NA = (A[:, 0] ** 2 - A[:, 0] * A[:, 1] + A[:, 1] ** 2).astype(float)
Ac = A[:, 0] + A[:, 1] * W3; alb = np.conj(Ac / np.abs(Ac))
aP = np.conj(td.G); aW = aP * W3 ** td.kl
top = NA > NMAX_DIRECT / 2
rho = 2 * math.pi / (9 * math.sqrt(3))
m_abs = np.sum(np.abs(aP[top])) / (rho * NMAX_DIRECT / 2)       # mean |a| per primary element (zeros included)
print("mean |a(A)| per primary element in (NMAX/2, NMAX]: %.4f" % m_abs)
nfail = 0
print("\n   w     combo   T direct (re)        T dual (re)          |diff|     E_trunc    TOL        |T| direct   |T| predicted   verdict")
for w in (0.1, 0.05):
    tail = quad(lambda t: math.exp(-(math.log(t / X)) ** 2 / (4 * w)), NMAX_DIRECT, 50 * NMAX_DIRECT, epsabs=0, epsrel=1e-10)[0]
    Etr = rho * m_abs * tail / math.sqrt(X)
    TOL = max(1e-9, 10 * Etr)
    for lab, a in (("P", aP), ("W", aW)):
        td_ = np.sum(a * alb * np.exp(-(np.log(NA / X)) ** 2 / (4 * w))) / math.sqrt(X)
        diff = abs(td_ - pred[(w, lab)]); ok = diff <= TOL
        if not ok: nfail += 1
        print("  %.2f    %s     % .10e   % .10e   %.1e   %.1e   %.1e   %.4e   %.4e     %s"
              % (w, lab, td_.real, pred[(w, lab)].real, diff, Etr, TOL, abs(td_), abs(pred[(w, lab)]), "PASS" if ok else "FAIL"))
        print("          (imaginary parts: direct %.1e, dual %.1e ; sum of |terms| of the direct sum: %.3e)"
              % (td_.imag, pred[(w, lab)].imag, np.sum(np.abs(a) * np.exp(-(np.log(NA / X)) ** 2 / (4 * w))) / math.sqrt(X)))
print("\nB7: %s (failures: %d)" % ("PASS" if nfail == 0 else "FAIL", nfail))
