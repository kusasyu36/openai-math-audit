# -*- coding: utf-8 -*-
"""
003dual stage 6 (lead, Claude Opus 5.5, 2026-11-12): anchors for the NORMALIZATION of the theta-coefficient tables.

Why: both sides of layers D/E use the same tables (kubota.py).  Codex (review_w9_003dual (5)) noted that a common
scale on both sides cannot hide a wrong C, but a wrong RELATIVE normalization between direct coefficients and cusp
coefficients could.  The relative normalization tau <-> d_+- is tested pointwise in layer B (automorphy).  What was
not tested against the paper's own statements is the ABSOLUTE normalization of tau.  This script ties it to two
statements of the paper itself, without using Patterson 1977 or a second transcription of DR:

  N1  eq:intro-theta-coefficients (paper2.tex:362-366): the coefficient of conj(theta) at lambda^-3 n b^3
      (n squarefree primary, (nb, 6) = 1) equals 3^{5/2} |b| conj(chi_n(lambda)^2) gamma_2(n).
      gamma_2 by brute force from :310-313 and the cubic symbol by Euler's criterion (own code, as in
      lead_gamma2_anchor.py).  Table: kubota.Kubota(...).theta.conj(), all entries with N(lambda^3 l) <= 700.
  N2  eq:cubic-theta-definition (:1650-1656): constant term 3^{5/2}/2 v^{2/3}.  Automorphy under
      E = [[0,-1],[1,0]] mixes the constant term with the coefficients, so it fixes their RATIO.  Test: with the
      table as is, theta(E w) = theta(w) at 6 points; with every coefficient scaled by 1.001 (constant kept),
      and with the constant scaled by 1.001 (coefficients kept), the identity must fail at some point.
      This shows the table's scale is the one forced by the paper's constant term, IF the automorphic form with
      this multiplier and this constant term is unique (assumed, not tested: one-dimensionality of the space).

CRITERIA (fixed before running):  N1 max |table - formula| <= 1e-10, over >= 100 entries, including n with 2 prime
  factors and some b != 1;  control (gamma_4 instead of gamma_2) differs by > 1e-3.
  N2 relative error <= 1e-9 at all 6 points for the table;  each scaled version has relative error > 1e-6 at >= 1 point.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy python -B -u audits/003dual/stage6_lead_normalization.py > audits/003dual/stage6_lead_normalization.out
"""
import sys, os, math, cmath, random
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "029reflect")); sys.path.insert(0, here)
import numpy as np
import kubota as kb
from lead_gamma2_anchor import cub_own, gamma_brute  # noqa: E402  (module runs its own check on import; output kept)

print(__doc__)
NMAX = 200000
K = kb.Kubota(NMAX, math.sqrt(NMAX / 27.0), os.path.join(here, "cache"))
d0 = K.theta.conj()
# ---------------- N1
idx = {}
for k in range(len(d0.num)):
    if d0.den == 3: idx[(int(d0.num[k, 0]), int(d0.num[k, 1]))] = complex(d0.coef[k])
worst = 0.0; ctl = 0.0; cnt = {}; nb = 0; ntot = 0
for A, c in list(idx.items()):
    NA = kb.norm(A)
    if NA > 700 or NA == 0 or not kb.is_primary(A) or (A[0] % 2 == 0 and A[1] % 2 == 0): continue
    if kb.exact_div(A, kb.LAM) is not None: continue
    u, k, fac = kb.factor(A)
    if any(e % 3 == 2 for _, e in fac): continue          # not of the form n b^3 with n squarefree
    primes = [pr.elt for pr, e in fac if e % 3 == 1]
    n = (1, 0)
    for p in primes: n = kb.mul(n, p)
    babs = (NA / kb.norm(n)) ** (1.0 / 6)
    chi2_lam = sum(cub_own(kb.LAM, p) for p in primes) % 3      # chi_n(lambda)^2 = (lambda/n)_3 = omega^k
    g2 = gamma_brute(n, primes, 1) if primes else 1.0
    g4 = gamma_brute(n, primes, 2) if primes else 1.0
    f = 3 ** 2.5 * babs * kb.W3 ** ((-chi2_lam) % 3) * g2
    worst = max(worst, abs(c - f)); ctl = max(ctl, abs(c - 3 ** 2.5 * babs * kb.W3 ** ((-chi2_lam) % 3) * g4))
    cnt[len(primes)] = cnt.get(len(primes), 0) + 1; nb += babs > 1.0001; ntot += 1
ok1 = worst <= 1e-10 and ntot >= 100 and cnt.get(2, 0) > 0 and nb > 0
print("N1   %s  entries %d (by number of prime factors of n: %s ; b != 1: %d) ; max |table - 3^{5/2}|b| conj(chi_n(lambda)^2) gamma_2(n)| = %.2e"
      % ("PASS" if ok1 else "FAIL", ntot, dict(sorted(cnt.items())), nb, worst))
print("NC1  %s  control gamma_4: max difference %.2e (must be > 1e-3)" % ("PASS" if ctl > 1e-3 else "FAIL", ctl))
# ---------------- N2
random.seed(2026111206)
th = K.theta
E = kb.mat_cx(kb.EMAT)
def series_scaled(sc_coef, sc_const):
    return th.twisted(np.full(len(th.coef), sc_coef), const=th.const * sc_const)
S0, S1, S2 = th, series_scaled(1.001, 1.0), series_scaled(1.0, 1.001)
errs = {0: [], 1: [], 2: []}
for i in range(6):
    z = complex(random.uniform(-0.4, 0.4), random.uniform(-0.4, 0.4)); v = random.uniform(0.75, 1.1)
    z2, v2 = kb.mobius(E, z, v)
    for key, S in ((0, S0), (1, S1), (2, S2)):
        a, b = S(z2, v2), S(z, v); errs[key].append(abs(a - b) / max(abs(a), abs(b)))
ok2 = max(errs[0]) <= 1e-9 and max(errs[1]) > 1e-6 and max(errs[2]) > 1e-6
print("N2   %s  theta(E w) = theta(w): max relative error %.2e at 6 points ; coefficients x1.001: max %.2e ; constant x1.001: max %.2e (both must be > 1e-6)"
      % ("PASS" if ok2 else "FAIL", max(errs[0]), max(errs[1]), max(errs[2])))
print("STAGE 6: %s" % ("all criteria met" if ok1 and ctl > 1e-3 and ok2 else "SOME CRITERIA NOT MET"))
