# -*- coding: utf-8 -*-
"""Gate controls G1-G5: does the tool kit reproduce the known automorphy of Patterson's theta?
Nothing downstream counts unless these pass.  Tolerance fixed in advance: relative error <= 1e-8."""
import sys, os, math, cmath, random
import numpy as np
from kubota import *

TOL = 1e-8
random.seed(20261111)
CACHE = os.environ.get("TMPDIR", "/tmp")
NMAX = 200000
K = Kubota(NMAX, math.sqrt(NMAX / 27.0), CACHE)
th, F10, F19 = K.theta, K.F10, K.F19
ck = K.td.checks
print("primary A = c d^3 with N(A) <= %d: %d ; primes with Gauss sums: %d" % (NMAX, len(K.td.A), ck["nprimes"]))
print("max | |g~(pi)| - 1 | = %.2e ;  max | g~(pi)^3 + pi/|pi| | = %.2e   (DR (1.8))" % (ck["abs"], ck["cube"]))
for c_elt, d in ck["mult"]:
    print("  twisted multiplicativity vs brute force, c = %s, N = %d: diff %.1e" % (c_elt, norm(c_elt), d))
print("theta series: %d terms, |mu| <= %.2f ; F10: %d terms, |mu| <= %.2f" % (len(th.mu), th.amu[-1], len(F10.mu), F10.amu[-1]))

nfail = 0
def report(name, lhs, rhs, expect_pass=True):
    global nfail
    err = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300)
    ok = (err <= TOL) == expect_pass
    if not ok: nfail += 1
    print("  %-58s lhs=% .12f%+.12fi  relerr=%.2e  %s" % (name, lhs.real, lhs.imag, err,
          ("OK" if ok else "FAIL") + ("" if expect_pass else " (negative control: must mismatch)")))
    return err

def rnd_w(vlo=0.5, vhi=1.2, zr=0.6):
    return complex(random.uniform(-zr, zr), random.uniform(-zr, zr)), random.uniform(vlo, vhi)

print("\nG1: SL2(Z) invariance and periods")
worst = 0
for i in range(4):
    z, v = rnd_w()
    worst = max(worst, report("theta(E w) = theta(w), w=(%.3f%+.3fi, %.3f)" % (z.real, z.imag, v), th(*mobius(mat_cx(EMAT), z, v)), th(z, v)))
    worst = max(worst, report("theta(w + 1) = theta(w)", th(z + 1, v), th(z, v)))
    worst = max(worst, report("theta(w + 3 omega) = theta(w)", th(z + 3 * W3, v), th(z, v)))
z, v = rnd_w()
report("theta(w + omega) = theta(w)   [not a period]", th(z + W3, v), th(z, v), expect_pass=False)
print("  G1 worst relative error: %.2e" % worst)

print("\nG2: theta(gamma w) = chi(gamma) theta(w) for gamma in Gamma_1(3), chi = (c/a)_3")
def rnd_small():
    return (random.randint(-2, 2), random.randint(-2, 2))
worst = 0; count = {0: 0, 1: 0, 2: 0}; tries = 0
while sum(count.values()) < 12 and tries < 20000:
    tries += 1
    s1, t1, t2 = rnd_small(), rnd_small(), rnd_small()
    if t1 == (0, 0): continue
    g = mat_mul(mat_mul(((ONE, (0, 0)), (mul((3, 0), t1), ONE)), Tmat(mul((3, 0), s1))), ((ONE, (0, 0)), (mul((3, 0), t2), ONE)))
    if g[1][0] == (0, 0): continue
    k = kubota(g)
    if count[k] >= 4: continue
    # choose w so that both w and gamma w have moderate height: w = gamma^{-1}(w2) with w2 near the "top" of the isometric sphere
    c = cx(g[1][0]); d = cx(g[1][1])
    z = -d / c + complex(random.uniform(-.3, .3), random.uniform(-.3, .3)) / abs(c); v = random.uniform(0.7, 1.3) / abs(c)
    z2, v2 = mobius(mat_cx(g), z, v)
    if min(v, v2) < 0.035: continue
    count[k] += 1
    e = report("gamma=%s chi=omega^%d v=%.3f v'=%.3f" % (str(g).replace(" ", ""), k, v, v2), th(z2, v2), W3 ** k * th(z, v))
    worst = max(worst, e)
    if k != 0 and count[k] == 1:
        report("   same gamma with conj(chi)", th(z2, v2), W3 ** (-k) * th(z, v), expect_pass=False)
print("  G2 worst relative error: %.2e  (chi values covered: %s)" % (worst, count))

print("\nG3: cusp expansions  theta(gamma_10 w) = F10*(w),  theta(gamma_19 w) = F19*(w)   (DR (5.13)-(5.15), App. A)")
worst = 0
for i in range(4):
    z, v = rnd_w(0.6, 1.4, 0.5)
    worst = max(worst, report("gamma_10, w=(%.3f%+.3fi, %.3f)" % (z.real, z.imag, v), th(*mobius(mat_cx(G10), z, v)), F10(z, v)))
    worst = max(worst, report("gamma_19", th(*mobius(mat_cx(G19), z, v)), F19(z, v)))
z, v = rnd_w(0.6, 1.4, 0.5)
report("gamma_10 against F19*   [swapped]", th(*mobius(mat_cx(G10), z, v)), F19(z, v), expect_pass=False)
print("  G3 worst relative error: %.2e" % worst)

print("\nG1': negative control, the prefactors e^{-2 pi i/9}, e^{+2 pi i/9} in (5.7) exchanged")
thc = build_series(K.td, math.sqrt(NMAX / 27.0), swap9=True)[0]
z, v = rnd_w()
report("swapped-prefactor series under E", thc(*mobius(mat_cx(EMAT), z, v)), thc(z, v), expect_pass=False)
print("  (note: replacing tau by conj(tau) is NOT a valid negative control: it gives conj(theta), also E-invariant)")

print("\nG5: translates  theta(z + a/c, v)  via coset decomposition  g = g1 * sigma * rep * T_t")
worst = 0; seen = {}
cs = [ONE, OMEGA, LAM, (3, 0), mul(LAM, (3, 0)), (9, 0), (-2, 0), (2, 3), (4, 3), mul(LAM, (2, 3)), mul((3, 0), (-2, 0)), (-5, 0), mul((2, 3), (4, 3))]
for c in cs:
    for rep in range(6):
        a = (random.randint(-6, 6), random.randint(-6, 6))
        if a == (0, 0): continue
        g0, _, _ = xgcd(a, c)
        if norm(g0) != 1: continue
        z = complex(random.uniform(-.4, .4), random.uniform(-.4, .4)) / abs(cx(c)); v = random.uniform(0.7, 1.1) / abs(cx(c))
        g, kexp, name, t = K.translate_parts(a, c)
        lhs = K.translate(a, c, z, v); rhs = th(z + cx(a) / cx(c), v)
        e = report("a/c=%s/%s rep=%s t=%s chi=w^%d" % (a, c, name, t, kexp), lhs, rhs)
        worst = max(worst, e); seen[name] = seen.get(name, 0) + 1
        if kexp != 0 and "neg" not in seen:
            seen["neg"] = 1
            K.conj_kappa = True
            report("   same with conj(chi)", K.translate(a, c, z, v), rhs, expect_pass=False)
            K.conj_kappa = False
print("  G5 worst relative error: %.2e ; coset types used: %s" % (worst, {k: v for k, v in seen.items() if k != "neg"}))

print("\nG4: X^{1/3} constant of the untwisted completed sum  X^{-1/2} sum_A g~(c)|d| U(N(A)/X)")
from scipy.special import gamma as Gam
from scipy.integrate import quad
def U(x):
    t = (np.asarray(x, dtype=float) - 1.25) / 0.75
    out = np.zeros_like(t); m = np.abs(t) < 1
    out[m] = np.exp(-1.0 / (1 - t[m] ** 2)); return out
Ut = quad(lambda x: float(U(np.array([x]))[0]) * x ** (-1.0 / 6), 0.5, 2.0, epsabs=1e-13)[0]
Rres = (2 * math.pi) ** (5.0 / 3) / (3 * Gam(2.0 / 3) * 27 ** (5.0 / 6))      # from DR (5.47) with r = 1
pred = Rres * Ut
A = K.td.A; NA = (A[:, 0] ** 2 - A[:, 0] * A[:, 1] + A[:, 1] ** 2).astype(float)
gt_d = np.conj(K.td.G)                                                          # g~(c)|d|
print("  predicted constant  Res * U~(5/6) = %.6f   (Res = %.6f, U~(5/6) = %.6f)" % (pred, Rres, Ut))
for X in (1000, 4000, 16000, 64000, 100000):
    S = np.sum(gt_d * U(NA / X)) / math.sqrt(X)
    print("  X=%6d   sum/X^(1/3) = %.6f%+.6fi   ratio to prediction = %.4f" % (X, (S / X ** (1 / 3)).real, (S / X ** (1 / 3)).imag, (S / X ** (1 / 3)).real / pred))

print("\nGATE: %s  (failures: %d, tolerance %.0e)" % ("PASS" if nfail == 0 else "FAIL", nfail, TOL))
