# -*- coding: utf-8 -*-
"""Task A (Q(omega) analogue): sextic-twisted cubic theta series  =  dual series with twist chi_p^{-j-2}.

Identity under test (my own derivation from DR's automorphy; the 003 analogue of 029 Lemma 6.2/6.3,
with no 3-adic projection, so that only the cusp at infinity occurs):

  Theta_{p,j}(z,v) := sum_{mu in lam^-3 O} conj(tau(mu)) chi_p^j(lam^3 mu) v K_{1/3}(4 pi |mu| v) e^(mu z)

  j in {1,2,3,5}:  Theta_{p,j}(w) = gamma_j gamma_{j+2} chi_p(lam)^{4j+4} * Theta_{p,-j-2}(w*)
  j = 4:           Theta_{p,4}(w) = gamma_4 chi_p(lam)^{-4} N^{-1/2} * [ sum conj(tau)(N 1_{p|n} - 1) ...(w*) + (N-1) sigma v*^{2/3} ]
  j = 0:           Theta_{p,0}(w) = (1 - 1/N) conj theta(w) - gamma_2 chi_p(lam)^4 N^{-1/2} * Theta_{p,-2}(w*)
                   (first run had a slip of mine here: I subtracted sigma v^{2/3} from the inactive term; the brute-force
                    layer agreed with the direct sum and exposed it.  Not an error of either paper.)
  w* = ( -conj(z) / (p^2 (|z|^2+v^2)),  v / (N(p)(|z|^2+v^2)) ),   n = lam^3 mu,
  gamma_j = N^{-1/2} sum_y chi_p(y)^j e^(y/(lam p)),  chi_p = sextic residue character.
Several primes r = prod p_i: the constant of p_i gets the extra factor chi_{p_i}(r/p_i)^{2 j_i + 2},  w* uses r.

Layer 1 (brute force): Theta = sum_h C_j(h) conj theta(z + lam^2 h/p), each translate through the coset machinery.
Layer 2 (closed form): the displayed formulas.
Tolerance fixed in advance: relative error <= 1e-8.
"""
import os, math, cmath, random, itertools
import numpy as np
from kubota import *

TOL = 1e-8
random.seed(2026111101)
CACHE = os.environ.get("TMPDIR", "/tmp")
NMAX = 200000
K = Kubota(NMAX, math.sqrt(NMAX / 27.0), CACHE)
thb = K.theta.conj()                     # series of conj(theta): coefficient conj(tau(-mu)) = conj(tau(mu))
nfail = 0; worst_pass = 0.0

def rel(lhs, rhs): return abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300)
def report(name, lhs, rhs, expect_pass=True):
    global nfail, worst_pass
    err = rel(lhs, rhs)
    ok = (err <= TOL) == expect_pass
    if not ok: nfail += 1
    if expect_pass: worst_pass = max(worst_pass, err)
    print("  %-46s lhs=% .10e%+.10ei  rhs=% .10e%+.10ei  relerr=%.1e  %s" % (name, lhs.real, lhs.imag, rhs.real, rhs.imag, err,
          ("OK" if ok else "FAIL") + ("" if expect_pass else " (negative control)")))
    return err

class Tw:
    """sextic twist data at a primary prime p"""
    def __init__(self, elt):
        self.pr = get_prime(elt); self.N = self.pr.N; self.c = cx(elt); self.P = self.pr.rat
        self.tab = {j: sextic_table(self.pr, j) for j in range(6)}
        self.gam = {j: gauss_gamma(self.pr, j) for j in range(1, 6)}
        self.div = np.zeros((self.P, self.P))
        for a in range(self.P):
            for b in range(self.P):
                if self.pr.sext((a, b)) is None: self.div[a, b] = 1
    def chi(self, x, j=1):
        k = self.pr.sext(x); return 0 if k is None else Z6 ** ((k * j) % 6)
    def on(self, num, j): return self.tab[j % 6][num[:, 0] % self.P, num[:, 1] % self.P]
    def divides(self, num): return self.div[num[:, 0] % self.P, num[:, 1] % self.P]
    def residues(self):
        if self.pr.inert: return [(a, b) for a in range(self.P) for b in range(self.P)]
        return [(a, 0) for a in range(self.P)]
    def ep(self, x):                       # e_p(x) = e^(x/(lam p))
        return cmath.exp(2j * cmath.pi * omcoef(mul(x, conj(self.pr.elt))) / self.N)
    def C(self, j, h):                     # finite Fourier coefficient of chi_p^j
        return sum(self.chi(y, j) * self.ep(neg(mul(h, y))) for y in self.residues() if self.chi(y, 0) != 0) / self.N

def wstar(z, v, relt):
    r = cx(relt); m = abs(z) ** 2 + v * v
    return -z.conjugate() / (r * r * m), v / (abs(r) ** 2 * m)

def local_factor(tw, j, rest, num):
    """per-prime dual factor: (constant, array over num).  rest = r/p (element)."""
    cross = tw.chi(rest, 2 * j + 2) if rest != ONE else 1.0
    if j % 6 in (1, 2, 3, 5):
        return tw.gam[j] * tw.gam[(j + 2) % 6] * tw.chi(LAM, 4 * j + 4) * cross, tw.on(num, -j - 2)
    if j % 6 == 4:
        # N^{-1/2} (N 1_{p|n} - 1); cross factor: chi_p(rest)^{-2} only (the n-dependence is unramified)
        cr = tw.chi(rest, -2) if rest != ONE else 1.0
        return tw.gam[4] * tw.chi(LAM, -4) * tw.N ** -0.5 * cr, tw.N * tw.divides(num) - 1.0
    cr = tw.chi(rest, 2) if rest != ONE else 1.0       # j = 0 active: chi(n rest^-2)^-2 chi(rest)^-2 = chi(n)^-2 chi(rest)^2
    return -tw.gam[2] * tw.N ** -0.5 * tw.chi(LAM, 4) * cr, tw.on(num, -2)

def direct(tws, js, z, v):
    f = np.ones(len(thb.num), dtype=complex)
    for tw, j in zip(tws, js): f = f * tw.on(thb.num, j)
    return thb.twisted(f, const=0.0)(z, v)

def closed(tws, js, z, v, dual_exps=None, scale_p_only=False):
    """closed-form dual side (all primes active; for a single prime with j=0 the inactive part is added by the caller)"""
    r = ONE
    for tw in tws: r = mul(r, tw.pr.elt)
    const = 1.0; f = np.ones(len(thb.num), dtype=complex); n0 = 1.0
    for i, (tw, j) in enumerate(zip(tws, js)):
        rest = exact_div(r, tw.pr.elt)
        c, arr = local_factor(tw, j, rest, thb.num)
        if dual_exps is not None: arr = tw.on(thb.num, dual_exps[i])
        const *= c; f = f * arr
        n0 *= (tw.N - 1.0) if (j % 6 == 4 and dual_exps is None) else 0.0      # value of the factor at n = 0
    zs, vs = wstar(z, v, r)
    if scale_p_only:
        m = abs(z) ** 2 + v * v; zs, vs = -z.conjugate() / (cx(r) * m), v / (abs(cx(r)) * m)
    return const * thb.twisted(f, const=SIGMA * n0)(zs, vs), const

def brute(tws, js, z, v):
    """sum over h of prod C_j(h_p) * conj theta(z + lam^2 sum h_p/p), via translates"""
    r = ONE
    for tw in tws: r = mul(r, tw.pr.elt)
    lam2 = mul(LAM, LAM); tot = 0; ntr = 0
    lists = []
    for tw, j in zip(tws, js):
        lists.append([(h, tw.C(j, h)) for h in tw.residues()])
    for combo in itertools.product(*lists):
        coef = 1.0; a = (0, 0)
        for (h, C), tw in zip(combo, tws):
            coef *= C; a = add(a, mul(h, exact_div(r, tw.pr.elt)))
        if abs(coef) < 1e-15: continue
        a = mul(lam2, a)
        ar, cr = reduce_frac(a, r)
        tot += coef * np.conj(K.translate(ar, cr, z, v)); ntr += 1
    return tot, ntr

def point(relt, scale=1.0):
    a = abs(cx(relt))
    return complex(random.uniform(-.45, .45), random.uniform(-.45, .45)) / a, random.uniform(0.55, 0.85) * scale / a

PRIMES = [(-2, -3), (1, 3), (4, 3), (-2, 3), (-5, 0), (-5, -6), (7, 3), (-11, 0)]
print("=== A1/A2: single prime, all j.  layer 2 = closed form, layer 1 = brute-force translates ===")
for elt in PRIMES:
    tw = Tw(elt)
    print("p = %s, N(p) = %d, chi_p(-1) = %+d, chi_p(lam) = zeta6^%d" % (elt, tw.N, round(tw.chi((-1, 0)).real), tw.pr.sext(LAM)))
    for j in range(6):
        z, v = point(elt)
        lhs = direct([tw], [j], z, v)
        rhs, const = closed([tw], [j], z, v)
        if j == 0: rhs = rhs + (1 - 1.0 / tw.N) * thb(z, v)       # inactive digit h = 0: the unreflected conj(theta), constant term included
        report("j=%d closed form, |const|=%.6f" % (j, abs(const)), lhs, rhs)
        if tw.N <= 31:
            b, ntr = brute([tw], [j], z, v)
            report("j=%d brute force (%d translates)" % (j, ntr), lhs, b)

print("\n=== negative controls (single prime, j = 1; dual twist must be chi^3) ===")
for elt in [(-2, -3), (4, 3), (-5, 0)]:
    tw = Tw(elt)
    pts = [point(elt) for _ in range(3)]
    for e in (3, 1, -1, 2, 4, 0):
        ratios = []
        for z, v in pts:
            lhs = direct([tw], [1], z, v); rhs, const = closed([tw], [1], z, v, dual_exps=[e])
            ratios.append(lhs / rhs)
        spread = max(abs(x - ratios[0]) for x in ratios) / abs(ratios[0])
        ok = (spread <= TOL and abs(ratios[0] - 1) <= TOL) == (e == 3)
        if not ok: nfail += 1
        print("  p=%s dual twist chi^%+d: lhs/rhs at 3 points = %s  spread=%.1e  %s" % (elt, e,
              " ".join("%.4f%+.4fi" % (x.real, x.imag) for x in ratios), spread,
              "OK (the correct twist)" if e == 3 and ok else ("OK (negative control)" if ok else "FAIL")))
    z, v = pts[0]
    report("p=%s inversion scale p instead of p^2" % (elt,), direct([tw], [1], z, v), closed([tw], [1], z, v, scale_p_only=True)[0], expect_pass=False)
    K.conj_kappa = True
    b, _ = brute([tw], [1], z, v)
    K.conj_kappa = False
    report("p=%s brute force with conj(chi_Kubota)" % (elt,), direct([tw], [1], z, v), b, expect_pass=False)

print("\n=== A3: two primes (cross terms) ===")
for e1, e2 in [((-2, -3), (4, 3)), ((-2, -3), (1, 3)), ((4, 3), (-2, 3)), ((-2, -3), (-5, 0))]:
    t1, t2 = Tw(e1), Tw(e2); r = mul(e1, e2)
    for js in [(1, 1), (1, 2), (5, 3), (4, 1), (2, 5), (0, 1)]:
        z, v = point(r)
        lhs = direct([t1, t2], js, z, v)
        if js[0] == 0:
            # inactive part at p1: (1 - 1/N1) * [twist by chi_{p2}^{j2} alone]
            rhs_act, const = closed([t1, t2], js, z, v)
            rhs_in, _ = closed([t2], [js[1]], z, v)
            rhs = rhs_act + (1 - 1.0 / t1.N) * rhs_in
        else:
            rhs, const = closed([t1, t2], js, z, v)
        report("r=%s*%s j=%s closed" % (e1, e2, js), lhs, rhs)
    if t1.N * t2.N <= 100:
        z, v = point(r)
        b, ntr = brute([t1, t2], (1, 1), z, v)
        report("r=%s*%s j=(1,1) brute (%d translates)" % (e1, e2, ntr), direct([t1, t2], (1, 1), z, v), b)
        # negative control: drop the cross term chi_p(r/p)^{2j+2}
        lhs = direct([t1, t2], (1, 1), z, v); rhs, const = closed([t1, t2], (1, 1), z, v)
        wrong = rhs / (t1.chi(e2, 4) * t2.chi(e1, 4))
        if abs(t1.chi(e2, 4) * t2.chi(e1, 4) - 1) > 1e-9:
            report("   without the cross factor", lhs, wrong, expect_pass=False)
        else:
            print("   (cross factor happens to be 1 for this pair; no negative control)")

print("\n=== A3': three primes, j = (1,1,1) ===")
tws = [Tw((-2, -3)), Tw((1, 3)), Tw((4, 3))]; r = mul(mul((-2, -3), (1, 3)), (4, 3))
z, v = 0.3 * cmath.exp(1j) / abs(cx(r)), 0.95 / abs(cx(r))
report("r=%s N=%d closed" % (r, norm(r)), direct(tws, (1, 1, 1), z, v), closed(tws, (1, 1, 1), z, v)[0])

print("\nworst relative error among identities expected to hold: %.2e" % worst_pass)
print("TASK A (Q(omega) analogue): %s  (failures: %d, tolerance %.0e)" % ("PASS" if nfail == 0 else "FAIL", nfail, TOL))
