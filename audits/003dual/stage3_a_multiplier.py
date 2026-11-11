# -*- coding: utf-8 -*-
"""
Stage 3, layer A (EXACT, no floating point): eq:ray-multiplier (paper2.tex:3175-3183) against the matrix g_1.

What is compared.  For a configuration (L, h_0, active primes r, representatives h_p = M^2 k_p) the
construction of :3039-3089 gives a, c, delta', b_g, H, g_1 = g H^-1.  Left: kappa(g_1) = (c_1/a_1)_3 computed
from the entries of g_1 (:3152-3154).  Right: the three cases of eq:ray-multiplier.  Both are exponents mod 3.
Both sides use the same routine jacobi3 (cubic reciprocity + supplementary laws calibrated by brute force);
jacobi3 itself is validated in J0 against kubota.cubic_symbol (factorisation + Euler criterion).  What is
NOT tested here: that kappa is the multiplier of theta (that is layer B, stage3_b) -- this layer only tests
the algebraic identity "(c_1/a_1)_3 of this particular g_1 equals the displayed product of symbols".

CRITERIA (fixed before the first run):
  J0  jacobi3(al, be) == kubota.cubic_symbol(al, be) for every primary non-unit be with N(be) <= 400 and
      every al in the box [-8, 8]^2 (None = not coprime must also agree).  0 mismatches.
  A1  kappa_literal == kappa_formula in every configuration.  0 mismatches.
  A2  kappa_0 (formula) and kappa_literal - (a/r)_3 are both constant when the k_p vary with (L, h_0, r) fixed. 0 violations.
  A3  (a/r)_3 == prod_p chi_p(a)^2 with chi_p from check_local.Local (another table). 0 mismatches.
  A4  in the case (c, lambda) = 1, replacing u_0 by u_0 + 3t (t in {1, omega, -1-omega}) does not change
      kappa_literal ("a fixed set of representatives modulo 3", :3087).  0 changes.
  COV all three cases of H occur, and in each case kappa takes at least two different values (otherwise
      the comparison is vacuous there).
  Configurations: L in {18, 36, 54}; every h_0 in O/(L); r in {1, p7, p13, -5, p19, p7*p13}; 4 tuples (k_p)
      for one prime, 3 for two primes (1 for r = 1).
  Negative controls (expected to FAIL; a control counts only if it is "live", i.e. changes something):
  NU  the other unit reading: c_0' := denominator of lambda^2 h_0 / L normalised by the rule of :3047
      applied to lambda^2 h_0 / L alone (so c = unit * c_0' * r).  Expected: mismatches > 0 among the
      configurations where c_0' != c/r.  (If it showed 0 mismatches the two readings would be
      indistinguishable at this layer; that would be reported, not hidden.)
  NU2 c_0' = c_0 * u for each unit u != 1: mismatches > 0 for each u != -1... (-1 is a cube: (-1/a)_3 = 1,
      so u = -1 is expected to give 0 mismatches in the case 3 | c; reported per case).
  NM  representatives h_p = k_p without the factor M^2: expected A2 violations > 0 (kappa_0 not fixed).
  ND  case v_lambda(c) = 1 with the first factor (-u_0/(a - u_0 b_g))_3 dropped: mismatches > 0.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage3_a_multiplier.py > audits/003dual/stage3_a_multiplier.out 2>/dev/null
"""
import sys, time
from stage2_lib import *

t0 = time.time()
print(__doc__)
RES = []
def report(name, ok, detail):
    RES.append((name, ok)); print("  %-4s %s  %s" % (name, "PASS" if ok else "FAIL", detail)); sys.stdout.flush()

# ---- J0
mis = tot = 0
for a in range(-20, 21):
    for b in range(-21, 22, 3):
        be = (a, b)
        if a % 3 != 1 or not (1 < norm(be) <= 400): continue
        for x in range(-8, 9):
            for y in range(-8, 9):
                tot += 1
                if jacobi3((x, y), be) != kb.cubic_symbol((x, y), be): mis += 1
report("J0", mis == 0, "%d mismatches / %d pairs (al, be)" % (mis, tot))

P7, P13, M5, P19 = (1, 3), (4, 3), (-5, 0), cl.split_prime(19)
assert norm(P19) == 19
RS = [[], [P7], [P13], [M5], [P19], [P7, P13]]
KT = {0: [[]], 1: [[(1, 0)], [(2, 0)], [(1, 1)], [(3, -2)]], 2: [[(1, 0), (1, 0)], [(2, 0), (1, 1)], [(0, 1), (3, 0)]]}

def c0_alt_unit(setup, cfg):
    """unit u with c_0' = c_0 u, c_0' = the denominator of lambda^2 h_0/L normalised by the rule of :3047"""
    alt = setup.build(cfg["h0"], [], [])
    q = exact_div(alt["c"], cfg["c0"])
    assert q is not None and norm(q) == 1
    return q

stats = {}
def bump(key, n=1): stats[key] = stats.get(key, 0) + n
kvals = {"3|c": set(), "v=1": set(), "coprime": set()}
for Lint in (18, 36, 54):
    setup = Setup((Lint, 0))
    for primes in RS:
        for x0 in range(Lint):
            for x1 in range(Lint):
                h0 = (x0, x1); k0s = []; k0lit = []; k0raw = []
                for ks in KT[len(primes)]:
                    cfg = setup.build(h0, primes, ks)
                    kl = kappa_literal(cfg); kf, k0 = kappa_formula(cfg)
                    bump("A1.tot"); bump("case." + cfg["case"]); kvals[cfg["case"]].add(kl)
                    if kl != kf: bump("A1.bad")
                    ar = jacobi3(cfg["a"], cfg["r"]); k0s.append(k0); k0lit.append((kl - ar) % 3)
                    # A3
                    e = 0
                    for p in primes: e += loc(p).k(cfg["a"])        # chi_p(a)^2 = zeta6^(2k) = omega^k
                    bump("A3.tot")
                    if e % 3 != ar: bump("A3.bad")
                    # A4
                    if cfg["case"] == "coprime":
                        for t in ((1, 0), (0, 1), (-1, -1)):
                            c2 = setup.build(h0, primes, ks, u0shift=t)
                            bump("A4.tot")
                            if kappa_literal(c2) != kl: bump("A4.bad")
                    # NU
                    u = c0_alt_unit(setup, cfg)
                    if u != ONE:
                        bump("NU.live")
                        c2 = setup.build(h0, primes, ks, c0_unit=u)
                        if kappa_formula(c2)[0] != kappa_literal(c2): bump("NU.bad")
                        bump("NU.live." + cfg["case"])
                        if kappa_formula(c2)[0] != kappa_literal(c2): bump("NU.bad." + cfg["case"])
                    # NU2
                    for ui, uu in enumerate(UNITS[1:], 1):
                        c2 = setup.build(h0, primes, ks, c0_unit=uu)
                        bump("NU2.tot.%d.%s" % (ui, cfg["case"]))
                        if kappa_formula(c2)[0] != kappa_literal(c2): bump("NU2.bad.%d.%s" % (ui, cfg["case"]))
                    # ND
                    if cfg["case"] == "v=1":
                        k_drop = (jacobi3(exact_div(cfg["c0"], cfg["u0"]), cfg["a"]) + ar) % 3
                        bump("ND.tot")
                        if k_drop != kl: bump("ND.bad")
                    # NM
                    if primes:
                        cr = setup.build(h0, primes, ks, raw=True)
                        k0raw.append((cr["c"], (kappa_literal(cr) - jacobi3(cr["a"], cr["r"])) % 3))
                if primes:
                    bump("A2.tot")
                    if len(set(k0s)) != 1 or len(set(k0lit)) != 1: bump("A2.bad")
                    bump("NM.tot")
                    if len(set(k0raw)) != 1: bump("NM.bad")
        print("    L = %d, r = %s done [%.0fs]" % (Lint, primes, time.time() - t0)); sys.stdout.flush()

g = lambda k: stats.get(k, 0)
print()
report("A1", g("A1.bad") == 0, "%d mismatches / %d configurations" % (g("A1.bad"), g("A1.tot")))
report("A2", g("A2.bad") == 0, "%d violations / %d (L, h_0, r)" % (g("A2.bad"), g("A2.tot")))
report("A3", g("A3.bad") == 0, "%d mismatches / %d" % (g("A3.bad"), g("A3.tot")))
report("A4", g("A4.bad") == 0, "%d changes / %d" % (g("A4.bad"), g("A4.tot")))
cov = all(g("case." + c) > 0 and len(kvals[c]) >= 2 for c in kvals)
report("COV", cov, "configurations per case: %s ; values of kappa (exponents) seen per case: %s"
       % ({c: g("case." + c) for c in kvals}, {c: sorted(kvals[c]) for c in kvals}))
print("  negative controls:")
report("NU", g("NU.bad") > 0, "alternative unit reading: live (c_0' != c/r) in %d configurations, mismatches in %d ; per case live/bad: %s"
       % (g("NU.live"), g("NU.bad"), {c: (g("NU.live." + c), g("NU.bad." + c)) for c in kvals}))
for ui in range(1, 6):
    print("       NU2 c_0 * zeta6^%d : mismatches/total per case: %s" % (ui, {c: (g("NU2.bad.%d.%s" % (ui, c)), g("NU2.tot.%d.%s" % (ui, c))) for c in kvals}))
report("NU2", all(g("NU2.bad.%d.3|c" % ui) > 0 for ui in (1, 2, 4, 5)), "units zeta6^k, k = 1,2,4,5, give mismatches in the case 3 | c (k = 3 is -1, a cube)")
report("NM", g("NM.bad") > 0, "without M^2: (c, kappa_0) not constant in %d / %d (L, h_0, r)" % (g("NM.bad"), g("NM.tot")))
report("ND", g("ND.bad") > 0, "case v=1 without the first factor: %d mismatches / %d" % (g("ND.bad"), g("ND.tot")))
print("\nLAYER A: %s   [%.0fs]" % ("all criteria met" if all(ok for _, ok in RES) else "SOME CRITERIA NOT MET: " + ",".join(n for n, ok in RES if not ok), time.time() - t0))
