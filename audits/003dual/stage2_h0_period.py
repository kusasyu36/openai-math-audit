# -*- coding: utf-8 -*-
"""
Counting check (written after all official runs; it changes no criterion).  Purpose: not to inflate counts.
For L = 18 the classes h_0, h_0 + 6, h_0 + 12 give shifts lambda^2 h_0/L differing by lambda^2*6/18 = -1 in Z,
a period of theta.  Checked here, exactly or to 1e-13:
  (1) phi-hat(h_0 + 6) = phi-hat(h_0) for all h_0 and every Psi_0 in A..F        (so the three classes are one translate
      with one coefficient; "324 classes" = 108 distinct translates, "36 classes in the support" = 12)
  (2) h_0 -> h_0 + 6*omega is NOT such a symmetry (shift -omega is not a period): number of h_0 with a different phi-hat
  (3) for the three classes the construction gives the same c, the same case and sigma, the same kappa_0, and
      a differing by a multiple of c (r in {1, p7}); delta' may differ
  (4) relation between Psi_0 = C and Psi_0 = F for P = {} (their |T| agree in stage3_d_mellin.out):
      is phi_F(n) = conj phi_C(n)?  is phi-hat_F(h) = conj phi-hat_C(-h)?
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage2_h0_period.py > audits/003dual/stage2_h0_period.out 2>/dev/null
"""
from stage2_lib import *
print(__doc__)
LINT = 18; setup = Setup((LINT, 0)); P7 = (1, 3)
TAB = {nm: phi_table(nm, LINT) for nm in "ABCDEF"}
PH = {nm: phi_hat(TAB[nm], LINT) for nm in "ABCDEF"}
for nm in "ABCDEF":
    d6 = max(abs(PH[nm][(x, y)] - PH[nm][((x + 6) % 18, y)]) for x in range(18) for y in range(18))
    nw = sum(1 for x in range(18) for y in range(18) if abs(PH[nm][(x, y)] - PH[nm][(x, (y + 6) % 18)]) > 1e-12)
    supp = [h for h, v in PH[nm].items() if abs(v) > 1e-13]
    orb = set((h[0] % 6, h[1]) for h in supp)
    print("  Psi_0 = %s: (1) max |phi-hat(h_0+6) - phi-hat(h_0)| = %.1e ; support %d classes = %d orbits under h_0 -> h_0+6 ; (2) h_0 with phi-hat(h_0+6w) != phi-hat(h_0): %d"
          % (nm, d6, len(supp), len(orb), nw))
bad = 0; tot = 0; dl_diff = 0
for primes, ks in (([], []), ([P7], [ONE])):
    for x in range(6):
        for y in range(18):
            cf = [setup.build((x + 6 * k, y), primes, ks) for k in range(3)]
            tot += 1
            same = all(c["c"] == cf[0]["c"] and c["case"] == cf[0]["case"] and c["sigma"] == cf[0]["sigma"]
                       and kappa_formula(c)[1] == kappa_formula(cf[0])[1] and divides(cf[0]["c"], sub(c["a"], cf[0]["a"])) for c in cf)
            bad += not same
            dl_diff += any(c["dl"] != cf[0]["dl"] for c in cf)
print("  (3) orbits {h_0, h_0+6, h_0+12} (r in {1, p7}): %d ; with different (c, case, sigma, kappa_0) or a not congruent mod c: %d ; with different delta': %d" % (tot, bad, dl_diff))
e1 = float(np.max(np.abs(TAB["F"] - np.conj(TAB["C"]))))
e2 = max(abs(PH["F"][h] - PH["C"][((-h[0]) % 18, (-h[1]) % 18)].conjugate()) for h in PH["F"])
print("  (4) max |phi_F - conj phi_C| = %.1e ; max |phi-hat_F(h) - conj phi-hat_C(-h)| = %.1e" % (e1, e2))
