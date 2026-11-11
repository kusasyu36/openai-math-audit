# -*- coding: utf-8 -*-
"""
Stage 5 (c): the (case, sigma, j) combinations that stage3_d_mellin did not run end to end and that the tables reach.
Same comparison, same weight family (log-Gaussian, NOT compactly supported), same bounds and same criteria as layer E:
  E1  |T_direct - T_dual| <= Btot (Btot = (T2) + (T3) + (Q) + rounding allowance (estimate), stage4_bounds.py);
      informative iff Btot <= 1e-7 |T_direct|.
  Controls NE1-NE7 on informative runs where live: |T_direct - T_ctl| > max(100 Btot, 1e-3 |T_direct|).
Runs (fixed before running; (w, X) are those already used for the same N(c) in stage3_d_mellin):
  3|c, sigma=0:      A {p7: 2}, {p7: 3}, {p7: 5}                       (w=0.05, X=64000)
  v=1, sigma=+:      C {p7: 0}, {p7: 3}, {p7: 4}, {p7: 5}              (w=0.1,  X=16000)
  v=1, sigma=-:      F {p7: 2}, {p7: 3}                                (w=0.1,  X=16000)
Still not reachable (and not run): two primes with N(c) >= 364.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage5_c_more.py > audits/003dual/stage5_c_more.out 2>/dev/null
"""
from stage5_lib import *
t0 = time.time()
print(__doc__)
P7 = (1, 3)
RUNS = [("A", {P7: j}, 0.05, 64000.0) for j in (2, 3, 5)] + [("C", {P7: j}, 0.1, 16000.0) for j in (0, 3, 4, 5)] + [("F", {P7: j}, 0.1, 16000.0) for j in (2, 3)]
XMAX = {0.05: 3.0e4, 0.1: 1.5e4}
E = EndToEnd()
nfail = ninf = 0; worst_abs = worst_rel = 0.0; ctl = {k: [0, 0, 1e9] for k in ("NE1", "NE2", "NE3", "NE4", "NE5", "NE6", "NE7")}
for nm, jd, w, X in RUNS:
    wt = WLogGauss(w, XMAX[w]); r = E.run(nm, jd, wt, X); Td = r["Td"]; diff = abs(Td - r["Tu"])
    inf_ = r["Btot"] <= 1e-7 * abs(Td); ok = diff <= r["Btot"]
    print("  Psi_0=%s P=%s w=%.2f X=%g #groups=%d N(c) in %s cases %s" % (nm, {norm(p): j for p, j in jd.items()}, w, X, len(r["G"]),
          sorted(set(g["key"][1] for g in r["G"])), sorted(set((g["cfg"]["case"], g["cfg"]["sigma"], "active" if g["cfg"]["primes"] else "inactive") for g in r["G"]))))
    print("      T_dir = % .12e %+.12ei   T_dual = % .12e %+.12ei" % (Td.real, Td.imag, r["Tu"].real, r["Tu"].imag))
    print("      |T_dir| = %.4e  |diff| = %.2e (abs)  %.2e (rel)   Btot = %.2e = direct tail %.1e + dual tail %.1e + quadrature %.1e + rounding allowance (estimate) %.1e   %s%s"
          % (abs(Td), diff, diff / abs(Td), r["Btot"], r["b_dir"], r["b_dual"], r["qb"], r["b_round"], "PASS" if ok else "FAIL", "" if inf_ else "  (NOT INFORMATIVE)"))
    if not inf_: ninf += 1; continue
    nfail += not ok; worst_abs = max(worst_abs, diff); worst_rel = max(worst_rel, diff / abs(Td))
    thr = max(100 * r["Btot"], 1e-3 * abs(Td)); line = []
    for k in sorted(ctl):
        Tc = E.control(k, nm, jd, wt, X, r["G"])
        if Tc is None: line.append("%s dead" % k); continue
        d = abs(Td - Tc); ctl[k][0] += 1; ctl[k][1] += d > thr; ctl[k][2] = min(ctl[k][2], d / abs(Td)); line.append("%s %.1e" % (k, d / abs(Td)))
    print("      controls, |T_dir - T_ctl|/|T_dir|: " + "  ".join(line)); sys.stdout.flush()
print("\n  E1 %s: %d failures / %d informative runs (%d not informative) ; worst |diff| = %.2e (ABSOLUTE), %.2e (RELATIVE)" % ("PASS" if nfail == 0 and ninf < len(RUNS) else "FAIL", nfail, len(RUNS) - ninf, ninf, worst_abs, worst_rel))
allok = nfail == 0
for k in sorted(ctl):
    ok = ctl[k][0] > 0 and ctl[k][1] == ctl[k][0]; allok &= ok
    print("  %s %s: live in %d runs, rejected in %d, smallest |T_dir - T_ctl|/|T_dir| = %.2e" % (k, "PASS" if ok else "FAIL", ctl[k][0], ctl[k][1], ctl[k][2]))
print("\nSTAGE 5 (c): %s   [%.0fs]" % ("all criteria met" if allok else "SOME CRITERIA NOT MET", time.time() - t0))
