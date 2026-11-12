# -*- coding: utf-8 -*-
"""
003dual stage 8b (lead, Claude Opus 5.5, 2026-11-12, POST HOC, after Codex w10 (5)): K3 compared each coefficient;
an odd perturbation changes individual coefficients, so "the odd direction does not affect K3" (REPORT 13.3) was a
non sequitur.  Here the stage-7 solutions S1, S2 (exact frozen settings, re-solved) are split into
even part (tau(mu)+tau(-mu))/2 and odd part (tau(mu)-tau(-mu))/2 over |mu| <= 3, and reported:
  max |even part - table|, max |odd part| (the table is even, so the odd part should be ~0).
Not pre-registered; reported against the K3 tolerances only for orientation.
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage7_lib as L
kb = L.kb
print(__doc__)
src = open(os.path.join(here, "stage7_hejhal.py")).read()
ns = {"__file__": os.path.join(here, "stage7_hejhal.py")}
exec(src.split("# ---------------------------------------------------------------- solves")[0].replace("print(__doc__, flush=True)", ""), ns)
T = ns["table"](9.0)
for name in ("S1", "S2"):
    nums, mu, A, b, sc, info = ns["build"](name)
    tau, sv, rank, resid = L.solve(A, b, sc)
    D = dict(zip(nums, tau)); a = dict(zip(nums, np.abs(mu)))
    ev, od, ratio_e = 0, 0, 0
    for m in nums:
        if a[m] > 3 * (1 + 1e-12): continue
        e = (D[m] + D[kb.neg(m)]) / 2; o = (D[m] - D[kb.neg(m)]) / 2
        de = abs(e - T.get(m, 0)); ev = max(ev, de); od = max(od, abs(o))
        ratio_e = max(ratio_e, de / ns["tol_of"](a[m]))
    print(f"{name}: rows {len(b)} (same as stage 7: {'yes' if (name=='S1' and len(b)==8690) or (name=='S2' and len(b)==10648) else 'NO'}); "
          f"|mu|<=3: max |even - table| {ev:.2e} (max/tol {ratio_e:.3f}), max |odd part| {od:.2e}")
