# -*- coding: utf-8 -*-
"""
003dual stage 7, POST-HOC diagnostic (lead, Claude Opus 5.5, 2026-11-12): the smallest singular value of the S1/S2
systems is isolated (S1: 7.6e-6 vs next 2.9e-3; S2: 8.0e-5 vs next 2.0e-2).  Where does the corresponding right
singular vector live (which |mu|), and how much does it move the compared coefficients (|mu| <= 3)?
Table-blind.  Not part of the pre-registered criteria.
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage7_lib as L
print(__doc__)
src = open(os.path.join(here, "stage7_hejhal.py")).read()
ns = {"__file__": os.path.join(here, "stage7_hejhal.py")}
exec(src.split("# ---------------------------------------------------------------- solves")[0].replace("print(__doc__, flush=True)", ""), ns)
for name in ("S1", "S2"):
    nums, mu, A, b, sc, info = ns["build"](name)
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    v = Vh[-1].conj()                       # scaled coordinates, unit norm
    a = np.abs(mu); w = np.abs(v) ** 2
    o = np.argsort(-w)[:8]
    print(f"\n{name}: sv min {s[-1]:.3e}, next {s[-2]:.3e}")
    print("  top weights (|mu|, weight):", ", ".join(f"({a[i]:.3f}, {w[i]:.3f})" for i in o))
    for lo, hi in ((0, 3), (3, 6), (6, 8), (8, 99)):
        print(f"  weight in {lo} < |mu| <= {hi}: {w[(a > lo) & (a <= hi)].sum():.3e}")
    # effect on the solution: component of the solution along this vector, in unscaled tau, over |mu| <= 3
    x = np.linalg.lstsq(A, b, rcond=1e-14)[0]
    coef = np.vdot(Vh[-1].conj(), x) if False else (Vh[-1] @ x)
    unscaled = np.abs(coef * v / sc)
    print(f"  solution component along it: {abs(coef):.3e}; its max contribution to |tau| over |mu| <= 3: {unscaled[a <= 3].max():.3e}")
    # perturbation sensitivity: residual-level noise 1e-13*|b| divided by s_min, mapped to |mu| <= 3
    eps = 1e-13 * np.linalg.norm(b) / s[-1]
    print(f"  noise amplification estimate along it (1e-13 |b| / s_min) -> max |tau| change over |mu| <= 3: {(eps * np.abs(v / sc))[a <= 3].max():.3e}")
