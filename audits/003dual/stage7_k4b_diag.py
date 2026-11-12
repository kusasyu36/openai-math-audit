# -*- coding: utf-8 -*-
"""
003dual stage 7, POST-HOC diagnostic of the failed control K4b (lead, Claude Opus 5.5, 2026-11-12).
K4b expected the system with the conjugate multiplier (chi -> conj chi, same reps and translations) to be
inconsistent; it was consistent (residual 1e-13).  Hypothesis written before running this script:
the solution is theta reflected in z, f(z, v) = theta(conj z, v), whose coefficient at mu is tau(conj mu),
because conj((c/a)_3) = (conj c / conj a)_3 and the reflection maps the group and the rep set {I, P, P2} to itself.
This is a diagnosis after the fact; it does not turn K4b into a pass.
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage7_lib as L
sys.argv = [sys.argv[0]]
print(__doc__)
# reuse build/table without rerunning the main script body: copy the two helpers
src = open(os.path.join(here, "stage7_hejhal.py")).read()
ns = {"__file__": os.path.join(here, "stage7_hejhal.py")}
exec(src.split("# ---------------------------------------------------------------- solves")[0].replace("print(__doc__, flush=True)", ""), ns)
T9 = ns["table"](9.0)
nums, mu, A, b, sc, info = ns["build"]("S1", conj_chi=True)
tau, sv, rank, resid = L.solve(A, b, sc)
kb = L.kb
def cj(m): return kb.conj(m)      # conj(lambda^3 mu) = conj(lambda)^3 conj(mu) = -lambda^3 conj(mu)  (conj lambda = -lambda)
best = {}
for lab, f in (("tau(mu)", lambda m: T9.get(m, 0)),
               ("tau(conj mu)", lambda m: T9.get(kb.neg(cj(m)), 0)),
               ("conj tau(conj mu)", lambda m: np.conj(T9.get(kb.neg(cj(m)), 0))),
               ("conj tau(-mu)", lambda m: np.conj(T9.get(kb.neg(m), 0)))):
    d = max(abs(t - f(m)) for m, t, a in zip(nums, tau, np.abs(mu)) if a <= 2.0)
    print(f"conj-multiplier solution vs {lab:18s}: max |diff| over |mu| <= 2 = {d:.2e}")
