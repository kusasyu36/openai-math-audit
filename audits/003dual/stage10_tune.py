# -*- coding: utf-8 -*-
"""
003dual stage 10, TABLE-BLIND tuning (lead, 2026-11-13).  Solves theta + F_10 + F_19 jointly (stage10_lib.py) in two
settings and prints only solved quantities: rank, singular values, residual, spread between the two solves.
The DR table is NOT read here.
"""
import sys, os, math, time
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage10_lib as T
print(__doc__, flush=True)
# run 1 (v0 0.45/0.42, cusp grid 22/24, per_point 4) gave a stable theta but F_10 unstable (spread 1e+07 at |mu| 2-3):
# the cusp grid was below the Nyquist spacing (needs spacing < 1/(4 M_cu)).  Run 2 (this setting) fixes the sampling.
CONF = {"S1": dict(v0=0.55, M_th=4.0, M_cu=3.5, n_th=18, n_cu=46, pp_cu=2, shift=(0.0, 0.0)),
        "S2": dict(v0=0.52, M_th=4.3, M_cu=3.7, n_th=19, n_cu=48, pp_cu=2, shift=(0.37, 0.61))}
sol = {}
for name, c in CONF.items():
    t0 = time.time()
    S = T.System(c["M_th"], c["M_cu"], c["v0"])
    A, b, info = S.build(c["n_th"], c["n_cu"], c["shift"], pp_cu=c["pp_cu"])
    print(f"  sampling: theta spacing {1/c['n_th']:.3f}, cusp spacing {3/c['n_cu']:.3f}, Nyquist 1/(4M) theta {1/(4*c['M_th']):.3f} cusp {1/(4*c['M_cu']):.3f}", flush=True)
    x, _, rank, sv = np.linalg.lstsq(A, b, rcond=1e-14)
    res = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
    print(f"{name} {c}: unknowns {S.n} {[len(S.nums[f]) for f in ('th','g10','g19')]}, rows {len(b)}, chi {info['chi']}, "
          f"skipped {info['skipped']}, rank {rank}, sv max {sv[0]:.3g}, smallest 6 {np.array2string(sv[-6:], precision=3)}, "
          f"rel residual {res:.2e}, {time.time()-t0:.0f}s", flush=True)
    sol[name] = (S, S.split(x)); del A
S1, X1 = sol["S1"]; S2, X2 = sol["S2"]
for fam, den in (("th", 3), ("g10", 6), ("g19", 6)):
    for lo, hi in ((0, 1), (1, 2), (2, 3)):
        ms = [m for m in X1[fam] if m in X2[fam] and lo < abs((m[0] + m[1] * T.W3) / T.kb.LAMC ** den) <= hi + 1e-12]
        if not ms: continue
        d = max(abs(X1[fam][m] - X2[fam][m]) for m in ms); s = max(abs(X1[fam][m]) for m in ms)
        print(f"spread {fam:4s} |mu| in ({lo},{hi}]: n={len(ms):4d} max|S1-S2| {d:.2e}  max|coef| {s:.3g}")
