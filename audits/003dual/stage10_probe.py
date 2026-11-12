# 003dual stage 10 PROBE (lead, 2026-11-13), run BEFORE the tuning and before the criteria of stage10_frozen.py were written.
# It plugged the DR table into a first joint system (M_th 6, M_cu 4, v0 0.30, small grids) as a code check.
# Output (copied from the session): chars (1, 1) (2, 2); unknowns th 1166, g10 1566, g19 1566; table entries outside
# index set 0/0/0; rows 1681, skipped 0; table residual 2.97e-07 (truncation level e^{-4 pi 4 0.3}).  Disclosed in REPORT 16.
import sys, os, time; sys.path.insert(0, os.path.expanduser('~/chiken/openaikoeru/audits/003dual'))
import numpy as np, stage10_lib as T
print("chars", T.char_of("g10"), T.char_of("g19"))
t0 = time.time()
S = T.System(6.0, 4.0, 0.30)
print({f: len(S.nums[f]) for f in S.nums}, "n", S.n)
xt, extra = T.table_vector(S); print("table entries outside index set (inside range):", extra)
A, b, info = S.build(10, 8)
print("rows", len(b), info, f"{time.time()-t0:.0f}s")
r = np.linalg.norm(A @ (xt * S.scale) - b) / np.linalg.norm(b)
print("table residual", r)
