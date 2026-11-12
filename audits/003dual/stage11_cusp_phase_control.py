# -*- coding: utf-8 -*-
"""
003dual stage 11 (lead, 2026-11-13): the cusp-side control proposed by Codex (w11 (6)), replacing the failed E4a.
Modified table: in d_10 ONLY, the additive phase e(mu) (kubota.build_series, ecek) is replaced by e(-mu);
indices, theta and d_19 unchanged.  Nothing is re-solved for the modified table.
CRITERIA (written before running):
  P1 the modified d_10 differs from the S1-solved d_10 (stage 10 frozen S1) by > 1 somewhere in |mu| <= 2;
  P2 in the S2 system of stage 10 (rows never used for S1): unmodified table relative residual <= 1e-8,
     modified table >= 1e-3.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy python -B -u audits/003dual/stage11_cusp_phase_control.py > audits/003dual/stage11_cusp_phase_control.out
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage10_lib as T
kb = T.kb
print(__doc__, flush=True)
CONF = {"S1": dict(v0=0.55, M_th=4.0, M_cu=3.5, n_th=18, n_cu=46, pp_cu=2, shift=(0.0, 0.0)),
        "S2": dict(v0=0.52, M_th=4.3, M_cu=3.7, n_th=19, n_cu=48, pp_cu=2, shift=(0.37, 0.61))}
def modify(S, xt):
    y = xt.copy(); o = S.off["g10"]
    for i, m6 in enumerate(S.nums["g10"]):
        m4 = kb.exact_div(m6, T.LAM2)
        if m4 is None: assert xt[o + i] == 0; continue
        y[o + i] *= np.exp(-4j * math.pi * (2 * m4[0] - m4[1]) / 9.0)    # e(mu) -> e(-mu)
    return y
c = CONF["S1"]; S = T.System(c["M_th"], c["M_cu"], c["v0"])
A, b, _ = S.build(c["n_th"], c["n_cu"], c["shift"], pp_cu=c["pp_cu"])
x = np.linalg.lstsq(A, b, rcond=1e-14)[0]; del A
xt, _ = T.table_vector(S); xm = modify(S, xt)
X = S.split(x)["g10"]; M = S.split(xm * S.scale)["g10"]
amu = lambda m: abs((m[0] + m[1] * T.W3) / kb.LAMC ** 6)
p1 = max(abs(X[m] - M[m]) for m in X if amu(m) <= 2)
print(f"P1 modified d_10 vs S1-solved, |mu| <= 2: max diff {p1:.3g}")
c = CONF["S2"]; S = T.System(c["M_th"], c["M_cu"], c["v0"])
A, b, _ = S.build(c["n_th"], c["n_cu"], c["shift"], pp_cu=c["pp_cu"])
xt, _ = T.table_vector(S); xm = modify(S, xt)
r0 = np.linalg.norm(A @ (xt * S.scale) - b) / np.linalg.norm(b); r1 = np.linalg.norm(A @ (xm * S.scale) - b) / np.linalg.norm(b)
print(f"P2 S2 system: unmodified table {r0:.2e}, modified table {r1:.2e}")
V = {"P1": p1 > 1, "P2": r0 <= 1e-8 and r1 >= 1e-3}
for k, v in V.items(): print(f"  {k}: {'PASS' if v else 'FAIL'}")
