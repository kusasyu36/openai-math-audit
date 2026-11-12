# -*- coding: utf-8 -*-
"""
003dual stage 10 POST-HOC (lead, 2026-11-13), after stage10_frozen.out.  Not part of the pre-registered criteria.
(1) Odd part (tau(mu) - tau(-mu))/2 of the theta coefficients from the JOINT solve (S1, S2 of stage 10), |mu| <= 3,
    to compare with the odd weak direction of stages 7-8 (stage8_evenodd_k3.out: odd part 8.3e-6 (S1), 3.9e-4 (S2)).
(2) Stage 9 on UNUSED points (Codex w10 (4)): tau solved in stage 7 S1 (v0 0.30, M 8) is fixed and put into the rows of
    stage 7 S2 (v0 0.26, other grid and shift; none of these rows were used in the S1 solve).
    Criteria written before running (2): H1 chi rows: relative residual <= 1e-3; H2 conj-chi rows: >= 1e-2;
    H3 trivial multiplier: >= 1e-2.
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage10_lib as T
import stage7_lib as L
kb = T.kb
print(__doc__, flush=True)
CONF = {"S1": dict(v0=0.55, M_th=4.0, M_cu=3.5, n_th=18, n_cu=46, pp_cu=2, shift=(0.0, 0.0)),
        "S2": dict(v0=0.52, M_th=4.3, M_cu=3.7, n_th=19, n_cu=48, pp_cu=2, shift=(0.37, 0.61))}
for name, c in CONF.items():
    S = T.System(c["M_th"], c["M_cu"], c["v0"])
    A, b, info = S.build(c["n_th"], c["n_cu"], c["shift"], pp_cu=c["pp_cu"])
    x = np.linalg.lstsq(A, b, rcond=1e-14)[0]; del A
    th = S.split(x)["th"]
    odd = max(abs(th[m] - th.get(kb.neg(m), 0)) / 2 for m in th if abs(L.mu_of([m])[0]) <= 3 + 1e-12)
    print(f"(1) joint {name}: max odd part of theta over |mu| <= 3: {odd:.2e}", flush=True)

# (2)
src = open(os.path.join(here, "stage7_hejhal.py")).read()
ns = {"__file__": os.path.join(here, "stage7_hejhal.py")}
exec(src.split("# ---------------------------------------------------------------- solves")[0].replace("print(__doc__, flush=True)", ""), ns)
n1, mu1, A1, b1, sc1, _ = ns["build"]("S1")
t1 = np.linalg.lstsq(A1, b1, rcond=1e-14)[0] / sc1; del A1
tau1 = dict(zip(n1, t1))
c = ns["CONF"]["S2"]
nums = L.dual_indices(c["M"]); mu = L.mu_of(nums)
sc = c["v0"] * L.kv(1.0 / 3, 4 * math.pi * np.abs(mu) * c["v0"])
P0, P1, KK, V2 = [], [], [], []
for z in L.sample_points(c["n_side"], c["v0"], c["shift"]):
    cand, _ = L.pullbacks(z, c["v0"], 40)
    cand = [q for q in cand if q[0] >= 1.05 * c["v0"]]
    if not cand: continue
    ch = [cand[0]]
    for q in cand[1:]:
        if len(ch) >= 4: break
        if all(q[1] != y[1] or q[2] != y[2] for y in ch): ch.append(q)
    for v2, k, rep, z2, h in ch:
        P0.append(L.phi_row(mu, z, c["v0"])); P1.append(L.phi_row(mu, z2, v2)); KK.append(k); V2.append(v2)
P0 = np.array(P0); P1 = np.array(P1); KK = np.array(KK); V2 = np.array(V2)
xv = np.array([tau1.get(m, 0) for m in nums])
def res(expo):
    f = L.W3 ** expo
    A = P0 - f[:, None] * P1; b = L.SIGMA * (f * V2 ** (2 / 3) - c["v0"] ** (2 / 3))
    return np.linalg.norm(A @ xv - b) / np.linalg.norm(b)
H = {"H1": res(KK) <= 1e-3, "H2": res((-KK) % 3) >= 1e-2, "H3": res(np.zeros_like(KK)) >= 1e-2}
print(f"(2) S1-solved tau in S2 rows ({len(KK)} rows): chi {res(KK):.2e}, conj chi {res((-KK) % 3):.2e}, trivial {res(np.zeros_like(KK)):.2e}")
for k, v in H.items(): print(f"  {k}: {'PASS' if v else 'FAIL'}")
