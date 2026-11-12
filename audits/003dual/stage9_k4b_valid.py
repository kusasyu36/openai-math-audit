# -*- coding: utf-8 -*-
"""
003dual stage 9 (lead, Claude Opus 5.5, 2026-11-13): a VALID multiplier control replacing K4b (stage 7).
K4b solved the conjugate-multiplier system and found it consistent: its solution is tau(conj mu)
(= theta reflected in z), see stage7_k4b_diag.out.  So "is the conj-chi system solvable" does not test the multiplier.
Valid version (Codex w10 (4)): FIX the coefficients and test the transformation rows themselves.
Same rows as stage 7 S1 (v0 0.30, M 8, grid 30, per_point 4, N(c) <= 40); only the factor f = omega^k changes.

CRITERIA (written before the first run of this script):
  V1 (code check) table tau, chi rows:                 relative residual <= 1e-10    (= K0)
  V2 table tau, conj-chi rows (f = omega^{-k}):         relative residual >= 1e-2
  V3 S1-solved tau (table never used), conj-chi rows:   relative residual >= 1e-2
  V4 table tau, trivial multiplier (f = 1 on all rows): relative residual >= 1e-2
  V5 the failure in V2 sits in the chi != 1 rows: on chi = 1 rows alone the residual is <= 1e-10,
     and on chi = omega, omega^2 rows it is >= 1e-2
  reported, no criterion: tau(conj mu) fixed in the conj-chi rows (should be small, consistent with the diagnosis).
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy python -B -u audits/003dual/stage9_k4b_valid.py > audits/003dual/stage9_k4b_valid.out
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage7_lib as L
kb = L.kb
print(__doc__, flush=True)
src = open(os.path.join(here, "stage7_hejhal.py")).read()
ns = {"__file__": os.path.join(here, "stage7_hejhal.py")}
exec(src.split("# ---------------------------------------------------------------- solves")[0].replace("print(__doc__, flush=True)", ""), ns)
T9 = ns["table"](9.0)
c = ns["CONF"]["S1"]
nums = L.dual_indices(c["M"]); mu = L.mu_of(nums)
sc = c["v0"] * L.kv(1.0 / 3, 4 * math.pi * np.abs(mu) * c["v0"])
pts = L.sample_points(c["n_side"], c["v0"], c["shift"])
P0, P1, KS = [], [], []          # phi at sample point, phi at pulled-back point, chi exponent, height
for z in pts:
    cand, _ = L.pullbacks(z, c["v0"], 40)
    cand = [x for x in cand if x[0] >= 1.05 * c["v0"]]
    if not cand: continue
    chosen = [cand[0]]
    for x in cand[1:]:
        if len(chosen) >= 4: break
        if all(x[1] != y[1] or x[2] != y[2] for y in chosen): chosen.append(x)
    for v2, k, rep, z2, h in chosen:
        P0.append(L.phi_row(mu, z, c["v0"]) / sc); P1.append(L.phi_row(mu, z2, v2) / sc); KS.append((k, v2))
P0 = np.array(P0); P1 = np.array(P1); kk = np.array([k for k, _ in KS]); v2s = np.array([v for _, v in KS])
print(f"rows {len(kk)}, chi exponents {[int((kk == j).sum()) for j in range(3)]}", flush=True)

def system(expo):
    f = L.W3 ** expo
    A = P0 - f[:, None] * P1
    b = L.SIGMA * (f * v2s ** (2.0 / 3) - c["v0"] ** (2.0 / 3))
    return A, b

def rres(A, b, x, mask=None):
    if mask is not None: A, b = A[mask], b[mask]
    return np.linalg.norm(A @ x - b) / np.linalg.norm(b)

x_tab = np.array([T9.get(m, 0) for m in nums]) * sc
x_ref = np.array([T9.get(kb.neg(kb.conj(m)), 0) for m in nums]) * sc     # tau(conj mu)
A_chi, b_chi = system(kk)
A_cj, b_cj = system((-kk) % 3)
A_one, b_one = system(np.zeros_like(kk))
x_s1, _, rank, _ = np.linalg.lstsq(A_chi, b_chi, rcond=1e-14)
print(f"S1 solve: rank {rank}/{len(nums)}, residual {rres(A_chi, b_chi, x_s1):.2e}")
V = {}
r = rres(A_chi, b_chi, x_tab); print(f"V1 table, chi rows: {r:.2e}"); V["V1"] = r <= 1e-10
r = rres(A_cj, b_cj, x_tab); print(f"V2 table, conj-chi rows: {r:.2e}"); V["V2"] = r >= 1e-2
r = rres(A_cj, b_cj, x_s1); print(f"V3 S1-solved, conj-chi rows: {r:.2e}"); V["V3"] = r >= 1e-2
r = rres(A_one, b_one, x_tab); print(f"V4 table, trivial multiplier: {r:.2e}"); V["V4"] = r >= 1e-2
r0 = rres(A_cj, b_cj, x_tab, kk == 0); r12 = rres(A_cj, b_cj, x_tab, kk != 0)
print(f"V5 table, conj-chi rows split: chi=1 rows {r0:.2e}, chi!=1 rows {r12:.2e}"); V["V5"] = r0 <= 1e-10 and r12 >= 1e-2
print(f"reported: tau(conj mu) in conj-chi rows {rres(A_cj, b_cj, x_ref):.2e}, in chi rows {rres(A_chi, b_chi, x_ref):.2e}")
print(f"reported: S1-solved in chi rows split: chi=1 {rres(A_chi, b_chi, x_s1, kk == 0):.2e}, chi!=1 {rres(A_chi, b_chi, x_s1, kk != 0):.2e}")
print("\nVERDICT:")
for k_, v in V.items(): print(f"  {k_}: {'PASS' if v else 'FAIL'}")
print("ALL PASS" if all(V.values()) else "NOT ALL PASS")
