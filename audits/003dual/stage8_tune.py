# -*- coding: utf-8 -*-
"""
003dual stage 8 tuning (lead, Claude Opus 5.5, 2026-11-12, TABLE-BLIND): does the isolated smallest singular value of
the stage-7 systems disappear when points near the cusps 10/19 are brought back?  Variant: for every point, also
allow pullbacks (rep I, P, P2) whose height is only >= g_low * v0 (g_low < 1), so that points whose best pullback was
g10/g19 still give equations (both sides expanded at infinity; truncation at the lower height is checked by M).
Prints: rows, skipped points, rank, the smallest 4 singular values, residual, and two-solve agreement.
"""
import sys, os, math, time
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage7_lib as L
print(__doc__)

def build(v0, M, n_side, shift, g_low, per_point=4):
    nums = L.dual_indices(M); mu = L.mu_of(nums)
    scale = v0 * L.kv(1.0 / 3, 4 * math.pi * np.abs(mu) * v0)
    rows, rhs = [], []; chi = [0, 0, 0]; skipped = 0; low = 0
    for z in L.sample_points(n_side, v0, shift):
        cand, _ = L.pullbacks(z, v0, 40, gain=g_low)
        hi = [x for x in cand if x[0] >= 1.05 * v0]
        lo = [x for x in cand if g_low * v0 <= x[0] < 1.05 * v0 and abs(x[0] - v0) > 1e-3 * v0]
        chosen = []
        for x in hi + lo:
            if len(chosen) >= per_point: break
            if all(x[1] != y[1] or x[2] != y[2] or abs(x[0] - y[0]) > 1e-6 for y in chosen): chosen.append(x)
        if not chosen: skipped += 1; continue
        for v2, k, rep, z2, h in chosen:
            if v2 < 1.05 * v0: low += 1
            f = L.W3 ** k
            rows.append((L.phi_row(mu, z, v0) - f * L.phi_row(mu, z2, v2)) / scale)
            rhs.append(L.SIGMA * (f * v2 ** (2.0 / 3) - v0 ** (2.0 / 3))); chi[k] += 1
    return nums, mu, np.array(rows), np.array(rhs), scale, chi, skipped, low

if __name__ == "__main__":
    g_low = float(sys.argv[1]) if len(sys.argv) > 1 else 0.85
    sols = []
    for (v0, M, ns, sh) in ((0.30, 9.0, 30, (0, 0)), (0.26, 10.0, 32, (0.37, 0.61))):
        t = time.time()
        nums, mu, A, b, sc, chi, sk, low = build(v0, M, ns, sh, g_low)
        tau, sv, rank, resid = L.solve(A, b, sc)
        print(f"g_low={g_low} v0={v0} M={M} unk={len(nums)} rows={len(b)} (low-height rows {low}) chi={chi} skipped={sk} "
              f"rank={rank} sv max {sv[0]:.3g} smallest4 {np.array2string(sv[-4:], precision=3)} resid={resid:.2e} {time.time()-t:.0f}s", flush=True)
        sols.append(dict(zip(nums, tau)))
    s1, s2 = sols
    for Mc in (1.0, 2.0, 3.0):
        d = [abs(s1[m] - s2[m]) for m in s1 if m in s2 and abs(L.mu_of([m])[0]) <= Mc * (1 + 1e-12)]
        print(f"  |mu|<={Mc}: n={len(d)} max|S1-S2|={max(d):.2e}")
