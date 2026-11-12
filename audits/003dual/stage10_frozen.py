# -*- coding: utf-8 -*-
"""
003dual stage 10 (lead, Claude Opus 5.5, 2026-11-13): theta AND the cusp functions F_10, F_19 solved jointly from
automorphy (stage10_lib.py), then compared ONCE with the DR table (kubota.build_series: (5.7), (5.13), (5.14), App. A).
This is the first check of d_10, d_19 (= d_+-) that does not use DR's cusp formulas, so it also tests the RELATIVE
normalisation of direct and cusp coefficients (Codex w9 (8), w10 (10)).

Tuning (stage10_tune.py/.out, table-blind; run 1 kept in stage10_tune_run1.out): run 1 had F_10 unstable because the
cusp grid was below the Nyquist spacing; run 2 fixed the sampling and gave sv_min/sv_max 0.32-0.34 and spreads
theta 5e-9 / 1.1e-5 / 8.2e-3, cusps 5.5e-9 / 5.9e-6 / 2.1e-2  (|mu| in (0,1] / (1,2] / (2,3]).
FROZEN = run 2:  S1 v0 0.55, M_th 4.0, M_cu 3.5, grid 18 / 46, per_point 4 (theta) / 2 (cusps), shift (0, 0)
                 S2 v0 0.52, M_th 4.3, M_cu 3.7, grid 19 / 48, same, shift (0.37, 0.61)
CRITERIA (written before the first comparison with the table):
  E0 (code check, not evidence): the DR table (theta, d_10, d_19) plugged into the S1 system: relative residual <= 1e-8
     (the truncation e^{-4 pi M_cu v0} ~ 3e-11 times the coefficient size).
  E1 (code check): every table entry of d_10, d_19 inside the range lies in the group-derived index set (stage10_lib).
  E2 (observed in tuning, so a property, not a test): both solves have full rank and sv_min / sv_max >= 1e-3.
  E3 comparison, all unknowns with |mu| <= 3 (table extended by 0 off its support), both solves, tolerances
     ~10 x the tuning spread:  theta  1e-7 / 1e-4 / 8e-2 ;  d_10, d_19  1e-7 / 6e-5 / 2e-1  for (0,1] / (1,2] / (2,3].
  E4 controls (power of the comparison, expected to pass trivially):
     (a) the swap9 table differs from S1 by > 1 somewhere in |mu| <= 2 for theta AND for d_10;
     (b) d_10, d_19 of the table multiplied by sqrt(3) (a wrong relative normalisation) differ from S1 by > 1.
  Reported, no criterion: the theta-only system of stage 7 type at the S1 setting (theta unknowns, theta sample points,
  rows only through reps I, P, P2): rank, skipped points, smallest singular values -- to see whether the isolated
  small singular value of stage 7 appears there and not in the joint system.
Tolerances are error estimates from the two-solve spread, not bounds.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy python -B -u audits/003dual/stage10_frozen.py > audits/003dual/stage10_frozen.out
"""
import sys, os, math, time
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage10_lib as T
kb = T.kb
print(__doc__, flush=True)
CONF = {"S1": dict(v0=0.55, M_th=4.0, M_cu=3.5, n_th=18, n_cu=46, pp_cu=2, shift=(0.0, 0.0)),
        "S2": dict(v0=0.52, M_th=4.3, M_cu=3.7, n_th=19, n_cu=48, pp_cu=2, shift=(0.37, 0.61))}
TOL = {"th": [(1, 1e-7), (2, 1e-4), (3, 8e-2)], "g10": [(1, 1e-7), (2, 6e-5), (3, 2e-1)], "g19": [(1, 1e-7), (2, 6e-5), (3, 2e-1)]}
DEN = {"th": 3, "g10": 6, "g19": 6}
V = {}; SOL = {}
def amu(fam, m): return abs((m[0] + m[1] * T.W3) / kb.LAMC ** DEN[fam])
def tol(fam, a):
    for e, t in TOL[fam]:
        if a <= e * (1 + 1e-12): return t

for name, c in CONF.items():
    t0 = time.time()
    S = T.System(c["M_th"], c["M_cu"], c["v0"])
    A, b, info = S.build(c["n_th"], c["n_cu"], c["shift"], pp_cu=c["pp_cu"])
    x, _, rank, sv = np.linalg.lstsq(A, b, rcond=1e-14)
    res = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
    print(f"\n== {name}: unknowns {S.n}, rows {len(b)}, chi {info['chi']}, rank {rank}, sv max {sv[0]:.3g} min {sv[-1]:.3g} "
          f"(ratio {sv[-1]/sv[0]:.2e}), rel residual {res:.2e}, {time.time()-t0:.0f}s", flush=True)
    print("   rows by (sample function, rep reached):", dict(sorted(info["rep"].items())))
    V[f"E2 {name}"] = rank == S.n and sv[-1] / sv[0] >= 1e-3
    xt, extra = T.table_vector(S)
    if name == "S1":
        r0 = np.linalg.norm(A @ (xt * S.scale) - b) / np.linalg.norm(b)
        print(f"E0 table in S1 system: {r0:.2e};  E1 table entries outside the index set: {extra}")
        V["E0"] = r0 <= 1e-8; V["E1"] = all(v == 0 for v in extra.values())
    del A
    X = S.split(x); TB = S.split(xt * S.scale); SOL[name] = X
    for fam in ("th", "g10", "g19"):
        worst = 0; by = {}
        for m, val in X[fam].items():
            a = amu(fam, m)
            if a > 3 * (1 + 1e-12): continue
            tv = TB[fam][m]; d = abs(val - tv); r = d / tol(fam, a); worst = max(worst, r)
            key = ("on" if tv != 0 else "off") + f" ({int(math.ceil(a - 1e-12)) - 1},{int(math.ceil(a - 1e-12))}]"
            e = by.setdefault(key, [0, 0.0, 0.0, 0.0]); e[0] += 1; e[1] = max(e[1], d); e[2] = max(e[2], r); e[3] = max(e[3], abs(tv))
        print(f"E3 {name} {fam}: worst diff/tol {worst:.3f}")
        for k in sorted(by):
            n, d, r, s = by[k]; print(f"     {k:14s} n={n:4d} max|solved-table| {d:.2e}  max diff/tol {r:.3f}  max|table| {s:.3g}")
        V[f"E3 {name} {fam}"] = worst <= 1
    if name == "S1":
        td = kb.ThetaData(200000, os.path.join(here, "cache"), verbose=False)
        thS, F10S, F19S = kb.build_series(td, 9.0, swap9=True)
        dsw = {}
        for fam, ser, conv in (("th", thS, lambda m: m), ("g10", F10S, lambda m: kb.mul(m, T.LAM2))):
            d = {}
            for m, cf in zip(map(tuple, ser.num), ser.coef): mm = conv(m); d[mm] = d.get(mm, 0) + cf
            dsw[fam] = max(abs(X[fam][m] - d.get(m, 0)) for m in X[fam] if amu(fam, m) <= 2)
        db = max(max(abs(X[f][m] - math.sqrt(3) * TB[f][m]) for m in X[f] if amu(f, m) <= 2) for f in ("g10", "g19"))
        print(f"E4a swap9 vs S1 (|mu|<=2): theta {dsw['th']:.3g}, d_10 {dsw['g10']:.3g};  E4b sqrt(3) x cusp table vs S1: {db:.3g}")
        V["E4a"] = dsw["th"] > 1 and dsw["g10"] > 1; V["E4b"] = db > 1

for fam in ("th", "g10", "g19"):
    ms = [m for m in SOL["S1"][fam] if m in SOL["S2"][fam] and amu(fam, m) <= 3]
    print(f"S1 vs S2 {fam} |mu|<=3: max diff {max(abs(SOL['S1'][fam][m]-SOL['S2'][fam][m]) for m in ms):.2e} ({len(ms)})")

# reported: theta-only, stage 7 type, at the S1 setting
c = CONF["S1"]; S = T.System(c["M_th"], c["M_cu"], c["v0"]); nth = len(S.nums["th"])
rows, rhs, skipped = [], [], 0
for z in T.sample_points("I", c["n_th"], c["shift"]):
    cands = []
    for g in T.candidates(z, c["v0"], 1.05 * c["v0"]):
        z2, v2 = kb.mobius(kb.mat_cx(g), z, c["v0"]); k, rep, t = T.dec(kb.mat_inv(g))
        if rep in ("I", "P", "P2"): cands.append((v2, k % 3, rep, T.S7.ROT[rep] * (z2 + kb.cx(t))))
    if not cands: skipped += 1; continue
    cands.sort(key=lambda q: -q[0]); ch = [cands[0]]
    for q in cands[1:]:
        if len(ch) >= 4: break
        if all(q[1] != y[1] or q[2] != y[2] for y in ch): ch.append(q)
    for v2, k, rep, z2 in ch:
        f = T.W3 ** k
        rows.append((S.block("th", z, c["v0"]) - f * S.block("th", z2, v2))[:nth] / S.scale[:nth])
        rhs.append(T.SIGMA * (f * v2 ** (2 / 3) - c["v0"] ** (2 / 3)))
A = np.array(rows); b = np.array(rhs)
x, _, rank, sv = np.linalg.lstsq(A, b, rcond=1e-14)
print(f"\nreported: theta-only (stage 7 type) at S1: rows {len(b)}, skipped points {skipped}, rank {rank}/{nth}, "
      f"smallest 5 sv {np.array2string(sv[-5:], precision=3)}, sv max {sv[0]:.3g}, rel residual {np.linalg.norm(A@x-b)/np.linalg.norm(b):.2e}")
print("\nVERDICT:")
for k, v in V.items(): print(f"  {k}: {'PASS' if v else 'FAIL'}")
print("ALL PASS" if all(V.values()) else "NOT ALL PASS")
