# -*- coding: utf-8 -*-
"""
003dual stage 7 (lead, Claude Opus 5.5, 2026-11-12): the theta coefficients ABOVE 3, by a route that does not use
DR (5.7)-(5.8).

Why: stage 6 (REPORT section 11) tied the coefficients prime to 3 to the paper's own formula, but the coefficients at
mu = +-lambda^{3n-3} c d^3 (n >= 1) and +-omega^k lambda^{3n-4} c d^3 (n >= 1) -- which the dual side of the
reflection uses through t_0 = tau (paper2.tex:3105-3115) -- were tied only indirectly (automorphy checks with the
DR table plugged in: 029reflect G1/G3/G5, layer B).
Here the coefficients are SOLVED from automorphy + constant term (Hejhal-type least squares, stage7_lib.py) and only
then compared with the DR table (kubota.build_series).  See stage7_lib.py for exactly which inputs are used.

Parameters were tuned beforehand on table-blind diagnostics only (rank, residual, agreement of two solves); the
tuning run (two solves, per_point = 4) gave max |S1 - S2| = 9.7e-7 (|mu| <= 1), 3.2e-5 (<= 2), 4.1e-4 (<= 3).
FROZEN:  S1: v0 = 0.30, M = 8.0, grid 30 x 90, shift (0, 0),      per_point 4, N(c) <= 40
         S2: v0 = 0.26, M = 9.0, grid 32 x 96, shift (0.37, 0.61), per_point 4, N(c) <= 40
CRITERIA (fixed before the first comparison with the table):
  K0 (code check, not evidence): the DR table plugged into the S1 system has relative residual <= 1e-10.
  K1 multiplier coverage: in each solve, rows with chi = omega and with chi = omega^2 are each >= 500.
  K2 no near-null direction: in each solve, full rank and sv_min / sv_max >= 1e-9.
  K3 comparison, for every mu in L* with |mu| <= 3 (table extended by 0 off its support), in both solves:
       |solved - table| <= 1e-5 (|mu| <= 1), 3e-4 (1 < |mu| <= 2), 4e-3 (2 < |mu| <= 3).
     Reported per class: (lambda-valuation of lambda^3 mu, unit class), off-support separately.
  K4 controls:  (a) the swap9 table (e^{+-2 pi i/9} exchanged) differs from S1 by > 1 somewhere in |mu| <= 2;
                (b) the system with the conjugate multiplier (chi -> conj chi) has relative residual >= 1e-6.
  Reported, no criterion: the system with chi = 1 rows only (rank, sv, max deviation from the table).
Tolerances are error estimates from the two-solve spread, not bounds.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy python -B -u audits/003dual/stage7_hejhal.py > audits/003dual/stage7_hejhal.out
"""
import sys, os, math, time
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import numpy as np
import stage7_lib as L
kb = L.kb

print(__doc__, flush=True)
TOL = [(1.0, 1e-5), (2.0, 3e-4), (3.0, 4e-3)]
MCMP = 3.0
CONF = {"S1": dict(v0=0.30, M=8.0, n_side=30, shift=(0.0, 0.0)),
        "S2": dict(v0=0.26, M=9.0, n_side=32, shift=(0.37, 0.61))}
verdict = {}


def tol_of(a):
    for edge, t in TOL:
        if a <= edge * (1 + 1e-12): return t
    return None


def build(name, conj_chi=False, chi1_only=False):
    c = CONF[name]
    nums = L.dual_indices(c["M"])
    mu = L.mu_of(nums)
    scale = c["v0"] * L.kv(1.0 / 3, 4 * math.pi * np.abs(mu) * c["v0"])
    pts = L.sample_points(c["n_side"], c["v0"], c["shift"])
    rows, rhs, ks = [], [], []
    info = {"chi": [0, 0, 0], "rep": {}, "skipped": 0, "skipped_rep": {}}
    for z in pts:
        cand, best_any = L.pullbacks(z, c["v0"], 40)
        cand = [x for x in cand if x[0] >= 1.05 * c["v0"]]
        if not cand:
            info["skipped"] += 1; info["skipped_rep"][best_any[1]] = info["skipped_rep"].get(best_any[1], 0) + 1
            continue
        chosen = [cand[0]]
        for x in cand[1:]:
            if len(chosen) >= 4: break
            if all(x[1] != y[1] or x[2] != y[2] for y in chosen): chosen.append(x)
        for v2, k, rep, z2, h in chosen:
            if chi1_only and k != 0: continue
            kk = (-k) % 3 if conj_chi else k
            f = L.W3 ** kk
            rows.append((L.phi_row(mu, z, c["v0"]) - f * L.phi_row(mu, z2, v2)) / scale)
            rhs.append(L.SIGMA * (f * v2 ** (2.0 / 3) - c["v0"] ** (2.0 / 3)))
            info["chi"][k] += 1; info["rep"][rep] = info["rep"].get(rep, 0) + 1
    return nums, mu, np.array(rows), np.array(rhs), scale, info


def table(Mmax, swap9=False):
    td = kb.ThetaData(200000, os.path.join(here, "cache"), verbose=False)
    th, _, _ = kb.build_series(td, Mmax, swap9=swap9)
    d = {}
    for m, cf in zip(map(tuple, th.num), th.coef):
        d[m] = d.get(m, 0) + cf
    return d


def klass(m):
    """(lambda-valuation of lambda^3 mu, unit index u in 0..5 of the cofactor mod lambda... ) as a label"""
    v = 0; x = m
    while True:
        q = kb.exact_div(x, kb.LAM)
        if q is None: break
        x = q; v += 1
    # cofactor x prime to lambda: x = +- omega^k * (primary)
    for k, u in enumerate(kb.UNITS):
        y = kb.mul(x, u)
        if kb.is_primary(y): break
    return v, k   # x * zeta6^k primary


# ---------------------------------------------------------------- solves
T9 = table(9.0)
for name in ("S1", "S2"):
    t0 = time.time()
    nums, mu, A, b, sc, info = build(name)
    tau, sv, rank, resid = L.solve(A, b, sc)
    print(f"\n== {name} {CONF[name]}: unknowns {len(nums)}, rows {len(b)}, chi exponents {info['chi']}, reps {info['rep']}, "
          f"points skipped {info['skipped']} {info['skipped_rep']}, rank {rank}, sv max {sv[0]:.3g} min {sv[-1]:.3g} "
          f"(ratio {sv[-1]/sv[0]:.2e}), smallest 5 sv {np.array2string(sv[-5:], precision=3)}, rel. residual {resid:.2e}, "
          f"{time.time()-t0:.0f}s", flush=True)
    verdict[f"K1 {name}"] = info["chi"][1] >= 500 and info["chi"][2] >= 500
    verdict[f"K2 {name}"] = rank == len(nums) and sv[-1] / sv[0] >= 1e-9
    if name == "S1":
        xt = np.array([T9.get(m, 0) for m in nums]) * sc
        r0 = np.linalg.norm(A @ xt - b) / np.linalg.norm(b)
        print(f"K0 code check: DR table plugged into the S1 system, relative residual {r0:.2e}")
        verdict["K0"] = r0 <= 1e-10
        S1 = dict(zip(nums, tau))
        # control (a): swap9
        T9s = table(9.0, swap9=True)
        dsw = max(abs(S1[m] - T9s.get(m, 0)) for m, a in zip(nums, np.abs(mu)) if a <= 2.0)
        print(f"K4a swap9 table vs S1, max |diff| over |mu| <= 2: {dsw:.3g}")
        verdict["K4a"] = dsw > 1
    else:
        S2 = dict(zip(nums, tau))
    # comparison
    sol = dict(zip(nums, tau))
    byclass = {}; worst_ratio = 0
    for m, a in zip(nums, np.abs(mu)):
        if a > MCMP * (1 + 1e-12): continue
        tv = T9.get(m, 0); dv = abs(sol[m] - tv); r = dv / tol_of(a)
        worst_ratio = max(worst_ratio, r)
        key = klass(m) if tv != 0 else "off-support"
        e = byclass.setdefault(key, [0, 0.0, 0.0, 9e9, 0.0])
        e[0] += 1; e[1] = max(e[1], dv); e[2] = max(e[2], r); e[3] = min(e[3], a); e[4] = max(e[4], a)
    print(f"K3 {name}: per class (lambda-valuation v of lambda^3 mu, unit k), count, max |solved - table|, max diff/tol, |mu| range")
    for key in sorted(byclass, key=str):
        n, dmax, rmax, amin, amax = byclass[key]
        lab = key if key == "off-support" else ("v=%d k=%d" % key)
        kind = ""
        if key != "off-support":
            v = key[0]
            kind = ("case 4, n=%d" % (v // 3)) if v % 3 == 0 else (("cases 1-3, n=%d" % ((v + 1) // 3)) if v % 3 == 2 else "??")
        print(f"   {lab:12s} {kind:16s} n={n:4d}  max diff {dmax:.2e}  max diff/tol {rmax:.3f}  |mu| in [{amin:.3f}, {amax:.3f}]")
    verdict[f"K3 {name}"] = worst_ratio <= 1
    if name == "S2":
        common = [m for m in S1 if m in S2 and abs(L.mu_of([m])[0]) <= MCMP * (1 + 1e-12)]
        print(f"S1 vs S2 over |mu| <= 3: max |diff| {max(abs(S1[m]-S2[m]) for m in common):.2e} ({len(common)} coefficients)")

# ---------------------------------------------------------------- controls on S1
nums, mu, A, b, sc, info = build("S1", conj_chi=True)
tau, sv, rank, resid = L.solve(A, b, sc)
print(f"\nK4b conjugate multiplier: rel. residual {resid:.2e} (rank {rank})")
verdict["K4b"] = resid >= 1e-6
nums, mu, A, b, sc, info = build("S1", chi1_only=True)
tau, sv, rank, resid = L.solve(A, b, sc)
dev = max(abs(t - T9.get(m, 0)) for m, t, a in zip(nums, tau, np.abs(mu)) if a <= 2.0)
print(f"chi = 1 rows only (report, no criterion): rows {len(b)}, rank {rank}/{len(nums)}, sv min/max {sv[-1]/sv[0]:.2e}, "
      f"rel. residual {resid:.2e}, max |solved - table| over |mu| <= 2: {dev:.2e}")

print("\nVERDICT:")
for k, v in verdict.items(): print(f"  {k}: {'PASS' if v else 'FAIL'}")
print("ALL PASS" if all(verdict.values()) else "NOT ALL PASS")
