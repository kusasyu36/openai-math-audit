# -*- coding: utf-8 -*-
"""
Stage 5 (a), POST-HOC re-evaluation with a tightened tail bound.  Written AFTER stage5_a_compact.out was seen: there the
one-prime runs B {p7: 1}, {p7: 4}, {p7: 0} missed the registered threshold "informative iff Btot <= 1e-2 |T|" by factors
1.2-3.2 (guaranteed digits 1.5-1.9), while the observed differences were 5e-10 to 2e-9 relative.  The registered verdict of
stage5_a_compact ("not informative") stands and is not replaced by this file.

What is tightened (two facts about the sequences d_+, d_- as implemented in kubota.build_series, both checked below on
the tables):  (i) d_+-(l) != 0 only for lambda^4 l in  -u (1 + 3O),  u in {1, omega, omega^2}: three translates of a
lattice of minimal distance 1/3 after division by 9 (the bound (L) of stage4_bounds used all of lambda^-4 O, minimal
distance 1/9: 3 times more points);  (ii) |d_+-(l)| <= 9 * 81^(1/6) N(l)^(1/6) = 18.7208 N(l)^(1/6)  (instead of 27).
For groups with sigma = +-:  tail <= |C| Bmax 18.7208 * min_A { (c X/N(c)^2)^-A Mt(A) * 3 S(1/3 + A, R; d = 1/3) }.
Groups with sigma = 0 keep the old bound.  Everything else (weight p = 4, a = 3, X = 49737, V^sharp, rounding allowance)
is as in stage5_a_compact.py.

CRITERIA (fixed before this script was run, but after stage5_a_compact):
  F1' |T_direct - T_dual| <= Btot' ; informative iff Btot' <= 1e-2 |T_direct| ; control NE1 (+i/81) rejected iff
      |T_direct - T_ctl| > 10 Btot' (the only control required at this precision, as registered in stage5_a_compact);
      NE2-NE5 are required if Btot' <= 1e-3 |T|.
Runs: B {p7: 1}, B {p7: 4}, B {p7: 0}, B {p7: 2}, B {p13: 1}, F {p7: 1}.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage5_a2_tight.py > audits/003dual/stage5_a2_tight.out 2>/dev/null
"""
from stage5_lib import *
import stage5_bump as s5
t0 = time.time()
print(__doc__)
P7, P13 = (1, 3), (4, 3); XB = 49737.0
RUNS = [("B", {P7: 1}), ("B", {P7: 4}), ("B", {P7: 0}), ("B", {P7: 2}), ("B", {P13: 1}), ("F", {P7: 1})]
bp = s5.Bump(3.0, 4); bp.setup_quadrature(); wt = WBump(bp, 3.0e7)
E = EndToEnd(); th = E.th
# checks of (i) and (ii) on the tables
for sg in ("+", "-"):
    s = th.d[sg]; nn = int(np.searchsorted(s.aell, th.complete[sg], side="right")); m = s.m[:nn]; ok = np.zeros(nn, dtype=bool)
    for u in ((1, 0), (0, 1), (-1, -1)):
        # m = -u A  <=>  A = -m conj(u)  (u a cube root of unity, u^-1 = conj u)
        ub = conj(u); A0 = -(m[:, 0] * ub[0] - m[:, 1] * ub[1]); A1 = -(m[:, 0] * ub[1] + m[:, 1] * ub[0] - m[:, 1] * ub[1])
        ok |= (A0 % 3 == 1) & (A1 % 3 == 0)
    okp = ok | False
    # the conjugate series lives on the negatives: accept either sign convention, but report which
    mm = -m; ok2 = np.zeros(nn, dtype=bool)
    for u in ((1, 0), (0, 1), (-1, -1)):
        ub = conj(u); A0 = -(mm[:, 0] * ub[0] - mm[:, 1] * ub[1]); A1 = -(mm[:, 0] * ub[1] + mm[:, 1] * ub[0] - mm[:, 1] * ub[1])
        ok2 |= (A0 % 3 == 1) & (A1 % 3 == 0)
    print("  d_%s: %d table entries ; in -u(1+3O): %d ; in +u(1+3O): %d ; max |d|/N^(1/6) = %.4f" % (sg, nn, int(okp.sum()), int(ok2.sum()), float(np.max(np.abs(s.coef[:nn]) / s.aell[:nn] ** (1.0 / 3)))))
    assert okp.all() or ok2.all(), "support claim (i) fails"
def tight_tail(X, Nc, R, Cabs, Bmax):
    best = float("inf")
    for A in s5.A_GRID:
        Mv = bp.Mt(A)
        if Mv == float("inf"): continue
        best = min(best, math.exp(-A * math.log(sb.CK * X / Nc ** 2)) * Mv * 3 * sb.S_pow(1.0 / 3 + A, R, 1.0 / 3))
    return Cabs * Bmax * 18.7208 * best
for nm, jd in RUNS:
    r = E.run(nm, jd, wt, XB); Td = r["Td"]; diff = abs(Td - r["Tu"])
    bd = 0.0
    for g in r["G"]:
        f = tight_tail if g["cfg"]["sigma"] != "0" else wt.dual_tail
        bd += f(XB, g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"])
    B2 = bd + r["qb"] + r["b_round"]; ratio = B2 / abs(Td); inf_ = ratio <= 1e-2
    print("  Psi_0=%s P=%s X=%g : |T_dir| = %.4e  |diff| = %.2e (abs) %.2e (rel) ; dual tail: registered %.2e -> tightened %.2e ; Btot' = %.2e ; digits guaranteed %.1f ; %s%s"
          % (nm, {norm(p): j for p, j in jd.items()}, XB, abs(Td), diff, diff / abs(Td), r["b_dual"], bd, B2, -math.log10(ratio), "within Btot'" if diff <= B2 else "EXCEEDS Btot'", "" if inf_ else " (NOT INFORMATIVE)"))
    if inf_:
        req = ["NE1"] + (["NE2", "NE3", "NE4", "NE5"] if ratio <= 1e-3 else [])
        line = []
        for k in req:
            Tc = E.control(k, nm, jd, wt, XB, r["G"])
            if Tc is None: line.append("%s dead" % k); continue
            d = abs(Td - Tc); line.append("%s: |T - T_ctl|/Btot' = %.1e (%s)" % (k, d / B2, "rejected" if d > 10 * B2 else "NOT rejected"))
        print("      " + " ; ".join(line))
    sys.stdout.flush()
print("\nSTAGE 5 (a), post-hoc tightened bound: done   [%.0fs]" % (time.time() - t0))
