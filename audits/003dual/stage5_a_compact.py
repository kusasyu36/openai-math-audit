# -*- coding: utf-8 -*-
"""
Stage 5 (a): the end-to-end comparison of layer E with a COMPACTLY SUPPORTED smooth weight, i.e. inside the hypothesis
of Proposition lem:reflection (W in C_c^infty((0, inf)), paper2.tex:670, :1293).

  W(y) = b(log(y)/a),  b(u) = exp(-(1-u^2)^-p) on |u| < 1,   p = 4, a = 3   (support [e^-3, e^3]).
  LEFT  T(X; Psi) of eq:T (:1287), a FINITE sum (W(N/X) = 0 for N >= X e^3; X e^3 <= 10^6 = table size, so no truncation).
  RIGHT sum over the groups of C * sum_l d(l) alpha(l) N(l)^-1/2 psi prod B V_*^sharp(N(l) X/N(c)^2)  (:1806, C of :3420,
        V^sharp of :1752), all table entries with N(l) <= R_g = min(table limit, 3e7 N(c)^2/X).
  The choice (p, a) and the runs were made from stage5_bump.out (bounds only).  p = 1 (the usual bump) was the first
  plan; its bounds were useless (>= 5e-2 per unit |C|), see stage5_bump.py.

CRITERIA (fixed before the first run of this script):
  F1  |T_direct - T_dual| <= Btot,  Btot = sum over groups of (T3') [upper bound; its ingredients Mt(A) are numerical
      quadratures of explicit envelopes] + (Q') quadrature of V^sharp [upper bound, same remark]
      + rounding allowance [ESTIMATE]: 1e-13 (sum|direct terms| + sum|dual terms|) + 2e-15 sum|dual coefficients|.
      A run is INFORMATIVE iff Btot <= 1e-2 |T_direct|.  For every run the number of digits log10(|T|/Btot) guaranteed by
      the bound is printed next to the observed |diff| (ABSOLUTE) and |diff|/|T| (RELATIVE).  The claim made for a run is
      "agreement within Btot", not within the observed difference.
  Controls (same definitions as NE1-NE7 of stage3_d_mellin), on informative runs where live; rejected iff
      |T_direct - T_ctl| > 10 Btot.   Required (chosen beforehand from the sizes seen with the log-Gaussian weight:
      NE1 2.0, NE2-NE5 >= 0.15, NE6/NE7 >= 0.058 relative):
        Btot <= 1e-2 |T|: NE1 ;   Btot <= 1e-3 |T|: also NE2, NE3, NE4, NE5 ;   Btot <= 1e-4 |T|: also NE6, NE7.
      Controls outside the required set are computed and printed but do not enter the verdict.
Runs (X = 49737 = the largest X with X e^3 <= 10^6, unless stated):
  expected informative from the bounds: B {} (X = 49737 and 20000), F {}, E {}, A {}, B {p7: 1}, B {p7: 4}, B {p7: 0}, B {p7: 2}
  expected NOT informative (run for the record; only the observed difference is reported): F {p7: 1}, A {p7: 1}, B {p13: 1}
  (C {} is the complex conjugate of F {} and is not run.)
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage5_a_compact.py > audits/003dual/stage5_a_compact.out 2>/dev/null
"""
from stage5_lib import *
import stage5_bump as s5
t0 = time.time()
print(__doc__)
P7, P13 = (1, 3), (4, 3)
XB = 49737.0
RUNS = [("B", {}, XB), ("B", {}, 20000.0), ("F", {}, XB), ("E", {}, XB), ("A", {}, XB),
        ("B", {P7: 1}, XB), ("B", {P7: 4}, XB), ("B", {P7: 0}, XB), ("B", {P7: 2}, XB),
        ("F", {P7: 1}, XB), ("A", {P7: 1}, XB), ("B", {P13: 1}, XB)]
bp = s5.Bump(3.0, 4); bp.setup_quadrature()
wt = WBump(bp, 3.0e7)
print("  weight: %s ; V^sharp quadrature: |tau| <= %.1f, step %.2f, %d nodes ; truncation bound %.1e ; max |V-hat|/envelope = %.3f  [%.0fs]"
      % (wt.name, bp.TM, bp.H, len(bp.tau), bp.trunc, bp.env_check, time.time() - t0))
E = EndToEnd()
REQ = [(1e-2, ["NE1"]), (1e-3, ["NE2", "NE3", "NE4", "NE5"]), (1e-4, ["NE6", "NE7"])]
ALLC = ["NE1", "NE2", "NE3", "NE4", "NE5", "NE6", "NE7"]
nfail = ninf = 0; ctl = {k: [0, 0, 1e9] for k in ALLC}; rows = []
for nm, jd, X in RUNS:
    t1 = time.time()
    r = E.run(nm, jd, wt, X); Td = r["Td"]; diff = abs(Td - r["Tu"])
    ratio = r["Btot"] / abs(Td); inf_ = ratio <= 1e-2; ok = diff <= r["Btot"]
    print("  Psi_0=%s P=%s X=%g #groups=%d N(c) in %s cases %s ; dual terms: %d" % (nm, {norm(p): j for p, j in jd.items()}, X, len(r["G"]),
          sorted(set(g["key"][1] for g in r["G"])), sorted(set((g["cfg"]["case"], g["cfg"]["sigma"]) for g in r["G"])), sum(g["key"][2] for g in r["G"])))
    print("      T_dir = % .12e %+.12ei   T_dual = % .12e %+.12ei" % (Td.real, Td.imag, r["Tu"].real, r["Tu"].imag))
    print("      |T_dir| = %.4e  |diff| = %.2e (abs)  %.2e (rel)   Btot = %.2e = dual tail %.1e + quadrature %.1e + rounding allowance (estimate) %.1e (direct tail %.1g) ; digits guaranteed by the bound: %.1f   %s%s  [%.0fs]"
          % (abs(Td), diff, diff / abs(Td), r["Btot"], r["b_dual"], r["qb"], r["b_round"], r["b_dir"], -math.log10(ratio), "within Btot" if ok else "EXCEEDS Btot", "" if inf_ else "  (NOT INFORMATIVE: Btot > 1e-2 |T|)", time.time() - t1))
    rows.append((nm, jd, X, abs(Td), diff, r["Btot"], inf_, ok))
    if not inf_: ninf += 1; continue
    nfail += not ok
    req = [k for lim, ks in REQ if ratio <= lim for k in ks]; line = []
    for k in ALLC:
        Tc = E.control(k, nm, jd, wt, X, r["G"])
        if Tc is None: line.append("%s dead" % k); continue
        d = abs(Td - Tc)
        if k in req: ctl[k][0] += 1; ctl[k][1] += d > 10 * r["Btot"]; ctl[k][2] = min(ctl[k][2], d / r["Btot"])
        line.append("%s%s %.1e" % (k, "*" if k in req else "", d / abs(Td)))
    print("      controls (* = required in this run), |T_dir - T_ctl|/|T_dir|: " + "  ".join(line)); sys.stdout.flush()
print("\n  F1 %s: %d runs exceed Btot / %d informative runs (%d of %d runs not informative)" % ("PASS" if nfail == 0 and ninf < len(RUNS) else "FAIL", nfail, len(RUNS) - ninf, ninf, len(RUNS)))
print("      non-informative runs, observed only: " + " ; ".join("%s %s: |diff|/|T| = %.1e, Btot/|T| = %.1e" % (nm, {norm(p): j for p, j in jd.items()}, d / T, B / T) for nm, jd, X, T, d, B, i, ok in rows if not i))
allok = nfail == 0
for k in ALLC:
    ok = ctl[k][1] == ctl[k][0]; allok &= ok
    print("  %s %s: required in %d runs, rejected (> 10 Btot) in %d, smallest |T_dir - T_ctl|/Btot = %.1e" % (k, ("PASS" if ok else "FAIL") if ctl[k][0] else "not required anywhere", ctl[k][0], ctl[k][1], ctl[k][2]))
print("\nSTAGE 5 (a): %s   [%.0fs]" % ("all criteria met" if allok else "SOME CRITERIA NOT MET", time.time() - t0))
