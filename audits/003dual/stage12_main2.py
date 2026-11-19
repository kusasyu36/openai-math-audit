# -*- coding: utf-8 -*-
"""
Stage 12, SECOND pre-registered comparison (main2).  Auditor: Claude Opus 5.5 (stage-12 subagent), 2026-11-19.
This docstring was written BEFORE the first run of this script.

HISTORY (read before the verdict).
  * stage12_main.out (first pre-registered run) printed "NOT MET" for both j = 1 and 3|c.  No run disagreed (every
    informative run was within Btot12); the cause was the control rule: NE2 was REQUIRED on the two A {} runs but its switch
    changes nothing there (alpha(c) = +-1 for every 3|c group, so conj(alpha)^2 = alpha^2), and the registered rule turned a
    dead required control into a failure of BOTH verdicts (condition (iv) was aggregated over all informative runs).
    That verdict stands and is NOT replaced by this file.
  * The rule below was corrected AFTER stage12_main.out was seen: the required controls are now LISTED PER RUN, from the
    liveness table of stage12_probe2.out (computed without any dual sum), and each verdict only looks at its own runs.
  * In stage12_main.out the six informative j = 1 main-frame runs already agreed (|diff| <= 8.2e-10, within Btot12).  With
    the same tables and the same formulas, the outcome of this run is therefore largely predictable; what this run adds is
    a clean pre-registered verdict on two weights never compared before.
  * Weights: p=3, a=2.5 and p=5, a=2.5, X = 82084 (planned in stage12_probe1.out, never compared; X e^a <= 10^6).
    p7bar is dropped (its run is exactly the complex conjugate of the p7 run, seen in stage12_main.out).

COMPARED, ERROR BUDGET: as in stage12_main.py, except that the structured tail is evaluated with struct_tail2
  (D at R (1 - 1e-9), closing a boundary gap of size ~1e-12 found after stage12_main: a table term with N(l) = R exactly
  could be in neither the table sum nor the bound).  Btot12 = structured tail (Mt by one-sided Riemann sum) + (Q') + rounding
  allowance (ESTIMATE).  Same premise and the same unbounded numerical steps as in stage12_main.py; the digit count is
  "digits of the evaluated error bound", not a guarantee.
CRITERIA.
  G0 as in stage12_main.py (gate on the derivation).
  F1: PASS iff |T_dir - T_dual| <= Btot12 ; INFORMATIVE iff Btot12 <= 1e-2 |T_dir| (unchanged 2-digit threshold).
  Controls: rejected iff |T_dir - T_ctl| > 10 Btot12.  REQUIRED controls, fixed here per run (from stage12_probe2.out;
  digits of the bound part: W3 3.61 / 2.69 / 3.81, W4 3.72 / 2.30 / 3.84, at least 0.16 digit from a tier boundary):
     B {p7:1}  (W3, W4): NE1, J, NE2, NE3, NE4, NE5
     B {p13:1} (W3, W4): NE1, J
     A {}      (W3, W4): NE1, NE3, NE4, NE5        [NE2 is dead by structure: alpha(c) = +-1 on all groups; not required]
  The script re-checks liveness in code before each control's dual sum; a REQUIRED control found dead = design failure
  (verdict NOT MET).  Non-required controls are computed and printed only.
  VERDICT j = 1 (runs B {p7:1}, B {p13:1} of W3 and W4 only): MET iff G0, at least one of these runs is informative, every
      informative one passes F1, and every required control of every informative one is live and rejected.
  VERDICT 3|c (runs A {} of W3 and W4 only): the same with these runs.
  SECONDARY (reference, registered here): the same two verdicts with the paper-only bound D_paper (stage12_lib; uses only
      the printed statement :1858-1867, 6 units, all m >= -4, |d| <= 27 3^{m/6}|b|) in place of the structured tail, with the
      same required controls.  Expected from probe 2: informative for B {p7:1} (2.94, 3.04 digits) and A {} (2.86, 2.89),
      B {p13:1} only for W3 (2.01, margin 0.01) and not for W4 (1.62).
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage12_main2.py > audits/003dual/stage12_main2.out
"""
from stage5_lib import *
import stage5_bump as s5
import stage12_lib as L
import stage12_ctl as K
t0 = time.time()
print(__doc__)
E = EndToEnd(); th = E.th; td = th.K.td
NA_all = E.NA; absG = np.abs(td.G)
P7, P13 = (1, 3), (4, 3)

# ---------------- G0 (same as stage12_main)
def table_sum(sg, s, R1, R2):
    S = th.d[sg]; N = S.aell ** 2; m = (N > R1) & (N <= R2)
    return float(np.sum(np.abs(S.coef[m]) * N[m] ** (-s)))
def formula_sum(sg, s, R1, R2):
    def win(scale, mult, units):
        Nl = NA_all * scale; m = (Nl > R1) & (Nl <= R2)
        return units * mult * float(np.sum(absG[m] * Nl[m] ** (-s)))
    if sg == "0":
        return sum(win(3.0 ** (3 * k - 3), 3 ** (k / 2 + 2.5), 2) for k in range(0, 12)) + sum(win(3.0 ** (3 * k - 4), 3 ** (k / 2 + 2), 6) for k in range(1, 12))
    return win(1 / 81.0, 9.0, 3)
g0a = 0.0
for sg in ("0", "+", "-"):
    lim = th.complete[sg] ** 2 * (1 - 1e-9)
    for (R1, R2) in ((0.0, lim / 100 * 1.000000001), (lim / 100 * 1.000000001, lim / 4 * 1.000000001), (lim / 4 * 1.000000001, lim)):
        for s in (0.0, 1.5, 3.0):
            a_ = table_sum(sg, s, R1, R2); b_ = formula_sum(sg, s, R1, R2); g0a = max(g0a, abs(a_ - b_) / a_)
g0b = True; pts = []
for s in (1.5, 2.5, 4.0):
    for M in (0.5, 10.0, 1e3, 1e5):
        msk = NA_all > M; part = float(np.sum(absG[msk] * NA_all[msk] ** (-s))); bd = L.P(s, M); g0b &= bd >= part; pts.append(bd / part)
G0 = g0a <= 1e-12 and g0b
print("G0 %s: (a) %.1e (limit 1e-12); (b) %s (smallest ratio %.3f)" % ("PASS" if G0 else "FAIL", g0a, g0b, min(pts)))
for p in (P7, P13): assert float(np.max(np.abs(loc(p).Btab(1)))) <= 1.0 + 1e-12
sys.stdout.flush()

REQUIRED = {("B", P7): ["NE1", "J", "NE2", "NE3", "NE4", "NE5"], ("B", P13): ["NE1", "J"], ("A", None): ["NE1", "NE3", "NE4", "NE5"]}
ALLC = ["NE1", "J", "NE2", "NE3", "NE4", "NE5", "NE6", "NE7"]
rows = []
for label, p, a in (("W3", 3, 2.5), ("W4", 5, 2.5)):
    E._VS.clear()
    X = 82084.0
    bp = s5.Bump(a, p); bp.setup_quadrature(); wt = WBump(bp, 3.0e7); mu = L.MtUp(bp)
    assert X * math.exp(a) <= 1.0e6
    print("\n=== frame %s: weight p=%d a=%g X=%d ; V^sharp quadrature |tau| <= %.1f, %d nodes, truncation bound %.1e, max |V-hat|/envelope = %.3f"
          % (label, p, a, X, bp.TM, len(bp.tau), bp.trunc, bp.env_check))
    for nm, jd in (("B", {P7: 1}), ("B", {P13: 1}), ("A", {})):
        t1 = time.time()
        key = (nm, next(iter(jd)) if jd else None); req = REQUIRED[key]
        r = E.run(nm, jd, wt, X); Td = r["Td"]; Tu = r["Tu"]; diff = abs(Td - Tu); G = r["G"]
        st = sum(L.struct_tail2(mu, X, g["cfg"]["sigma"], g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"])[0] for g in G)
        sp = sum(L.struct_tail2(mu, X, g["cfg"]["sigma"], g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"], paper=True)[0] for g in G)
        B12 = st + r["qb"] + r["b_round"]; Bp = sp + r["qb"] + r["b_round"]
        inf_ = B12 <= 1e-2 * abs(Td); ok = diff <= B12; infp = Bp <= 1e-2 * abs(Td); okp = diff <= Bp
        print("  [%s] Psi_0=%s P=%s : groups %d, N(c) in %s, dual terms %d" % (label, nm, jd, len(G), sorted(set(g["key"][1] for g in G)), sum(g["key"][2] for g in G)))
        print("      T_dir = % .13e %+.13ei   T_dual = % .13e %+.13ei" % (Td.real, Td.imag, Tu.real, Tu.imag))
        print("      |T_dir| = %.4e  |diff| = %.2e (abs) %.2e (rel) ; Btot12 = %.2e (structured tail %.2e + quadrature %.1e + rounding allowance (estimate) %.1e), digits %.2f, %s%s ;"
              " paper-only: Btot = %.2e, digits %.2f, %s%s  [%.0fs]"
              % (abs(Td), diff, diff / abs(Td), B12, st, r["qb"], r["b_round"], math.log10(abs(Td) / B12), "within" if ok else "EXCEEDS", "" if inf_ else " (NOT INFORMATIVE)",
                 Bp, math.log10(abs(Td) / Bp), "within" if okp else "EXCEEDS", "" if infp else " (NOT INFORMATIVE)", time.time() - t1))
        ctl = {}; line = []
        if inf_ or infp:
            for k in ALLC:
                if k == "J":
                    if not any(j == 1 for j in jd.values()): continue
                    lv, Gj = K.liveness(E, "J", nm, jd, wt, X, G)
                else:
                    lv = K.liveness(E, k, nm, jd, wt, X, G); Gj = None
                if not (lv > 0):
                    ctl[k] = (None, k in req); line.append("%s%s DEAD" % (k, "*" if k in req else "")); continue
                if k == "J": Tc = E.dual(Gj, wt, X)[0]; del Gj
                elif k == "NE5": Tc = E.dual(E.groups(nm, jd, wt, X, E.PHno2[nm]), wt, X)[0]
                else: Tc = E.dual(G, wt, X, variant=k)[0]
                d = abs(Td - Tc); ctl[k] = (d, k in req)
                line.append("%s%s live(%.1e) |T-Tc|/|T|=%.1e /Btot12=%.1e /Btot_paper=%.1e" % (k, "*" if k in req else "", lv, d / abs(Td), d / B12, d / Bp))
            print("      controls (* = required): " + " ; ".join(line))
        rows.append(dict(label=label, nm=nm, jd=jd, T=abs(Td), diff=diff, B=B12, Bp=Bp, inf=inf_, ok=ok, infp=infp, okp=okp, ctl=ctl, req=req))
        del G, r; sys.stdout.flush()

def verdict(sel, paper=False):
    rs = [r for r in rows if sel(r)]
    Bk, ik, ok_ = ("Bp", "infp", "okp") if paper else ("B", "inf", "ok")
    inf_ = [r for r in rs if r[ik]]
    c2 = len(inf_) > 0; c3 = all(r[ok_] for r in inf_)
    c4 = all((k in r["ctl"]) and r["ctl"][k][0] is not None and r["ctl"][k][0] > 10 * r[Bk] for r in inf_ for k in r["req"])
    return G0 and c2 and c3 and c4, len(inf_), len(rs), c3, c4
print("\n=== summary")
for r in rows:
    print("  %s %s %-12s |T|=%.3e |diff|=%.2e | Btot12=%.2e digits %.2f %s %s | paper-only Btot=%.2e digits %.2f %s %s | required: %s"
          % (r["label"], r["nm"], r["jd"], r["T"], r["diff"], r["B"], math.log10(r["T"] / r["B"]), "inf" if r["inf"] else "NOT-inf", "within" if r["ok"] else "EXCEEDS",
             r["Bp"], math.log10(r["T"] / r["Bp"]), "inf" if r["infp"] else "NOT-inf", "within" if r["okp"] else "EXCEEDS",
             ", ".join("%s:%s" % (k, ("DEAD" if r["ctl"].get(k, (None,))[0] is None else "%.1e*Btot12" % (r["ctl"][k][0] / r["B"]))) for k in r["req"])))
for name, sel in (("j=1", lambda r: r["nm"] == "B"), ("3|c", lambda r: r["nm"] == "A")):
    for paper in (False, True):
        v, ni, nt, c3, c4 = verdict(sel, paper)
        print("%s VERDICT %s : %s   (G0 %s, informative %d of %d, all informative within %s, all required controls live and rejected %s)"
              % ("SECONDARY (paper-only bound)" if paper else "MAIN2", name, "MET" if v else "NOT MET", G0, ni, nt, c3, c4))
print("[%.0fs]" % (time.time() - t0))
