# -*- coding: utf-8 -*-
"""
Stage 12 (main): the end-to-end comparison of layer E with a COMPACTLY SUPPORTED weight (inside the hypothesis of
Proposition lem:reflection, W in C_c^infty((0, inf)), paper2.tex :670, :1293), aimed at the one-prime case j = 1 and at 3|c.
Auditor: Claude Opus 5.5 (stage-12 subagent), 2026-11-19.  This docstring was written BEFORE the first run of this script.

HISTORY THAT MATTERS FOR READING THE VERDICT.
  * The tail bound used here (stage12_lib: structure of the coefficients, |d| = 9|G(A)|, |G(A)| = N(b)^{1/2}, summed as a
    double sum over (n, b)) was designed AFTER the stage-5 results had been seen (stage5_a_compact.out: B {p7:1} at
    p=4, a=3, X=49737 missed the threshold with 1.8 digits; observed |diff| 6.9e-12).
  * The weights and X of the MAIN FRAME below were chosen from stage12_probe1.out, which contains only bounds and the
    direct sum |T_dir| (no dual sum, no comparison).  No comparison has ever been run with these (weight, X).
  * The KNOWN FRAME re-uses the stage-5 setting (p=4, a=3, X=49737) whose differences were already seen; it is reported
    separately and is labelled post hoc.  It does not enter the main verdict.

WEIGHT.  W(y) = b(log(y)/a), b(u) = exp(-(1-u^2)^-p) on |u| < 1 (stage5_bump.Bump).  X = floor(10^6 e^-a), so X e^a <= 10^6 and
  the direct sum is a FINITE sum over the table (no truncation).
COMPARED.  T_dir = T(X; Psi) of eq:T (:1287)   vs   T_dual = sum over groups of C sum_l d alpha N^-1/2 psi prod B V^sharp
  (:1806, C of :3420, V^sharp of :1752), table entries with N(l) <= R_g (stage5_lib.EndToEnd, unchanged).
ERROR BUDGET  Btot12 = (sum over groups of the STRUCTURED tail bound of stage12_lib, with Mt evaluated by the one-sided
  Riemann sum Mt_upper) + (Q') quadrature bound of V^sharp (stage5_bump, unchanged) + rounding allowance (ESTIMATE, as in
  stage 5: 1e-13 (sum|direct terms| + sum|dual terms|) + 2e-15 sum|dual coefficients|).
  Status of Btot12: the formula is an upper bound under the PREMISE that the coefficient formulas of kubota.build_series
  (DR (5.7), (5.13), (5.14)) hold beyond the table; its EVALUATION contains unbounded numerical steps (the 800-node
  Gauss-Legendre integral inside the envelope, the tau-integral beyond the grid end, floating-point rounding).  Following
  REPORT section 10 the resulting digit count is called "digits of the evaluated error bound", not "guaranteed digits".

CRITERIA (the stage-5 criteria, unchanged in substance):
  G0 (gate on the derivation, not on the comparison): (a) the family decomposition of stage12_lib reproduces
      sum |d_sigma(l)| N(l)^-s on windows of the table (window ends moved off the lattice by a factor 1 - 1e-9; probe 1
      showed a boundary artefact at exactly N(l) = 10^6/81) to relative 1e-12, for sigma in {0,+,-} and s in {0, 1.5, 3};
      (b) the P-bound is >= the table partial sum at 12 points.  If G0 fails, no verdict is drawn.
  F1  per run: PASS iff |T_dir - T_dual| <= Btot12.   INFORMATIVE iff Btot12 <= 1e-2 |T_dir|  (the stage-5 threshold of
      2 digits; unchanged).  Digits = log10(|T_dir|/Btot12).  The claim for a run is "agreement within Btot12".
  Controls, on informative runs only; a control is REJECTED iff |T_dir - T_ctl| > 10 Btot12.
      Required (stage-5 rule):  Btot12 <= 1e-2|T|: NE1 and J ;  <= 1e-3|T|: also NE2, NE3, NE4, NE5 ;  <= 1e-4|T|: also NE6, NE7.
      NE1-NE7 as in stage3_d_mellin / stage5_lib.  NEW control J (dual side only): the factor B_{p,1} of the active prime
      is replaced by B_{p,5} in the dual sum (direct side and the constant C unchanged).  Expected: |T_dir - T_J| of the
      order of |T_dir| (a different character sum), i.e. far above 10 Btot12 when Btot12 <= 1e-2 |T|.
      LIVENESS (lesson 8): before the dual sum of a control is computed, the script checks in code that the switch
      changes its ingredient (max |changed - original| > 0 over the groups: C for NE1/NE2/NE6, arr for NE3/J, the phi-hat
      table for NE5, V^sharp at sample points for NE4/NE7) and prints the size.  A control whose switch changes nothing is
      reported DEAD and is not counted; if a REQUIRED control is dead, that is reported as a design failure.
  MAIN VERDICT for j = 1: "met" iff (i) G0 holds, (ii) at least one MAIN-FRAME run with a j = 1 prime is informative,
      (iii) every informative main-frame run passes F1, (iv) every required control on the informative main-frame runs
      is live and rejected.  The same statement is made separately for 3|c (runs with Psi_0 = A).
RUNS.
  MAIN FRAME (never compared before):
    W1: p=4, a=2.5, X=82084 :  B {p7:1}, B {p7bar:1}, B {p13:1}, A {}, B {p7:2} ; for the record (expected NOT informative
        from probe 1): F {p7:1}, A {p7:1}
    W2: p=4, a=2,   X=135335:  B {p7:1}, B {p7bar:1}, B {p13:1}, A {} ; for the record: F {p7:1}
    (p7 = (1,3), p7bar = (-2,-3), p13 = (4,3) as integer pairs a + b omega.)
    Expected digits from probe 1: W1: 4.0, 4.0, 2.8, 3.3, 3.9 (F 1.6, A{p7:1} -0.3); W2: 3.9, 3.9, 2.8, 3.5 (F 1.1).
  KNOWN FRAME (post hoc, stage-5 setting p=4, a=3, X=49737): B {p7:1}, B {p7:0}, B {p7:4}, A {}.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage12_main.py > audits/003dual/stage12_main.out
SMOKE TEST (machinery only, on the KNOWN setting B {p7:1}, p=4 a=3 X=49737, whose result was already seen in stage 5):
  STAGE12_SMOKE=1 ... stage12_main.py > audits/003dual/stage12_smoke.out   (run before the main run; not part of the verdict)
"""
from stage5_lib import *
import stage5_bump as s5
import stage12_lib as L
t0 = time.time()
print(__doc__)
E = EndToEnd(); th = E.th; td = th.K.td
NA_all = E.NA; absG = np.abs(td.G)
P7, P7b, P13 = (1, 3), (-2, -3), (4, 3)

# ------------------------------------------------------------------ G0
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
print("G0 %s: (a) worst relative difference table vs family decomposition = %.1e (limit 1e-12); (b) P-bound >= partial sums at 12 points: %s (smallest ratio %.3f)"
      % ("PASS" if G0 else "FAIL", g0a, g0b, min(pts)))
for p in (P7, P7b, P13): assert float(np.max(np.abs(loc(p).Btab(1)))) <= 1.0 + 1e-12
sys.stdout.flush()

# ------------------------------------------------------------------ controls: liveness and J
def groups_jdual(nm, jd, jd_dual, wt, X):
    """E.groups with B_{p, jd_dual[p]} in the dual coefficients, everything else (C included) from jd"""
    G = E.groups(nm, jd, wt, X)
    for g in G:
        cfg = g["cfg"]; sg, Nc, n = g["key"]; m = th.d[sg].m[:n]
        arr = g["arr"].copy()
        for p in jd:
            if p not in [q for q in cfg["primes"]]: continue
            lp = loc(p); i = lp.idx_arr(m)
            old = lp.Btab(jd[p])[i]; new = lp.Btab(jd_dual[p])[i]
            nz = old != 0
            arr[nz] = arr[nz] / old[nz] * new[nz]
            assert np.all(new[~nz] == 0)
        g["arr_orig"] = g["arr"]; g["arr"] = arr
    return G

def liveness(k, nm, jd, wt, X, G):
    """size of the change made by the switch of control k to its own ingredient (0 means dead)"""
    if k == "NE1": return max(float(np.max(np.abs(2 * (-1j / 81) * np.conj(g["ac"]) ** 2 * g["Cp"]))) for g in G)
    if k == "NE2": return max(abs(-1j / 81 * (g["ac"] ** 2 - np.conj(g["ac"]) ** 2) * g["Cp"]) for g in G)
    if k == "NE6": return max(abs(-1j / 81 * np.conj(g["ac"]) ** 2 * g["Cp"] * (W3 ** ((2 * g["k0"]) % 3) - 1)) for g in G)
    if k == "NE3":
        return max(float(np.max(np.abs(g["arr"] * (np.conj(th.d[g["key"][0]].ell[:g["key"][2]]) / th.d[g["key"][0]].ell[:g["key"][2]] - 1)))) if g["key"][2] else 0.0 for g in G)
    if k == "NE5":
        if nm not in E.PHno2: return 0.0
        a_, b_ = E.PH[nm], E.PHno2[nm]
        return max(abs(a_.get(h, 0) - b_.get(h, 0)) for h in set(a_) | set(b_))
    if k in ("NE4", "NE7"):
        best = 0.0
        for g in G:
            sg, Nc, n = g["key"]
            if n == 0: continue
            mm = th.d[sg].m[:min(n, 50)]; Nm = (mm[:, 0] ** 2 - mm[:, 0] * mm[:, 1] + mm[:, 1] ** 2).astype(float)
            x0 = Nm / 81.0 * X / (Nc * Nc)
            v0 = wt.vsharp(x0)
            if k == "NE4": v1 = wt.vsharp(Nm / 81.0 * X / Nc)
            else: v1 = wt.vsharp(x0, invert=True)
            best = max(best, float(np.max(np.abs(v1 - v0))))
        return best
    if k == "J":
        Gj = groups_jdual(nm, jd, {p: 5 for p in jd}, wt, X)
        val = max(float(np.max(np.abs(g["arr"] - g["arr_orig"]))) if g["key"][2] else 0.0 for g in Gj)
        return val, Gj
    raise KeyError(k)

REQ = [(1e-2, ["NE1", "J"]), (1e-3, ["NE2", "NE3", "NE4", "NE5"]), (1e-4, ["NE6", "NE7"])]
ALLC = ["NE1", "J", "NE2", "NE3", "NE4", "NE5", "NE6", "NE7"]

def run_frame(label, p, a, X, runs, frame):
    E._VS.clear()
    bp = s5.Bump(a, p); bp.setup_quadrature(); wt = WBump(bp, 3.0e7); mu = L.MtUp(bp)
    assert X * math.exp(a) <= 1.0e6
    print("\n=== %s frame %s: weight p=%d a=%g X=%d ; V^sharp quadrature |tau| <= %.1f, %d nodes, truncation bound %.1e, max |V-hat|/envelope = %.3f"
          % (frame, label, p, a, X, bp.TM, len(bp.tau), bp.trunc, bp.env_check))
    res = []
    for nm, jd in runs:
        t1 = time.time()
        r = E.run(nm, jd, wt, X); Td = r["Td"]; Tu = r["Tu"]; diff = abs(Td - Tu); G = r["G"]
        st = 0.0
        for g in G: st += L.struct_tail(mu, X, g["cfg"]["sigma"], g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"])[0]
        B12 = st + r["qb"] + r["b_round"]; ratio = B12 / abs(Td); inf_ = ratio <= 1e-2; ok = diff <= B12
        j1 = any(j == 1 for j in jd.values())
        print("  [%s] Psi_0=%s P=%s : groups %d, N(c) in %s, cases %s, dual terms %d"
              % (label, nm, {q: j for q, j in jd.items()}, len(G), sorted(set(g["key"][1] for g in G)), sorted(set((g["cfg"]["case"], g["cfg"]["sigma"]) for g in G)), sum(g["key"][2] for g in G)))
        print("      T_dir = % .13e %+.13ei   T_dual = % .13e %+.13ei" % (Td.real, Td.imag, Tu.real, Tu.imag))
        print("      |T_dir| = %.4e  |diff| = %.2e (abs) %.2e (rel) ; Btot12 = %.2e = structured tail %.2e + quadrature %.1e + rounding allowance (estimate) %.1e"
              " [stage-5 registered tail %.2e] ; digits of the evaluated bound %.2f ; %s%s  [%.0fs]"
              % (abs(Td), diff, diff / abs(Td), B12, st, r["qb"], r["b_round"], r["b_dual"], -math.log10(ratio), "within Btot12" if ok else "EXCEEDS Btot12",
                 "" if inf_ else "  (NOT INFORMATIVE: Btot12 > 1e-2 |T|)", time.time() - t1))
        row = dict(frame=frame, label=label, nm=nm, jd=jd, j1=j1, T=abs(Td), diff=diff, B=B12, inf=inf_, ok=ok, ctl={})
        if inf_:
            req = [k for lim, ks in REQ if ratio <= lim for k in ks]; line = []
            for k in ALLC:
                if k == "J":
                    if not any(j == 1 for j in jd.values()): continue
                    lv, Gj = liveness("J", nm, jd, wt, X, G)
                else:
                    lv = liveness(k, nm, jd, wt, X, G); Gj = None
                if not (lv > 0):
                    row["ctl"][k] = ("dead", k in req); line.append("%s%s DEAD(switch changes nothing)" % (k, "*" if k in req else "")); continue
                if k == "J": Tc = E.dual(Gj, wt, X)[0]; del Gj
                elif k == "NE5": Tc = E.dual(E.groups(nm, jd, wt, X, E.PHno2[nm]), wt, X)[0]
                else: Tc = E.dual(G, wt, X, variant=k)[0]
                d = abs(Td - Tc); rej = d > 10 * B12
                row["ctl"][k] = ("rejected" if rej else "NOT rejected", k in req)
                line.append("%s%s live(%.1e) |T-Tc|/|T|=%.1e /Btot12=%.1e %s" % (k, "*" if k in req else "", lv, d / abs(Td), d / B12, "rej" if rej else "NOT-rej"))
            print("      controls (* = required): " + " ; ".join(line))
        res.append(row); del G, r
        sys.stdout.flush()
    return res

W1 = [("B", {P7: 1}), ("B", {P7b: 1}), ("B", {P13: 1}), ("A", {}), ("B", {P7: 2}), ("F", {P7: 1}), ("A", {P7: 1})]
W2 = [("B", {P7: 1}), ("B", {P7b: 1}), ("B", {P13: 1}), ("A", {}), ("F", {P7: 1})]
KN = [("B", {P7: 1}), ("B", {P7: 0}), ("B", {P7: 4}), ("A", {})]
rows = []
if os.environ.get("STAGE12_SMOKE"):
    # smoke test of the machinery on the KNOWN (already seen) setting only; output kept in stage12_smoke.out
    run_frame("K-smoke", 4, 3.0, 49737.0, KN[:1], "KNOWN(post hoc, smoke test)")
    sys.exit(0)
rows += run_frame("W1", 4, 2.5, 82084.0, W1, "MAIN")
rows += run_frame("W2", 4, 2.0, 135335.0, W2, "MAIN")
rows += run_frame("K", 4, 3.0, 49737.0, KN, "KNOWN(post hoc)")

print("\n=== summary")
for r in rows:
    print("  %-16s %-3s %s %-14s |T|=%.3e |diff|=%.2e Btot12=%.2e digits %.2f %s %s ; controls: %s"
          % (r["frame"], r["label"], r["nm"], {q: j for q, j in r["jd"].items()}, r["T"], r["diff"], r["B"], math.log10(r["T"] / r["B"]),
             "informative" if r["inf"] else "not-informative", "within" if r["ok"] else "EXCEEDS",
             ", ".join("%s%s:%s" % (k, "*" if v[1] else "", v[0]) for k, v in r["ctl"].items()) or "-"))
def verdict(sel):
    main = [r for r in rows if r["frame"] == "MAIN" and sel(r)]
    inf_ = [r for r in main if r["inf"]]
    c2 = len(inf_) > 0
    c3 = all(r["ok"] for r in [r for r in rows if r["frame"] == "MAIN" and r["inf"]])
    c4 = all(v[0] == "rejected" for r in rows if r["frame"] == "MAIN" and r["inf"] for k, v in r["ctl"].items() if v[1])
    return G0 and c2 and c3 and c4, (G0, c2, c3, c4, len(inf_), len(main))
v1, d1 = verdict(lambda r: r["j1"]); v3, d3 = verdict(lambda r: r["nm"] == "A")
print("MAIN VERDICT j=1 : %s   (G0 %s, informative j=1 main runs %d of %d, all informative main runs within Btot12 %s, all required controls live and rejected %s)"
      % ("MET" if v1 else "NOT MET", d1[0], d1[4], d1[5], d1[2], d1[3]))
print("MAIN VERDICT 3|c : %s   (G0 %s, informative Psi_0=A main runs %d of %d, all informative main runs within Btot12 %s, all required controls live and rejected %s)"
      % ("MET" if v3 else "NOT MET", d3[0], d3[4], d3[5], d3[2], d3[3]))
print("[%.0fs]" % (time.time() - t0))
