# -*- coding: utf-8 -*-
"""
Stage 12, probe 2 (PLANNING ONLY for stage12_main2; written after stage12_main.out was seen).  E.run and E.dual are NOT
called: no dual sum, no comparison.  For the two weights planned in probe 1 but never compared (p=3, a=2.5 and p=5, a=2.5;
X = 82084) and the runs B {p7:1}, B {p13:1}, A {}:
  - per group: sigma, N(c), alpha(c), |alpha(c)^4 - 1|, k0 (inputs of NE2 / NE6);
  - the bound part of Btot12 with struct_tail2 (D evaluated at R (1 - 1e-9)), digits, the control tier it implies and
    its distance to the tier boundaries; the paper-only bound D_paper as a reference column;
  - the liveness (size of the change of its own ingredient) of every control NE1-NE7 and J.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage12_probe2.py > audits/003dual/stage12_probe2.out
"""
from stage5_lib import *
import stage5_bump as s5
import stage12_lib as L
import stage12_ctl as K
t0 = time.time()
print(__doc__)
E = EndToEnd()
P7, P13 = (1, 3), (4, 3)
RUNS = [("B", {P7: 1}), ("B", {P13: 1}), ("A", {})]
X = 82084.0
for (p, a) in ((3, 2.5), (5, 2.5)):
    bp = s5.Bump(a, p); bp.setup_quadrature(); wt = WBump(bp, 3.0e7); mu = L.MtUp(bp)
    assert X * math.exp(a) <= 1.0e6
    print("\n=== weight p=%d a=%g X=%d" % (p, a, X))
    for nm, jd in RUNS:
        Td, Sd = E.direct(nm, jd, wt, X)
        G = E.groups(nm, jd, wt, X)
        info = sorted(set((g["cfg"]["sigma"], g["key"][1], complex(np.round(g["ac"], 6)), round(abs(g["ac"] ** 4 - 1), 9), g["k0"]) for g in G), key=str)
        print("  %s %s : |T_dir| = %.4e ; distinct (sigma, N(c), alpha(c), |alpha^4-1|, k0): %s" % (nm, jd, abs(Td), info))
        st = sum(L.struct_tail2(mu, X, g["cfg"]["sigma"], g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"])[0] for g in G)
        sp = sum(L.struct_tail2(mu, X, g["cfg"]["sigma"], g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"], paper=True)[0] for g in G)
        dg = math.log10(abs(Td) / st)
        tier = "NE1,J(if j=1)" if dg < 3 else ("+NE2-NE5" if dg < 4 else "+NE6,NE7")
        print("      bound part of Btot12 = %.3e -> digits %.3f (tier: %s ; distance to nearest tier boundary %.2f digits) ; paper-only bound %.3e -> digits %.2f"
              % (st, dg, tier, min(abs(dg - 2), abs(dg - 3), abs(dg - 4)), sp, math.log10(abs(Td) / sp)))
        line = []
        for k in ("NE1", "J", "NE2", "NE3", "NE4", "NE5", "NE6", "NE7"):
            if k == "J":
                if not any(j == 1 for j in jd.values()): continue
                lv, Gj = K.liveness(E, "J", nm, jd, wt, X, G); del Gj
            else: lv = K.liveness(E, k, nm, jd, wt, X, G)
            line.append("%s %s(%.1e)" % (k, "live" if lv > 0 else "DEAD", lv))
        print("      liveness: " + " ; ".join(line))
        del G; sys.stdout.flush()
print("\nprobe 2 done [%.0fs]" % (time.time() - t0))
