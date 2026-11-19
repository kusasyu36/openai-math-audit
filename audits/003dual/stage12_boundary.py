# -*- coding: utf-8 -*-
"""
Stage 12: size of the boundary gap of stage12_main.out (written after stage12_main.out was seen; bounds only, no dual sum).
For every run of the W1/W2 frames and the K frame of stage12_main.py: the structured tail with D evaluated at R (struct_tail,
as used in stage12_main) and at R (1 - 1e-9) (struct_tail2, which also counts every l with N(l) = R).  The difference bounds
what a term with N(l) = R exactly, missed by both the table sum and the bound, can contribute.  Also the number of table
entries with |N(l) - R| <= 1e-9 R in the groups (where the floating-point comparison decides).
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage12_boundary.py > audits/003dual/stage12_boundary.out
"""
from stage5_lib import *
import stage5_bump as s5
import stage12_lib as L
print(__doc__)
E = EndToEnd(); th = E.th
P7, P7b, P13 = (1, 3), (-2, -3), (4, 3)
FR = [("W1", 4, 2.5, 82084.0, [("B", {P7: 1}), ("B", {P7b: 1}), ("B", {P13: 1}), ("A", {}), ("B", {P7: 2}), ("F", {P7: 1}), ("A", {P7: 1})]),
      ("W2", 4, 2.0, 135335.0, [("B", {P7: 1}), ("B", {P7b: 1}), ("B", {P13: 1}), ("A", {}), ("F", {P7: 1})]),
      ("K", 4, 3.0, 49737.0, [("B", {P7: 1}), ("B", {P7: 0}), ("B", {P7: 4}), ("A", {})])]
for label, p, a, X, runs in FR:
    bp = s5.Bump(a, p); bp.setup_quadrature(); wt = WBump(bp, 3.0e7); mu = L.MtUp(bp)
    for nm, jd in runs:
        G = E.groups(nm, jd, wt, X)
        t1 = t2 = 0.0; nb = 0
        for g in G:
            args = (X, g["cfg"]["sigma"], g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"])
            t1 += L.struct_tail(mu, *args)[0]; t2 += L.struct_tail2(mu, *args)[0]
            N = th.d[g["key"][0]].aell ** 2; nb += int(np.sum(np.abs(N - g["R"]) <= 1e-9 * g["R"]))
        print("  %s %s %-14s struct_tail %.6e  struct_tail2 %.6e  difference %.2e (relative %.1e) ; table entries within 1e-9 R of R over the groups: %d"
              % (label, nm, jd, t1, t2, t2 - t1, (t2 - t1) / t1, nb))
        del G; sys.stdout.flush()
