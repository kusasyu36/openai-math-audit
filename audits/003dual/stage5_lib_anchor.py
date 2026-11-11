# -*- coding: utf-8 -*-
"""
Stage 5: anchor for stage5_lib.py (a copy of the formulas of stage3_d_mellin.py with the weight as a parameter).
Recomputes two runs of stage3_d_mellin.out and prints T_dir, T_dual to be compared digit by digit with that file:
  A {p7: 1}, w = 0.05, X = 64000   and   B {p7: 1}, w = 0.1, X = 2000.
Criterion (fixed before running): both numbers agree with the lines of stage3_d_mellin.out to all 13 printed digits.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage5_lib_anchor.py > audits/003dual/stage5_lib_anchor.out 2>/dev/null
"""
from stage5_lib import *
print(__doc__)
E = EndToEnd(); P7 = (1, 3)
for nm, w, X, xm in (("A", 0.05, 64000.0, 3.0e4), ("B", 0.1, 2000.0, 1.5e4)):
    r = E.run(nm, {P7: 1}, WLogGauss(w, xm), X)
    print("  Psi_0=%s P={7: 1} w=%.2f X=%g" % (nm, w, X))
    print("      T_dir = % .12e %+.12ei   T_dual = % .12e %+.12ei" % (r["Td"].real, r["Td"].imag, r["Tu"].real, r["Tu"].imag))
