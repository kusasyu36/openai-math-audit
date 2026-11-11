# -*- coding: utf-8 -*-
"""
Stage 5, note on the control NB3 of stage3_b_translate (delta' replaced by delta' + 1): it was not rejected in 3 of the
1557 well-conditioned configurations (registered requirement: >= 95%).  Explanation examined here: replacing delta' by
delta' + 1 shifts z' by -1/c; for c a unit and sigma = 0 this is a shift by a unit... only c = 1 gives the shift -1,
a period of theta, so the control cannot see it.  Counted exactly: configurations (all 324 h_0, r in {1,p7,p13,-5,p7*p13})
with c = 1 and sigma = 0; and, for these, whether e~(-l) = 1 on the first 5000 support points of d_0.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage5_nb3_check.py > audits/003dual/stage5_nb3_check.out 2>/dev/null
"""
from stage2_lib import *
print(__doc__)
th = Theta(200000); setup = Setup((18, 0))
n1 = []; nunit = 0
for primes in ([], [(1, 3)], [(4, 3)], [(-5, 0)], [(1, 3), (4, 3)]):
    for x0 in range(18):
        for x1 in range(18):
            cfg = setup.build((x0, x1), primes, [ONE] * len(primes))
            if norm(cfg["c"]) == 1:
                nunit += 1
                if cfg["sigma"] == "0": n1.append(((x0, x1), cfg["c"]))
m = th.d["0"].m[:5000]
dev = float(np.max(np.abs(AddChar(neg(ONE), LAM3)(m) - 1)))       # e~(-l) = e(-m/lambda^3)
print("configurations with N(c) = 1: %d ; of these with sigma = 0: %d : %s" % (nunit, len(n1), n1))
print("max |e~(-l) - 1| over the first 5000 support points of d_0: %.1e  (0 means the shift by -1 is invisible)" % dev)
