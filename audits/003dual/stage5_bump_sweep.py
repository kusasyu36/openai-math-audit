# -*- coding: utf-8 -*-
"""
Stage 5 (a), second planning table (bounds only, no comparison), written after stage5_a_compact.out showed that with
(p, a) = (4, 3) the one-prime runs B {p7: j} are guaranteed to only 1.5-2.2 digits.  Question: is there a (p, a) in the
family of stage5_bump.py whose (T3') bound for B {p7: 1} (N(c) = 28 and 7, cusp tables, R = 111.11^2) is at least 10 times
smaller relative to X-independent scale?  Printed: the bound per unit |C| Bmax at the largest admissible X = 0.999e6 e^-a.
Decision rule (fixed before running): rerun the one-prime runs with another (p, a) only if its bound for N(c) = 28 is
below 1.2e-3 (10 times better than (4, 3): 1.2e-2).  Otherwise nothing is rerun and the digits of stage5_a_compact stand.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage5_bump_sweep.py > audits/003dual/stage5_bump_sweep.out 2>/dev/null
"""
import math, sys
import stage5_bump as s5
print(__doc__)
R = 111.11 ** 2
for p in (3, 4, 5, 6, 8):
    for a in (1.0, 1.5, 2.0, 2.5, 3.0):
        bp = s5.Bump(a, p); X = math.floor(1e6 * math.exp(-a) * 0.999)
        print("  p = %d  a = %.1f  X = %6d : N(c) = 28: %.2e   N(c) = 7: %.2e   N(c) = 84: %.2e   N(c) = 36 (theta table): %.2e"
              % (p, a, X, bp.dual_tail_bound(X, 28, R, 1, 1), bp.dual_tail_bound(X, 7, R, 1, 1), bp.dual_tail_bound(X, 84, R, 1, 1), bp.dual_tail_bound(X, 36, 192.45 ** 2, 1, 1)))
        sys.stdout.flush()
