"""Exact check of the overlap-rank fraction bound in
17-initial-list-gates-and-resampling-estimates.tex lines 196-207:
  fraction of slot-index choices with overlap rank j <= (2 n T_s)^{3j} L^{-j},
where Q (<= n T_s) occurrences each pick a uniform index in [L] (worst case: all in
one cell), one index (0) is globally pinned, and
  j = #occurrences - #distinct unpinned indices visited.
Also checks '#occurrences on repeated-or-pinned slots <= 2j'.
Exact Fractions, full enumeration."""
from fractions import Fraction
from itertools import product
from collections import Counter
worst = None
for Q in range(1, 6):
    for L in range(2, 8):
        cnt = Counter(); bad2j = 0
        for idx in product(range(L), repeat=Q):
            c = Counter(idx)
            distinct_unpinned = sum(1 for k in c if k != 0)
            j = Q - distinct_unpinned
            used = sum(m for k, m in c.items() if k == 0 or m >= 2)
            if used > 2 * j: bad2j += 1
            cnt[j] += 1
        tot = L ** Q
        for j, m in cnt.items():
            frac = Fraction(m, tot)
            bound = Fraction((2 * Q) ** (3 * j), L ** j)
            ratio = frac / bound
            if worst is None or ratio > worst[0]:
                worst = (ratio, Q, L, j)
            assert frac <= bound, (Q, L, j, frac, bound)
        assert bad2j == 0
print('bound holds in all cases; max frac/bound =', worst)
