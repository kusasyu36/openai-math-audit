"""Independent exact re-implementation of the two finite certificates of the companion paper
"Exact Fourier certificates for complex Hadamard matrices of order six" (sections 02, 04, 05).

Written from the TeX specification (05-verification.tex).  Where the specification leaves a
convention open (canonical representatives, key order, pivot priority lists I), the SAME
conventions as the published verifier are adopted so that the three printed lines can be
compared number-for-number.  No code from verify.py is imported or executed.

Mathematics implemented:
  * injective moments m_H(A) as formal symbols; symmetries: row/col permutations, negation, transpose;
  * Lemma contraction (eq:injective-contraction) as the recursive reduction U_0 (eq:recursive-pair-reduction);
  * seven-row relations (eq:seven-row-contraction) enumerated per 05-verification.tex:139-175;
  * exact substitution lemma (Lemma exact-elimination) via modular echelon + CRT + rational
    reconstruction + exact integer re-check  W P = 0;
  * Lemma partitions (eq:verified-partition-expansion);
  * mixed moments, Lemma mixed-reductions (i)-(v), the finite mixed table (05-verification.tex:250-307);
  * Young-symmetrized forms (eq:Young-functions), Gram entries (Lemma gram), Schur positivity;
  * final constant / residual coefficient check (eq:exact-residual-coefficients).
"""
import itertools, sys, ast, time
from fractions import Fraction as Fr
from functools import lru_cache
from math import gcd
from collections import Counter

t_start = time.time()
def log(*a):
    print(f"[{time.time()-t_start:7.1f}s]", *a, flush=True)

# ---------------- basic charge / matrix utilities ----------------
def order(a): return sum(x for x in a if x > 0)
def negv(a): return tuple(-x for x in a)
def addv(a, b): return tuple(x + y for x, y in zip(a, b))
def transpose(A): return tuple(zip(*A)) if A else ()
def negm(A): return tuple(negv(r) for r in A)
def strip_sort_rows(A):
    """remove all-zero rows and sort the rows"""
    return tuple(sorted(r for r in A if any(r)))
def total_weight(A): return sum(order(r) for r in A)

def canonical_cols_then_rows(C):
    """For a matrix C (rows x cols), permute columns so that column orders are nondecreasing
    (all such permutations) and return (sorted column orders, min over them of sorted rows)."""
    cols = transpose(C)
    co = tuple(order(c) for c in cols)
    target = tuple(sorted(co))
    best = None
    for p in itertools.permutations(range(len(cols))):
        if tuple(co[i] for i in p) != target: continue
        M = tuple(sorted(transpose(tuple(cols[i] for i in p))))
        if best is None or M < best: best = M
    return target, best

@lru_cache(maxsize=None)
def canon_pair(A):
    """canonical representative of an exponent matrix under row perms, col perms, negation, transpose;
    zero rows and columns removed (05-verification.tex:111-120)."""
    A = strip_sort_rows(transpose(strip_sort_rows(transpose(A)))) if A else ()
    if not A: return ()
    cands = []
    for C in (A, negm(A), transpose(A), negm(transpose(A))):
        if len(C) >= len(C[0]):
            cands.append(canonical_cols_then_rows(C))
    return min(cands)[1]

# ---------------- sparse linear expressions ----------------
def acc(d, key, coef):
    v = d.get(key, 0) + coef
    if v: d[key] = v
    else: d.pop(key, None)
def acc_all(d, g, scale=1):
    for k, x in g.items(): acc(d, k, x * scale)
def apply_linear(g, F):
    r = {}
    for k, x in g.items(): acc_all(r, F(k), x)
    return r

# ---------------- Lemma contraction: U_0 reduction ----------------
def contractions_of_row(B, i):
    """unsigned sum over l != i of U_0( canonical( B with row i added into row l and row i deleted ) )"""
    v = {}
    a = B[i]; rest = B[:i] + B[i+1:]
    for j, b in enumerate(rest):
        M = rest[:j] + (addv(a, b),) + rest[j+1:]
        acc_all(v, U0(canon_pair(M)))
    return v

@lru_cache(maxsize=None)
def U0(A):
    """recursive reduction to moments all of whose nonzero rows/cols have order >= 2 (eq:recursive-pair-reduction)"""
    for B in (A, transpose(A)):
        assert len(B) <= 6 and all(sum(r) == 0 for r in B), ("invalid exponent matrix", A)
    for B in (A, transpose(A)):
        for i, r in enumerate(B):
            if order(r) == 1:
                v = {}
                acc_all(v, contractions_of_row(B, i), -Fr(1, 7 - len(B)))
                return v
    return {A: Fr(1)}

# ---------------- enumeration helpers ----------------
@lru_cache(maxsize=None)
def nondecreasing_parts(k, n, lo, hi):
    """nondecreasing sequences of n integers in [lo,hi] with sum k"""
    if n == 0: return ((),) if k == 0 else ()
    out = []
    for i in range(lo, min(k // n, hi) + 1):
        for rest in nondecreasing_parts(k - i, n - 1, i, hi): out.append((i,) + rest)
    return tuple(out)
@lru_cache(maxsize=None)
def capped_compositions(k, n, cap):
    """compositions of k into n nonnegative parts with part i <= cap[i]"""
    if n == 0: return ((),) if k == 0 else ()
    out = []
    lo = max(0, k - sum(cap[1:])); hi = min(k, cap[0])
    for i in range(lo, hi + 1):
        for rest in capped_compositions(k - i, n - 1, cap[1:]): out.append((i,) + rest)
    return tuple(out)
def balanced_vectors(k, n):
    """all balanced integer vectors of length n and order k"""
    for u in capped_compositions(k, n, (k,) * n):
        for v in capped_compositions(k, n, tuple(k if ui == 0 else 0 for ui in u)):
            yield addv(u, negv(v))

def seven_row_arrays(row_orders, col_orders, emit):
    """all 7 x q integer arrays with zero row/col sums, prescribed row orders (nondecreasing) and
    column orders; rows of equal order in lexicographic nondecreasing order; deduplicated by canon_pair."""
    q = len(col_orders); seen = set()
    def rec(i, pos_cap, neg_cap, prev, rows):
        if i == 7:
            C = tuple(rows)
            key = canon_pair(C)
            if key not in seen:
                seen.add(key); emit(C)
            return
        k = row_orders[i]; rest = sum(row_orders[i+1:])
        for u in capped_compositions(k, q, pos_cap):
            pc = addv(pos_cap, negv(u))
            for v in capped_compositions(k, q, tuple(c if ui == 0 else 0 for ui, c in zip(u, neg_cap))):
                a = addv(u, negv(v))
                if prev is not None and a < prev: continue
                nc = addv(neg_cap, negv(v))
                if max(addv(pc, nc)) > rest: continue
                rec(i + 1, pc, nc, a if (i + 1 < 7 and row_orders[i+1] == k) else None, rows + [a])
    rec(0, tuple(col_orders), tuple(col_orders), None, [])

# ---------------- exact substitution (Lemma exact-elimination) ----------------
def clear_denominators(row):
    d = 1
    for x in row.values(): d = d * x.denominator // gcd(d, x.denominator)
    return d, {k: int(x * d) for k, x in row.items()}

def modular_echelon(E, p):
    """E: list of dict col->int.  Pivot = largest column index of the reduced row.  Returns
    (rel: pivot -> normalized reduced row dict, witnesses: list of original rows giving pivots)."""
    rel = {}; wit = []
    for vec in E:
        s = {i: x % p for i, x in vec.items() if x % p}
        while s:
            k = max(s); b = s[k]
            if k not in rel:
                binv = pow(b, -1, p)
                rel[k] = {i: z * binv % p for i, z in s.items()}; wit.append(vec); break
            # subtract b * rel[k]
            for i, z in rel[k].items():
                h = (s.get(i, 0) - z * b) % p
                if h: s[i] = h
                else: s.pop(i, None)
    # back substitution so that each pivot row has no other pivot column
    for k in sorted(rel):
        s = rel[k]
        for i, z in list(s.items()):
            if i != k and i in rel:
                for j, y in rel[i].items():
                    h = (s.get(j, 0) - y * z) % p
                    if h: s[j] = h
                    else: s.pop(j, None)
    return rel, wit

def rational_reconstruct(z, M):
    """short rational with numerator*denominator small, congruent to z mod M (half-extended Euclid)"""
    a0, a1 = M, z; u0, u1 = 0, 1
    while a1 * a1 > M // 2:
        qq = a0 // a1
        a0, a1 = a1, a0 - qq * a1
        u0, u1 = u1, u0 - qq * u1
    return Fr(a1, u1)

def exact_substitutions(rows, keys):
    """returns dict key -> {key: Fraction} giving the certified projection P (rows of Lemma exact-elimination)
    and prints the rank table (nonzero equations, coordinates, independent witness rows)."""
    E0 = []
    for row in rows:
        g = clear_denominators(row)[1]
        if g: E0.append(g)
    num = {k: i for i, k in enumerate(keys)}
    assert len(num) == len(keys)
    E = [{num[k]: x for k, x in g.items()} for g in E0]
    p0 = 1000000007
    R, witnesses = modular_echelon(E, p0)
    M = p0; proj = None
    for q in (1000000009, 998244353, 1000003, 2013265921, 754974721, 32452843, 86028121):
        S, _ = modular_echelon(witnesses, q)
        if S.keys() != R.keys(): continue
        iq = pow(M, -1, q)
        for k, v in R.items():
            u = S[k]
            R[k] = {i: v.get(i, 0) + ((u.get(i, 0) - v.get(i, 0)) * iq % q) * M for i in v.keys() | u.keys()}
        M *= q
        proj = {}; scaled_rules = []
        ok = True
        for i, k in enumerate(keys):
            if i in R:
                f = {j: -rational_reconstruct(z, M) for j, z in R[i].items() if j != i}
                big = max((max(abs(x.numerator), x.denominator) for x in f.values()), default=0)
                if big * big > M // 200: ok = False; break
            else:
                f = {i: Fr(1)}
            proj[k] = {keys[j]: x for j, x in f.items()}
            scaled_rules.append(clear_denominators(f))
        if ok: break
        proj = None
    assert proj is not None and len(proj) == len(keys) and len(R) == len(witnesses)
    # exact check W P = 0 on every retained (witness) equation
    for g in witnesses:
        m = 1
        for i in g: m = m * scaled_rules[i][0] // gcd(m, scaled_rules[i][0])
        s = {}
        for i, x in g.items():
            d, f = scaled_rules[i]; acc_all(s, f, (m // d) * x)
        assert not s, "exact substitution check failed"
    print("rank table:", len(E0), len(keys), len(R), flush=True)
    return proj

# ---------------- pair table ----------------
def build_pair_table():
    rows = []
    def emit(A):
        assert order(A[0]) == 1 and len(A) == 7
        assert all(sum(r) == 0 for r in A) and all(sum(c) == 0 for c in transpose(A))
        rows.append(contractions_of_row(A, 0))          # (eq:pair-table-equations)
    for k in range(7, 12):
        for q in range(2, min(4, k // 2) + 1):
            for c in nondecreasing_parts(k, q, 2, k // 2):
                for r in nondecreasing_parts(k, 7, 1, k // 2):
                    seven_row_arrays(r, c, emit)
    log("pair equations generated:", len(rows))
    keys = sorted({k for g in rows for k in g}, key=lambda A: (total_weight(A), A))
    I = [0, 1, 2, 5, 4, 6, 7, 8, 9, 55, 80, 228]
    keys = [keys[i] for i in I] + [k for i, k in enumerate(keys) if i not in I]
    return exact_substitutions(rows, keys)

# ---------------- Lemma partitions: expansion of a product of characters ----------------
@lru_cache(maxsize=None)
def set_partitions(n):
    if n == 0: return ((),)
    return tuple(a + (j,) for a in set_partitions(n - 1) for j in range(max(a, default=-1) + 2))
def falling(n, k):
    r = 1
    for i in range(k): r *= (n - i)
    return r
PAIR_PROJ = None
@lru_cache(maxsize=None)
def expand_pair(A):
    """A: rows = coordinates, columns = charges.  Returns the projected linear expression
    sum_S (6)_|S| / 6^n U(A_S)  (eq:verified-partition-expansion)."""
    if not A:
        return apply_linear(U0(A), PAIR_PROJ.__getitem__)
    n = len(A[0]); v = {}
    for s in set_partitions(n):
        k = 1 + max(s)
        B = tuple(tuple(sum(x for i, x in zip(s, row) if i == j) for j in range(k)) for row in A)
        acc_all(v, U0(canon_pair(B)), Fr(falling(6, k), 6 ** n))
    return apply_linear(v, PAIR_PROJ.__getitem__)
def pair_of_charges(*charges):
    A = transpose(charges)                                 # rows = coordinates
    A = strip_sort_rows(transpose(strip_sort_rows(transpose(A)))) if A else ()
    return expand_pair(A)

ALPHA = (1, 1, 1, -1, -1, -1)
def check_pair_certificate():
    w = (-3, -2, 0, 1, 2, 2)
    gam = [(-1, -1, 0, 0, 1, 1), (-2, 0, 0, 0, 1, 1), (-1, -1, 0, 1, 0, 1), (-1, -1, 0, 1, 1, 0)]
    poly = {(w,): -1}
    for c, g in zip((2, 1, 1, 1), gam): poly[(g, addv(w, negv(g)))] = c
    f = {}
    for a, x in poly.items():
        for b, y in poly.items():
            acc_all(f, pair_of_charges(*(a + tuple(negv(t) for t in b))), 3 * x * y)
    results = {}
    for name, q in (("q_cross", (2, 0, 0, -2, 0, 0)), ("q_same", (2, -2, 0, 0, 0, 0))):
        eq = dict(f); acc_all(eq, pair_of_charges(ALPHA, q, negv(ALPHA), negv(q)))
        results[name] = (len(eq), max((abs(x) for x in eq.values()), default=0))
        log(f"pair certificate with {name}: nonzero coefficients after substitution = {len(eq)}")
    return results

# ---------------- mixed moments ----------------
def normalize_labels(c):
    new = {}
    return tuple(new.setdefault(s, len(new)) for s in c)
@lru_cache(maxsize=None)
def canon_mixed(c, A):
    """A: rows = coordinates (sorted, nonzero), columns = factors; c: labels per factor.
    Canonical form under simultaneous factor/label permutation, label renaming by first occurrence,
    coordinate-row sorting, global negation, removal of zero factors (Lemma mixed-reductions (i))."""
    B = transpose(A)                                       # factors
    keep = tuple(i for i, b in enumerate(B) if any(b))
    best = None
    for p in itertools.permutations(keep):
        T = transpose(tuple(B[i] for i in p))
        key = (normalize_labels(tuple(c[i] for i in p)), min(strip_sort_rows(T), strip_sort_rows(negm(T))))
        if best is None or key < best: best = key
    return best
def mixed_key(c, *charges):
    return canon_mixed(c, strip_sort_rows(transpose(charges)))
def is_alpha_type(a): return sorted(a) == sorted(ALPHA)
ZERO_PAIR_ROWS = {}
@lru_cache(maxsize=None)
def expand_mixed(K):
    """Lemma mixed-reductions: (ii) zero rules, (iv) single-label -> pair moments, (iii) singleton order-two label."""
    c, A = K
    assert len(A) <= 6 and all(sum(r) == 0 for r in A) and all(sum(col) == 0 for col in transpose(A))
    B = transpose(A)
    zero = any(order(b) == 1 or is_alpha_type(b) for b in B)
    types = sorted(set(c))
    if len(types) <= 1:
        v = expand_pair(A)
        v = {((), a): x for a, x in v.items()}
        if zero:
            ZERO_PAIR_ROWS[K] = v; return {}
        return v
    if zero: return {}
    for i, b in enumerate(B):
        if order(b) == 2 and c.count(c[i]) == 1:
            v = {}
            for t in types:
                if t != c[i]:
                    acc_all(v, expand_mixed(canon_mixed(c[:i] + (t,) + c[i+1:], A)), -Fr(1, 7 - len(types)))
            return v
    return {K: Fr(1)}

def shift(a, i, j): return tuple(x + (t == i) - (t == j) for t, x in enumerate(a))
def mixed_equations_for(A, rows, used, limit):
    def val(labels, charges):
        v = expand_mixed(mixed_key(labels, *charges)); used.update(v); return v
    for s in set_partitions(3): val(s, transpose(A))          # registers zero conditions / coordinates
    for a, b, c in sorted(set(itertools.permutations(transpose(A)))):
        if a > b: continue
        combs = {(i, j): (shift(a, i, j), shift(b, j, i), c) for i in range(6) for j in range(6)}
        if max(total_weight(t) for t in combs.values()) <= limit:
            label_sets = ((0, 1, 2), (0, 1, 0), (0, 1, 1)) if any(c) else ((0, 1, 0),)
            for s in label_sets:                                 # (eq:different-matrices) x third factor
                v = {}
                for (i, j), t in combs.items():
                    if i != j: acc_all(v, val(s, t))
                if v: rows.append(v)
        for s in (((0, 0, 0), (0, 0, 1)) if any(c) else ((0, 0, 0),)):   # (eq:same-matrix) x third factor
            for i in range(6):
                for reversed_shift in (0, 1):
                    comb = [combs[(j, i) if reversed_shift else (i, j)] for j in range(6)]
                    if max(total_weight(t) for t in comb) <= limit:
                        v = {}
                        for t in comb: acc_all(v, val(s, t))
                        acc_all(v, val(s[-2:], (negv(c), c)), -1)
                        if v: rows.append(v)
def charge_rep(a): return min(tuple(sorted(a)), tuple(sorted(negv(a))))
def charge_triples(limit):
    g = [list(balanced_vectors(k, 6)) for k in range(limit // 2 + 1)]
    for n, ga in enumerate(g):
        for a in sorted(set(map(charge_rep, ga))):
            s = set()
            for gb in g[:n + 1]:
                for b in gb:
                    c = negv(addv(a, b))
                    if order(c) <= order(b) and total_weight((a, b, c)) <= limit:
                        s.add(tuple(sorted(zip(a, b, c))))
            yield from sorted(s)
MIXED_I = [0, 1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 16, 20, 24, 25, 27, 31, 35, 36, 42, 44, 45, 46, 47, 57, 59,
        60, 61, 62, 63, 64, 66, 67, 75, 76, 77, 78, 82, 85, 92, 104, 105, 123, 125, 126, 127, 128, 129,
        132, 133, 141, 143, 144, 145, 147, 151, 153, 158, 159, 161, 162, 165, 166, 167, 172, 173, 174,
        193, 201, 203, 209, 219, 223, 226, 229, 231, 234, 240, 241, 249, 258, 259, 265, 281, 282, 293,
        304, 311, 316, 317, 323, 340, 341, 342, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 420,
        421, 422, 423, 425, 429, 442, 446, 447, 449, 450, 452, 457, 458, 461, 465, 467, 470, 472, 473,
        474, 476, 484, 499, 500, 501, 502, 503, 506, 507, 514, 524, 526, 528, 535, 539, 542, 545, 550,
        551, 552, 554, 555, 558, 559, 560, 561, 567, 575, 582, 583, 588, 591, 597, 598, 600, 601, 602,
        603, 604, 605, 606, 607, 608, 609, 610, 613, 616, 631, 632, 633, 634, 640, 649, 650, 654, 657,
        658, 660, 672, 674, 680, 685, 686, 693, 695, 696, 763, 764, 778, 780, 784, 795, 801, 811, 825,
        832, 840, 847, 870, 880, 893, 897, 900, 904, 924, 948, 949, 967, 1011, 1012, 1013, 1015, 1018,
        1019, 1021, 1022, 1026, 1028, 1031, 1042, 1060, 1075, 1078, 1080, 1083, 1087, 1105, 1120, 1122,
        1140, 1218, 1219, 1220, 1221, 1224, 1226, 1235, 1237, 1238, 1239, 1240, 1241, 1242, 1243, 1244,
        1246, 1247, 1248, 1250, 1251, 1252, 1253, 1254, 1255, 1258, 1259, 1260, 1261, 1263, 1268, 1271,
        1272, 1275, 1278, 1279, 1283, 1287, 1288, 1290, 1291, 1294, 1295, 1297, 1299, 1301, 1303, 1305,
        1309, 1311, 1318, 1319, 1320, 1323, 1326, 1327, 1330, 1338, 1339, 1341, 1344, 1346, 1350, 1352,
        1355, 1376, 1380, 1381, 1383, 1384, 1386, 1387, 1396, 1398, 1399, 1402, 1404, 1408, 1410, 1412,
        1417, 1421, 1423, 1426, 1427, 1430, 1431, 1432, 1441, 1444, 1446, 1447, 1448, 1452, 1458, 1459,
        1461, 1462, 1468, 1471, 1475, 1479, 1481, 1482, 1485, 1489, 1491, 1493, 1494, 1510, 1525, 1537,
        1539, 1542, 1549, 1554, 1565, 1567, 1573, 1575, 1577, 1585, 1588, 1591, 1593, 1595, 1597, 1599,
        1600, 1601, 1606, 1609, 1612, 1623, 1625, 1627, 1670, 1671, 1674, 1678, 1679, 1683, 1684, 1685,
        1687, 1693, 1694, 1695, 1696, 1702, 1706, 1715, 1728, 1731, 1738, 1745, 1749, 1751, 1754, 1756,
        1762, 1767, 1776, 1798, 1799, 1808, 1809, 1816, 1817, 1819, 1821, 1823, 1826, 1827, 1836, 1844,
        1845, 1849, 1859, 1868, 1882, 1892, 1895, 1896, 1899, 1901, 1908, 1909, 1921, 1946, 1960, 1961,
        1965, 1970, 1971, 1972, 1976, 1981, 1982, 1983, 1997, 2003, 2004, 2008, 2011, 2014, 2015, 2018,
        2019, 2024, 2025, 2039, 2040, 2076, 2077, 2096, 2099, 2100, 2105, 2107, 2108, 2109, 2110, 2183,
        2206, 2212, 2235, 2239, 2241, 2259, 2272, 2280, 2296, 2321, 2370, 2375, 2377, 2379, 2380, 2395,
        2396, 2403, 2479, 2483, 2508, 2608, 2610, 2655, 2745, 2831, 2832, 2833, 2834, 2841, 2845, 2849,
        2850, 2851, 2856, 2864, 2868, 2917, 2922, 2927, 3036, 3037, 3038, 3039, 3040, 3041, 3043, 3044,
        3045, 3046, 3047, 3048, 3049, 3050, 3051, 3052, 3053, 3054, 3055, 3056, 3057, 3058, 3059, 3061,
        3062, 3063, 3064, 3066, 3067, 3068, 3069, 3070, 3072, 3074, 3075, 3076, 3077, 3082, 3089, 3093,
        3095, 3096, 3098, 3099, 3101, 3102, 3103, 3104, 3105, 3106, 3107, 3109, 3110, 3111, 3112, 3113,
        3114, 3115, 3119, 3120, 3121, 3122, 3125, 3126, 3128, 3130, 3132, 3133, 3134, 3141, 3142, 3143,
        3144, 3145, 3147, 3148, 3149, 3150, 3152, 3154, 3158, 3160, 3162, 3163, 3165, 3168, 3169, 3170,
        3172, 3173, 3175, 3176, 3177, 3178, 3179, 3182, 3186, 3187, 3188, 3189, 3190, 3191, 3192, 3193,
        3194, 3195, 3197, 3198, 3199, 3200, 3201, 3204, 3211, 3212, 3221, 3231, 3232, 3233, 3235, 3236,
        3237, 3238, 3239, 3240, 3241, 3246, 3251, 3252, 3261, 3271, 3273, 3275, 3282, 3283, 3284, 3285,
        3292, 3293, 3295, 3296, 3299, 3300, 3301, 3303, 3311, 3316, 3317, 3318, 3320, 3321, 3325, 3327,
        3329, 3338, 3341, 3346, 3350, 3359, 3365, 3370, 3393, 3405, 3407, 3408]
def build_mixed_table(limit=11):
    rows = []; keys = set()
    ntr = 0
    for A in charge_triples(limit):
        mixed_equations_for(A, rows, keys, limit); ntr += 1
    log("charge triples:", ntr, " mixed equations:", len(rows), " zero-pair rows:", len(ZERO_PAIR_ROWS))
    for g in ZERO_PAIR_ROWS.values():
        keys.update(g); rows.append(g)
    keys = sorted(keys, key=lambda t: (bool(t[0]), total_weight(t[1]), t))
    keys = [keys[i] for i in MIXED_I] + [k for i, k in enumerate(keys) if i not in MIXED_I]
    return exact_substitutions(rows, keys)

# ---------------- Young-symmetrized forms and the Gram certificate ----------------
def block_permutations_with_sign(blocks):
    """all permutations of {0..5} fixing 0 that permute within each block, with sign"""
    for choice in itertools.product(*[itertools.permutations(b) for b in blocks]):
        p = list(range(6)); inv = 0
        for blk, img in zip(blocks, choice):
            for i, z in enumerate(img):
                p[blk[i]] = z; inv += sum(x > z for x in img[:i])
        yield tuple(p), (-1) ** inv
def young_form(lam, t, v):
    rows_ = []; h = 1
    for r in lam: rows_.append(tuple(range(h, h + r))); h += r
    cols_ = [tuple(s[j] for s in rows_ if len(s) > j) for j in range(max(lam))]
    f = {}
    for q, sgn in block_permutations_with_sign(cols_):
        for p, _ in block_permutations_with_sign(rows_):
            vq = tuple(v[i] for i in q); vqp = tuple(vq[i] for i in p)
            acc(f, (t, vqp), sgn)
    return f
def schur_positive(Z):
    B = [[Fr(x) for x in row] for row in Z]; pivots = []
    while B:
        q = B.pop(); p = q.pop()
        assert p > 40000, "pivot too small"
        pivots.append(p)
        B = [[x - y * a / p for x, a in zip(b, q)] for b, y in zip(B, q)]
    return pivots
def make_Z(m, a):
    assert len(a) * 2 == m * (m + 1)
    it = iter(a); Z = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1): Z[i][j] = Z[j][i] = next(it)
    return Z
@lru_cache(maxsize=None)
def gram_keys(X, Y):
    (s, a), (t, v) = X, Y
    b = negv(v); args = (a, b, negv(addv(a, b)))          # (a, -v, v-a)  Lemma gram
    K = mixed_key((s, t, 0), *args)
    if not (s and t): return K, None
    return K, mixed_key((1, 2, 0), *args)
def gram_certificate(data, MIXED_PROJ):
    E = max(e for lam, ws, e, a in data)
    h = Counter(); npiv = 0
    for lam, ws, e, a in data:
        assert sum(lam) == 5 and all(t in (0, 1) and len(v) == 6 and sum(v) == 0 for t, v in ws)
        forms = [young_form(lam, t, v) for t, v in ws]
        m = len(forms); Z = make_Z(m, a); npiv += len(schur_positive(Z))
        for i in range(m):
            for j in range(i + 1):
                z = Z[i][j] * (1 + (i != j)) * 10 ** (E - e)
                for X, p in forms[i].items():
                    for Y, q in forms[j].items():
                        K, L = gram_keys(X, Y); s = z * p * q
                        if L is None: h[K] += 5 * s
                        else: h[K] += s; h[L] += 4 * s
    log("positive pivots verified:", npiv)
    ans = apply_linear(apply_linear(h, expand_mixed), MIXED_PROJ.__getitem__)
    f = ans.pop(((), ()))
    v = sum(abs(y) for y in ans.values())
    print(f, v, 5 * 10 ** E, flush=True)
    return f, v

if __name__ == "__main__":
    src = open("data/openai_math/preprints/Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026/verification/verify.py").read()
    i = src.index("\ndata=["); j = src.index("\nEpair=pairtable()")
    data = ast.literal_eval(src[i + len("\ndata="):j].strip())     # integer data only, parsed as a literal
    log("building pair table")
    PAIR_PROJ = build_pair_table()
    res = check_pair_certificate()
    log("pair certificate results (nonzero coefficient counts):", res)
    assert all(n == 0 for n, _ in res.values()), "pair certificate FAILED"
    log("building mixed table")
    MIXED_PROJ = build_mixed_table(11)
    f, v = gram_certificate(data, MIXED_PROJ)
    log("constant f =", f, " residual sum =", v, " f+v =", f + v)
    print("6480 f =", 6480 * f, "(paper: -613302797399911)")
    print("6480 v =", 6480 * v, "(paper: 179721388988719)")
    print("f + v < 0:", f + v < 0, " f+v =", f + v, "(paper: -2007321335237/30)")
    log("done")
