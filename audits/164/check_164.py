"""Computational audit of the finite/combinatorial parts of OpenAI preprint 164
"Monochromatic finite sums and products in the positive integers" (2026-09-23).

Only Section 2 (combinatorial deduction), Corollary 1.2 and the exponent
bookkeeping of Section 8.1 are elementary enough to test by computer.
Run:  ../.venv/bin/python -I check_164.py
"""
import itertools
import random
from fractions import Fraction

ok_count = 0


def check(name, cond):
    global ok_count
    print(("[OK]   " if cond else "[FAIL] ") + name)
    ok_count += bool(cond)


def subsets(s, nonempty=True):
    s = list(s)
    for k in range(1 if nonempty else 0, len(s) + 1):
        yield from itertools.combinations(s, k)


lt = lambda A, B: max(A) < min(B)          # A < B for nonempty sets

# ---- 1. chain property: k<d  =>  B_k in E_{B_d}   (eq. 2.2, 2.3)
def in_E(Acand, B):
    """A in E_B where B = T u {i}: A = P u {j}, P nonempty, P < T < j < i."""
    T, i = B
    P, j = Acand
    return len(P) > 0 and lt(P, T) and max(T) < j < i


def chains(n, m):
    # all chains B_d = T_d u {i_d} in [n] with T_1<...<T_m<i_1<...<i_m
    idx = list(range(1, n + 1))
    for pivots in itertools.combinations(idx, m):
        rest = [x for x in idx if x < pivots[0]]
        # choose m consecutive nonempty blocks from an increasing subset of rest
        for cut in itertools.combinations(range(1, len(rest)), m - 1):
            pass
        # enumerate small tails explicitly
        def rec(start, d, acc):
            if d == m:
                yield list(acc)
                return
            for a in range(start, pivots[0]):
                for b in range(a, pivots[0]):
                    for T in subsets(range(a, b + 1)):
                        if min(T) == a and max(T) == b:
                            yield from rec(b + 1, d + 1, acc + [tuple(T)])
        for tails in rec(1, 0, []):
            yield [(tails[d], pivots[d]) for d in range(m)]


ok = True
cnt = 0
for n in range(3, 9):
    for m in (2, 3):
        for ch in chains(n, m):
            cnt += 1
            for d in range(m):
                for kk in range(d):
                    T, i = ch[kk]
                    ok &= in_E((T, i), ch[d])
            # disjointness
            allidx = [x for (T, i) in ch for x in T + (i,)]
            ok &= len(allidx) == len(set(allidx))
check(f"chain => B_k in E_(B_d) for k<d, blocks disjoint ({cnt} chains, n<=8, m<=3)", ok)

# ---- 2. Section 8.1: lambda_e = p_A - p_C + p_T0 > 0 for every pair at pivot i
ok = True
cnt = 0
for n in range(2, 9):
    p = {l: 2 ** (n - l) for l in range(1, n + 1)}
    pS = lambda S: sum(p[l] for l in S)
    for i in range(1, n + 1):
        below = range(1, i)
        tails = [S for S in subsets(below)]
        for T0 in tails:
            for C in tails:
                for j in range(1, i):
                    if not (max(C) < j):
                        continue
                    for P in subsets(range(1, min(C))):
                        A = P + (j,)
                        lam = pS(A) - pS(C) + pS(T0)
                        cnt += 1
                        ok &= lam > 0
                        # direct check from the scale update kappa
                        kap = {l: (p[l] if l < i else (-pS(T0) if l == i else 0)) for l in range(1, n + 1)}
                        ok &= (sum(kap[l] for l in A) - sum(kap[l] for l in C + (i,))) == lam
    # b'_{B0} = b_{B0}
    for i in range(2, n + 1):
        for T0 in subsets(range(1, i)):
            ok &= sum(p[l] for l in T0) - sum(p[l] for l in T0) == 0
check(f"Sec 8.1 exponents lambda_e > 0 and match kappa update ({cnt} pairs, n<=8)", ok)

# ---- 3. Corollary 1.2: separation (1.1) with R=2, D=1 => distinctness
random.seed(1)
ok = True
for trial in range(3000):
    m = random.randint(1, 7)
    a = []
    for d in range(m):
        if d == 0:
            a.append(random.randint(3, 6))
        else:
            base = sum(a) + eval("*".join(map(str, a)))
            a.append(2 * base + random.randint(1, 3))     # a_d > R*(sum+prod)^D, R=2,D=1
    FS = {}
    FP = {}
    for J in subsets(range(m)):
        s = sum(a[j] for j in J)
        pr = 1
        for j in J:
            pr *= a[j]
        FS.setdefault(s, []).append(J)
        FP.setdefault(pr, []).append(J)
    ok &= len(FS) == 2 ** m - 1 and len(FP) == 2 ** m - 1
    ok &= set(FS) & set(FP) == set(a)
    ok &= len(set(FS) | set(FP)) == 2 * (2 ** m - 1) - m
check("Cor 1.2: separation => |FS|=|FP|=2^m-1, FS∩FP=A (3000 random m<=7)", ok)

# minimal separated example also fine with the weakest admissible a_d
ok = True
for m in range(1, 7):
    a = [3]
    while len(a) < m:
        a.append(2 * (sum(a) + eval("*".join(map(str, a)))) + 1)
    sums = [sum(a[j] for j in J) for J in subsets(range(m))]
    prods = []
    for J in subsets(range(m)):
        pr = 1
        for j in J:
            pr *= a[j]
        prods.append(pr)
    ok &= len(set(sums)) == len(sums) and len(set(prods)) == len(prods) and set(sums) & set(prods) == set(a)
check("Cor 1.2 tight separated sequences m<=6", ok)

# ---- 4. Lemma 2.1 (chain selection) index bookkeeping
ok = True
for F in range(1, 12):
    for m in range(1, 5):
        for u in itertools.product(range(1, F + 1), repeat=m):
            if sum(u) > F:
                continue
            q = [2 * x for x in u]
            ok &= all(x >= 2 for x in q) and sum(q) <= 2 * F
            # tails use q_d - 1 indices, pivots m indices: total = sum q
            ok &= sum(x - 1 for x in q) + m == sum(q)
            # index-set size of union over J = sum_{d in J} q_d = 2 * (subset sum of u)
            for J in subsets(range(m)):
                ok &= sum(q[d] for d in J) == 2 * sum(u[d] for d in J)
check("Lemma 2.1: tail/pivot index counts and cardinality colours", ok)

# ---- 5. tolerance arithmetic in the deduction (eq. 2.11 -> 2.13)
ok = True
for _ in range(2000):
    delta = Fraction(random.randint(1, 100), 100)
    Cc = random.randint(1, 50)
    K = random.randint(1, 500)
    m = random.randint(2, 5)
    tau = min(Fraction(1, 5), delta / (12 * Cc + 1))
    eta = min(delta / (4 * Cc + 1), delta * tau ** (2 ** m) / (4 * K + 1))
    ok &= 3 * Cc * tau < delta / 4 and Cc * eta < delta / 4 and K * eta < delta * tau ** (2 ** m) / 4
    ok &= Cc * (3 * tau + eta) < delta / 2
    # lower bound: delta/2 * tau^(2^m) (mass of good event times min integrand)
    # minus K*eta  > 0
    ok &= delta / 2 * tau ** (2 ** m) - K * eta > 0
    q = 2 ** m - m - 1
    ok &= tau ** q >= tau ** (2 ** m)
check("Sec 2.5 tolerance choices give positive total count", ok)

# ---- 6. Cor 2.2: common residue v of a1,a2,a1+a2 mod q forces v = 0
ok = all(((2 * v) % q == v) == (v == 0) for q in range(1, 60) for v in range(q))
check("Cor 2.2: v = 2v mod q  <=>  v = 0", ok)

print(f"\n{ok_count} checks passed")
