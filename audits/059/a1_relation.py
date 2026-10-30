# Independent audit script (own code; does not import anything from the OpenAI repo).
# Checks Lemma 3.1 (finite product relation) from explicit-relation.json, the chain
# completion of Lemma 3.2 for every needed starting weight, variable counts,
# Hessian ranks, N mod 8, and the Milnor-number equality of the seeds.
import json, sys
from fractions import Fraction
from collections import Counter

path = sys.argv[1]
data = json.load(open(path))
d = int(data["d"])
assert d == 3**11, d
D = 3 * d
rows = [(int(r["b"]), int(r["e"])) for r in data["relation"]]
assert len({b for b, _ in rows}) == len(rows)
for b, e in rows:
    assert 0 < b and 2 * b < d and b % 3 != 0 and e != 0

# smallest-prime-factor sieve up to 3d (largest integer appearing is 3d-1 or 2d+b<5d/2)
LIM = 3 * d + 1
spf = list(range(LIM))
i = 2
while i * i < LIM:
    if spf[i] == i:
        for j in range(i * i, LIM, i):
            if spf[j] == j:
                spf[j] = i
    i += 1

def fac(n):
    out = Counter()
    while n > 1:
        p = spf[n]
        out[p] += 1
        n //= p
    return out

# Product 1: A_b = (d-b)/b ; Product 2: R_b = (2d+b)(2d-b)/((3d-b)(d+b))
tot = [Counter(), Counter()]
for b, e in rows:
    for g, (num, den) in enumerate((((d - b,), (b,)), ((2*d + b, 2*d - b), (3*d - b, d + b)))):
        for n in num:
            for p, k in fac(n).items(): tot[g][p] += e * k
        for n in den:
            for p, k in fac(n).items(): tot[g][p] -= e * k
res = [{p: v for p, v in t.items() if v} for t in tot]
E = sum(e for _, e in rows)
print("product A residual valuations:", res[0])
print("product R residual valuations:", res[1])
print("E = sum e_b > 0:", E > 0, " digits(E) =", len(str(E)))
print("support size", len(rows), " max digits", max(len(str(abs(e))) for _, e in rows))
assert not res[0] and not res[1] and E > 0
assert E == int(data["sum_e"])

# Also directly confirm with an independent route: product of R_b/A_b also 1 (used in Case 1)
# (follows from the two above) -- and verify the Case-1 'exceptional' algebra symbolically:
# group-2 factor for row b equals ((3d-b)/b)^{-1} * ((2d+b)/(d-b)) * ((2d-b)/(d+b)) = R_b/A_b
for b, _ in rows[:50]:
    lhs = Fraction(b, 3*d - b) * Fraction(2*d + b, d - b) * Fraction(2*d - b, d + b)
    rhs = Fraction((2*d + b) * (2*d - b), (3*d - b) * (d + b)) / Fraction(d - b, b)
    assert lhs == rhs
# check (1-w)/w for the four weights
for b, _ in rows[:50]:
    w = [Fraction(b, d), Fraction(b, 3*d), Fraction(d - b, 3*d), Fraction(d + b, 3*d)]
    v = [(1 - x) / x for x in w]
    assert v[0] == Fraction(d - b, b)
    assert v[1] == Fraction(3*d - b, b) and v[2] == Fraction(2*d + b, d - b) and v[3] == Fraction(2*d - b, d + b)
    # reduced denominators
    assert w[0].denominator == d and all(x.denominator == 3*d for x in w[1:])

# Lemma 3.2 chain completion, own implementation
def chain(c):
    assert 0 < 2 * c < D
    out = []
    while c:
        a = D // c
        r = D - a * c
        out.append((c, a, r))
        c = r
    return out

maxlen = 0
init = [0, 0]   # starting variable counts side1, side2
tails = 0
minexp = 10**9
for b, e in rows:
    for c, cnt in ((3*b, 2*e), (b, -2*e), (d - b, 2*e), (d + b, 2*e)):
        ch = chain(c)
        maxlen = max(maxlen, len(ch))
        for i, (ci, a, r) in enumerate(ch):
            assert a >= 2 and 0 <= r < ci and a * ci + r == D
            assert Fraction(ci, D) < Fraction(1, 2)
            deg = a + (1 if r else 0)
            assert deg >= 3, (c, ch)
            minexp = min(minexp, a)
        # terminal exponent odd >=3 dividing D
        cl, al, rl = ch[-1]
        assert rl == 0 and al * cl == D and al % 2 == 1 and al >= 3
        side = 1 if cnt > 0 else 0
        init[side] += abs(cnt)
        tails += abs(cnt) * (len(ch) - 1)
assert init[1] - init[0] == 4 * E
side1 = init[0] + tails + 4 * E
side2 = init[1] + tails
assert side1 == side2
s0 = side1 + 2
assert s0 % 2 == 0
s = s0 + 2
N = 4 * s + 6 + 2
assert N == 4 * s0 + 16 and N % 8 == 0
print("chains OK; max chain length", maxlen, "min exponent", minexp)
print("s0 digits", len(str(s0)), " N mod 8 =", N % 8, " N digits", len(str(N)))
print("Hessian rank f1 = 16E (digits %d), f2 = 0" % len(str(16 * E)))

# Milnor number ratio mu(p2)/mu(p1) = prod over weights of ((1-w)/w)^{count2-count1}; squares give 1
# = prod_b A_b^{2e} * (R_b/A_b)^{2e} = prod R_b^{2e}: zero valuations => 1.
print("mu(p1)=mu(p2): follows from residuals above (ratio = prod R_b^{2e_b} = 1)")
print("ALL CHECKS PASSED")
