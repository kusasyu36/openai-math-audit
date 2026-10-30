# Independent re-implementation of the Appendix A rank certificate (Lemma 3.1, first proof):
# cutoff 4450, iterative deletion of singly-supported columns, rank over F_1009 of the
# valuation matrix augmented by a row of ones. Expected: 1118/1100, 896/867, rank 868.
from collections import Counter
d = 3**11; K = 4450; P = 1009
LIM = 3*d + 1
spf = list(range(LIM))
i = 2
while i*i < LIM:
    if spf[i] == i:
        for j in range(i*i, LIM, i):
            if spf[j] == j: spf[j] = i
    i += 1
def fac(n):
    c = Counter()
    while n > 1:
        p = spf[n]; c[p] += 1; n //= p
    return c
cols = {}
for b in range(1, (d+1)//2):
    if b % 3 == 0: continue
    ok = all(max(fac(n), default=1) <= K for n in (d-b, b, 2*d+b, 2*d-b, 3*d-b, d+b))
    if not ok: continue
    col = Counter()
    for p, k in fac(d-b).items(): col[('A', p)] += k
    for p, k in fac(b).items(): col[('A', p)] -= k
    for n in (2*d+b, 2*d-b):
        for p, k in fac(n).items(): col[('R', p)] += k
    for n in (3*d-b, d+b):
        for p, k in fac(n).items(): col[('R', p)] -= k
    cols[b] = {r: v for r, v in col.items() if v}
def rows_of(cs):
    c = Counter()
    for col in cs.values():
        for r in col: c[r] += 1
    return c
print("after cutoff: columns", len(cols), "nonzero rows", len(rows_of(cols)))
while True:
    cnt = rows_of(cols)
    single = {r for r, v in cnt.items() if v == 1}
    if not single: break
    cols = {b: c for b, c in cols.items() if not (set(c) & single)}
cnt = rows_of(cols)
print("after deletions: columns", len(cols), "nonzero rows", len(cnt))
# rank over F_P of augmented matrix: rows = valuation rows + ones row; columns = b
rowkeys = sorted(cnt, key=str) + ['ONE']
idx = {r: i for i, r in enumerate(rowkeys)}
vecs = []
for b, c in cols.items():
    v = [0]*len(rowkeys)
    for r, x in c.items(): v[idx[r]] = x % P
    v[idx['ONE']] = 1
    vecs.append(v)
# rank of column set = rank of matrix
basis = {}
rank = 0
for v in vecs:
    v = v[:]
    for piv in sorted(basis):
        if v[piv]:
            f = v[piv]; bv = basis[piv]
            v = [(x - f*y) % P for x, y in zip(v, bv)]
    nz = next((i for i, x in enumerate(v) if x), None)
    if nz is None: continue
    inv = pow(v[nz], P-2, P)
    v = [x*inv % P for x in v]
    # reduce existing basis by new vector to keep it reduced
    for piv in list(basis):
        if basis[piv][nz]:
            f = basis[piv][nz]
            basis[piv] = [(x - f*y) % P for x, y in zip(basis[piv], v)]
    basis[nz] = v
    rank += 1
print("rank over F_1009 of augmented matrix:", rank)
