"""Audit 201: Prop 10.1-type solvability with a nontrivial coefficient algebra (split form).
Coefficient factor: 2x2 matrices I, J with I^2 = a, J^2 = 1, IJ = -JI (split quaternion basis).
w = I(z1 z2) + J z1 + I(z2 z1 z1) + J z2 : all coefficients in span{I, J}, so conjugation by IJ
sends W -> -W and eigenvalues come in linked pairs (lambda, -lambda) (nontrivial linkage partition).
Input factor: m = 11 blocks of size b = 1 (inputs are scalars per block), output size 2*11.
W(Z) = sum_omega a_omega (x) Z^omega.  Same unknown pattern and test as check_prop10_1.py.
Usage: python -I check_linkage_d2.py LMAX SEED [noreverse]
"""
import sys, itertools
import numpy as np

p = 1048573
LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
CONTROL = len(sys.argv) > 3 and sys.argv[3] == "noreverse"
rng = np.random.default_rng(SEED)
m = 11
aq = 7
I2 = np.array([[0, 1], [aq, 0]], dtype=np.int64)
J2 = np.array([[1, 0], [0, p - 1]], dtype=np.int64)
assert np.all((I2 @ J2 + J2 @ I2) % p == 0)
ONE = np.eye(2, dtype=np.int64)

words = [()]
for l in range(1, LMAX + 1):
    words += list(itertools.product((1, 2), repeat=l))
widx = {w: k for k, w in enumerate(words)}
NW = len(words)
pu, pv, po = [], [], []
for u in words:
    for v in words:
        if len(u) + len(v) <= LMAX:
            pu.append(widx[u]); pv.append(widx[v]); po.append(widx[u + v])
pu, pv, po = map(np.array, (pu, pv, po))

def smul(A, B):
    prod = np.matmul(A[pu], B[pv]) % p
    out = np.zeros((NW,) + A.shape[1:], dtype=np.int64)
    np.add.at(out, po, prod)
    return out % p

def rand(*s):
    return rng.integers(0, p, size=s, dtype=np.int64)

Z = [np.zeros((NW, m, m), dtype=np.int64) for _ in range(2)]
for c in range(2):
    for i in range(m):
        Z[c][0][i, i] = rand()
for letter, offs in ((1, (1, 2)), (2, (3, 4))):
    for r in offs:
        for j in range(m):
            i = (j + r) % m
            for c in range(2):
                Z[c][widx[(letter,)]][i, j] = rand()

def kron_series(a, S):
    return np.stack([np.kron(a, S[k]) % p for k in range(NW)])

def wser(Z1, Z2):
    t1 = kron_series(I2, smul(Z1, Z2))
    t2 = kron_series(J2, Z1)
    t3 = kron_series(I2, smul(smul(Z2, Z1), Z1))
    t4 = kron_series(J2, Z2)
    return (t1 + t2 + t3 + t4) % p

N = 2 * m
W0 = wser(*Z)[0]
def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return out
q = [1]
for i in range(m):
    idx = [i, m + i]
    B = W0[np.ix_(idx, idx)]
    # off-block entries must vanish at the constant term
    tr = int(B[0, 0] + B[1, 1]) % p
    det = int(B[0, 0] * B[1, 1] - B[0, 1] * B[1, 0]) % p
    q = polymul(q, [det, (-tr) % p, 1])
print("block traces at constant term (expect 0, symmetric spectrum):",
      sorted({int((W0[i, i] + W0[m + i, m + i]) % p) for i in range(m)}))

def qW(Z1, Z2):
    W = wser(Z1, Z2)
    acc = np.zeros((NW, N, N), dtype=np.int64); acc[0] = np.eye(N, dtype=np.int64) * q[-1]
    for coef in reversed(q[:-1]):
        acc = smul(acc, W)
        acc[0] = (acc[0] + np.eye(N, dtype=np.int64) * coef) % p
    return acc

def dsl(l):
    return [widx[w] for w in words if len(w) == l]

def rank_mod(M):
    M = M.copy() % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        nzr = np.nonzero(M[r:, c])[0]
        if len(nzr) == 0:
            continue
        piv = r + nzr[0]
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows:
            break
    return r, M

QW = qW(*Z)
print("initial nonzeros by degree 0..2:", [int(np.count_nonzero(QW[dsl(l)])) for l in range(3)])
for l in range(3, LMAX + 1):
    unk = []
    if not CONTROL:
        for r in (1, 2, 3, 4):
            for i in range(m):
                j = (i - r) % m
                for w in words:
                    if len(w) == l - 1:
                        for c in range(2):
                            unk.append((c, widx[w], j, i))
    for i in range(m):
        for w in words:
            if len(w) == l:
                for c in range(2):
                    unk.append((c, widx[w], i, i))
    sl = dsl(l)
    base = qW(*Z)[sl].reshape(-1)
    cols = []
    for (c, wi, a, bb) in unk:
        Zp = [Z[0].copy(), Z[1].copy()]
        Zp[c][wi][a, bb] = (Zp[c][wi][a, bb] + 1) % p
        cols.append((qW(*Zp)[sl].reshape(-1) - base) % p)
    Af = np.stack(cols, axis=1)
    Rp = rand(600, base.size)
    A = (Rp @ Af) % p
    rhs = (-(Rp @ base)) % p
    rA, _ = rank_mod(A)
    rAug, Red = rank_mod(np.concatenate([A, rhs[:, None]], axis=1))
    ok = rA == rAug
    print(f"[d=2 linkage seed={SEED}{' CONTROL-no-reverse' if CONTROL else ''}] l={l}: unknowns={len(unk)}, "
          f"rank A={rA}, rank[A|b]={rAug}, consistent={ok}  (#eigen-equations={22*len(sl)})")
    if not ok:
        break
    x = np.zeros(len(unk), dtype=np.int64)
    for i in range(rAug):
        nz = np.nonzero(Red[i, :-1])[0]
        if len(nz):
            x[nz[0]] = Red[i, -1]
    for k, (c, wi, a, bb) in enumerate(unk):
        Z[c][wi][a, bb] = (Z[c][wi][a, bb] + x[k]) % p
    QW = qW(*Z)
    print("   after correction nonzeros by degree:", [int(np.count_nonzero(QW[dsl(ll)])) for ll in range(l + 1)])
