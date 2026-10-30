"""Audit 201: numerical (mod p) test of Proposition 10.1 / Section 11.3 linear systems
in the simplest split case d = 1 (A = k), block size b = 2, m = 11 blocks.

Z = (Z1, Z2) truncated free series in s1, s2 with N x N matrix coefficients (N = m*b).
Constant term R = diag(P_i), linear terms = random arrows on offsets {1,2} (s1) and {3,4} (s2).
q = product of the characteristic polynomials of the diagonal blocks of w(R) (distinct roots generically).
At step l >= 3 unknowns: degree-(l-1) coefficients on reverse edges (offsets -1..-4) and
degree-l coefficients on diagonal blocks. Equation: [q(w(Z))]_omega = 0 for |omega| = l.
We test that the affine system is consistent (rank A == rank [A|b]), solve, and continue.
Negative control: drop the reverse-edge unknowns.

All arithmetic exact mod a prime p (stand-in for a generic char-0 point).
Usage: python -I check_prop10_1.py WCHOICE LMAX SEED
"""
import sys, itertools
import numpy as np

p = 1048573  # prime < 2^20
WCH = sys.argv[1] if len(sys.argv) > 1 else "comm"
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 1
CONTROL = len(sys.argv) > 4 and sys.argv[4] == "noreverse"
rng = np.random.default_rng(SEED)
m, b = 11, 2
N = m * b

words = [()]
for l in range(1, LMAX + 1):
    words += list(itertools.product((1, 2), repeat=l))
widx = {w: k for k, w in enumerate(words)}
NW = len(words)
pairs_u, pairs_v, pairs_o = [], [], []
for u in words:
    for v in words:
        if len(u) + len(v) <= LMAX:
            pairs_u.append(widx[u]); pairs_v.append(widx[v]); pairs_o.append(widx[u + v])
pairs_u = np.array(pairs_u); pairs_v = np.array(pairs_v); pairs_o = np.array(pairs_o)

def smul(A, B):
    prod = np.matmul(A[pairs_u], B[pairs_v]) % p
    out = np.zeros_like(A)
    np.add.at(out, pairs_o, prod)
    return out % p

def sconst(M):
    S = np.zeros((NW, N, N), dtype=np.int64); S[0] = M % p; return S

def blk(i, j):
    return slice(b * (i % m), b * (i % m) + b), slice(b * (j % m), b * (j % m) + b)

def rand(*shape):
    return rng.integers(0, p, size=shape, dtype=np.int64)

# constant term and linear arrows
Z = [np.zeros((NW, N, N), dtype=np.int64) for _ in range(2)]
for c in range(2):
    for i in range(m):
        Z[c][0][blk(i, i)] = rand(b, b)
for letter, offs in ((1, (1, 2)), (2, (3, 4))):
    for r in offs:
        for j in range(m):
            i = (j + r) % m
            for c in range(2):
                Z[c][widx[(letter,)]][blk(i, j)] = rand(b, b)

def wser(Z1, Z2):
    if WCH == "comm":   # w = z1 z2 - z2 z1   (trace identically 0 -> g = 1 type obstruction)
        return (smul(Z1, Z2) - smul(Z2, Z1)) % p
    if WCH == "gen":    # w = z1 z2 z1 + z2^2 + z1
        return (smul(smul(Z1, Z2), Z1) + smul(Z2, Z2) + Z1) % p
    if WCH == "comm2":  # w = [z1,z2] + [z1,[z1,z2]]  (trace 0)
        C = (smul(Z1, Z2) - smul(Z2, Z1)) % p
        return (C + smul(Z1, C) - smul(C, Z1)) % p
    raise ValueError

def polymul(a, bb):
    out = [0] * (len(a) + len(bb) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(bb):
            out[i + j] = (out[i + j] + x * y) % p
    return out

W0 = wser(Z[0], Z[1])[0]
q = [1]  # coefficients low->high
for i in range(m):
    B = W0[blk(i, i)]
    tr = int(B[0, 0] + B[1, 1]) % p
    det = int(B[0, 0] * B[1, 1] - B[0, 1] * B[1, 0]) % p
    q = polymul(q, [det, (-tr) % p, 1])

def qW(Z1, Z2):
    W = wser(Z1, Z2)
    acc = sconst(np.eye(N, dtype=np.int64) * q[-1])
    for coef in reversed(q[:-1]):
        acc = (smul(acc, W) + sconst(np.eye(N, dtype=np.int64) * coef)) % p
    return acc

def degree_slice(l):
    return [widx[w] for w in words if len(w) == l]

def rank_mod(M):
    M = M.copy() % p
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i; break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        nz = np.nonzero(M[:, c])[0]
        for i in nz:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows:
            break
    return r, M

def solve_mod(A, rhs):
    aug = np.concatenate([A, rhs[:, None]], axis=1) % p
    r_aug, R = rank_mod(aug)
    r_A, _ = rank_mod(A)
    if r_A != r_aug:
        return r_A, r_aug, None
    x = np.zeros(A.shape[1], dtype=np.int64)
    for i in range(r_aug):
        nzc = np.nonzero(R[i, :-1])[0]
        if len(nzc) == 0:
            continue
        x[nzc[0]] = R[i, -1]  # free vars = 0
    return r_A, r_aug, x

# check degrees 1, 2 vanish for the initial substitution
QW = qW(*Z)
for l in (0, 1, 2):
    nz = int(np.count_nonzero(QW[degree_slice(l)]))
    print(f"[{WCH} seed={SEED}] initial: nonzero entries of q(w(Z)) in degree {l}: {nz}")

for l in range(3, LMAX + 1):
    unk = []  # (component, word index, block (row,col), entry)
    if not CONTROL:
        for r in (1, 2, 3, 4):
            for i in range(m):
                j = (i - r) % m  # reverse edge j <- i
                for w in words:
                    if len(w) == l - 1:
                        for c in range(2):
                            for e in range(b * b):
                                unk.append((c, widx[w], j, i, e))
    for i in range(m):
        for w in words:
            if len(w) == l:
                for c in range(2):
                    for e in range(b * b):
                        unk.append((c, widx[w], i, i, e))
    sl = degree_slice(l)
    base = qW(*Z)[sl].reshape(-1)
    nproj = 700
    Rp = rand(nproj, base.size)
    cols = []
    for (c, wi, bi, bj, e) in unk:
        Zp = [Z[0].copy(), Z[1].copy()]
        rs, cs = blk(bi, bj)
        Zp[c][wi][rs.start + e // b, cs.start + e % b] = 1
        col = (qW(*Zp)[sl].reshape(-1) - base) % p
        cols.append(col)
    Afull = np.stack(cols, axis=1)
    A = (Rp @ Afull) % p
    rhs = (-(Rp @ base)) % p
    rA, rAug, x = solve_mod(A, rhs)
    rFull, _ = rank_mod(Afull[np.random.default_rng(SEED+7).permutation(Afull.shape[0])[:0]] if False else A)
    print(f"[{WCH} seed={SEED}{' CONTROL-no-reverse' if CONTROL else ''}] l={l}: unknowns={len(unk)}, "
          f"base error nonzeros={int(np.count_nonzero(base))}, rank A={rA}, rank [A|b]={rAug}, consistent={x is not None}")
    if x is None:
        break
    for k, (c, wi, bi, bj, e) in enumerate(unk):
        if x[k]:
            rs, cs = blk(bi, bj)
            Z[c][wi][rs.start + e // b, cs.start + e % b] = (Z[c][wi][rs.start + e // b, cs.start + e % b] + x[k]) % p
    QW = qW(*Z)
    print(f"   after correction: nonzeros of q(w(Z)) in degrees <= {l}:",
          [int(np.count_nonzero(QW[degree_slice(ll)])) for ll in range(l + 1)])
