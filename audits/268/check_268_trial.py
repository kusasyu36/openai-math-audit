"""Exact-integer recomputation of the periodic matrix-product trial energies
(trials.tex, Lemma trial-contraction, eq. trial-integer-results) for audit 268.

The 26x26 integer matrices A_s are rebuilt from the tex entry rules alone
(trials.tex:72-137), then
  V_N = Tr(T^N),  U_N = Tr(T^{N-2} O),  T = sum_s A_s (x) A_s,
are computed with Python integers via the difference blocks d-e, for N = 72, 120.
Also: block dimensions, comparison with the JSON matrices, and a brute-force
check of Lemma trial-contraction for N = 4, 6 (explicit psi_N and explicit H_N).
"""
import sys, json, time
import numpy as np
from fractions import Fraction as Fr

ell = [1, 1, 1, 1, 3, 3, 3, 5]
Dpar = {1: -289, 2: -84, 3: -9, 4: 33, 5: -22, 6: 27, 7: 71, 8: -11}
Ppar = {(1,5): -11, (1,6): -35, (1,7): -28, (2,5): -23, (2,6): 16, (2,7): -49, (3,5): -6, (3,6): -25,
        (3,7): -20, (4,5): -34, (4,6): 9, (4,7): -32, (5,8): 1, (6,8): 2, (7,8): -16}
labels = [(al, d) for al in range(1, 9) for d in range(-ell[al-1], ell[al-1]+1, 2)]
assert len(labels) == 26
idx = {lab: i for i, lab in enumerate(labels)}

def entry(al, d, be, e, s):
    l, k = ell[al-1], ell[be-1]
    if e != d + 2*s or abs(k-l) > 2 or (k == l and al != be):
        return 0
    r = (d + l)//2
    if k == l + 2:
        p = Ppar[(al, be)]
        fac = {-1: 2*(l+1-r)*(l+2-r), 0: 2*(r+1)*(l-r+1), 1: (r+1)*(r+2)}[s]
    elif k == l:
        p = Dpar[al]
        fac = {-1: 2*(l-r+1), 0: d, 1: -(r+1)}[s]
    else:
        p = -Ppar[(be, al)]
        fac = {-1: 2, 0: -2, 1: 1}[s]
    return p*fac

A = {}
for s in (-1, 0, 1):
    M = np.zeros((26, 26), dtype=object)
    for (al, d) in labels:
        for (be, e) in labels:
            M[idx[(al, d)], idx[(be, e)]] = entry(al, d, be, e, s)
    A[s] = M

# compare with the paper's JSON (read only)
jpath = sys.argv[1] if len(sys.argv) > 1 else None
if jpath:
    with open(jpath) as fh:
        J = json.load(fh)
    jl = [tuple(x) if len(x)==2 else tuple(x) for x in J['labels']]
    print("JSON labels (first 5):", J['labels'][:5], " tex order (first 5):", labels[:5])
    same = all(np.array_equal(np.array(J['matrices'][str(s)], dtype=object), A[s]) for s in (-1, 0, 1))
    print("A_s rebuilt from tex rules == matrix-data.json:", same)
    if not same:
        for s in (-1,0,1):
            Mj = np.array(J['matrices'][str(s)], dtype=object)
            diff = [(labels[i], labels[j], A[s][i,j], Mj[i,j]) for i in range(26) for j in range(26) if A[s][i,j] != Mj[i,j]]
            print(f"  s={s}: {len(diff)} differing entries; first:", diff[:5])

def kron(X, Y):
    return np.kron(X, Y)  # object dtype OK
T = sum(kron(A[s], A[s]) for s in (-1, 0, 1))
Bm = {(s, t): A[s].dot(A[t]) for s in (-1,0,1) for t in (-1,0,1)}
O = sum(s*t*kron(Bm[(s,t)], Bm[(s,t)]) for s in (-1,0,1) for t in (-1,0,1))
for s in (-1, 0):
    for t in (0, 1):
        O = O + kron(Bm[(s,t)], Bm[(s+1,t-1)]) + kron(Bm[(s+1,t-1)], Bm[(s,t)])

# difference blocks
pairs = [(i, j) for i in range(26) for j in range(26)]
delta = {pr: labels[pr[0]][1] - labels[pr[1]][1] for pr in pairs}
blocks = {}
for n, pr in enumerate(pairs):
    blocks.setdefault(delta[pr], []).append(n)
dims = {dl: len(v) for dl, v in sorted(blocks.items())}
print("block dimensions by d-e:", dims)
expected = {-10:1,-8:8,-6:32,-4:80,-2:136,0:162,2:136,4:80,6:32,8:8,10:1}
print("match eq:trial-blocks:", dims == expected)
# check T and O are block diagonal
Tnz = [(pairs[i], pairs[j]) for i in range(676) for j in range(676) if T[i,j] != 0 and delta[pairs[i]] != delta[pairs[j]]]
Onz = [(pairs[i], pairs[j]) for i in range(676) for j in range(676) if O[i,j] != 0 and delta[pairs[i]] != delta[pairs[j]]]
print("off-block nonzeros in T, O:", len(Tnz), len(Onz))

def matpow(M, n):
    R = None; P = M
    while n:
        if n & 1:
            R = P if R is None else R.dot(P)
        n >>= 1
        if n: P = P.dot(P)
    return R

def contractions(N):
    U = 0; V = 0
    for dl, ids in blocks.items():
        Td = T[np.ix_(ids, ids)]; Od = O[np.ix_(ids, ids)]
        P = matpow(Td, N-2)
        U += int(np.sum(P * Od.T)); V += int(np.sum(P * Td.dot(Td).T))
    return U, V

# brute-force check of Lemma trial-contraction for small N
Sz = np.diag([1.0, 0.0, -1.0]); Sp = np.zeros((3,3)); Sp[0,1] = Sp[1,2] = np.sqrt(2); Sm = Sp.T
def Hper(N):
    import functools
    def op(o, j):
        mats = [np.eye(3)]*N; mats[j] = o
        return functools.reduce(np.kron, mats)
    H = np.zeros((3**N, 3**N))
    for j in range(N):
        k = (j+1) % N
        H += op(Sz,j)@op(Sz,k) + 0.5*(op(Sp,j)@op(Sm,k) + op(Sm,j)@op(Sp,k))
    return H
Af = {s: A[s].astype(float) for s in A}
for N in (4, 6):
    U, V = contractions(N)
    # explicit psi: basis index order m=+1,0,-1 per site (matches Sz = diag(1,0,-1))
    psi = np.zeros(3**N)
    import itertools
    for conf in itertools.product((1, 0, -1), repeat=N):
        M = np.eye(26)
        for s in conf: M = M @ Af[s]
        ind = 0
        for s in conf: ind = 3*ind + (1 - s)
        psi[ind] = np.trace(M)
    H = Hper(N)
    nrm = psi @ psi; en = psi @ H @ psi
    E0 = np.linalg.eigvalsh(H)[0]
    print(f"N={N}: V_N={V} vs |psi|^2={nrm:.6e} (rel diff {abs(V-nrm)/nrm:.1e});  N U_N={N*U} vs <psi|H|psi>={en:.6e} (rel diff {abs(N*U-en)/abs(en):.1e});"
          f"  energy/site = {N*U/V/N:.6f} >= E0/N = {E0/N:.6f}: {N*U/V >= E0 - 1e-9}")

a = Fr(700741, 500000)
cN = {72: 1401482112601, 120: 1401482105174}
for N in (72, 120):
    t0 = time.time()
    U, V = contractions(N)
    ok1 = V > 0
    ok2 = -cN[N]*V <= 10**12*U < (-cN[N]+1)*V
    ratio = Fr(U, V)
    print(f"\nN={N}: V>0: {ok1};  -c_N V <= 1e12 U < (-c_N+1) V: {ok2};  U/V = {float(ratio):.15f};  "
          f"U/V < -a: {ratio < -a};  margin per site -a - U/V = {float(-a-ratio):.3e};  shifted energy bound N(U/V + a) = {float(N*(ratio+a)):.3e}  [{time.time()-t0:.0f}s]")
    print(f"   digits: V has {len(str(V))} digits, U has {len(str(U))} digits")
    if jpath:
        import os
        cert = os.path.join(os.path.dirname(jpath), f"certificate{N}.json")
        with open(cert) as fh: c = json.load(fh)
        from math import gcd
        if 'norm' in c:
            cV, cU = int(c['norm']), int(c['one_bond_energy_numerator']); cfl = int(c['scaled_floor'])
        else:
            cV, cU = int(c['V']), int(c['U']); cfl = int(c['floor_1e12_U_over_V'])
        g = gcd(U, V)
        print(f"   JSON certificate{N}: raw integers match: V {cV == V}, U {cU == U};  reduced ratio U/V matches: {Fr(cU, cV) == ratio};"
              f"  floor(1e12 U/V) match: {cfl == (10**12*U)//V};  gcd(U,V) has {len(str(g))} digits, JSON V has {len(str(cV))} digits, own V {len(str(V))} digits")

# MPS energy per bond as a function of N (float, scaled), to see the N -> infinity limit
Tf = T.astype(float); Of = O.astype(float)
w = np.linalg.eigvals(Tf); lam = np.max(np.abs(w))
Ts = Tf/lam; Os = Of/lam**2
print("\nMPS energy per bond e(N) = U_N/V_N (float, normalized transfer matrix):")
for N in (8, 12, 16, 24, 36, 48, 72, 120, 200, 400):
    P = np.linalg.matrix_power(Ts, N-2)
    print(f"   N={N:4d}: e(N) = {np.trace(P@Os)/np.trace(P@Ts@Ts):.12f}")
print(f"   -a = {float(-a):.12f};  White-Huse e_0(inf) = -1.401484039")
