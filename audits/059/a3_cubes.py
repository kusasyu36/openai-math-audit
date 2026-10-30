# Independent audit script: Lemma 5.3 (six-cube factor), signature table,
# plus elementary checks of the Phi_d(1) identity and the symmetrization algebra.
import cmath, math, itertools
from fractions import Fraction

def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def tr(a): return [list(r) for r in zip(*a)]
def kron(a, b):
    return [[a[i][j]*b[k][l] for j in range(len(a[0])) for l in range(len(b[0]))] for i in range(len(a)) for k in range(len(b))]
def inv2(a):
    det = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[Fraction(a[1][1], det), Fraction(-a[0][1], det)], [Fraction(-a[1][0], det), Fraction(a[0][0], det)]]
def neg(a): return [[-x for x in r] for r in a]

A = [[1, -1], [0, 1]]
Sc = neg(A)
t = neg(mm(inv2(A), tr(A)))           # -A^{-1}A^t
assert t == [[0, -1], [1, -1]], t
# also equals (-1)^1 S^{-1} S^t with S=-A
t2 = neg(mm(inv2(Sc), tr(Sc)))
assert t2 == t
I2 = [[1, 0], [0, 1]]
assert mm(mm(t, t), t) == I2
J = mm(A, mm(t, t))
assert J == [[0, 1], [-1, 0]], J
assert mm(mm(tr(t), J), t) == J      # t-invariant
print("cube: t, t^3=1, J=A t^2, t^T J t = J : OK")

# six-fold
def kpow(m, n):
    r = m
    for _ in range(n-1): r = kron(r, m)
    return r
T6 = kpow(t, 6); J6 = kpow(J, 6); A6 = kpow(A, 6)
Su = neg(A6)                            # paper: S_u = (-1)^15 (-A)^{x6} = -A^{x6}
# check symmetric form B_u = S_u T6^2 equals -J6 and is symmetric, T6-invariant
T6sq = mm(T6, T6)
Bu = mm(Su, T6sq)
assert Bu == neg(J6)
assert Bu == tr(Bu)
assert mm(mm(tr(T6), Bu), T6) == Bu
# and T6 = S_u^{-1} S_u^t  <=>  S_u T6 = S_u^t
assert mm(Su, T6) == tr(Su)
print("six cubes: B_u = S_u T^2 = -J^{x6}, symmetric, invariant, and S_u T = S_u^t : OK")

# signatures of x -> J6(x, conj y) on eigenspaces of T6 (numerically, via projectors)
w = cmath.exp(2j*math.pi/3)
n = 64
T6c = [[complex(x) for x in r] for r in T6]
def mpow(m, k):
    r = [[1.0+0j if i == j else 0j for j in range(n)] for i in range(n)]
    for _ in range(k): r = mm(r, m)
    return r
Tp = [mpow(T6c, k) for k in range(3)]
def inertia_hermitian(G):
    # Gaussian elimination congruence (Sylvester): count signs of pivots w/ pivoting
    G = [r[:] for r in G]; m = len(G); pos = negc = zero = 0
    idx = list(range(m))
    while idx:
        # pick largest |diag|
        k = max(idx, key=lambda i: abs(G[i][i]))
        if abs(G[k][k]) < 1e-9:
            # try to create a nonzero diagonal by adding a row/col
            found = False
            for i in idx:
                for j in idx:
                    if i != j and abs(G[i][j]) > 1e-9:
                        # replace basis vector i by e_i + c e_j to get nonzero diag
                        c = 1 if (G[i][j]).real != 0 or True else 1
                        for cc in (1, 1j):
                            val = G[i][i] + G[j][j]*abs(cc)**2 + 2*(cc.conjugate()*G[i][j] if isinstance(cc, complex) else G[i][j]).real
                            if abs(val) > 1e-9:
                                for r in range(m): G[r][i] += G[r][j]*(cc.conjugate() if isinstance(cc, complex) else cc)
                                for r in range(m): G[i][r] += G[j][r]*cc
                                found = True; break
                        if found: break
                if found: break
            if not found:
                zero += len(idx); break
            continue
        p = G[k][k].real
        if p > 0: pos += 1
        else: negc += 1
        for i in idx:
            if i != k:
                f = G[i][k]/G[k][k]
                for j in idx: G[i][j] -= f*G[k][j]
        idx.remove(k)
        for i in idx: G[i][k] = 0; G[k][i] = 0
    return pos, negc, zero

for e, lam in ((0, 1), (1, w), (2, w*w)):
    # projector onto lambda-eigenspace: (1/3) sum lam^{-k} T^k
    P = [[sum(lam**(-k)*Tp[k][i][j] for k in range(3))/3 for j in range(n)] for i in range(n)]
    # column basis by greedy rank selection
    cols = []
    basis = []
    for j in range(n):
        v = [P[i][j] for i in range(n)]
        # Gram-Schmidt
        u = v[:]
        for b in basis:
            c = sum(b[i].conjugate()*u[i] for i in range(n))
            u = [u[i]-c*b[i] for i in range(n)]
        nr = math.sqrt(sum(abs(x)**2 for x in u))
        if nr > 1e-8:
            basis.append([x/nr for x in u]); cols.append(v)
    dim = len(cols)
    # hermitian form h(x,y) = J6(x, conj y) = x^T J6 conj(y)
    G = [[sum(cols[a][i]*J6[i][j]*cols[b][j].conjugate() for i in range(n) for j in range(n) if J6[i][j]) for b in range(dim)] for a in range(dim)]
    # hermitian check
    assert all(abs(G[a][b]-G[b][a].conjugate()) < 1e-8 for a in range(dim) for b in range(dim))
    print(f"eigenvalue omega^{e}: dim {dim}, inertia of J^(x6)(x, conj y) =", inertia_hermitian(G))

# Phi_d(1) identity: compute cyclotomic polys by division for d up to 300
def polydiv(a, b):
    a = a[:]; q = [0]*(len(a)-len(b)+1)
    for i in range(len(a)-len(b), -1, -1):
        c = a[i+len(b)-1] // b[-1]; q[i] = c
        for j, y in enumerate(b): a[i+j] -= c*y
    assert not any(a); return q
phi = {}
for dd in range(1, 450):
    p = [-1] + [0]*(dd-1) + [1]
    for e2 in range(1, dd):
        if dd % e2 == 0: p = polydiv(p, phi[e2])
    phi[dd] = p
def pf(x):
    s = set(); q = 2
    while q*q <= x:
        while x % q == 0: s.add(q); x //= q
        q += 1
    if x > 1: s.add(x)
    return s
for dd in range(2, 450):
    v = sum(phi[dd]); ps = pf(dd)
    assert v == (list(ps)[0] if len(ps) == 1 else 1)
assert sum(phi[143]) == 1 and sum(phi[429]) == 1
print("Phi_d(1) = p for prime powers, 1 otherwise (d<450): OK")
