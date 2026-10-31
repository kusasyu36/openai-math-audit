"""
109 監査: 05-layers.tex §5.1 と 03-motifs.tex Prop complex-motif-interface の位相フレーム恒等式を厳密計算で確かめる。
 - C = H S H = aI + bX, H_0 = b S C S, C^2 = X, C^4 = I, C^{-1} = bI + aX = -iZCZ, C = iZC^{-1}Z
 - C_U = H^{⊗m} diag(i^{q_U(x)}) H^{⊗m} (q_U(x) = wt(P_U x) mod 4, P_U は F_2 上の射影)
 - U ⊂ V, 残差の正規直交基底 {v} に対し C_V C_U^{-1} = Π (aI + bX_v)^{±1}
 - X_a -> Y_a 経路: C_{u^perp} C_{<u>}^{-1} = supp(u) 上 C^{-1}, 他は C; Z 補正で C^{⊗m}
 - H_0^{⊗D} = b^D S_D C^{⊗D} S_D
 - 残差ベクトルの重みが奇数、i^{ε[v·x]} = (1+i^ε)/2 + (1-i^ε)/2 (-1)^{v·x}
実行: uv run --with sympy --with mpmath --with numpy python check_phase.py
"""
import itertools
import sympy as sp
from sympy import I, Rational, sqrt, Matrix, eye, zeros, kronecker_product as kp

ok_all = True
def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK' if cond else 'NG'}] {name} {detail}")

a = (1 + I) / 2
b = (1 - I) / 2
H = Matrix([[1, 1], [1, -1]]) / sqrt(2)
H0 = Matrix([[1, 1], [1, -1]]) / 2
S = sp.diag(1, I)
Z = sp.diag(1, -1)
X = Matrix([[0, 1], [1, 0]])
C = sp.simplify(H * S * H)
check("C = HSH = aI + bX", sp.simplify(C - (a * eye(2) + b * X)) == zeros(2, 2))
check("H_0 = b S C S", sp.simplify(H0 - b * S * C * S) == zeros(2, 2))
check("C^2 = X", sp.simplify(C * C - X) == zeros(2, 2))
check("C^4 = I", sp.simplify(C**4 - eye(2)) == zeros(2, 2))
Cinv = C.inv()
check("C^-1 = bI + aX", sp.simplify(Cinv - (b * eye(2) + a * X)) == zeros(2, 2))
check("C^-1 = -i Z C Z", sp.simplify(Cinv + I * Z * C * Z) == zeros(2, 2))
check("C = i Z C^-1 Z", sp.simplify(C - I * Z * Cinv * Z) == zeros(2, 2))
check("b^2 = -i/2", sp.simplify(b**2 + I / 2) == 0)
check("H Z H = X", sp.simplify(H * Z * H - X) == zeros(2, 2))

def tensor(ms):
    out = ms[0]
    for M in ms[1:]:
        out = kp(out, M)
    return out

def vecs(m):
    return list(itertools.product([0, 1], repeat=m))

def dot(x, y): return sum(xi * yi for xi, yi in zip(x, y)) % 2

def gf2_proj(basis, x):
    """projection onto span(basis) along its orthogonal complement over F_2; basis need not be orthonormal.
    Solve: y = sum c_i b_i with (x - y) . b_j = 0 for all j."""
    k = len(basis)
    if k == 0: return tuple([0] * len(x))
    G = Matrix([[dot(bi, bj) for bj in basis] for bi in basis])
    rhs = Matrix([dot(x, bj) for bj in basis])
    # solve G c = rhs over F_2 by brute force (k small)
    for cvec in itertools.product([0, 1], repeat=k):
        if all((sum(G[j, i] * cvec[i] for i in range(k)) - rhs[j]) % 2 == 0 for j in range(k)):
            y = [0] * len(x)
            for i in range(k):
                if cvec[i]:
                    y = [(yi + bi) % 2 for yi, bi in zip(y, basis[i])]
            return tuple(y)
    raise ValueError("degenerate")

def C_U(m, basis):
    Hm = tensor([H] * m)
    D = sp.diag(*[I ** (sum(gf2_proj(basis, x)) % 4) for x in vecs(m)])
    return sp.simplify(Hm * D * Hm)

def X_v(m, v):
    """translation (X_v f)(x) = f(x+v): matrix entries [x, x+v]"""
    vs = vecs(m); idx = {x: i for i, x in enumerate(vs)}
    M = zeros(2**m, 2**m)
    for x in vs:
        y = tuple((xi + vi) % 2 for xi, vi in zip(x, v))
        M[idx[x], idx[y]] = 1
    return M

def kernel(m, v, eps):
    K = a * eye(2**m) + b * X_v(m, v)
    return K if eps == 1 else K.inv()

m = 3
# --- C_F = C^{⊗m}, C_0 = I ---
check("C_0 = I", sp.simplify(C_U(m, []) - eye(2**m)) == zeros(2**m, 2**m))
check("C_F = C^{⊗m}", sp.simplify(C_U(m, [tuple(int(i == j) for i in range(m)) for j in range(m)]) - tensor([C] * m)) == zeros(2**m, 2**m))

# --- generic nested nondegenerate labels U ⊂ V with orthonormal residual basis ---
# U = <(1,1,1)> (norm 1), V = <(1,1,1),(1,1,0)+? > need orthonormal residual: pick v = (0,1,1)? (0,1,1).(1,1,1) = 0 but norm 0 -> alternating.
# choose u = (1,0,0), residual basis v1 = (0,1,0), v2 = (0,0,1): all orthonormal
u = (1, 0, 0); v1 = (0, 1, 0); v2 = (0, 0, 1)
CU = C_U(m, [u]); CV = C_U(m, [u, v1, v2])
prod = kernel(m, v1, 1) * kernel(m, v2, 1)
check("C_V C_U^-1 = C_{v1} C_{v2} (increasing edge, standard basis residual)", sp.simplify(CV * CU.inv() - prod) == zeros(8, 8))
check("decreasing edge C_U C_V^-1 = C_{v1}^-1 C_{v2}^-1", sp.simplify(CU * CV.inv() - kernel(m, v1, -1) * kernel(m, v2, -1)) == zeros(8, 8))
# non-standard orthonormal residual: U = <(1,1,1)>, V = <(1,1,1),(1,1,0)+(0,1,1)=(1,0,1)>? (1,0,1).(1,1,1)=0 norm 0 no.
# m=4: u=(1,1,1,0) norm 1; v=(0,0,0,1) orthonormal; w=(1,1,0,0)? norm 0. try v=(1,0,0,1): .u = 1 no. v = (1,1,1,1)? .u = 1 no.
# m=4: u=(1,0,0,0), v=(0,1,1,1) (weight 3, norm 1, u.v=0): check with weight-3 residual vector (coefficient 3 ≡ -1 mod 4)
m4 = 4
u4 = (1, 0, 0, 0); v4 = (0, 1, 1, 1)
CU4 = C_U(m4, [u4]); CV4 = C_U(m4, [u4, v4])
check("weight-3 residual vector: C_V C_U^-1 = C_v^{-1} (coefficient wt(v)=3 ≡ -1 mod 4)",
      sp.simplify(CV4 * CU4.inv() - kernel(m4, v4, -1)) == zeros(16, 16))
check("weight-3 residual: it is NOT C_v", sp.simplify(CV4 * CU4.inv() - kernel(m4, v4, 1)) != zeros(16, 16))
# phase identity i^{eps [v.x]} = (1+i^eps)/2 + (1-i^eps)/2 (-1)^{v.x}
for eps in (1, -1):
    for bit in (0, 1):
        lhs = I ** (eps * bit); rhs = (1 + I**eps) / 2 + (1 - I**eps) / 2 * (-1) ** bit
        check(f"phase identity eps={eps} bit={bit}", sp.simplify(lhs - rhs) == 0)

# --- X_a -> Y_a route: C_{u^perp} C_{<u>}^{-1} = C^{-1} on supp(u), C elsewhere; Z-corrections ---
m = 3
u = (1, 1, 1)   # weight 3 (odd), norm 1 -> models the 27-support vector u_a
basis_uperp = [(1, 1, 0), (0, 1, 1)]
Cu = C_U(m, [u]); Cup = C_U(m, basis_uperp)
route = sp.simplify(Cup * Cu.inv())
expected = tensor([Cinv if ui else C for ui in u])
check("C_{u^perp} C_u^{-1} = ⊗(C^{-1} on supp u, C elsewhere)", sp.simplify(route - expected) == zeros(8, 8))
Zu = tensor([Z if ui else eye(2) for ui in u])
corrected = sp.simplify(I ** sum(u) * Zu * route * Zu)
check("i^{wt u} Z_u (route) Z_u = C^{⊗m}", sp.simplify(corrected - tensor([C] * m)) == zeros(8, 8))
# with k columns: tensor power; sign i^{27k} applied once
k = 2
routek = tensor([route] * k); Zuk = tensor([Zu] * k)
check("k=2 columns: i^{wt(u) k} Z^{⊗k} route^{⊗k} Z^{⊗k} = C^{⊗mk}", sp.simplify(I ** (sum(u) * k) * Zuk * routek * Zuk - tensor([C] * (m * k))) == zeros(64, 64))
# q_{U^perp}(x) - q_U(x) = wt(x) + 2[u.x] mod 4
okq = True
for x in vecs(m):
    lhs = (sum(gf2_proj(basis_uperp, x)) - sum(gf2_proj([u], x))) % 4
    rhs = (sum(x) + 2 * dot(u, x)) % 4
    okq &= lhs == rhs
check("q_{u^perp}(x) - q_u(x) = wt(x) + 2[u.x] mod 4", okq)

# --- H_0^{⊗D} = b^D S_D C^{⊗D} S_D ---
for D in (1, 2, 3):
    SD = tensor([S] * D)
    check(f"H_0^{{⊗{D}}} = b^{D} S_D C^{{⊗{D}}} S_D", sp.simplify(tensor([H0] * D) - b**D * SD * tensor([C] * D) * SD) == zeros(2**D, 2**D))
# (C^{-1})^{⊗f} = (-i)^f Z^{⊗f} C^{⊗f} Z^{⊗f}
for f in (1, 2, 3):
    check(f"(C^-1)^{{⊗{f}}} = (-i)^{f} Z^{{⊗{f}}} C^{{⊗{f}}} Z^{{⊗{f}}}", sp.simplify(tensor([Cinv] * f) - (-I) ** f * tensor([Z] * f) * tensor([C] * f) * tensor([Z] * f)) == zeros(2**f, 2**f))
# conjugation P_M X_v P_M^-1 = X_{Mv}
Mmat = Matrix([[1, 1, 0], [0, 1, 0], [1, 0, 1]])  # in GL(3,F_2)
vs = vecs(3); idx = {x: i for i, x in enumerate(vs)}
PM = zeros(8, 8)
for x in vs:
    y = tuple(int(sum(Mmat[i, j] * x[j] for j in range(3)) % 2) for i in range(3))
    PM[idx[y], idx[x]] = 1   # moves value at x to Mx
v = (1, 1, 0)
Mv = tuple(int(sum(Mmat[i, j] * v[j] for j in range(3)) % 2) for i in range(3))
check("P_M X_v P_M^-1 = X_{Mv}", sp.simplify(PM * X_v(3, v) * PM.inv() - X_v(3, Mv)) == zeros(8, 8))

# --- entries of C_U are Gaussian integers / 2^m ---
CUtest = C_U(3, [(1, 0, 0), (0, 1, 1)]) if False else C_U(3, [(1, 0, 0), (0, 1, 0)])  # nondegenerate label
ints = all(sp.simplify(e * 8).is_integer or (sp.re(sp.simplify(e * 8)).is_integer and sp.im(sp.simplify(e * 8)).is_integer) for e in CUtest)
check("entries of C_U are in Z[i]/2^m", ints)

print("\nALL OK" if ok_all else "\nSOME CHECKS FAILED")
