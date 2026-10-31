"""
109 監査: 04-swap.tex の Lemma lower-lower（下三角 x 部分置換 x 下三角の分解）と Lemma matrix-shear の
3 ステップ更新、Prop power-interchange の再帰 (eq:swap-recurrence) と Lemma chunk-swap の
行数条件・幾何和を自作実装で確かめる。
実行: uv run --with sympy --with mpmath --with numpy python check_shear.py
"""
import random
from fractions import Fraction as Fr
from math import log, ceil, comb
import sympy as sp

ok_all = True
def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK' if cond else 'NG'}] {name} {detail}")

def lower_lower(A):
    """Lemma lower-lower の消去をそのまま実装: A = E1 Π E2, E1,E2 下三角可逆, Π 部分置換 (rank(A) 個の 1)"""
    A = sp.Matrix(A); n = A.rows
    Lm = sp.eye(n); Rm = sp.eye(n)   # maintain  L A R = current
    cur = A.copy()
    active_rows = list(range(n)); active_cols = list(range(n))
    while True:
        # topmost nonzero active row
        piv_row = None
        for i in active_rows:
            if any(cur[i, j] != 0 for j in active_cols):
                piv_row = i; break
        if piv_row is None: break
        piv_col = max(j for j in active_cols if cur[piv_row, j] != 0)   # rightmost nonzero entry
        # clear entries to the left in pivot row by right-multiplying lower triangular elementary matrices
        # column op: col_j <- col_j - (cur[piv_row,j]/pivot) col_{piv_col}, j < piv_col  => matrix E with E[piv_col, j] = -c: lower triangular (row index > col index)
        for j in active_cols:
            if j < piv_col and cur[piv_row, j] != 0:
                cfac = cur[piv_row, j] / cur[piv_row, piv_col]
                E = sp.eye(n); E[piv_col, j] = -cfac
                cur = cur * E; Rm = Rm * E
        # clear entries below the pivot by left-multiplying lower triangular elementary matrices
        for i in active_rows:
            if i > piv_row and cur[i, piv_col] != 0:
                cfac = cur[i, piv_col] / cur[piv_row, piv_col]
                E = sp.eye(n); E[i, piv_row] = -cfac   # row_i <- row_i - c row_piv : lower triangular since i > piv_row
                cur = E * cur; Lm = E * Lm
        # scale pivot to one
        E = sp.eye(n); E[piv_row, piv_row] = 1 / cur[piv_row, piv_col]
        cur = E * cur; Lm = E * Lm
        active_rows.remove(piv_row); active_cols.remove(piv_col)
    Pi = cur
    return Lm.inv(), Pi, Rm.inv()

def is_lower(M): return all(M[i, j] == 0 for i in range(M.rows) for j in range(M.cols) if j > i)

rng = random.Random(7)
allok = True
for trial in range(40):
    n = rng.choice([3, 4, 5, 6])
    rank = rng.randint(0, n)
    # random rational matrix of given rank
    U = sp.Matrix(n, rank, lambda i, j: sp.Rational(rng.randint(-3, 3), rng.randint(1, 3)))
    V = sp.Matrix(rank, n, lambda i, j: sp.Rational(rng.randint(-3, 3), rng.randint(1, 3)))
    A = U * V if rank else sp.zeros(n, n)
    E1, Pi, E2 = lower_lower(A)
    ones = sum(1 for i in range(n) for j in range(n) if Pi[i, j] == 1)
    zeros_ok = all(Pi[i, j] in (0, 1) for i in range(n) for j in range(n))
    rowsum = all(sum(Pi[i, :]) <= 1 for i in range(n)); colsum = all(sum(Pi[:, j]) <= 1 for j in range(n))
    good = (E1 * Pi * E2 == A) and is_lower(E1) and is_lower(E2) and E1.det() != 0 and E2.det() != 0 and zeros_ok and rowsum and colsum and ones == A.rank()
    allok &= good
    if not good:
        print("   failure:", A, Pi)
check("Lemma lower-lower: A = E1 Π E2 with lower-triangular invertible E1,E2 and rank(A) ones (40 random rational matrices)", allok)

# Lemma matrix-shear: pivot (i,j): D_j <- D_j + H_i ; swap(H_i, D_j) ; D_j <- H_i - D_j  gives (h+d, d)
okp = True
for h in range(-5, 6):
    for d in range(-5, 6):
        H, D = h, d
        D = D + H; H, D = D, H; D = H - D
        okp &= (H, D) == (h + d, d)
check("pivot update (D+=H; swap; D=H-D) yields (h+d, d)", okp)
# full shear: D <- E2 D ; H <- E1^{-1} H ; H <- H + Π D ; H <- E1 H ; D <- E2^{-1} D  == (H + A D, D)  mod q^b
q, bdig = 5, 3; Q = q**bdig
n = 4
A = sp.Matrix(n, n, lambda i, j: sp.Rational(rng.randint(-3, 3), rng.choice([1, 2, 3])))
E1, Pi, E2 = lower_lower(A)
def modmat(M):
    return sp.Matrix(M.rows, M.cols, lambda i, j: (int(M[i, j].p) * pow(int(M[i, j].q), -1, Q)) % Q)
# q must avoid denominators and numerators of diagonal units: check & choose q accordingly
dens = {int(x.q) for M in (A, E1, E2, E1.inv(), E2.inv()) for x in M} | {abs(int(x.p)) for M in (E1, E2) for i in range(n) for x in [M[i, i]]}
q = next(p for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113) if all(dd % p for dd in dens if dd))
Q = q**bdig
Hv = sp.Matrix([rng.randrange(Q) for _ in range(n)]); Dv = sp.Matrix([rng.randrange(Q) for _ in range(n)])
def mm(M, v): return sp.Matrix([int(x) % Q for x in (M * v)])
D1 = mm(modmat(E2), Dv); H1 = mm(modmat(E1.inv()), Hv)
H2 = sp.Matrix([(H1[i] + sum(int(Pi[i, j]) * D1[j] for j in range(n))) % Q for i in range(n)])
H3 = mm(modmat(E1), H2); D3 = mm(modmat(E2.inv()), D1)
target = mm(modmat(A), Dv)
check(f"shear sequence equals (H + A D, D) over Z/{q}^{bdig}", H3 == sp.Matrix([(Hv[i] + target[i]) % Q for i in range(n)]) and D3 == Dv)
# three-step interchange: D<-D-H ; H<-H+D ; D<-H-D  sends (H,D) to (D,H)
H, D = 11, 29
D = D - H; H = H + D; D = H - D
check("D-=H; H+=D; D=H-D sends (H,D) to (D,H)", (H, D) == (29, 11))

# Φ_{M'} Φ_M = Φ_{M'+M} : address map (H,D)->(H+MD,D)
M1 = sp.Matrix(2, 2, lambda i, j: rng.randint(-2, 2)); M2 = sp.Matrix(2, 2, lambda i, j: rng.randint(-2, 2))
Hv = sp.Matrix([3, 4]); Dv = sp.Matrix([1, 2])
step = Hv + M1 * Dv; step = step + M2 * Dv
check("Φ_{M'}Φ_M = Φ_{M'+M}", step == Hv + (M1 + M2) * Dv)

# recurrence (eq:swap-recurrence): F_k <= a F_{k-1} + C, a = s/W ; F_k = O(m^{k tau})
m = 10**6
s_b, W_b = 177176569088785861287000000, 177176569091445000000
a = Fr(s_b, W_b)
tau = 1 - Fr(1, 2**50)
Fk = Fr(1)
ok_rec = True
for k in range(1, 8):
    Fk = a * Fk + 1
    ratio = float(Fk) / float(sp.N(sp.Integer(m)**(k * tau), 30))
    ok_rec &= ratio < 2  # bounded by a constant (here <2) as k grows since a < m^tau
    print(f"   k={k}: F_k / m^(k tau) = {ratio:.6f}")
check("F_k / m^{k tau} stays bounded (a = s_b/W_b < m^tau)", ok_rec)
# geometric sum bound sum_j a^j = O(m^{k tau}) requires a <= m^tau: a/m^tau
print(f"   a/m^tau = {float(a) / float(sp.N(sp.Integer(m)**tau, 30)):.15f}")

# Lemma chunk-swap: rho = ceil(k log W / (2 log q)), R = q^{2 rho} >= W^k ; rho < e and rho = O(log e) for large e
q = 3; Wb = W_b
for e in (10**7, 10**9, 10**12, 10**20):
    k = ceil(log(e) / log(m))
    rho = ceil(k * log(Wb) / (2 * log(q)))
    print(f"   e={e:.0e}: k={k}, rho={rho}, q^(2rho)>=W^k: {2*rho*log(q) >= k*log(Wb)}, rho<e: {rho < e}")
# sum_j a_j m^{j tau} <= (m-1) sum_{j<=k} m^{j tau} = O(m^{k tau}) = O(e^tau) since m^k < m e
e = 10**9; k = ceil(log(e) / log(m))
ssum = sum((m - 1) * m**(j * float(tau)) for j in range(k + 1))
check("sum_j a_j m^{j tau} = O(e^tau) (ratio bounded)", ssum / e**float(tau) < m * (m - 1) * 2, f"{ssum / e**float(tau):.3e}")

# descriptor bound 04-swap.tex:246-256: V_j >= q^{2 e_0}, j log W <= (log_m e_0) log W = O(e_0)
e0 = 10**6; V0 = q**(2 * e0 + 5)
k0 = ceil(log(e0) / log(m))
check("j log W <= (log_m e0) log W = O(e0)", k0 * log(Wb) < e0)

print("\nALL OK" if ok_all else "\nSOME CHECKS FAILED")
