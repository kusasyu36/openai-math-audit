"""Audit 321, part 1: Fox matrix of the presentation (04-construction.tex l.14-21),
the character rho (l.92-94), the matrix D (l.137-144), the minor (l.161-163),
the rows of t_1, t_2, the Laurent matrix D(z) of Lemma 2.1 applied to the
infinite cyclic cover of B given by s, the abelianised matrices, and a
module-theoretic test of Lemma 4.1 (K_ab = 0  <=>  Fox rows of t_1,t_2 lie in
the left ZG-span of the Fox rows of r_1..r_4), tested after mapping ZG to
  (a) Q[G_ab] = Q[Z x Z/6], component-wise over the four fields (PID test), and
  (b) Q[A], Q[B] for the finite quotients A = Z/15 x| Z/2, B = Z/26 x| Z/3.

Run:  uv run --with sympy python check_fox.py
Note: radsimp prints some exact values as floats (e.g. -6.0 - 2.0*sqrt(3)*I = -4(1-zeta^2));
all equality checks are performed symbolically and exactly.
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from fox import gen, mul, inv, power, fox, evaluate, scalar_rep, fox_matrix

ok_all = True
def check(name, cond):
    global ok_all
    ok_all &= bool(cond)
    print(("OK   " if cond else "FAIL ") + name)

# generators: 0=x1, 1=a1, 2=x2, 3=a2, 4=s   (paper's order, l.135)
x1, a1, x2, a2, s = (gen(i) for i in range(5))
r1 = mul(x1, a1, inv(x1), power(a1, -4))      # x1 a1 x1^-1 a1^-4
r2 = mul(x2, a2, inv(x2), power(a2, -3))      # x2 a2 x2^-1 a2^-3
r3 = mul(x1, power(a2, -13))                  # x1 a2^-13
r4 = mul(s, x2, inv(s), power(a1, -5))        # s x2 s^-1 a1^-5
t1 = power(x1, 2)
t2 = power(x2, 3)
R = [r1, r2, r3, r4]

# exact roots of unity
zeta = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2      # primitive cube root
w6 = sp.Rational(1, 2) + sp.sqrt(3) * sp.I / 2         # primitive sixth root
simp = lambda e: sp.radsimp(sp.expand(e))

check("zeta^3 = 1, zeta != 1", simp(zeta**3 - 1) == 0 and simp(zeta - 1) != 0)

# ---- character rho, l.92-94 ----
rho = scalar_rep([-1, zeta, zeta**2, -1, 1])
for lab, r in zip(["r1", "r2", "r3", "r4", "t1", "t2"], R + [t1, t2]):
    check(f"rho({lab}) = 1", simp(rho(r) - 1) == 0)

# ---- Fox matrix D (l.137-144) ----
D = sp.Matrix(fox_matrix(R, 5, rho)).applyfunc(simp)
D_paper = sp.Matrix([
    [1 - zeta, -2, 0, 0, 0],
    [0, 0, 2, zeta**2 - 1, 0],
    [1, 0, 0, -1, 0],
    [0, -(1 + zeta), 1, 0, 1 - zeta**2],
]).applyfunc(simp)
print("D computed (zeta = exp(2 pi i/3)):")
sp.pprint(D)
check("D equals the matrix in eq:fox-matrix (l.137-144)", (D - D_paper).applyfunc(simp) == sp.zeros(4, 5))
minor = simp(D[:, [0, 1, 2, 4]].det())
check(f"minor (cols x1,a1,x2,s) = -4(1-zeta^2)  [{minor}]", simp(minor + 4 * (1 - zeta**2)) == 0)
check("rank D = 4  => H_2(B;L_rho) = 0 (l.164-167)", D.rank(simplify=simp) == 4)
check("chi(B) = 1-5+4 = 0, h0 = 0, h2 = 0 => h_1(B;L_rho) = h_1(G;C_rho) = 0", 1 - 5 + 4 == 0)

# ---- rows of t1, t2 under rho vanish => H_2(B';L_rho) = C^2 ----
rows_t = sp.Matrix(fox_matrix([t1, t2], 5, rho)).applyfunc(simp)
check("Fox rows of t1 = x1^2 and t2 = x2^3 vanish at rho", rows_t == sp.zeros(2, 5))
check("hence h_2(B';L_rho) = 2 = chi(B') and the two 3-cells kill exactly this", 1 - 5 + 6 == 2)

# ---- Lemma 2.1 on the infinite cyclic cover of B given by s -> z ----
z = sp.symbols("z")
rho_z = scalar_rep([-1, zeta, zeta**2, -1, z])
Dz = sp.Matrix(fox_matrix(R, 5, rho_z)).applyfunc(simp)
print("D(z) for the Z-cover (s -> z):")
sp.pprint(Dz)
check("D(1) = D", (Dz.subs(z, 1) - D).applyfunc(simp) == sp.zeros(4, 5))
minors_z = [sp.factor(simp(Dz[:, list(c)].det())) for c in itertools.combinations(range(5), 4)]
print("4x4 minors of D(z):", minors_z)
check("some 4x4 minor of D(z) is a nonzero Laurent polynomial => H_2(B-hat;p^*L) = 0",
      any(m != 0 for m in minors_z))

# ---- abelianisation: integer Fox matrices (trivial rep) ----
triv = scalar_rep([1, 1, 1, 1, 1])
Jab = sp.Matrix(fox_matrix(R + [t1, t2], 5, triv))
print("abelianised 6x5 matrix (rows r1..r4,t1,t2):"); sp.pprint(Jab)
from sympy.matrices.normalforms import smith_normal_form
snf = smith_normal_form(Jab, domain=sp.ZZ)
diag = [abs(snf[i, i]) for i in range(5)]
print("Smith normal form diagonal:", diag)
check("H_1(G;Z) = Z + Z/6 (SNF diag has three 1's, one 6, one 0)", sorted(diag) == [0, 1, 1, 1, 6])
check("abelianised rows r1..r4 have rank 4 over Q (H_1(P;Q) = Q)", Jab[:4, :].rank() == 4)

# =====================================================================
# Module membership test of Lemma 4.1 (K_ab = 0)
# =====================================================================
# (a) over Q[G_ab] = Q[Z x Z/6] = product over the 4 fields of K[z,z^-1].
#     G_ab: s free; a1 -> w^2 (order 3), a2 -> w^3 (order 2), x1 = a2 -> w^3, x2 = 2 a1 -> w^4.
#     Membership in the K[z^-1,z]-span <=> z^k v in the K[z]-span for some k.

def poly(e, ext):
    return sp.Poly(sp.expand(e), z, extension=ext) if ext is not None else sp.Poly(sp.expand(e), z, domain="QQ")

def echelon_over_polys(rows, ext):
    """Row echelon form over K[z] (Euclidean domain) via gcd-type row operations.
    rows: list of lists of sympy expressions polynomial in z."""
    M = [[poly(e, ext) for e in row] for row in rows]
    nrows, ncols = len(M), len(M[0])
    pivots = []
    r = 0
    for c in range(ncols):
        # bring a row with minimal-degree nonzero entry in column c to position r, eliminate others
        while True:
            cand = [i for i in range(r, nrows) if not M[i][c].is_zero]
            if not cand:
                break
            i0 = min(cand, key=lambda i: M[i][c].degree())
            M[r], M[i0] = M[i0], M[r]
            done = True
            for i in range(r + 1, nrows):
                if not M[i][c].is_zero:
                    q, rem = M[i][c].div(M[r][c])
                    M[i] = [M[i][j] - q * M[r][j] for j in range(ncols)]
                    if not rem.is_zero:
                        done = False
            if done:
                pivots.append((r, c))
                r += 1
                break
        if r == nrows:
            break
    return M, pivots

def in_span_laurent(rows, v, ext, maxk=12):
    """Is v in the K[z,z^-1]-span of rows?  rows, v: lists of Laurent-polynomial exprs."""
    # clear negative powers
    N = 0
    for e in list(v) + [e for row in rows for e in row]:
        num, den = sp.fraction(sp.cancel(sp.together(sp.expand(e))))
        dp = sp.Poly(den, z)
        N = max(N, dp.degree())
    rows_p = [[sp.expand(e * z**N) for e in row] for row in rows]
    E, piv = echelon_over_polys(rows_p, ext)
    for k in range(maxk + 1):
        vec = [poly(sp.expand(e * z**(N + k)), ext) for e in v]
        okk = True
        for (ri, ci) in piv:
            if vec[ci].is_zero:
                continue
            q, rem = vec[ci].div(E[ri][ci])
            if not rem.is_zero:
                okk = False
                break
            vec = [vec[j] - q * E[ri][j] for j in range(len(vec))]
        if okk and all(p.is_zero for p in vec):
            return True, k
    return False, None

print("\n(a) membership over Q[G_ab] components")
for label, w, ext in [("trivial", 1, None), ("sign", -1, None), ("zeta3", w6**2, sp.sqrt(-3)), ("zeta6", w6, sp.sqrt(-3))]:
    rep = scalar_rep([w**3, w**2, w**4, w**3, z])
    M = [[simp(e) for e in row] for row in fox_matrix(R, 5, rep)]
    T = [[simp(e) for e in row] for row in fox_matrix([t1, t2], 5, rep)]
    for ti, trow in enumerate(T):
        res, k = in_span_laurent(M, trow, ext)
        check(f"  [{label:7s}] Fox row of t{ti+1} lies in the K[z^+-]-span of rows r1..r4 (z^k needed: {k})", res)

# (b) over Q[A] and Q[B] for the finite quotients (left ZG-module generated by the rows).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_group_lib import SemiDirect, word_images

def membership_finite(Gfin, imgs, label):
    els = Gfin.elements()
    idx = {g: i for i, g in enumerate(els)}
    n = len(els)
    def rep(word):
        # returns the group element (as index) -- we build vectors in Q[G]^5 manually
        r = Gfin.e()
        for g, e in word:
            r = Gfin.mul(r, Gfin.pw(imgs[g], e))
        return r
    def row_vec(relator, left=None):
        v = [0] * (5 * n)
        for j in range(5):
            d = fox(relator, j)
            for word, c in d.items():
                g = rep(word)
                if left is not None:
                    g = Gfin.mul(left, g)
                v[j * n + idx[g]] += c
        return v
    span_rows = [row_vec(r, left=g) for r in R for g in els]
    S = sp.Matrix(span_rows)
    rk = S.rank()
    for ti, t in enumerate([t1, t2]):
        tv = sp.Matrix([row_vec(t)])
        rk2 = sp.Matrix.vstack(S, tv).rank()
        check(f"  [{label}] Fox row of t{ti+1} in left Q[{label}]-span of rows r1..r4 (rank {rk} -> {rk2})", rk2 == rk)
    return rk

print("\n(b) membership over the group rings of the finite quotients A (order 30) and B (order 78)")
A = SemiDirect(15, 2, 4); B = SemiDirect(26, 3, 3)
imgA = {0: (0, 1), 1: (1, 0), 2: A.pw((1, 0), 5), 3: (0, 1), 4: A.e()}
imgB = {0: B.pw((1, 0), 13), 1: B.pw((0, 1), 2), 2: (0, 1), 3: (1, 0), 4: B.e()}
membership_finite(A, imgA, "A")
membership_finite(B, imgB, "B")

print("\nALL OK" if ok_all else "\nSOME CHECKS FAILED")
