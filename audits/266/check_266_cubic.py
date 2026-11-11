"""Companion paper, Section 03-complete-family (cubic obstruction): exact recomputation.
Values in Q(w), w = exp(2 pi i/3), are represented as pairs (x, y) = x + y*w with Fractions,
using w^2 = -1 - w.  All arithmetic is exact."""
import itertools
from fractions import Fraction as Fr

def t_exp(i, j):
    if i == 0 or j == 0 or i == j: return 0
    return 1 if (i - j) % 5 in (1, 4) else -1       # adjacent on the 5-cycle 1-2-3-4-5-1
tE = [[t_exp(i, j) % 3 for j in range(6)] for i in range(6)]

class Qw:
    def __init__(s, x, y=0): s.x = Fr(x); s.y = Fr(y)
    def __add__(s, o): return Qw(s.x+o.x, s.y+o.y)
    def __sub__(s, o): return Qw(s.x-o.x, s.y-o.y)
    def __mul__(s, o):
        if not isinstance(o, Qw): o = Qw(o)
        return Qw(s.x*o.x - s.y*o.y, s.x*o.y + s.y*o.x - s.y*o.y)
    __rmul__ = __mul__
    def conj(s): return Qw(s.x - s.y, -s.y)          # conj(w) = w^2 = -1-w
    def abs2(s): return s.x*s.x - s.x*s.y + s.y*s.y
    def __eq__(s, o): return s.x == o.x and s.y == o.y
    def iszero(s): return s.x == 0 and s.y == 0
    def __repr__(s): return f"({s.x}+{s.y}w)" if s.y else f"{s.x}"
W = [Qw(1,0), Qw(0,1), Qw(-1,-1)]        # w^0, w^1, w^2
def cube_sum(N, denom):                  # (N0 + N1 w + N2 w^2)/denom
    return Qw(Fr(N[0]-N[2], denom), Fr(N[1]-N[2], denom))

# T and T∘2 Hadamard, exactly: inner product of rows i,j = sum_k w^{t_ik - t_jk}
def row_ip(E, i, j):
    N = [0,0,0]
    for k in range(6): N[(E[i][k]-E[j][k]) % 3] += 1
    return cube_sum(N, 1)
def is_hadamard(E):
    return all(row_ip(E, i, j).iszero() for i in range(6) for j in range(6) if i != j) and all(row_ip(E, i, i) == Qw(6) for i in range(6))
print("T Hadamard (exact):", is_hadamard(tE))
tE2 = [[(2*x) % 3 for x in row] for row in tE]
print("T∘2 Hadamard (exact):", is_hadamard(tE2), "; T∘2 == conj(T):", all(tE2[i][j] == (-tE[i][j]) % 3 for i in range(6) for j in range(6)))
plus = sorted((i, j) for i in range(1,6) for j in range(i+1,6) if t_exp(i,j) == 1)
print("plus-graph edges (5-cycle 1-2-3-4-5-1 expected):", plus)
S6 = [[0,0,0,0,0,0],[0,0,1,1,2,2],[0,1,0,2,2,1],[0,1,2,0,1,2],[0,2,2,1,0,1],[0,2,1,2,1,0]]  # exponents of S_6 as recalled from Tadej-Zyczkowski (93)-(94)
pm = [0,1,2,5,4,3]
print("T_ij == S6[p(i),p(j)], p=(0,1,2,5,4,3), with recalled S_6 table:", all(tE[i][j] == S6[pm[i]][pm[j]] for i in range(6) for j in range(6)), "(S_6 table recalled from memory, see report)")

def v(b, k): return sum(b[i]*tE[i][k] for i in range(6)) % 3
def Gcounts(b):
    N = [0,0,0]
    for k in range(6): N[v(b,k)] += 1
    return tuple(N)
def G(b): return cube_sum(Gcounts(b), 6)
e = lambda i: tuple(1 if t == i else 0 for t in range(6))
def add(a, b): return tuple(x+y for x,y in zip(a,b))
def neg(a): return tuple(-x for x in a)
def order(a): return sum(x for x in a if x > 0)
D = [add(e(i), neg(e(j))) for i in range(6) for j in range(6) if i != j]; Dset = set(D)

reps = [(-3,0,0,0,0,3),(-3,0,0,0,1,2),(-3,0,0,1,1,1),(-2,-1,0,0,1,2),(-2,-1,0,1,1,1),(-2,0,0,0,0,2),(-2,0,0,0,1,1),(-1,-1,-1,1,1,1),(-1,-1,0,0,1,1),(-1,0,0,0,0,1),(0,0,0,0,0,0)]
tableN = [(6,0,0),(2,2,2),(3,0,3),(4,1,1),(2,2,2),(2,2,2),(3,3,0),(3,3,0),(4,1,1),(2,2,2),(6,0,0)]
tableU = [(0,0,0),(0,0,0),(0,0,0),(60,42,42),(144,144,144),(0,0,0),(180,180,72),(540,540,216),(480,336,336),(576,576,576),(0,0,0)]
def orbit_rep(b): return min(tuple(sorted(b)), tuple(sorted(neg(b))))
B3 = {b for b in itertools.product(range(-3,4), repeat=6) if sum(b) == 0 and order(b) <= 3}
reps_found = sorted({orbit_rep(b) for b in B3}, key=lambda r: (order(r), r))
print("orbit reps of balanced charges of order<=3:", len(reps_found), "; equals table rows:", sorted(reps) == sorted(reps_found))

a0 = (-1,-1,0,0,1,1); c0 = (0,0,0,0,-1,1)
O = {(tuple(a0[p[i]] for i in range(6)), tuple(c0[p[i]] for i in range(6))) for p in itertools.permutations(range(6))}
print("|O| =", len(O), "(paper 180)")
Oparam = set()
for p_ in range(6):
    for q_ in range(6):
        if p_ == q_: continue
        for h in range(6):
            for i in range(h+1, 6):
                if h in (p_,q_) or i in (p_,q_): continue
                Oparam.add((add(add(e(p_),e(q_)), neg(add(e(h),e(i)))), add(e(q_), neg(e(p_)))))
print("parametrization A=e_p+e_q-e_h-e_i, C=e_q-e_p recovers O exactly:", Oparam == O)

def R_direct(b):   # (eq:cubic-R) evaluated exactly in Q(w)
    tot = Qw(0)
    for (A, C) in O:
        for eps in (-1, 1):
            eb = tuple(eps*x for x in b)
            if add(eb, neg(A)) in Dset:
                prod = G(add(A, neg(C))) * G(add(add(eb, neg(A)), C))
                tot = tot + (prod if eps == 1 else prod.conj())
    return tot * Fr(1,2)
def U_loops(b):    # the integer loop rule (03-complete-family.tex:309-318)
    U = [0,0,0]
    for p_ in range(6):
        for q_ in range(6):
            if p_ == q_: continue
            for h in range(6):
                for i in range(h+1, 6):
                    if h in (p_,q_) or i in (p_,q_): continue
                    A = add(add(e(p_),e(q_)), neg(add(e(h),e(i)))); u = add(add(e(p_),e(p_)), neg(add(e(h),e(i))))
                    for s in (-1,1):
                        sb = tuple(s*x for x in b)
                        if add(sb, neg(A)) in Dset:
                            for k in range(6):
                                for l in range(6):
                                    U[(s*(v(u,k) + v(add(sb,neg(u)),l))) % 3] += 1
    return tuple(U)
okN = okU = okR = True
lam = {}
for j, (rep, Nt, Ut) in enumerate(zip(reps, tableN, tableU)):
    Nc = Gcounts(rep); Uc = U_loops(rep); Rd = R_direct(rep); Gv = G(rep)
    okN &= Nc == Nt; okU &= Uc == Ut; okR &= (cube_sum(Uc, 72) == Rd)
    if not Gv.iszero():
        l = Rd.x/Gv.x if Gv.x != 0 else Rd.y/Gv.y
        lam[j] = l if (Gv*l) == Rd else "not proportional"
    else:
        lam[j] = None
    print(f"j={j} w={rep}: N={Nc} table{Nt} | U={Uc} table{Ut} | G={Gv} R={Rd} R/G={lam[j]}")
print("Table cubic-counts: N column reproduced:", okN, "; U column reproduced:", okU, "; loop rule == direct formula (eq:cubic-R):", okR)
print("multipliers lambda_3,6,7,8 =", [lam[j] for j in (3,6,7,8)], "(paper 1/2, 3, 9, 4); R=0 elsewhere:", all(R_direct(reps[j]).iszero() for j in (0,1,2,4,5,9,10)))
print("G=0 on orbits 1,4,5,9:", all(G(reps[j]).iszero() for j in (1,4,5,9)), "; |G|^2 at types 6,8:", G(reps[6]).abs2(), G(reps[8]).abs2(), "(paper 1/4)")
print("|G(b)| constant on each orbit (all order<=3 charges):", all(G(b).abs2() == G(orbit_rep(b)).abs2() for b in B3))
lam_tab = {3: Fr(1,2), 6: Fr(3), 7: Fr(9), 8: Fr(4)}
rep_index = {r: j for j, r in enumerate(reps)}
print("R(b) = lambda_j G(b) pointwise for ALL order<=3 charges (not only reps):", all(R_direct(b) == G(b)*lam_tab.get(rep_index[orbit_rep(b)], Fr(0)) for b in B3))
def orbit_size(rep):
    s = set()
    for p in itertools.permutations(range(6)):
        bb = tuple(rep[p[i]] for i in range(6)); s.add(bb); s.add(neg(bb))
    return len(s)
sizes = {j: orbit_size(reps[j]) for j in (3,6,7,8)}
coef = {j: sizes[j]*lam_tab[j] for j in sizes}
print("orbit sizes:", sizes, "(paper 360,120,20,90); coefficients:", coef, "(paper 180,360,180,360); /180:", {j: c/180 for j,c in coef.items()})
def classify(wj):
    cnt = {}
    for d in D:
        j = rep_index[orbit_rep(add(wj,d))]; cnt[j] = cnt.get(j,0)+1
    return dict(sorted(cnt.items()))
c5 = classify(reps[5]); c8 = classify(reps[8])
print("w5+D by orbit:", c5, "(paper {0:1,1:8,3:12,6:8,9:1})")
print("w8+D by orbit:", c8, "(paper {3:4,4:8,6:4,7:2,8:8,9:4})")
# Linear system in S_j (orbit averages of Re f conj G).  Known: S1=S4=S5=S9=0 (G=0), S6=S8=-(1/5)|G|^2=-1/20.
known = {1:Fr(0),4:Fr(0),5:Fr(0),9:Fr(0),6:Fr(-1,20),8:Fr(-1,20)}
import sympy as sp
S0,S3,S7 = sp.symbols('S0 S3 S7'); Ssym = {0:S0,3:S3,7:S7}
def lin(cnt): return sum(c*(Ssym[j] if j in Ssym else sp.Rational(known[j])) for j,c in cnt.items())
eq1 = sum(sp.Rational(coef[j]/180)*(Ssym[j] if j in Ssym else sp.Rational(known[j])) for j in coef)
eq5 = lin(c5); eq8 = lin(c8)
print("eq (first relation):", sp.expand(eq1), "= 0   (paper S3+2S6+S7+2S8=0)")
print("eq from a=c=w5 :", sp.expand(eq5), "= 0   (paper: S0+12S3 = 2/5)")
print("eq from a=c=w8 :", sp.expand(eq8), "= 0   (paper: 4S3+2S7 = 3/5)")
sol = sp.solve([eq1, eq5, eq8], [S0,S3,S7], dict=True)[0]
print("solution:", sol, "(paper S3=S7=1/10, S0=-4/5)")
print("E|sum z_i^3|^2 = 6 + 30*S0 =", 6+30*sol[S0], "(paper -18); contradiction:", 6+30*sol[S0] < 0)
print("orbit of 3(e5-e0) has size", orbit_size(reps[0]), "= |D| = 30 ->", orbit_size(reps[0]) == 30)
# Lemma cubic-symmetry (03-complete-family.tex:167-172): count normalized symmetries of T by brute force on row permutations:
# for each row permutation sigma, does there exist a column permutation tau, conjugation flag, and phases making T invariant?
# Dephased form: check whether (conj?)(T[sigma(i), tau(j)]) rephased equals T: equivalently the dephased matrix of the permuted T equals T.
def dephase(E):  # exponent matrix dephased so that row 0 and column 0 are zero
    return [[(E[i][j]-E[i][0]-E[0][j]+E[0][0]) % 3 for j in range(6)] for i in range(6)]
count = 0; rows_with_sym = set()
for sg in itertools.permutations(range(6)):
    found = False
    for tau in itertools.permutations(range(6)):
        for cj in (1, 2):
            E2 = dephase([[ (cj*tE[sg[i]][tau[j]]) % 3 for j in range(6)] for i in range(6)])
            if E2 == tE:
                count += 1; found = True
    if found: rows_with_sym.add(sg)
print("normalized symmetries (row perm, col perm, conj) found by brute force:", count, "(paper: 720); row permutations realized:", len(rows_with_sym), "of 720")
