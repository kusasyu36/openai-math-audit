"""
109 監査: 03-motifs.tex の有限ネットワークを小さい h で自作実装して確かめる。
 (1) Lemma scalar-motifs: 8 行のスケジュール x 3 段で (X,Y) -> (-Y,X)、scratch 復元（bit: F_2, complex: Q(i)）
 (2) Lemma motif-residuals: ラベル部分空間の入れ子・非退化・残差次元の表（Q 上の形式 I-J/9、F_2 上の内積）
 (3) 次元の総和 = Wm - 2N + 2L（記号計算、一般の h）
 (4) Prop bit-motif-interface: rank(P_V - P_U) = dim V - dim U, 2P_{U_a} の階数 1, P_{U_a^perp}+P_{U_a} = I
実行: uv run --with sympy --with mpmath --with numpy python check_network.py
"""
import random, sys
SKIP_Q = "--skip-q" in sys.argv
from itertools import combinations, product
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

ok_all = True
def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK' if cond else 'NG'}] {name} {detail}")

# ======================= (1) scalar network simulation =======================
def simulate(h, mode, seed=0):
    """mode='bit' (F_2) or 'complex' (Fraction with Gaussian rationals represented as (re,im) Fractions)."""
    rng = random.Random(seed)
    T = list(combinations(range(h), 3))
    v = len(T)
    idx = {t: i for i, t in enumerate(T)}
    if mode == 'bit':
        def rnd(): return rng.randrange(2)
        def add(a, b): return (a + b) % 2
        def neg(a): return a
        def scale(c, a): return (c * a) % 2 if isinstance(c, int) else None
        def neighbors(S, U): return S != U and len(set(S) & set(U)) == 1
        def side_coef(S, U): return 1
        def scatter_coef(): return None
    else:
        # Gaussian rationals as complex Fractions: pairs
        def rnd(): return (Fr(rng.randrange(-8, 9), 2 ** rng.randrange(0, 3)), Fr(rng.randrange(-8, 9), 2 ** rng.randrange(0, 3)))
        def add(a, b): return (a[0] + b[0], a[1] + b[1])
        def neg(a): return (-a[0], -a[1])
        def scale(c, a): return (c * a[0], c * a[1])
        def neighbors(S, U): return S != U and len(set(S) & set(U)) % 2 == 0
        def side_coef(S, U): return -Fr(len(set(S) & set(U)) - 1, 2)
    zero = 0 if mode == 'bit' else (Fr(0), Fr(0))

    # data banks indexed by a in T^3
    X = {a: rnd() for a in product(T, repeat=3)}
    Y = {a: rnd() for a in product(T, repeat=3)}
    X0, Y0 = dict(X), dict(Y)

    def invocation(j, fixed, src, tgt, inverse):
        """one invocation at stage j (coordinate j varies), fixed = other two coords.
        src,tgt are dicts (the logical source/target banks). inverse: run the inverse schedule."""
        def key(t):
            a = list(fixed)
            a.insert(j, t)
            return tuple(a)
        # scratch wires with arbitrary initial values
        A = {(S, U): rnd() for S in T for U in T if neighbors(S, U)}
        C = {i: rnd() for i in range(h)}
        Cs = rnd()
        A0, C0, Cs0 = dict(A), dict(C), Cs
        sgn = -1 if inverse else 1

        def side_injection(s):   # y_S += s * sum coef * A_{SU}
            for S in T:
                acc = zero
                for U in T:
                    if neighbors(S, U):
                        term = A[(S, U)] if mode == 'bit' else scale(side_coef(S, U), A[(S, U)])
                        acc = add(acc, term)
                tgt[key(S)] = add(tgt[key(S)], acc if s == 1 else neg(acc))
        def scatter(s):
            for S in T:
                if mode == 'bit':
                    acc = zero
                    for i in S: acc = add(acc, C[i])
                else:
                    acc = zero
                    for i in S: acc = add(acc, C[i])
                    acc = scale(Fr(1, 2), add(acc, neg(Cs)))
                tgt[key(S)] = add(tgt[key(S)], acc if s == 1 else neg(acc))
        def copy(s):
            for (S, U) in A:
                x = src[key(U)]
                A[(S, U)] = add(A[(S, U)], x if s == 1 else neg(x))
        def gather(s):
            nonlocal Cs
            for i in range(h):
                acc = zero
                for U in T:
                    if i in U: acc = add(acc, src[key(U)])
                C[i] = add(C[i], acc if s == 1 else neg(acc))
            if mode == 'complex':
                acc = zero
                for U in T: acc = add(acc, src[key(U)])
                Cs = add(Cs, acc if s == 1 else neg(acc))
        rows = [lambda: side_injection(-1), lambda: scatter(-1), lambda: copy(1), lambda: gather(1),
                lambda: scatter(1), lambda: side_injection(1), lambda: gather(-1), lambda: copy(-1)]
        inv_rows = [lambda: side_injection(1), lambda: scatter(1), lambda: copy(-1), lambda: gather(-1),
                    lambda: scatter(-1), lambda: side_injection(-1), lambda: gather(1), lambda: copy(1)]
        if not inverse:
            for r in rows: r()
        else:
            for r in range(7, -1, -1): inv_rows[r]()
        restored = (A == A0 and C == C0 and Cs == Cs0)
        return restored

    restored_all = True
    for stage in (1, 2, 3):
        j = stage - 1
        src, tgt, inverse = (X, Y, False) if stage != 2 else (Y, X, True)
        others = [k for k in range(3) if k != j]
        for fixed in product(T, repeat=2):
            restored_all &= invocation(j, fixed, src, tgt, inverse)
    # expected: (X,Y) -> (-Y, X)
    expX = {a: neg(Y0[a]) for a in X0}
    expY = {a: X0[a] for a in X0}
    return restored_all, X == expX, Y == expY

for mode, h in (('bit', 5), ('bit', 6), ('complex', 5)):
    restored, okX, okY = simulate(h, mode)
    check(f"scalar network {mode} h={h}: scratch restored", restored)
    check(f"scalar network {mode} h={h}: (X,Y)->(-Y,X)", okX and okY)

# ======================= (3) symbolic dimension sum =======================
hh, j = sp.symbols('h j', positive=True)
a = hh**(j - 1); f = hh**(3 - j)
X_wire = (a - 1) * (hh - 1) + (hh - 1)
Y_wire = (a - 1) * (hh - 1) + (hh - 1)
side_wire = (a - 1) + ((a - 1) * (hh - 1) + 1) + (hh - 2) + 1 + a * hh * (f - 1)
center_wire = (a - 1) * hh + hh + hh + hh + a * hh * (f - 1)
m = hh**3
check("X wire total over 3 stages = m-1", sp.simplify(sum(X_wire.subs(j, k) for k in (1, 2, 3)) - (m - 1)) == 0)
check("Y wire total over 3 stages = m-1", sp.simplify(sum(Y_wire.subs(j, k) for k in (1, 2, 3)) - (m - 1)) == 0)
check("side wire total = m", sp.simplify(side_wire - m) == 0)
check("center wire total = m + 2h", sp.simplify(center_wire - (m + 2 * hh)) == 0)
vv, I_, z, c, Nn = sp.symbols('v I z c N', positive=True)
Wsym = 2 * Nn + I_ * (vv * z + c)
total = 2 * Nn * (m - 1) + I_ * vv * z * m + I_ * c * (m + 2 * hh)
check("sum |Δdim| = Wm - 2N + 2L with L = I c h", sp.simplify(total - (Wsym * m - 2 * Nn + 2 * I_ * c * hh)) == 0)

# ======================= (2) explicit label subspaces, small h =======================
def tensor3(*vs):
    out = vs[0]
    for w in vs[1:]:
        out = np.kron(out, w)
    return out

def indicator(h, t):
    e = np.zeros(h, dtype=int); e[list(t)] = 1; return e

# ---- Q-version with form G = I - J/9 (bit network labels) ----
def q_checks(h, a_triple):
    G1 = sp.eye(h) - sp.ones(h, h) / 9
    Gm = sp.kronecker_product(G1, G1, G1)  # tensor bilinear form on F^{⊗3}
    mdim = h**3
    T = list(combinations(range(h), 3))
    def col(vec): return sp.Matrix(vec.tolist())
    def span(cols):  # returns Matrix with independent columns
        if not cols: return sp.zeros(mdim, 0)
        M = sp.Matrix.hstack(*cols)
        return M.columnspace() and sp.Matrix.hstack(*M.columnspace())
    def dim(B): return B.shape[1] if B.shape[1] else 0
    def contains(small, big):
        if dim(small) == 0: return True
        return sp.Matrix.hstack(big, small).rank() == big.rank()
    def nondeg(B):
        if dim(B) == 0: return True
        return (B.T * Gm * B).det() != 0
    def perp_in(small, big):  # orthogonal complement of small inside big
        if dim(small) == 0: return big
        if dim(big) == 0: return big
        # vectors big*c with small^T G big c = 0
        K = (small.T * Gm * big).nullspace()
        if not K: return sp.zeros(mdim, 0)
        return sp.Matrix.hstack(*[big * k for k in K])
    def ortho_sum(B1, B2):
        if dim(B1) == 0: return B2
        if dim(B2) == 0: return B1
        return span([B1[:, i] for i in range(dim(B1))] + [B2[:, i] for i in range(dim(B2))])
    def tensor_spaces(*spaces):  # spaces: list of (basis matrix over its own dimension)
        cols = []
        for combo in product(*[[S[:, i] for i in range(S.shape[1])] for S in spaces]):
            vec = combo[0]
            for w in combo[1:]:
                vec = sp.kronecker_product(vec, w)
            cols.append(vec)
        return span(cols)
    Fh = sp.eye(h)
    full = sp.eye(mdim)
    def line(t): return sp.Matrix(indicator(h, t).tolist())
    def perp_h(B):  # complement in F = Q^h
        K = (B.T * G1).nullspace()
        return sp.Matrix.hstack(*K)
    results = []
    for stage in (1, 2, 3):
        jj = stage - 1
        A = tensor_spaces(*[Fh] * jj) if jj > 0 else sp.Matrix([[1]])
        P = tensor_spaces(*[line(a_triple[k]) for k in range(jj)]) if jj > 0 else sp.Matrix([[1]])
        Q = tensor_spaces(*[line(a_triple[k]) for k in range(jj + 1, 3)]) if jj < 2 else sp.Matrix([[1]])
        GA = sp.kronecker_product(*([G1] * jj)) if jj > 0 else sp.Matrix([[1]])
        # B = P^perp in A
        if jj > 0:
            K = (P.T * GA * A).nullspace()
            B = sp.Matrix.hstack(*[A * k for k in K])
        else:
            B = sp.zeros(1, 0)
        tX = line(a_triple[jj])
        # a neighboring triple for the side wire: |S∩T|=1 in the bit construction
        tY = next(line(u) for u in T if u != a_triple[jj] and len(set(u) & set(a_triple[jj])) == 1)
        tXp, tYp = perp_h(tX), perp_h(tY)
        def lab(*parts):  # parts: list of (spaces to tensor) to be orthogonally summed, then tensor Q
            summ = None
            for sp_list in parts:
                t = tensor_spaces(*sp_list)
                summ = t if summ is None else ortho_sum(summ, t)
            return tensor_spaces(summ, Q)
        # gate labels (03-motifs.tex:329-342), all tensored by Q
        L0 = lab([B, tY]) if dim(B) else sp.zeros(mdim, 0)
        L1 = lab([B, Fh]) if dim(B) else sp.zeros(mdim, 0)
        L2 = lab([B, Fh], [P, tX]) if dim(B) else lab([P, tX])
        L3 = lab([A, Fh])
        L4 = L1
        L5 = lab([B, Fh], [P, tYp]) if dim(B) else lab([P, tYp])
        L6 = L3
        L7 = L3
        L2X = lab([B, Fh], [P, tX]) if dim(B) else lab([P, tX])   # X's own gate at time 2 uses t_X
        L5Y = lab([B, Fh], [P, tXp]) if dim(B) else lab([P, tXp])  # Y's gate at time 5 uses its own t (=tX here as the Y wire with triple a)
        # incoming/outgoing data labels
        X_in = lab([A, tX]); X_out = lab([A, Fh])
        Y_in = lab([B, tX]) if dim(B) else sp.zeros(mdim, 0)
        Y_out = lab([B, Fh], [P, tXp]) if dim(B) else lab([P, tXp])
        aa, ff = h**jj, h**(2 - jj)
        edges = {
            "X in->2": (X_in, L2X, (aa - 1) * (h - 1)),
            "X 2->3": (L2X, L3, h - 1),
            "X 3->6": (L3, L6, 0), "X 6->7": (L6, L7, 0), "X 7->out": (L7, X_out, 0),
            "Y in->0(own t)": (Y_in, (lab([B, tX]) if dim(B) else sp.zeros(mdim, 0)), 0),
            "Y 0->1": ((lab([B, tX]) if dim(B) else sp.zeros(mdim, 0)), L1, (aa - 1) * (h - 1)),
            "Y 1->4": (L1, L4, 0), "Y 4->5": (L4, L5Y, h - 1), "Y 5->out": (L5Y, Y_out, 0),
            "center src->1": (sp.zeros(mdim, 0), L1, (aa - 1) * h),
            "center 1->3": (L1, L3, h), "center 3->4 (decreasing)": (L4, L3, h), "center 4->6": (L4, L6, h),
            "center 6->sink": (L6, full, aa * h * (ff - 1)),
            "side src->0": (sp.zeros(mdim, 0), L0, aa - 1),
            "side 0->2": (L0, L2, (aa - 1) * (h - 1) + 1),
            "side 2->5": (L2, L5, h - 2), "side 5->7": (L5, L7, 1),
            "side 7->sink": (L7, full, aa * h * (ff - 1)),
        }
        for name, (small, big, expected) in edges.items():
            inc = contains(small, big)
            nd = nondeg(small) and nondeg(big)
            res = perp_in(small, big)
            d_res = dim(res)
            ok = inc and nd and d_res == expected and d_res == dim(big) - dim(small)
            results.append(ok)
            if not ok:
                print(f"   stage {stage} {name}: contains={inc} nondeg={nd} dim(res)={d_res} expected={expected} dims={dim(small)},{dim(big)}")
        # stage boundary: X_out == next X_in, Y_out == next Y_in
        if stage < 3:
            jn = jj + 1
            An = tensor_spaces(*[Fh] * jn)
            Pn = tensor_spaces(*[line(a_triple[k]) for k in range(jn)])
            Qn = tensor_spaces(*[line(a_triple[k]) for k in range(jn + 1, 3)]) if jn < 2 else sp.Matrix([[1]])
            GAn = sp.kronecker_product(*([G1] * jn))
            Kn = (Pn.T * GAn * An).nullspace()
            Bn = sp.Matrix.hstack(*[An * k for k in Kn])
            tXn = line(a_triple[jn])
            X_in_next = tensor_spaces(tensor_spaces(An, tXn), Qn)
            Y_in_next = tensor_spaces(tensor_spaces(Bn, tXn), Qn)
            results.append(contains(X_out, X_in_next) and contains(X_in_next, X_out))
            results.append(contains(Y_out, Y_in_next) and contains(Y_in_next, Y_out))
        else:
            # final labels: X_out = full, Y_out = U_a^perp
            u = tensor3(indicator(h, a_triple[0]), indicator(h, a_triple[1]), indicator(h, a_triple[2]))
            U = sp.Matrix(u.tolist())
            Uperp = sp.Matrix.hstack(*(U.T * Gm).nullspace())
            results.append(dim(X_out) == mdim)
            results.append(contains(Y_out, Uperp) and contains(Uperp, Y_out))
            # projections: P_V - P_U rank, 2P_U rank 1, P_{U^perp} + P_U = I
            def proj(Bm):
                if dim(Bm) == 0: return sp.zeros(mdim, mdim)
                return Bm * (Bm.T * Gm * Bm).inv() * Bm.T * Gm
            PU, PUp = proj(U), proj(Uperp)
            results.append((PU + PUp - sp.eye(mdim)).is_zero_matrix)
            results.append((2 * PU).rank() == 1)
            results.append((proj(L3) - proj(L1)).rank() == dim(L3) - dim(L1))
    return all(results)

h_small = 5   # m = 125 (h != 9 keeps I-J/9 nondegenerate; h=4 has no |S∩T|=1 neighbors)
a_triple = ((0, 1, 2), (0, 3, 4), (1, 2, 4))
if not SKIP_Q:
    check(f"Q labels (form I-J/9), h={h_small}: nesting, nondegeneracy, residual dims, projections", q_checks(h_small, a_triple))
else:
    print("   (Q label check skipped; see check_network.out from the full run: OK)")

# ---- F_2-version (complex network labels): nesting, residual dims, nondegenerate + nonalternating ----
def gf2_rank(M):
    M = M.copy() % 2; r = 0; rows, cols = M.shape
    for c_ in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c_]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c_]: M[i] ^= M[r]
        r += 1
        if r == rows: break
    return r
def gf2_nullspace(M):
    M = M.copy() % 2; rows, cols = M.shape; r = 0; pivots = []
    for c_ in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c_]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c_]: M[i] ^= M[r]
        pivots.append(c_); r += 1
        if r == rows: break
    free = [c_ for c_ in range(cols) if c_ not in pivots]
    basis = []
    for fcol in free:
        vec = np.zeros(cols, dtype=np.uint8); vec[fcol] = 1
        for i, pc in enumerate(pivots):
            if M[i, fcol]: vec[pc] = 1
        basis.append(vec)
    return np.array(basis, dtype=np.uint8).T if basis else np.zeros((cols, 0), dtype=np.uint8)
def gf2_colspace(M):
    if M.shape[1] == 0: return M
    Mt = M.T.copy() % 2; rows, cols = Mt.shape; r = 0
    for c_ in range(cols):
        piv = next((i for i in range(r, rows) if Mt[i, c_]), None)
        if piv is None: continue
        Mt[[r, piv]] = Mt[[piv, r]]
        for i in range(rows):
            if i != r and Mt[i, c_]: Mt[i] ^= Mt[r]
        r += 1
        if r == rows: break
    return Mt[:r].T.copy()

def f2_checks(h, a_triple):
    mdim = h**3
    T = list(combinations(range(h), 3))
    def line(t): return indicator(h, t).astype(np.uint8).reshape(-1, 1)
    def ten(*Ss):
        out = Ss[0]
        for S in Ss[1:]:
            if out.shape[1] == 0 or S.shape[1] == 0:
                out = np.zeros((out.shape[0] * S.shape[0], 0), dtype=np.uint8); continue
            out = np.kron(out, S) % 2
        return gf2_colspace(out.astype(np.uint8))
    def contains(small, big):
        if small.shape[1] == 0: return True
        return gf2_rank(np.hstack([big, small])) == gf2_rank(big)
    def perp_in(small, big):
        if small.shape[1] == 0 or big.shape[1] == 0: return big
        K = gf2_nullspace((small.T @ big) % 2)
        return gf2_colspace((big @ K) % 2) if K.shape[1] else np.zeros((mdim, 0), dtype=np.uint8)
    def perp_h(B):
        return gf2_nullspace(B.T % 2)
    def nondeg(B):
        if B.shape[1] == 0: return True
        return gf2_rank((B.T @ B) % 2) == B.shape[1]
    def nonalt(B):  # contains a vector of norm 1  <=> some basis vector has odd weight
        return B.shape[1] == 0 or any(int(B[:, i].sum()) % 2 == 1 for i in range(B.shape[1]))
    def osum(B1, B2):
        if B1.shape[1] == 0: return B2
        if B2.shape[1] == 0: return B1
        return gf2_colspace(np.hstack([B1, B2]))
    Fh = np.eye(h, dtype=np.uint8); one = np.array([[1]], dtype=np.uint8)
    full = np.eye(mdim, dtype=np.uint8)
    res_ok = True
    for stage in (1, 2, 3):
        jj = stage - 1
        A = ten(*([Fh] * jj)) if jj else one
        P = ten(*[line(a_triple[k]) for k in range(jj)]) if jj else one
        Q = ten(*[line(a_triple[k]) for k in range(jj + 1, 3)]) if jj < 2 else one
        B = (gf2_colspace((A @ gf2_nullspace((P.T @ A) % 2)) % 2) if jj else np.zeros((1, 0), dtype=np.uint8))
        tX = line(a_triple[jj])
        tY = next(line(u) for u in T if u != a_triple[jj] and len(set(u) & set(a_triple[jj])) % 2 == 0)
        tXp, tYp = perp_h(tX), perp_h(tY)
        def lab(*parts):
            summ = None
            for sl in parts:
                t = ten(*sl)
                summ = t if summ is None else osum(summ, t)
            return ten(summ, Q)
        Z = np.zeros((mdim, 0), dtype=np.uint8)
        L0 = lab([B, tY]) if B.shape[1] else Z
        L1 = lab([B, Fh]) if B.shape[1] else Z
        L2 = lab([B, Fh], [P, tX]) if B.shape[1] else lab([P, tX])
        L3 = lab([A, Fh]); L4 = L1
        L5 = lab([B, Fh], [P, tYp]) if B.shape[1] else lab([P, tYp])
        L5Y = lab([B, Fh], [P, tXp]) if B.shape[1] else lab([P, tXp])
        L6 = L7 = L3
        X_in = lab([A, tX]); X_out = lab([A, Fh])
        Y_in = lab([B, tX]) if B.shape[1] else Z
        Y_out = L5Y
        aa, ff = h**jj, h**(2 - jj)
        edges = {
            "X in->2": (X_in, L2, (aa - 1) * (h - 1)), "X 2->3": (L2, L3, h - 1),
            "Y 0->1": (Y_in, L1, (aa - 1) * (h - 1)), "Y 4->5": (L4, L5Y, h - 1),
            "center src->1": (Z, L1, (aa - 1) * h), "center 1->3": (L1, L3, h),
            "center 3->4": (L4, L3, h), "center 6->sink": (L6, full, aa * h * (ff - 1)),
            "side src->0": (Z, L0, aa - 1), "side 0->2": (L0, L2, (aa - 1) * (h - 1) + 1),
            "side 2->5": (L2, L5, h - 2), "side 5->7": (L5, L7, 1), "side 7->sink": (L7, full, aa * h * (ff - 1)),
        }
        for name, (small, big, expected) in edges.items():
            inc = contains(small, big)
            res = perp_in(small, big)
            ok = inc and nondeg(small) and nondeg(big) and res.shape[1] == expected and nondeg(res) and nonalt(res)
            if not ok:
                print(f"   F2 stage {stage} {name}: contains={inc} dim(res)={res.shape[1]} expected={expected} nondeg={nondeg(res)} nonalt={nonalt(res)}")
            res_ok &= ok
        if stage == 3:
            u = tensor3(*[indicator(h, t) for t in a_triple]).astype(np.uint8).reshape(-1, 1)
            Up = gf2_nullspace(u.T)
            res_ok &= contains(Y_out, Up) and contains(Up, Y_out) and X_out.shape[1] == mdim
            res_ok &= int(u.sum()) % 2 == 1  # norm one
    return res_ok

check(f"F_2 labels, h=5: nesting, residual dims, nondegenerate & nonalternating residuals", f2_checks(5, ((0, 1, 2), (0, 1, 3), (1, 2, 4))))
# h=6 では X∪Y=[6] となる隣接対があり <t_X,t_Y>^perp が交代形式になる（論文の議論は h>6 を使う; 03-motifs.tex:433-435）。h=7 で検査。
check(f"F_2 labels, h=7: same", f2_checks(7, ((0, 1, 2), (3, 4, 5), (0, 2, 6))))

print("\nALL OK" if ok_all else "\nSOME CHECKS FAILED")
