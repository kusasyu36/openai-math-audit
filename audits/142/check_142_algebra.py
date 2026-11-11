# -*- coding: utf-8 -*-
r"""
監査 142: 代数的な補題を小さい例で検算する。
 A1  Z[xi] で q = unit * lambda^{q-1}, R/(lambda) = F_q (05-geometry.tex:88-97)。
 A2  格子 Lambda での障害: Ne ∈ lambda^{m-1}Λ \ lambda^mΛ, m = 1+(q-1)v_q(N) (05-geometry.tex:449-457),
     W_t の係数増大 ≤ N(2q)^m (07-odd-split.tex:231-237)。
 A3  Frobenius 濾過の補題 (05-geometry.tex:305-327) の抽象版: (R/λ^m)^n 上で Π = 1 + λA なら Π^{q^r} - 1 ≡ 0 (mod λ^{1+(q-1)r})。
 A4  偶数次分割 Prop 2.3 (02-finite-fields.tex:224-272) を N=4,6 で実際の F_p 演算で再現。
 A5  Thm 6.1 の局所行列模型 (03-norms.tex:169-192) の関係式を q=3 で記号的に確認。
 A6  Lemma 6.4 (03-norms.tex:321-352): ブロック上三角代数の共役からトレース形式の核で階数 1 冪等元を取り出す。
 A7  Prop 7.1 (06-division.tex:101-111): 部分和の最大 V が (1-σ)V = E を満たす。
r"""
import random
from fractions import Fraction
import itertools
from sympy import symbols, Matrix, expand, Poly, rem, eye, zeros

# ---------------- Z[xi] (xi: q 次原始根) ----------------
def cyc_reduce(c, q):
    """c: 長さ任意の整数係数 (T^i). T^q = 1 で畳み, Phi_q で T^{q-1} を消す. 返り値: 長さ q-1."""
    r = [0] * q
    for i, x in enumerate(c):
        r[i % q] += x
    top = r[q - 1]
    return [r[i] - top for i in range(q - 1)]

def cyc_mul(a, b, q):
    prod = [0] * (2 * q)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            prod[i + j] += x * y
    return cyc_reduce(prod, q)

def cyc_add(a, b): return [x + y for x, y in zip(a, b)]
def cyc_scale(a, s): return [x * s for x in a]
def cyc_one(q): return [1] + [0] * (q - 2)
def cyc_lam(q):  # 1 - xi
    v = [0] * (q - 1); v[0] = 1; v[1] -= 1
    return v

def cyc_pow(a, e, q):
    r = cyc_one(q)
    for _ in range(e):
        r = cyc_mul(r, a, q)
    return r

def cyc_eval1_mod_q(a, q):
    """R/(lambda) = F_q への像: xi -> 1."""
    return sum(a) % q

def cyc_divide_by_lambda(a, q):
    """a が (1-xi) で割り切れるなら a/(1-xi) を返す, さもなくば None."""
    if sum(a) % q:
        return None
    # 代表多項式 X(T) (deg <= q-2) に c*Phi_q を足して X(1)=0 にし, (1-T) で割る
    c = -sum(a) // q
    X = [x + c for x in a] + [c]  # 長さ q, Phi_q = 1+T+...+T^{q-1}
    assert sum(X) == 0
    # (1 - T) * Y = X  -> Y の係数: Y_0 = X_0, Y_i = X_i + Y_{i-1}
    Y = []
    acc = 0
    for i in range(q - 1):
        acc += X[i]
        Y.append(acc)
    # 検算: (1-T)Y = X  (最後の項 -Y_{q-2} T^{q-1} が X_{q-1} と一致するはず)
    assert -Y[-1] == X[q - 1]
    return cyc_reduce(Y, q)

def cyc_norm(a, q):
    """Z[xi] のノルム: 行列表現の行列式 (Fraction で厳密)."""
    n = q - 1
    cols = []
    basis = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    M = [[0] * n for _ in range(n)]
    for j in range(n):
        img = cyc_mul(a, basis[j], q)
        for i in range(n):
            M[i][j] = Fraction(img[i])
    # 行列式
    det = Fraction(1)
    M = [row[:] for row in M]
    for c in range(n):
        pr = next((i for i in range(c, n) if M[i][c] != 0), None)
        if pr is None:
            return 0
        if pr != c:
            M[c], M[pr] = M[pr], M[c]; det = -det
        det *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return det

def test_A1():
    print("== A1  q = u * lambda^{q-1}, u 単数 (05-geometry.tex:88-97)")
    ok = True
    for q in [3, 5, 7, 11]:
        lam = cyc_lam(q)
        lq = cyc_pow(lam, q - 1, q)
        # q / lambda^{q-1}: q を lambda で q-1 回割る
        a = cyc_scale(cyc_one(q), q)
        for _ in range(q - 1):
            a = cyc_divide_by_lambda(a, q)
            assert a is not None
        u = a
        back = cyc_mul(u, lq, q)
        res = (back == cyc_scale(cyc_one(q), q)) and abs(cyc_norm(u, q)) == 1 and cyc_divide_by_lambda(u, q) is None
        ok &= res
        print(f"  q={q}: q = u*lambda^(q-1) 再構成 {back == cyc_scale(cyc_one(q), q)}, N(u) = {cyc_norm(u, q)}, u ∉ (lambda): {cyc_divide_by_lambda(u, q) is None} -> {'OK' if res else 'NG'}")
    print("  A1:", "OK" if ok else "NG")
    return ok

# ---------------- 無限遠格子 Lambda (巡回座標) ----------------
def lat_apply_1_minus_sigma(W):
    """sigma ∞_j = ∞_{j+1}. ((1-σ)W)_j = W_j - W_{j-1}."""
    q = len(W)
    return [W[j] - W[j - 1] for j in range(q)]

def lat_divide(W):
    """(1-σ)W' = W の Λ 内の解 (一意) を返す. 存在しなければ None."""
    q = len(W)
    if sum(W) != 0:
        return None
    partial = []
    acc = 0
    for j in range(q):
        acc += W[j]
        partial.append(acc)   # W'_j = W'_0 + sum_{i=1..j} W_i  ただし j=0 の項は W_0 ... 調整
    # W'_j - W'_{j-1} = W_j for j>=1, and W'_0 - W'_{q-1} = W_0.
    # 取り方: W'_j = c + sum_{i=1}^{j} W_i  (j>=0), 条件 j=0: c - W'_{q-1} = W_0 <=> -sum_{i=1}^{q-1} W_i = W_0 <=> sum W = 0 OK
    S = [sum(W[1:j + 1]) for j in range(q)]
    tot = sum(S)
    if tot % q:
        return None
    c = -tot // q
    Wp = [c + s for s in S]
    assert lat_apply_1_minus_sigma(Wp) == W and sum(Wp) == 0
    return Wp

def vq(n, q):
    v = 0
    while n % q == 0:
        n //= q; v += 1
    return v

def test_A2():
    print("== A2  Ne の lambda 可除性 (05-geometry.tex:449-457) と W_t の係数増大 (07-odd-split.tex:231-237)")
    ok = True
    for q in [3, 5, 7]:
        for N in [q, 2 * q, q * q, 4 * q, 2 * q * q, q ** 3, 7 * q * q]:
            r = vq(N, q); m = 1 + (q - 1) * r
            e = [0] * q; e[0] = -1; e[1] = 1       # e = ∞_1 - ∞_0
            W = [N * x for x in e]
            depth = 0
            maxcoef = [max(abs(x) for x in W)]
            cur = W
            while True:
                nxt = lat_divide(cur)
                if nxt is None:
                    break
                cur = nxt; depth += 1
                maxcoef.append(max(abs(x) for x in cur))
            # W_t = lambda^{-(t-1)} N e for t=1..m exists, W_m ∉ lambda Λ  <=> depth == m-1
            bound_ok = all(maxcoef[t] <= N * (2 * q) ** (t + 1) for t in range(len(maxcoef)))
            growth_ok = all(maxcoef[t + 1] <= 2 * q * maxcoef[t] for t in range(len(maxcoef) - 1))
            res = (depth == m - 1) and bound_ok and growth_ok
            ok &= res
            print(f"  q={q} N={N:4d} r={r} m={m:2d}: 可除回数={depth} (期待 {m-1}), max|W_t|={maxcoef}, 1 ステップ ≤2q 倍: {growth_ok} -> {'OK' if res else 'NG'}")
    print("  A2:", "OK" if ok else "NG")
    return ok

# ---------------- 濾過の補題 (抽象版) ----------------
def mat_mul(A, B, q):
    n = len(A)
    C = [[[0] * (q - 1) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            acc = [0] * (q - 1)
            for k in range(n):
                acc = cyc_add(acc, cyc_mul(A[i][k], B[k][j], q))
            C[i][j] = acc
    return C

def lam_valuation(a, q, cap):
    """a の lambda 進付値 (cap で打ち切り). a=0 なら cap."""
    if all(x == 0 for x in a):
        return cap
    v = 0
    while v < cap:
        b = cyc_divide_by_lambda(a, q)
        if b is None:
            return v
        a = b; v += 1
    return cap

def test_A3():
    print("== A3  Π = 1 + λA (A 任意) ⇒ Π^{q^r} - 1 ≡ 0 mod λ^{1+(q-1)r} (05-geometry.tex:305-327 の抽象版)")
    ok = True
    random.seed(3)
    for q, r, n in [(3, 1, 2), (3, 2, 2), (5, 1, 2), (3, 3, 1), (7, 1, 2), (5, 2, 1)]:
        m = 1 + (q - 1) * r
        lam = cyc_lam(q)
        A = [[[random.randint(-3, 3) for _ in range(q - 1)] for _ in range(n)] for _ in range(n)]
        Pi = [[cyc_add(cyc_mul(lam, A[i][j], q), cyc_one(q) if i == j else [0] * (q - 1)) for j in range(n)] for i in range(n)]
        P = Pi
        for _ in range(q ** r - 1):
            P = mat_mul(P, Pi, q)
        # P - I
        minv = 10 ** 9
        for i in range(n):
            for j in range(n):
                d = cyc_add(P[i][j], cyc_scale(cyc_one(q), -1) if i == j else [0] * (q - 1))
                minv = min(minv, lam_valuation(d, q, m + 5))
        # 鋭さの対照: Π^{q^{r-1}} - 1 は一般に λ^{m} で割れない
        P2 = Pi
        for _ in range(q ** (r - 1) - 1):
            P2 = mat_mul(P2, Pi, q)
        minv2 = 10 ** 9
        for i in range(n):
            for j in range(n):
                d = cyc_add(P2[i][j], cyc_scale(cyc_one(q), -1) if i == j else [0] * (q - 1))
                minv2 = min(minv2, lam_valuation(d, q, m + 5))
        res = minv >= m
        ok &= res
        print(f"  q={q} r={r} n={n} m={m}: v_λ(Π^(q^r)-1) ≥ {minv}  (必要 ≥ {m}); 対照 v_λ(Π^(q^(r-1))-1) = {minv2} -> {'OK' if res else 'NG'}")
    print("  A3:", "OK" if ok else "NG")
    return ok

# ---------------- 偶数次分割 Prop 2.3 ----------------
def test_A4():
    print("== A4  Prop 2.3 偶数次分割 (02-finite-fields.tex:224-272): 向き付けと行カウント")
    ok = True
    random.seed(11)
    trials = 0; fails = 0
    for p in [97, 113, 193, 257, 1009]:
        s = vq(p - 1, 2); t = (p - 1) >> s
        # 2-primary 生成元 u: 非平方 c から u = c^t
        c = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1)
        u = pow(c, t, p)
        assert pow(u, 2 ** (s - 1), p) == p - 1
        # log to base u in the 2-group (位数 2^s なので表で十分)
        logtab = {pow(u, k, p): k for k in range(2 ** s)}
        for N in [2, 4, 6, 8]:
            for _ in range(30):
                roots = random.sample(range(p), N)
                # 向き付け: ordered pair (i,j), Δ = x_i - x_j, log(Δ^t) ∈ [0, 2^{s-1}) なら選択
                sel = {}
                for i in range(N):
                    for j in range(N):
                        if i == j: continue
                        d = pow((roots[i] - roots[j]) % p, t, p)
                        sel[(i, j)] = logtab[d] < 2 ** (s - 1)
                # 各非順序対でちょうど一方が選ばれる
                exactly_one = all(sel[(i, j)] != sel[(j, i)] for i in range(N) for j in range(i + 1, N))
                rowcounts = [sum(1 for j in range(N) if j != i and sel[(i, j)]) for i in range(N)]
                not_all_equal = len(set(rowcounts)) > 1
                total_ok = sum(rowcounts) == N * (N - 1) // 2
                # gcd テスト: c - a が消える成分の集合 (a = rowcount 値) は真部分集合
                a = rowcounts[0]
                subset = [i for i in range(N) if rowcounts[i] == a]
                proper = 0 < len(subset) < N
                trials += 1
                if not (exactly_one and not_all_equal and total_ok and proper):
                    fails += 1
        print(f"  p={p} (p-1 = 2^{s}·{t}): 試行累計 {trials}, 失敗 {fails}")
    ok = fails == 0
    print("  A4:", "OK" if ok else "NG")
    return ok

# ---------------- 局所行列模型 (q=3, 記号) ----------------
def test_A5():
    print("== A5  局所行列模型の関係式 (03-norms.tex:169-192), q=3, 記号計算 (z^2+z+1=0)")
    z, F, R, G, S, u = symbols('z F R G S u')
    q = 3
    def red(M):
        return M.applyfunc(lambda e: rem(expand(e), z**2 + z + 1, z))
    zinv = z**2  # z^{-1} = z^2
    # 分岐点の模型: V' e_j = z^j R e_j, Y e_j = e_{j+1}, Y e_{q-1} = F e_0
    Vp = zeros(q, q); Y = zeros(q, q)
    for j in range(q):
        Vp[j, j] = z**j * R
        if j < q - 1:
            Y[j + 1, j] = 1
        else:
            Y[0, j] = F
    rel1 = red(Vp * Y - z * Y * Vp)
    rel2 = red(Y**3 - F * eye(q))
    rel3 = red(Vp**3 - R**3 * eye(q))
    ram_ok = rel1.is_zero_matrix and rel2.is_zero_matrix and rel3.is_zero_matrix
    # 非分岐の模型: Y e_j = z^{-j} S e_j, V e_j = e_{j+1}, V e_{q-1} = G e_0
    Y2 = zeros(q, q); V2 = zeros(q, q)
    for j in range(q):
        Y2[j, j] = zinv**j * S
        if j < q - 1:
            V2[j + 1, j] = 1
        else:
            V2[0, j] = G
    rel4 = red(V2 * Y2 - z * Y2 * V2)
    rel5 = red(V2**3 - G * eye(q))
    rel6 = red(Y2**3 - S**3 * eye(q))
    unr_ok = rel4.is_zero_matrix and rel5.is_zero_matrix and rel6.is_zero_matrix
    # 可換子因子: (Y^{-w} V)^q = zeta^{-w q(q-1)/2} F^{-w} G, q 奇数なので zeta 因子は 1 (03-norms.tex:155-160)
    comm_ok = all((w * q * (q - 1) // 2) % q == 0 for w in range(-10, 11)) and (q * (q - 1) // 2) % 2 == 1  # 2 では成立しない対照
    comm_ok_q2 = ((2 * (2 - 1) // 2) % 2 == 0)  # q=2 なら因子 zeta^{w} が残る
    ok = ram_ok and unr_ok and comm_ok and (not comm_ok_q2)
    print(f"  分岐模型: V'Y = ζYV' {rel1.is_zero_matrix}, Y^3 = F {rel2.is_zero_matrix}, V'^3 = R^3 {rel3.is_zero_matrix}")
    print(f"  非分岐模型: VY = ζYV {rel4.is_zero_matrix}, V^3 = G {rel5.is_zero_matrix}, Y^3 = S^3 {rel6.is_zero_matrix}")
    print(f"  ζ^(w q(q-1)/2) = 1 (q 奇数) : {comm_ok}; q=2 では消えない対照: {not comm_ok_q2}")
    # トレース座標 (03-norms.tex:133-137): 非単位基底単項式のトレースが 0 — 分岐模型で確認
    traces = {}
    for a in range(q):
        for b in range(q):
            Wab = red(Y**a * Vp**b)
            traces[(a, b)] = rem(expand(Wab.trace()), z**2 + z + 1, z)
    tr_ok = all((traces[k] == 0) for k in traces if k != (0, 0)) and traces[(0, 0)] == q
    print(f"  tr(W_ab) = 0 for (a,b) ≠ (0,0), tr(1) = q : {tr_ok}")
    ok &= tr_ok
    print("  A5:", "OK" if ok else "NG")
    return ok

# ---------------- Lemma 6.4 ----------------
def rank_mod(M, p):
    M = [row[:] for row in M]; n = len(M); m = len(M[0]); r = 0
    for c in range(m):
        pr = next((i for i in range(r, n) if M[i][c] % p), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(n):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
        if r == n: break
    return r

def nullspace_mod(M, p):
    M = [row[:] for row in M]; rows = len(M); n = len(M[0]); piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, rows) if M[i][c] % p), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == rows: break
    free = [c for c in range(n) if c not in piv]
    out = []
    for fc in free:
        v = [0] * n; v[fc] = 1
        for i, pc in enumerate(piv): v[pc] = (-M[i][fc]) % p
        out.append(v)
    return out

def mm(A, B, p):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) % p for j in range(m)] for i in range(n)]

def inv_mod(A, p):
    n = len(A)
    M = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        pr = next((i for i in range(c, n) if M[i][c] % p), None)
        if pr is None: return None
        M[c], M[pr] = M[pr], M[c]
        inv = pow(M[c][c], -1, p)
        M[c] = [(x * inv) % p for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[c])]
    return [row[n:] for row in M]

def test_A6():
    print("== A6  Lemma 6.4 (03-norms.tex:321-352): トレース形式の核 → 共通核 U → 階数 1 冪等元 vw ∈ B")
    p = 101
    random.seed(5)
    ok = True
    for blocks in [(2, 1, 2), (1, 1, 1), (3, 2), (5,), (1, 2, 2)]:
        n = sum(blocks)
        # ブロック上三角代数の基底 (行列単位 E_ij で i のブロック ≤ j のブロック)
        bidx = []
        for bi, sz in enumerate(blocks):
            bidx += [bi] * sz
        basis0 = []
        for i in range(n):
            for j in range(n):
                if bidx[i] <= bidx[j]:
                    E = [[0] * n for _ in range(n)]; E[i][j] = 1
                    basis0.append(E)
        # ランダム共役
        while True:
            Pm = [[random.randrange(p) for _ in range(n)] for _ in range(n)]
            Pinv = inv_mod(Pm, p)
            if Pinv is not None: break
        basis = [mm(mm(Pm, E, p), Pinv, p) for E in basis0]
        d = len(basis)
        flat = lambda M: [M[i][j] for i in range(n) for j in range(n)]
        # トレース形式 Gram 行列
        Gm = [[sum(mm(basis[a], basis[b], p)[i][i] for i in range(n)) % p for b in range(d)] for a in range(d)]
        ker = nullspace_mod(Gm, p)   # N の係数ベクトル
        Nmats = []
        for kv in ker:
            M = [[0] * n for _ in range(n)]
            for c, Bm in zip(kv, basis):
                if c:
                    for i in range(n):
                        for j in range(n):
                            M[i][j] = (M[i][j] + c * Bm[i][j]) % p
            Nmats.append(M)
        expected_dimN = sum(blocks[i] * blocks[j] for i in range(len(blocks)) for j in range(len(blocks)) if i < j)
        # 共通核 U
        if Nmats:
            stacked = [row for M in Nmats for row in M]
            U = nullspace_mod(stacked, p)
        else:
            U = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
        dimU_ok = len(U) == blocks[0]
        v = U[0]
        # w: wv = 1 となる行ベクトル
        i0 = next(i for i in range(n) if v[i] % p)
        w = [0] * n; w[i0] = pow(v[i0], -1, p)
        a0 = [[(v[i] * w[j]) % p for j in range(n)] for i in range(n)]
        idem = mm(a0, a0, p) == a0
        rank1 = rank_mod(a0, p) == 1
        # a0 ∈ span(basis)?
        A = [[flat(Bm)[k] for Bm in basis] + [flat(a0)[k]] for k in range(n * n)]
        member = rank_mod(A, p) == rank_mod([row[:-1] for row in A], p)
        res = (len(ker) == expected_dimN) and dimU_ok and idem and rank1 and member
        ok &= res
        print(f"  blocks={blocks}: dim N={len(ker)} (期待 {expected_dimN}), dim U={len(U)} (期待 {blocks[0]}), a0^2=a0 {idem}, rank 1 {rank1}, a0 ∈ B {member} -> {'OK' if res else 'NG'}")
    print("  A6:", "OK" if ok else "NG")
    return ok

# ---------------- Prop 7.1 の max トリック ----------------
def test_A7():
    print("== A7  Prop 7.1 (06-division.tex:101-111): V = max_j S_j は (1-σ)V = E, V ≥ 0")
    ok = True
    random.seed(9)
    for q in [3, 5, 7]:
        for _ in range(200):
            # 点集合: 複数の σ-軌道 (長さ q) と σ-固定点 (分岐点, 長さ 1)
            norb = random.randint(1, 3); nfix = random.randint(0, 2)
            # E: 各軌道上の係数ベクトル, 条件 Σ_a σ^a E = 0 ⇔ 各軌道で係数和 0, 固定点で係数 0
            orbits = []
            for _ in range(norb):
                c = [random.randint(-4, 4) for _ in range(q)]
                c[-1] -= sum(c)
                orbits.append(c)
            fixed = [0] * nfix
            def sigma(E):   # σ: 軌道内で j -> j+1 (点の移動). 関数への作用は逆だが divisor への作用として定義
                return [[o[(j - 1) % q] for j in range(q)] for o in E[0]], E[1][:]
            def add(E1, E2):
                return [[a + b for a, b in zip(o1, o2)] for o1, o2 in zip(E1[0], E2[0])], [a + b for a, b in zip(E1[1], E2[1])]
            def vmax(Es):
                return [[max(o[j] for o in [E[0][k] for E in Es]) for j in range(q)] for k in range(norb)], [max(E[1][k] for E in Es) for k in range(nfix)]
            E = (orbits, fixed)
            S = [([[0] * q for _ in range(norb)], [0] * nfix)]
            cur = E
            for j in range(1, q):
                S.append(add(S[-1], cur) if j > 1 else cur)
                cur = sigma(cur)
            # S_j = sum_{a<j} σ^a E ; 上の構成を正確に: S_1 = E, S_2 = E + σE, ...
            S = [([[0] * q for _ in range(norb)], [0] * nfix)]
            term = E
            for j in range(1, q):
                S.append(add(S[-1], term))
                term = sigma(term)
            # S_q = 0 の確認
            Sq = add(S[-1], term)
            zero_ok = all(x == 0 for o in Sq[0] for x in o) and all(x == 0 for x in Sq[1])
            V = vmax(S)
            sV = sigma(V)
            lhs = ([[a - b for a, b in zip(o1, o2)] for o1, o2 in zip(V[0], sV[0])], [a - b for a, b in zip(V[1], sV[1])])
            eff = all(x >= 0 for o in V[0] for x in o) and all(x >= 0 for x in V[1])
            res = zero_ok and eff and lhs[0] == E[0] and lhs[1] == E[1]
            ok &= res
            if not res:
                print("   NG example", q, E, V)
        print(f"  q={q}: 200 件 {'OK' if ok else 'NG'}")
    print("  A7:", "OK" if ok else "NG")
    return ok

if __name__ == "__main__":
    results = [test_A1(), test_A2(), test_A3(), test_A4(), test_A5(), test_A6(), test_A7()]
    print("ALL:", "OK" if all(results) else "NG")
