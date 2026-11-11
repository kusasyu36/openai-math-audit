# -*- coding: utf-8 -*-
"""
監査 142: §3 (01-table.tex) の補助素数と表の構成、eq:auxiliary の性格を小さい例で検算する。
 T1  Lemma 3.1: U_q = R_ell^{H_q} は F_p 上 q 次の体 (条件 p^{(ell-1)/q} != 1 mod ell が要)。
     負の対照: p^{(ell-1)/q} == 1 mod ell なら体にならない (q 個の体の積)。
 T2  q=2: v^p = -v の解 v から c = v^2 が非平方 (01-table.tex:100-117)。
 T3  eq:auxiliary の q=2 版は「p を法とする最小の素数平方非剰余 (ただし ell ≡ 1 mod 24)」と同値
     (相互法則)。小さい p で全数確認し、最小素数非剰余との比を記録。
 T4  LTE 恒等式 v_q(Q^d-1) = v_q(Q-1) + v_q(d) (02-finite-fields.tex:139-148, 01-table.tex:153-156)。
"""
import sys
from sympy import Poly, symbols, isprime, primerange, factorint, GF
from sympy import Matrix
import random

def vq(n, q):
    v = 0
    while n % q == 0:
        n //= q; v += 1
    return v

def is_aux(ell, p, q):
    return (isprime(ell) and ell not in (2, 3, p, q) and ell % (12 * q) == 1
            and pow(p, (ell - 1) // q, ell) != 1)

def least_aux(p, q, limit=10**7):
    ell = 12 * q + 1
    while ell < limit:
        if is_aux(ell, p, q):
            return ell
        ell += 12 * q
    return None

# ---------- 線形代数 mod p (素朴な実装) ----------
def nullspace_mod(M, p):
    """M: 行のリスト (各行は長さ n のリスト). 返り値: 核の基底 (列ベクトルのリスト)."""
    M = [row[:] for row in M]
    rows, n = len(M), len(M[0])
    piv_cols = []
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, rows):
            if M[i][c] % p:
                pr = i; break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        piv_cols.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(n) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [0] * n
        v[fc] = 1
        for i, pc in enumerate(piv_cols):
            v[pc] = (-M[i][fc]) % p
        basis.append(v)
    return basis

def solve_coords(basis, vec, p):
    """vec を basis (列ベクトル) の一次結合として表す係数. 解けなければ None."""
    n = len(vec); m = len(basis)
    # 拡大係数行列 [basis | vec]
    A = [[basis[j][i] for j in range(m)] + [vec[i]] for i in range(n)]
    # 消去
    r = 0; piv = []
    for c in range(m):
        pr = None
        for i in range(r, n):
            if A[i][c] % p:
                pr = i; break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][c], -1, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] % p:
                f = A[i][c]
                A[i] = [(a - f * b) % p for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    for i in range(r, n):
        if A[i][m] % p:
            return None
    x = [0] * m
    for i, c in enumerate(piv):
        x[c] = A[i][m]
    return x

# ---------- R_ell = F_p[T]/Phi_ell ----------
def reduce_mod_phi(coeffs, ell, p):
    """coeffs: T^i の係数 (i < ell). Phi_ell = 1+T+...+T^{ell-1}, 基底 T^0..T^{ell-2}."""
    c = [x % p for x in coeffs] + [0] * (ell - len(coeffs))
    top = c[ell - 1]
    # T^{ell-1} = -(1 + T + ... + T^{ell-2})
    out = [(c[i] - top) % p for i in range(ell - 1)]
    return out

def tau_matrix(a, ell, p):
    """tau_a: T -> T^a の行列 (列 = 基底 T^i の像)."""
    n = ell - 1
    cols = []
    for i in range(n):
        img = [0] * ell
        img[(a * i) % ell] = 1
        cols.append(reduce_mod_phi(img, ell, p))
    # 行列 M[row][col]
    return [[cols[j][i] for j in range(n)] for i in range(n)]

def mul_R(x, y, ell, p):
    n = ell - 1
    prod = [0] * (2 * n - 1)
    for i, xi in enumerate(x):
        if xi:
            for j, yj in enumerate(y):
                if yj:
                    prod[i + j] = (prod[i + j] + xi * yj) % p
    # 次数 >= ell を T^ell = 1 で落とす
    red = [0] * ell
    for d, c in enumerate(prod):
        red[d % ell] = (red[d % ell] + c) % p
    return reduce_mod_phi(red, ell, p)

def invariant_algebra(p, q, ell):
    n = ell - 1
    Hq = sorted({pow(a, q, ell) for a in range(1, ell)})
    # H_q は巡回群 (Z/ell)^* の部分群なので巡回. 生成元 g を一つ取り, ker(tau_g - I) が H_q の固定空間.
    # (全 a について連立させると行数 ell^2/q になり遅いので, 生成元 1 個で済ませる. 数学的には同じ空間.)
    g = None
    for cand in Hq:
        if cand == 1:
            continue
        # cand の位数が |H_q| か
        o, x = 1, cand
        while x != 1:
            x = (x * cand) % ell; o += 1
        if o == len(Hq):
            g = cand; break
    assert g is not None
    M = tau_matrix(g, ell, p)
    rows = []
    for i in range(n):
        row = M[i][:]
        row[i] = (row[i] - 1) % p
        rows.append(row)
    basis = nullspace_mod(rows, p)
    # 念のため: 基底が H_q の別の元 (ランダムに 3 個) でも固定されるか確認
    import random as _r
    for a in _r.sample(Hq, min(3, len(Hq))):
        Ma = tau_matrix(a, ell, p)
        for b in basis:
            img = [sum(Ma[i][j] * b[j] for j in range(n)) % p for i in range(n)]
            assert img == b, ("fixed-space check failed", p, q, ell, a)
    return Hq, basis

def minpoly_degree_irreducible(elem, basis, ell, p):
    """elem の U_q 内での最小多項式を求め、(次数, 既約か) を返す."""
    T = symbols('T')
    powers = [[1] + [0] * (ell - 2)]  # 1
    cur = powers[0]
    dim = len(basis)
    for k in range(1, dim + 1):
        cur = mul_R(cur, elem, ell, p)
        powers.append(cur)
    # 最小の k で 1, elem, ..., elem^k が一次従属
    for k in range(1, dim + 1):
        # 解く: elem^k = sum c_i elem^i (i<k)
        coords = solve_coords(powers[:k], powers[k], p)
        if coords is not None:
            poly_coeffs = [1] + [(-c) % p for c in reversed(coords)]
            P = Poly(poly_coeffs, T, modulus=p)
            return k, P.is_irreducible, P
    return None, None, None

def test_T1():
    print("== T1  Lemma 3.1 (01-table.tex:64-91): U_q は q 次の体か")
    ok = True
    cases = []
    for (p, q) in [(7, 3), (13, 3), (5, 3), (11, 5), (7, 5), (19, 3), (3, 5)]:
        ell = least_aux(p, q)
        cases.append((p, q, ell, True))
    # 負の対照: ell ≡ 1 mod 12q だが p が q 乗剰余
    for (p, q) in [(7, 3), (13, 3), (11, 5)]:
        ell = 12 * q + 1
        while True:
            if isprime(ell) and ell not in (2, 3, p, q) and pow(p, (ell - 1) // q, ell) == 1:
                break
            ell += 12 * q
        cases.append((p, q, ell, False))
    for (p, q, ell, aux) in cases:
        Hq, basis = invariant_algebra(p, q, ell)
        dim = len(basis)
        # 基底の第一元 (1 でないもの) の最小多項式
        one = [1] + [0] * (ell - 2)
        elem = None
        for b in basis:
            if b != one:
                elem = b; break
        # ランダム結合で試す (体なら 1 以外の任意元の最小多項式は次数 q で既約、または F_p の元)
        random.seed(1)
        found_field = False
        split_evidence = None
        for _ in range(6):
            coeffs = [random.randrange(p) for _ in basis]
            e = [0] * (ell - 1)
            for c, b in zip(coeffs, basis):
                e = [(x + c * y) % p for x, y in zip(e, b)]
            deg, irr, P = minpoly_degree_irreducible(e, basis, ell, p)
            if deg == dim and irr:
                found_field = True
                break
            if deg is not None and not irr:
                split_evidence = P
        verdict = (dim == q) and (found_field == aux)
        ok &= verdict
        print(f"  p={p:3d} q={q} ell={ell:4d} aux={aux!s:5} |H_q|={len(Hq):4d} dim U_q={dim} "
              f"field_found={found_field} -> {'OK' if verdict else 'NG'}"
              + ("" if aux else f"  (非体の証拠: minpoly {split_evidence.as_expr() if split_evidence is not None else '?'} が可約)"))
    print("  T1:", "OK" if ok else "NG")
    return ok

def test_T2():
    print("== T2  q=2 の非平方構成 (01-table.tex:100-117)")
    ok = True
    for p in [5, 7, 11, 13, 17, 29, 41, 97, 101, 113]:
        q = 2
        ell = least_aux(p, q)
        Hq, basis = invariant_algebra(p, q, ell)
        assert len(basis) == 2, (p, ell, len(basis))
        # Frobenius v -> v^p を U_2 上で: v^p を計算 (繰り返し二乗)
        def pw(x, e):
            r = [1] + [0] * (ell - 2)
            b = x
            while e:
                if e & 1: r = mul_R(r, b, ell, p)
                b = mul_R(b, b, ell, p); e >>= 1
            return r
        # 線形方程式 v^p + v = 0 を U_2 の基底で解く
        imgs = [[(a + b) % p for a, b in zip(pw(bv, p), bv)] for bv in basis]
        # 2 本の列ベクトル imgs[0], imgs[1] の一次従属関係
        # 解く c0*imgs[0] + c1*imgs[1] = 0
        rows = [[imgs[0][i], imgs[1][i]] for i in range(ell - 1)]
        ker = nullspace_mod(rows, p)
        assert len(ker) == 1, (p, ell, len(ker))
        c0, c1 = ker[0]
        v = [(c0 * a + c1 * b) % p for a, b in zip(basis[0], basis[1])]
        v2 = mul_R(v, v, ell, p)
        # v^2 は F_p の元 (定数) か
        is_const = all(x == 0 for x in v2[1:])
        c = v2[0]
        nonsq = pow(c, (p - 1) // 2, p) == p - 1
        s = vq(p - 1, 2)
        omega = pow(c, (p - 1) >> s, p)
        # omega の位数が 2^s か
        order_ok = pow(omega, 2 ** (s - 1), p) == p - 1
        res = is_const and nonsq and order_ok
        ok &= res
        print(f"  p={p:4d} ell_2={ell:4d} c=v^2={c:4d} const={is_const} nonsquare={nonsq} ord(omega_2)=2^{s}: {order_ok} -> {'OK' if res else 'NG'}")
    print("  T2:", "OK" if ok else "NG")
    return ok

def legendre(a, p):
    return pow(a % p, (p - 1) // 2, p)

def test_T3():
    print("== T3  eq:auxiliary (q=2) と最小素数平方非剰余の関係 (00-introduction.tex:52-55)")
    ok = True
    worst_ratio = 0
    rows = []
    for p in primerange(5, 3000):
        ell = least_aux(p, 2)
        # 相互法則で言い換えた条件: 素数 ell ≡ 1 (mod 24), (ell/p) = -1
        ell2 = 25
        while True:
            if isprime(ell2) and ell2 != p and legendre(ell2, p) == p - 1:
                break
            ell2 += 24
        # 最小素数非剰余 (制約なし)
        n_p = 2
        while legendre(n_p, p) != p - 1:
            n_p += 1
            while not isprime(n_p):
                n_p += 1
        ok &= (ell == ell2)
        rows.append((p, ell, n_p))
        worst_ratio = max(worst_ratio, ell / n_p)
    print(f"  素数 p < 3000 の全件で  最小 ell_2 == 最小の素数 ell ≡ 1 (24), (ell/p) = -1 : {'OK' if ok else 'NG'}")
    big = sorted(rows, key=lambda r: -r[1])[:5]
    print("  ell_2 が大きい例 (p, ell_2, 最小素数非剰余 n_p):", big)
    print(f"  max ell_2 / n_p = {worst_ratio:.1f}")
    print("  注: ell_2 は p を法とする素数平方非剰余であり、その最小値を polylog(p) で抑えるのは Ankeny 型 (GRH) の命題。")
    print("  T3:", "OK" if ok else "NG")
    return ok

def test_T4():
    print("== T4  LTE: v_q(Q^d - 1) = v_q(Q-1) + v_q(d), q 奇素数, q | Q-1")
    ok = True
    random.seed(7)
    cnt = 0
    for q in [3, 5, 7, 11]:
        for _ in range(200):
            s = random.randint(1, 3)
            t = random.randint(1, 50)
            while t % q == 0:
                t = random.randint(1, 50)
            Q = q ** s * t + 1
            d = random.randint(1, 60)
            lhs = vq(Q ** d - 1, q)
            rhs = vq(Q - 1, q) + vq(d, q)
            ok &= (lhs == rhs); cnt += 1
    # q=2 では一般に不成立 (Q ≡ 3 mod 4, d=2) の対照
    bad = vq(3 ** 2 - 1, 2) == vq(3 - 1, 2) + vq(2, 2)
    print(f"  奇素数 q, {cnt} 件: {'OK' if ok else 'NG'};  q=2 の対照 (Q=3,d=2) で等式が{'成立' if bad else '不成立'} (論文も q 奇数を仮定: 02-finite-fields.tex:137)")
    print("  T4:", "OK" if ok else "NG")
    return ok

if __name__ == "__main__":
    results = [test_T1(), test_T2(), test_T3(), test_T4()]
    print("ALL:", "OK" if all(results) else "NG")
