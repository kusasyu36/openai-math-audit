# -*- coding: utf-8 -*-
r"""
監査 142: §4 (05-geometry.tex) の核心 Lemma 4.3 / Cor 4.5 を、最小の実例 (q=3, N=3, 種数 1) で総当たり検証する。
  C: Y^3 = F(X) = (X-x1)(X-x2)(X-x3) over K = F_p (p ≡ 1 mod 3), k = F_{p^3} ([k:K] = q^r, r = v_3(3) = 1, m = 3).
  平面 3 次曲線なので J = Pic^0(C) ≅ C 自身 (原点 O = ∞_0). 弦接線法で群演算を実装し、
  σ_J(P) = σ(P) ⊖ ∞_1 ([σP − σO] = [σP − ∞_1]) で σ を J に誘導, λ = 1 − σ_J.
検証項目:
  C1  群法則の健全性 (結合律のランダム検査, 1+σ+σ^2 = 0 on J(k)).
  C2  div(Y) の関係 Σ[P_i − O] = [∞_1 + ∞_2 − 2O] (05-geometry.tex:204-214): P1⊕P2⊕P3 = ∞_1⊕∞_2.
  C3  Lemma 4.3 の帰結: J[λ^t] ⊂ J(k) for t ≤ m = 3 (|J(k)[λ^t]| = 3^t), しかし J(K) では t = 1 までしか保証されない.
  C4  Cor 4.5: [P_i − O] の λ-割り算の鎖 (t ≤ m) は全て J(k) 内に存在する: |λ^{-1}(P_i) ∩ J(k)| = 3, |λ^{-2}(P_i) ∩ J(k)| = 9.
      対照: 深さ 3 (λ^3 d = [P_i − O], t = 4 > m) では J(k) 内の解が 27 個に達しないことを観察 (保証なし).
  C5  Prop 4.6 の格子側: 3e = u λ^2 e, W_1 = 3e, W_2, W_3 ∉ λΛ (A2 と同じだが q=3,N=3 の値を明示).
"""
import itertools, random, sys

p = 7
# k = F_{p^3} = F_p[t]/(t^3 - 2)  (2 は mod 7 で非立方)
def kadd(a, b): return tuple((x + y) % p for x, y in zip(a, b))
def ksub(a, b): return tuple((x - y) % p for x, y in zip(a, b))
def kneg(a): return tuple((-x) % p for x in a)
def kmul(a, b):
    c = [0] * 5
    for i in range(3):
        for j in range(3):
            c[i + j] = (c[i + j] + a[i] * b[j]) % p
    # t^3 = 2, t^4 = 2t
    return ((c[0] + 2 * c[3]) % p, (c[1] + 2 * c[4]) % p, c[2] % p)
def kconst(c): return (c % p, 0, 0)
ZERO = kconst(0); ONE = kconst(1)
def kpow(a, e):
    r = ONE
    while e:
        if e & 1: r = kmul(r, a)
        a = kmul(a, a); e >>= 1
    return r
def kinv(a):
    assert a != ZERO
    return kpow(a, p ** 3 - 2)
def kis_zero(a): return a == ZERO
def kscal(a, c): return tuple((x * c) % p for x in a)

K_elems = [kconst(c) for c in range(p)]
k_elems = [(a, b, c) for a in range(p) for b in range(p) for c in range(p)]
assert kpow((0, 1, 0), p ** 3 - 1) == ONE  # t^(343-1)=1 (体であることの弱い検査)
# ζ: 3 次原始根 in F_7: 2 (2^3 = 8 ≡ 1)
zeta = kconst(2); zeta2 = kmul(zeta, zeta)
assert kmul(zeta2, zeta) == ONE and zeta != ONE

roots = [0, 1, 3]
rk = [kconst(x) for x in roots]
q = 3; N = 3

# 射影 3 次曲線 G(X,Y,Z) = Y^3 − Π (X − x_i Z)
def G(P):
    X, Y, Z = P
    prod = ONE
    for xi in rk:
        prod = kmul(prod, ksub(X, kmul(xi, Z)))
    return ksub(kpow(Y, 3), prod)
def gradG(P):
    X, Y, Z = P
    lin = [ksub(X, kmul(xi, Z)) for xi in rk]
    gX = ZERO; gZ = ZERO
    for i in range(3):
        others = ONE
        for j in range(3):
            if j != i: others = kmul(others, lin[j])
        gX = ksub(gX, others)
        gZ = kadd(gZ, kmul(rk[i], others))
    gY = kscal(kmul(Y, Y), 3)
    return (gX, gY, gZ)

def normalize(P):
    X, Y, Z = P
    if not kis_zero(Z):
        inv = kinv(Z); return (kmul(X, inv), kmul(Y, inv), ONE)
    if not kis_zero(X):
        inv = kinv(X); return (ONE, kmul(Y, inv), ZERO)
    inv = kinv(Y); return (ZERO, ONE, ZERO)

def padd(P, Q): return tuple(kadd(a, b) for a, b in zip(P, Q))
def psub(P, Q): return tuple(ksub(a, b) for a, b in zip(P, Q))
def pscal(P, c): return tuple(kmul(a, c) for a in P)

def cross_zero(P, Q):
    X1, Y1, Z1 = P; X2, Y2, Z2 = Q
    return (kis_zero(ksub(kmul(Y1, Z2), kmul(Z1, Y2))) and kis_zero(ksub(kmul(Z1, X2), kmul(X1, Z2)))
            and kis_zero(ksub(kmul(X1, Y2), kmul(Y1, X2))))

half = kconst(pow(2, -1, p))

def third_point(P, Q):
    """弦接線法の P*Q: P,Q を通る直線 (P=Q なら接線) と C の第三交点."""
    if cross_zero(P, Q):
        # 接線: ∇G(P)·v = 0 を満たす P と独立な v
        a, b, c = gradG(P)
        if not kis_zero(a):
            ia = kinv(a); cands = [(kneg(kmul(b, ia)), ONE, ZERO), (kneg(kmul(c, ia)), ZERO, ONE)]
        elif not kis_zero(b):
            ib = kinv(b); cands = [(ONE, ZERO, ZERO), (ZERO, kneg(kmul(c, ib)), ONE)]
        else:
            cands = [(ONE, ZERO, ZERO), (ZERO, ONE, ZERO)]
        Qp = cands[0] if not cross_zero(P, cands[0]) else cands[1]
        assert not cross_zero(P, Qp)
        # G(sP + tQ') = t^2 (γ s + δ t); γ = (G(P+Q')+G(P−Q'))/2, δ = (G(P+Q')−G(P−Q'))/2 ; R = δP − γQ'
        gp = G(padd(P, Qp)); gm = G(psub(P, Qp))
        gam = kmul(kadd(gp, gm), half); dl = kmul(ksub(gp, gm), half)
        R = psub(pscal(P, dl), pscal(Qp, gam))
    else:
        # G(sP + tQ) = s t (α s + β t); α = (G(P+Q) − G(P−Q))/2, β = (G(P+Q) + G(P−Q))/2 ; R = βP − αQ
        gp = G(padd(P, Q)); gm = G(psub(P, Q))
        al = kmul(ksub(gp, gm), half); be = kmul(kadd(gp, gm), half)
        R = psub(pscal(P, be), pscal(Q, al))
    assert not all(kis_zero(c) for c in R), ("degenerate", P, Q)
    R = normalize(R)
    assert kis_zero(G(R)), ("third point not on curve", P, Q, R)
    return R

O = normalize((ONE, ONE, ZERO))          # ∞_0 : Y/X = 1
inf1 = normalize((ONE, zeta, ZERO))      # ∞_1 : Y/X = ζ
inf2 = normalize((ONE, zeta2, ZERO))
OO = third_point(O, O)

def gadd(P, Q): return third_point(O, third_point(P, Q))
def gneg(P): return third_point(OO, P)
def gsub(P, Q): return gadd(P, gneg(Q))
def gmul(n, P):
    R = O; B = P
    if n < 0: B = gneg(B); n = -n
    while n:
        if n & 1: R = gadd(R, B)
        B = gadd(B, B); n >>= 1
    return R

def points_over(field):
    pts = set()
    cubes = {}
    for y in field:
        cubes.setdefault(kpow(y, 3), []).append(y)
    for x in field:
        Fx = ONE
        for xi in rk: Fx = kmul(Fx, ksub(x, xi))
        for y in cubes.get(Fx, []):
            pts.add(normalize((x, y, ONE)))
    # 無限遠点 (Y/X)^3 = 1
    for z3 in [ONE, zeta, zeta2]:
        if z3 in set(field):  # 定義体に ζ があるか (F_7 にはある)
            pts.add(normalize((ONE, z3, ZERO)))
    return pts

def sigma_pt(P): return normalize((P[0], kmul(zeta, P[1]), P[2]))
def sigma_J(P): return gsub(sigma_pt(P), inf1)
def lam(P): return gsub(P, sigma_J(P))

def main():
    ok = True
    Ck = points_over(k_elems); CK = points_over(K_elems)
    print(f"p={p}, F = (X-0)(X-1)(X-3), q=3, N=3, g=1. |C(F_7)| = {len(CK)}, |C(F_343)| = {len(Ck)}")
    # C1 群法則の健全性
    random.seed(1)
    L = list(Ck)
    assoc = all(gadd(gadd(A, B), C) == gadd(A, gadd(B, C)) for A, B, C in [tuple(random.choice(L) for _ in range(3)) for _ in range(300)])
    ident = all(gadd(A, O) == A for A in L)
    inv_ok = all(gadd(A, gneg(A)) == O for A in L)
    n_pts = len(Ck)
    order_ok = all(gmul(n_pts, A) == O for A in random.sample(L, 40))
    orbit_ok = all(gadd(gadd(A, sigma_J(A)), sigma_J(sigma_J(A))) == O for A in L)
    sigma_hom = all(sigma_J(gadd(A, B)) == gadd(sigma_J(A), sigma_J(B)) for A, B in [tuple(random.choice(L) for _ in range(2)) for _ in range(300)])
    c1 = assoc and ident and inv_ok and order_ok and orbit_ok and sigma_hom
    ok &= c1
    print(f"C1 群法則: 結合律 {assoc}, 単位元 {ident}, 逆元 {inv_ok}, |J(k)|·P = O {order_ok}; σ_J は準同型 {sigma_hom}; 1+σ+σ^2 = 0 on J(k): {orbit_ok} -> {'OK' if c1 else 'NG'}")

    # C2 div(Y) 関係
    P = [normalize((x, ZERO, ONE)) for x in rk]
    lhs = gadd(gadd(P[0], P[1]), P[2]); rhs = gadd(inf1, inf2)
    c2 = lhs == rhs
    ok &= c2
    print(f"C2 Σ[P_i − O] = [∞_1 + ∞_2 − 2O]: {c2}")
    # ε_i の関係 λ[P_i − O] = [e], e = ∞_1 − ∞_0
    e_cls = gsub(inf1, O)
    c2b = all(lam(Pi) == e_cls for Pi in P)
    ok &= c2b
    print(f"   λ[P_i − O] = [∞_1 − ∞_0] for all i: {c2b}")

    # C3 J[λ^t] ∩ J(field)
    def lam_kernel_sizes(field_pts, tmax):
        sizes = []
        ker = {A for A in field_pts if lam(A) == O}
        sizes.append(len(ker))
        for t in range(2, tmax + 1):
            ker = {A for A in field_pts if lam(A) in ker}
            sizes.append(len(ker))
        return sizes
    sk = lam_kernel_sizes(Ck, 5); sK = lam_kernel_sizes(CK, 3)
    # 幾何的に |J[λ^t]| = 3^t (N-2 = 1 なので (R/λ^t)^1)
    c3 = sk[:3] == [3, 9, 27] and sK[0] == 3
    ok &= c3
    print(f"C3 |J(k)[λ^t]| (t=1..5) = {sk} (Lemma 4.3 の予想: t ≤ m=3 で 3^t = 3, 9, 27 が全て k 有理) -> {'OK' if sk[:3]==[3,9,27] else 'NG'}")
    print(f"   |J(K)[λ^t]| (t=1..3) = {sK} (K=F_7, r=0, m=1: t=1 のみ保証 (=3))")

    # C4 lift 鎖
    def preimages(target_set, field_pts):
        return {A for A in field_pts if lam(A) in target_set}
    c4 = True
    for i, Pi in enumerate(P):
        lv1 = preimages({Pi}, Ck)          # λ d = [P_i − O]  (t=2 ≤ m)
        lv2 = preimages(lv1, Ck)           # λ^2 d = [P_i − O] (t=3 ≤ m)
        lv3 = preimages(lv2, Ck)           # λ^3 d = [P_i − O] (t=4 > m, 保証なし)
        lv1K = preimages({Pi}, CK)
        good = len(lv1) == 3 and len(lv2) == 9
        c4 &= good
        print(f"C4 P_{i+1}: |λ^-1| = {len(lv1)} (期待 3), |λ^-2| = {len(lv2)} (期待 9), |λ^-3| = {len(lv3)} (保証なし; 幾何的には 27); K=F_7 上の |λ^-1| = {len(lv1K)} -> {'OK' if good else 'NG'}")
    ok &= c4
    # 全ての鎖の元 d_{i,3} について λ^3 d = [e]
    dm_all = set()
    for Pi in P:
        lv1 = preimages({Pi}, Ck); lv2 = preimages(lv1, Ck)
        dm_all |= lv2
    c4b = all(lam(lam(lam(d))) == e_cls for d in dm_all)
    ok &= c4b
    print(f"   全 d_(i,3) で λ^3 d = [e]: {c4b}")

    # C4' Prop 4.6 の前提: d_3 = Σ_i d_{i,3} (任意の選択) は λ^3 d_3 = [3e]. 全ての選択で確認.
    chains = []
    for Pi in P:
        lv1 = preimages({Pi}, Ck); lv2 = preimages(lv1, Ck)
        chains.append(list(lv2))
    cnt = 0; bad = 0
    for trip in itertools.product(*chains):
        d3 = gadd(gadd(trip[0], trip[1]), trip[2])
        if lam(lam(lam(d3))) != gmul(3, e_cls): bad += 1
        cnt += 1
    c4c = bad == 0
    ok &= c4c
    print(f"   Σ_i d_(i,3) の全 {cnt} 通りの選択で λ^3 d_3 = [3e]: {c4c}")
    # 3e の格子的性質: 3e = u λ^2 e. [3e] = λ^2[u e] in J でもあるので d_3 − [u e]... は J[λ^3] の元.
    # J 側では λ^3 は J(k̄) 上全射で d_3 は常に存在する. 格子側では W_3 = λ^{-2}(3e) ∉ λΛ. (A2 で確認済み)

    # C5 格子 (q=3,N=3): W_1 = 3e = (−3,3,0), W_2, W_3
    def lat_div(W):
        qq = len(W)
        S = [sum(W[1:j + 1]) for j in range(qq)]
        tot = sum(S)
        if tot % qq: return None
        c = -tot // qq
        return [c + s for s in S]
    W = [-3, 3, 0]; Ws = [W]
    while True:
        nxt = lat_div(Ws[-1])
        if nxt is None: break
        Ws.append(nxt)
    c5 = len(Ws) == 3
    ok &= c5
    print(f"C5 W_t = λ^-(t-1)(3e): {Ws}; W_3 ∉ λΛ (割れない): {c5}")
    print("ALL:", "OK" if ok else "NG")

if __name__ == "__main__":
    main()
