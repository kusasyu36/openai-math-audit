"""
成果178 監査: foundations.tex 補題 4.1 (configuration pencil) の恒等式を小さな部分ペアリングで検算。
- Q_j の漸化式 Q_j = A_e Q_{j-1} - (D_e - I) Q_{j-2} (j>=3), Q_1=A_e, Q_2=A_e^2-D_e   (foundations.tex:110-111)
- W = Σ (σh)^j Q_j = (1-h^2)[I - σhA_e + h^2(D_e-I)]^{-1} = M(h)^{-1}/m_h          (foundations.tex:113-114, 36-38)
- 空ペアリングで W = I on 1^⊥                                                        (foundations.tex:41-42)
- 1 組を解決したときの更新: 𝒟' = 𝒟 - UBU^T, A_e' = A_e + UFBU^T, M' = M + U(hI-σF)BU^T (foundations.tex:84-86)
- Woodbury: W' = W - W U H U^T W, H = (H_0^{-1} + D_p)^{-1}, [m_h(hI-σF)B]^{-1} + G_0 = H_0^{-1} (foundations.tex:88-93, 120-123)
- 行和: Q_j の各行和 = d q^{j-1} (j>=1)                                               (foundations.tex:105)
- 縮約恒等式 (4.2)                                                                    (foundations.tex:178-185)
"""
import sys
import numpy as np
np.set_printoptions(precision=5, suppress=True, linewidth=130)
out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

rng = np.random.default_rng(7)

def build(n, d, npairs):
    """n 頂点 d 本ずつのスタブ. npairs 組をランダムに解決(単純性は要求しない; 恒等式の検証なので多重辺も許す)."""
    stubs = [(v, k) for v in range(n) for k in range(d)]
    N = len(stubs)
    perm = rng.permutation(N)
    partner = -np.ones(N, dtype=int)
    for i in range(npairs):
        a, b = perm[2*i], perm[2*i+1]
        partner[a] = b; partner[b] = a
    return stubs, partner

def matrices(n, d, stubs, partner, h, sigma):
    N = len(stubs)
    S = [i for i in range(N) if partner[i] < 0]
    l = len(S)
    C = np.zeros((n, N))             # 全スタブの接続行列
    for i, (v, k) in enumerate(stubs): C[v, i] = 1
    CS = C[:, S]                     # 未解決スタブ
    ones_l = np.ones(l)
    Pi_l = np.eye(l) - np.outer(ones_l, ones_l)/l
    Xi = CS @ Pi_l
    Dcal = Xi @ Xi.T
    t = CS @ ones_l
    A = np.zeros((n, n))
    for i in range(N):
        j = partner[i]
        if j > i:
            A[stubs[i][0], stubs[j][0]] += 1; A[stubs[j][0], stubs[i][0]] += 1
    Ae = A + np.outer(t, t)/l
    # P_*: 解決スタブは相手と交換, 未解決スタブは平均化射影
    Pst = np.zeros((N, N))
    for i in range(N):
        j = partner[i]
        if j >= 0: Pst[i, j] = 1
    for i in S:
        for j in S: Pst[i, j] = 1/l
    De = C @ Pst @ Pst @ C.T
    assert np.allclose(De, d*np.eye(n) - Dcal), "D_e = dI - 𝒟 failed"
    q = d - 1
    z = 1/h + q*h; mh = h/(1-h**2)
    M = z*np.eye(n) - sigma*Ae - h*Dcal
    return dict(S=S, l=l, C=C, Xi=Xi, Dcal=Dcal, Ae=Ae, Pst=Pst, De=De, M=M, mh=mh, z=z, A=A, t=t)

def perp(n):
    ones = np.ones(n)/np.sqrt(n)
    Q, _ = np.linalg.qr(np.column_stack([ones, rng.standard_normal((n, n-1))]))
    return Q[:, 1:]   # 1^⊥ の正規直交基底

for (n, d, npairs) in [(8, 3, 5), (10, 4, 9), (12, 3, 10), (9, 5, 12)]:
    stubs, partner = build(n, d, npairs)
    h = 0.8/np.sqrt(d-1); sigma = rng.choice([1, -1])
    m = matrices(n, d, stubs, partner, h, sigma)
    Bp = perp(n)
    N = len(stubs); C = m['C']; Pst = m['Pst']; Ae = m['Ae']; De = m['De']
    # Q_j
    Q = [np.eye(n)]
    X = C @ Pst
    Y = (C.T @ C - np.eye(N)) @ Pst
    for j in range(1, 40):
        Q.append(X @ np.linalg.matrix_power(Y, j-1) @ C.T)
    ok1 = np.allclose(Q[1], Ae) and np.allclose(Q[2], Ae @ Ae - De)
    ok2 = all(np.allclose(Q[j], Ae @ Q[j-1] - (De - np.eye(n)) @ Q[j-2]) for j in range(3, 40))
    rows = all(np.allclose(Q[j].sum(1), d*(d-1)**(j-1)) for j in range(1, 10))
    P(f"[n={n} d={d} pairs={npairs} l={m['l']} σ={sigma}] Q1=A_e,Q2=A_e^2-D_e: {ok1}; recurrence j=3..39: {ok2}; row sums d q^(j-1): {rows}")
    # W の 3 つの表現 (1^⊥ 上)
    Wseries = sum((sigma*h)**j * Q[j] for j in range(40))
    Wpencil = (1-h**2)*np.linalg.inv(np.eye(n) - sigma*h*Ae + h**2*(De - np.eye(n)))
    Minv = np.linalg.inv(Bp.T @ m['M'] @ Bp)/m['mh']
    e_series = np.abs(Bp.T @ Wseries @ Bp - Bp.T @ Wpencil @ Bp).max()
    e_pencil = np.abs(Bp.T @ Wpencil @ Bp - Minv).max()
    P(f"    W: series(39 terms) vs pencil on 1^⊥: {e_series:.2e};  pencil vs M^-1/m_h on 1^⊥: {e_pencil:.2e}")
    # 露出更新
    S = m['S']; l = m['l']
    if l > 2:
        ia, ib = rng.choice(len(S), 2, replace=False)
        u, v = S[ia], S[ib]
        partner2 = partner.copy(); partner2[u] = v; partner2[v] = u
        m2 = matrices(n, d, stubs, partner2, h, sigma)
        Xi = m['Xi']; U = Xi[:, [ia, ib]]
        J2 = np.ones((2, 2)); G0 = np.eye(2) - J2/l; B = np.eye(2) + J2/(l-2)
        F = np.array([[0., 1.], [1., 0.]])
        assert np.allclose(np.linalg.inv(G0), B)
        okD = np.allclose(m2['Dcal'], m['Dcal'] - U @ B @ U.T)
        okA = np.allclose(m2['Ae'], m['Ae'] + U @ F @ B @ U.T)
        okM = np.allclose(m2['M'], m['M'] + U @ (h*np.eye(2) - sigma*F) @ B @ U.T)
        P(f"    exposure: 𝒟'={okD}, A_e'={okA}, M'={okM}")
        # Woodbury
        Wb = np.linalg.inv(Bp.T @ m['M'] @ Bp)/m['mh']      # 1^⊥ 座標
        W2b = np.linalg.inv(Bp.T @ m2['M'] @ Bp)/m2['mh']
        Ub = Bp.T @ U
        H0 = -sigma*h*F @ B
        K = (h*np.eye(2) - sigma*F) @ B
        okH0 = np.allclose(np.linalg.inv(m['mh']*K) + G0, np.linalg.inv(H0))
        G = Ub.T @ Wb @ Ub; Dp = G - G0
        H = np.linalg.inv(np.linalg.inv(H0) + Dp)
        okW = np.allclose(W2b, Wb - Wb @ Ub @ H @ Ub.T @ Wb)
        P(f"    [m_h K]^-1 + G_0 = H_0^-1: {okH0};  Woodbury W' = W - WUHU^TW: {okW};  ||D_p|| = {np.linalg.norm(Dp):.3f}")
        # 慣性保存: M>0 なら M'>0 (||D_p|| small の場合)
        ev1 = np.linalg.eigvalsh(Bp.T @ m['M'] @ Bp).min(); ev2 = np.linalg.eigvalsh(Bp.T @ m2['M'] @ Bp).min()
        P(f"    min eig M on 1^⊥ = {ev1:.4f}, after exposure = {ev2:.4f}")

# 空ペアリング: W = I on 1^⊥
n, d = 10, 3
stubs, partner = build(n, d, 0)
h = 0.7/np.sqrt(d-1)
m = matrices(n, d, stubs, partner, h, 1)
Bp = perp(n)
Wb = np.linalg.inv(Bp.T @ m['M'] @ Bp)/m['mh']
P(f"[empty pairing n={n} d={d}] max|W - I| on 1^⊥ = {np.abs(Wb - np.eye(n-1)).max():.2e}")

# 縮約恒等式 (4.2): I ⊆ S, |I| = l-k, a,c 任意配列
for trial in range(3):
    l = int(rng.integers(6, 12)); k = int(rng.integers(0, 3)); sigma = rng.choice([1, -1]); h = 0.5
    I = list(range(l - k))
    a = rng.standard_normal(l); c = rng.standard_normal(l)
    F = np.array([[0., 1.], [1., 0.]]); J2 = np.ones((2, 2)); B = np.eye(2) + J2/(l-2); H0 = -sigma*h*F @ B
    tot = 0.0; cnt = 0
    for u in I:
        for v in I:
            if u == v: continue
            ap = np.array([a[u], a[v]]); cp = np.array([c[u], c[v]])
            tot += ap @ H0 @ cp; cnt += 1
    lhs = (tot/cnt)/(-2*sigma*h)
    sa = a[I].sum(); sc = c[I].sum(); sac = (a[I]*c[I]).sum(); m_ = len(I)
    rhs = (1 + 1/(l-2))*sa*sc/(m_*(m_-1)) - k*sac/((l-2)*m_*(m_-1))
    P(f"[contraction (4.2)] l={l} k={k}: lhs={lhs:.6f} rhs={rhs:.6f} diff={abs(lhs-rhs):.2e}")

open(sys.argv[0].replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
