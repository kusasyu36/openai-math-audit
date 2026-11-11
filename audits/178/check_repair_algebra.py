"""
成果178 監査: repair.tex (定理 8.3) の 2x2 代数と期待値公式の検算。
- σF + 2P_tree = (1/(d-2))[[s, σd],[σd, s]], det = -1, 逆行列 = K_σ  (repair.tex:380-388)
- Tr(K_σ M_i) = -2 x_i^σ/(d-2)  (repair.tex:447-449)
- Tr(K_σ M_i K_σ M_i) > 0 for |λ| < d  (repair.tex:503-504)
- switch: A - A' = U F U^T  (repair.tex:359)
- E_{ν0}[z_i z_k^T] = δ_ik M_i / n  (repair.tex:439)  小グラフで総当り
- 傾斜後の平均 (8.68) を総当りで検証 (repair.tex:483-489)
- Schur/慣性 (8.66): 近傍 Schur 補行列 = diag(x) + U_E S_B^{-1} U_E^T  (repair.tex:396-413)
"""
import itertools, sys
import numpy as np
import sympy as sp
import networkx as nx

np.set_printoptions(precision=6, suppress=True, linewidth=120)
out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

# ---------- 1. 記号計算 ----------
d, lam, sig = sp.symbols('d lambda sigma', real=True)
q = d - 1
s = 2*sp.sqrt(q)
hc = 1/sp.sqrt(q)
mc = hc/(1 - hc**2)
F = sp.Matrix([[0, 1], [1, 0]])
Ptree = mc*sp.Matrix([[1, sig*hc], [sig*hc, 1]])
ideal = sig*F + 2*Ptree
ideal_claim = sp.Rational(1, 1)/(d-2)*sp.Matrix([[s, sig*d], [sig*d, s]])
K_claim = d/(d-2)*sp.Matrix([[-s/d, sig], [sig, -s/d]])
P("[1] m_c = sqrt(q)/(d-2) ?", sp.simplify(mc - sp.sqrt(q)/(d-2)) == 0)
P("[1] σF+2P_tree == (1/(d-2))[[s,σd],[σd,s]] ?", sp.simplify((ideal - ideal_claim).subs(sig**2, 1)) == sp.zeros(2, 2))
P("[1] det(ideal) == -1 ?", sp.simplify(ideal_claim.det().subs(sig**2, 1)) == -1)
prod = sp.simplify((ideal_claim*K_claim).subs(sig**2, 1))
P("[1] ideal * K_σ == I ?", prod == sp.eye(2))
Mi = sp.Matrix([[1, lam/d], [lam/d, 1]])
x = s - sig*lam
tr1 = sp.simplify((K_claim*Mi).trace())
P("[1] Tr(K_σ M_i) + 2x/(d-2) == 0 ?", sp.simplify(tr1 + 2*x/(d-2)) == 0)
tr2 = sp.simplify(((K_claim*Mi)**2).trace().subs(sig**2, 1))
P("[1] Tr(K M K M) =", sp.factor(tr2))
for dv in [3, 4, 5, 8]:
    vals = [tr2.subs({d: dv, lam: lv, sig: 1}) for lv in np.linspace(-dv+1e-6, dv-1e-6, 7)]
    P(f"    d={dv}: min Tr(KMKM) over |λ|<d =", float(min(vals)))

# ---------- 2. 小グラフでの総当り ----------
rng = np.random.default_rng(1)
for dv, n in [(3, 20), (4, 18), (5, 16)]:
    G = nx.random_regular_graph(dv, n, seed=int(rng.integers(1 << 30)))
    A = nx.to_numpy_array(G, nodelist=range(n))
    w, V = np.linalg.eigh(A)
    # 定数ベクトルを除く固有基底
    idx = [i for i in range(n) if abs(w[i] - dv) > 1e-8]
    assert len(idx) == n - 1, "connected?"
    lam_v = w[idx]; Vp = V[:, idx]
    sv = 2*np.sqrt(dv-1)
    oriented = [(u, v) for u in range(n) for v in G.neighbors(u)]
    assert len(oriented) == dv*n
    # E[z_i z_k^T]
    err = 0.0
    for i in range(n-1):
        for k in range(n-1):
            Z = np.zeros((2, 2))
            for (u, v) in oriented:
                zi = np.array([Vp[u, i], Vp[v, i]]); zk = np.array([Vp[u, k], Vp[v, k]])
                Z += np.outer(zi, zk)
            Z /= len(oriented)
            target = (np.array([[1, lam_v[i]/dv], [lam_v[i]/dv, 1]])/n) if i == k else np.zeros((2, 2))
            err = max(err, np.abs(Z - target).max())
    P(f"[2] d={dv} n={n}: max|E[z_i z_k^T] - δ M_i/n| =", f"{err:.2e}")
    # switch identity A - A' = U F U^T
    Fn = np.array([[0., 1.], [1., 0.]])
    ok = True
    for _ in range(20):
        (u, v), (xx, yy) = oriented[rng.integers(len(oriented))], oriented[rng.integers(len(oriented))]
        if len({u, v, xx, yy}) < 4 or A[u, yy] or A[xx, v]:
            continue
        A2 = A.copy()
        A2[u, v] = A2[v, u] = 0; A2[xx, yy] = A2[yy, xx] = 0
        A2[u, yy] = A2[yy, u] = 1; A2[xx, v] = A2[v, xx] = 1
        U = np.zeros((n, 2)); U[u, 0] = 1; U[xx, 0] = -1; U[v, 1] = 1; U[yy, 1] = -1
        ok &= np.allclose(A - A2, U @ Fn @ U.T)
    P(f"[2] d={dv}: A - A' = U F U^T for valid switches ?", ok)

    # (8.68) 傾斜平均の総当り検証 (σ=+1,-1; T_τ, K_τ を任意に与えて恒等式として確認)
    for sgm in [1, -1]:
        xs = sv - sgm*lam_v
        Ksig = dv/(dv-2)*np.array([[-sv/dv, sgm], [sgm, -sv/dv]])
        # near subspace E0: 固有値の slack が小さい方から 4 本 (恒等式検証なので任意でよい)
        E0 = list(np.argsort(xs)[:4])
        # T_τ: 各符号で slack 最小の 2 本
        T = {}
        for tau in [1, -1]:
            xt = sv - tau*lam_v
            T[tau] = list(np.argsort(xt)[:2])
        Bstar = 50.0; JT = 3.0
        pairs = [(g, j) for g in oriented for j in oriented]
        nu0 = 1.0/len(pairs)
        def Ktau(tau):
            return dv/(dv-2)*np.array([[-sv/dv, tau], [tau, -sv/dv]])
        dens = []
        for (g, j) in pairs:
            br = 0.0
            for tau in [1, -1]:
                for i in T[tau]:
                    zg = np.array([Vp[g[0], i], Vp[g[1], i]]); zj = np.array([Vp[j[0], i], Vp[j[1], i]])
                    br += zg @ Ktau(tau) @ zj
            dens.append(nu0*(1 + (1/Bstar)*(-(2*n/JT))*br))
        dens = np.array(dens)
        P(f"[2] d={dv} σ={sgm}: tilt density sums to 1 ? {dens.sum():.6f}, min density {dens.min():.3e} (>0 ?)")
        # E_ν H^σ (E0 x E0)
        H = np.zeros((len(E0), len(E0)))
        for (g, j), pr in zip(pairs, dens):
            zg = Vp[[g[0], g[1]], :][:, E0]; zj = Vp[[j[0], j[1]], :][:, E0]
            dz = zg - zj   # 2 x |E0|
            H += pr * (dz.T @ Ksig @ dz)
        # 式 (8.68)
        H68 = np.zeros_like(H)
        for a_, i in enumerate(E0):
            Mi_ = np.array([[1, lam_v[i]/dv], [lam_v[i]/dv, 1]])
            val = -4*xs[i]/((dv-2)*n)
            for tau in [1, -1]:
                if i in T[tau]:
                    val += 4/(Bstar*JT*n)*np.trace(Ksig @ Mi_ @ Ktau(tau) @ Mi_)
            H68[a_, a_] = val
        P(f"[2] d={dv} σ={sgm}: max|E_ν H - (8.68)| = {np.abs(H - H68).max():.2e}  (|E_ν H| max {np.abs(H).max():.2e})")

    # (8.66) Schur 補行列の恒等式: 2 本の頂点非交差 switch
    sgm = 1
    xs = sv - sgm*lam_v
    found = None
    for _ in range(2000):
        g1, j1, g2, j2 = (oriented[rng.integers(len(oriented))] for _ in range(4))
        verts = [*g1, *j1, *g2, *j2]
        if len(set(verts)) < 8: continue
        if A[g1[0], j1[1]] or A[j1[0], g1[1]] or A[g2[0], j2[1]] or A[j2[0], g2[1]]: continue
        found = (g1, j1, g2, j2); break
    (g1, j1, g2, j2) = found
    U = np.zeros((n, 4))
    for c, (g, j) in enumerate([(g1, j1), (g2, j2)]):
        U[g[0], 2*c] = 1; U[j[0], 2*c] = -1; U[g[1], 2*c+1] = 1; U[j[1], 2*c+1] = -1
    C = np.kron(np.eye(2), sgm*Fn)
    A2 = A - U @ C @ U.T
    assert np.allclose(A2, A2.T) and set(np.unique(A2)) <= {0.0, 1.0} and np.allclose(A2.sum(1), dv)
    # 固有基底座標
    Ue = Vp.T @ U   # (n-1) x 4
    E0 = list(np.argsort(xs)[:3]); Bk = [i for i in range(n-1) if i not in E0]
    Dfull = np.diag(xs) + Ue @ C @ Ue.T          # sI - σA' を固有基底で
    full_eig = np.linalg.eigvalsh(Dfull)
    direct = np.linalg.eigvalsh(sv*np.eye(n) - sgm*A2)
    kc = int(np.argmin(np.abs(direct - (sv - sgm*dv))))   # 定数ベクトルの固有値 s-σd を除く
    direct_np = np.delete(direct, kc)
    P(f"[2] d={dv}: spec(sI-σA') via eigenbasis vs direct, max diff = {np.abs(np.sort(full_eig) - np.sort(direct_np)).max():.2e} (constant mode {direct[kc]:.3f} removed)")
    DB = np.diag(xs[Bk]); UB = Ue[Bk, :]; UE = Ue[E0, :]
    RB = np.diag(1/xs[Bk])
    SB = np.linalg.inv(C) + UB.T @ RB @ UB
    bulk_new = DB + UB @ C @ UB.T
    schur_direct = np.diag(xs[E0]) + UE @ C @ UE.T - (UE @ C @ UB.T) @ np.linalg.solve(bulk_new, UB @ C @ UE.T)
    schur_66 = np.diag(xs[E0]) + UE @ np.linalg.solve(SB, UE.T)
    P(f"[2] d={dv}: near Schur complement (8.66) identity, max diff = {np.abs(schur_direct - schur_66).max():.2e}")
    def inertia(Mx):
        ev = np.linalg.eigvalsh(Mx); return (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()), int((abs(ev) <= 1e-9).sum()))
    P(f"[2] d={dv}: In(D_B)={inertia(DB)} In(-S_B)={inertia(-SB)} In(-C^-1)={inertia(-np.linalg.inv(C))} In(bulk_new)={inertia(bulk_new)}")
    lhs = tuple(np.add(inertia(DB), inertia(-SB))); rhs = tuple(np.add(inertia(-np.linalg.inv(C)), inertia(bulk_new)))
    P(f"[2] d={dv}: Haynsworth inertia identity holds ? {lhs == rhs}")

open(sys.argv[0].replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
