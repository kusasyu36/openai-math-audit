"""Numerical test of the completion-stage inequalities (completion.tex, phase-cover.tex)
on EXACT complex Hadamard matrices of order 6 with seed centers perturbed inside the
stated seed radius.  Every inequality below is claimed to be *necessary* for an exact
completion; so it must hold for these exact matrices.  Independent numpy code."""
import numpy as np, itertools
rng = np.random.default_rng(2026)
TWO_PI = 2*np.pi

def fourier_family(a, b):
    F = np.exp(2j*np.pi*np.outer(np.arange(6), np.arange(6))/6)
    R = np.zeros((6,6)); 
    for r in (1,3,5): R[r] = [0, a, b, 0, a, b]
    return F*np.exp(1j*R)
def dita(c, d=None):
    # Dita-type construction [[P, Q],[P, -Q]] with P = F_3 and Q = diag(1,e^{ic},e^{id}) F_3 (Hadamard for all c,d)
    if d is None: d = c*0.618033988749895 + 0.3
    w3 = np.exp(2j*np.pi/3); F3 = np.array([[w3**(j*k) for k in range(3)] for j in range(3)])
    Q = np.diag([1, np.exp(1j*c), np.exp(1j*d)])@F3
    return np.block([[F3, Q],[F3, -Q]])
def bjorck():
    d = (1-np.sqrt(3))/2 + 1j*np.sqrt(np.sqrt(3)/2)
    first = np.array([1, 1j*d, -d, -1j, -np.conj(d), 1j*np.conj(d)])
    return np.array([np.roll(first, k) for k in range(6)])
def cubic():
    w = np.exp(2j*np.pi/3)
    def t(i,j):
        if i==0 or j==0 or i==j: return 0
        return 1 if (i-j)%5 in (1,4) else -1
    return np.array([[w**t(i,j) for j in range(6)] for i in range(6)])
def is_hadamard(H, tol=1e-9):
    return np.allclose(np.abs(H), 1, atol=tol) and np.allclose(H@H.conj().T, 6*np.eye(6), atol=tol)
def dephase(H):
    H = H / H[:, [0]]; H = H / H[[0], :]
    return H
def random_equivalent(H):
    P = rng.permutation(6); Q = rng.permutation(6)
    H = H[P][:, Q]
    if rng.random() < .5: H = H.T
    if rng.random() < .5: H = H.conj()
    return dephase(H)

def generators():
    yield "F6(a,b)", lambda: fourier_family(rng.uniform(0, TWO_PI), rng.uniform(0, TWO_PI))
    yield "D6(c)", lambda: dita(rng.uniform(0, TWO_PI))
    yield "C6 Bjorck", bjorck
    yield "S6 cubic", cubic
    yield "F6^T(a,b)", lambda: fourier_family(rng.uniform(0, TWO_PI), rng.uniform(0, TWO_PI)).T

# ---------- completion.tex machinery (ideal real arithmetic in double precision) ----------
r_s = .0320001; R_s = TWO_PI*r_s
def seed_blocks(Hc):
    a = Hc[1,1]; b = Hc[1,2:]; c = Hc[2:,1]
    E0 = np.array([[1,1],[1,a]]); B0 = np.vstack([np.ones(4), b]); C0 = np.column_stack([np.ones(4), c])
    return a, b, c, E0, B0, C0
def frames(B0, C0):
    GB = B0.conj().T @ np.linalg.inv(B0 @ B0.conj().T)        # G_B = B0^*(B0B0^*)^{-1}, 4x2
    GC = C0 @ np.linalg.inv(C0.conj().T @ C0)                 # G_C = C0 (C0^*C0)^{-1}, 4x2
    # orthonormal bases of ker B0 (4x2) and ker C0^* (4x2)
    def kernel(M):
        u, s, vh = np.linalg.svd(M); return vh[2:].conj().T
    KB = kernel(B0); KC = kernel(C0.conj().T)
    mB = np.linalg.svd(B0, compute_uv=False).min(); mC = np.linalg.svd(C0, compute_uv=False).min()
    return GB, GC, KB, KC, mB, mC
def polar(Y):
    u, s, vh = np.linalg.svd(Y); return u@vh
def u2_params(V):
    det = np.linalg.det(V); alpha = np.angle(det)/2
    if alpha < 0: alpha += np.pi          # alpha in [0, pi) radians  (cycles [0,1/2])
    W = np.exp(-1j*alpha)*V
    z, w_ = W[0,0], W[0,1]
    T = np.arctan2(abs(w_), abs(z))        # in [0, pi/2]  (t in [0,1/4] cycles)
    beta = np.angle(z) if abs(z) > 1e-14 else 0.0
    gamma = np.angle(w_) if abs(w_) > 1e-14 else 0.0
    return alpha, T, beta % TWO_PI, gamma % TWO_PI
def Umat(al, T, be, ga):
    return np.exp(1j*al)*np.array([[np.exp(1j*be)*np.cos(T), np.exp(1j*ga)*np.sin(T)],[-np.exp(-1j*ga)*np.sin(T), np.exp(-1j*be)*np.cos(T)]])
def dU(al,T,be,ga,th):
    x,y,z,wv = th; U = Umat(al,T,be,ga)
    dA = 1j*U
    dT = np.exp(1j*al)*np.array([[-np.exp(1j*be)*np.sin(T), np.exp(1j*ga)*np.cos(T)],[-np.exp(-1j*ga)*np.cos(T), -np.exp(-1j*be)*np.sin(T)]])
    dB = np.array([[1j*U[0,0],0],[0,-1j*U[1,1]]]); dG = np.array([[0,1j*U[0,1]],[-1j*U[1,0],0]])
    return x*dA + y*dT + z*dB + wv*dG

results = {}
def rec(name, ok, val=None):
    results.setdefault(name, [0,0,[]]); results[name][0] += ok; results[name][1] += 1
    if val is not None: results[name][2].append(val)

ntests = 0
for fam, gen in generators():
    for trial in range(60 if "a,b" in fam or "(c)" in fam else 25):
        H = gen()
        assert is_hadamard(H), fam
        Hx = random_equivalent(H)                 # exact dephased Hadamard H'
        if not is_hadamard(Hx): continue
        # exact seed phases (cycles) and a perturbed center inside the seed ball
        q_exact = np.angle(Hx)/TWO_PI
        noise = rng.standard_normal(9); noise *= rng.uniform(0, 1)*r_s/np.linalg.norm(noise)
        Hc = Hx.copy()
        pos = [(1,1),(1,2),(1,3),(1,4),(1,5),(2,1),(3,1),(4,1),(5,1)]
        for (i,j), dn in zip(pos, noise): Hc[i,j] = np.exp(TWO_PI*1j*(q_exact[i,j] + dn))
        a, b, c, E0, B0, C0 = seed_blocks(Hc)
        GB, GC, KB, KC, mB, mC = frames(B0, C0)
        m = min(mB, mC)
        if m <= 1.24: rec("frame check m>1.24 (rounding.tex:374)", False); continue
        rec("frame check m>1.24 (rounding.tex:374)", True)
        # seed displacement s in radians (shortest angles)
        s = np.array([np.angle(Hx[i,j]/Hc[i,j]) for (i,j) in pos])
        assert np.linalg.norm(s) <= R_s + 1e-12
        E_, B_, C_, D_ = Hx[:2,:2], Hx[:2,2:], Hx[2:,:2], Hx[2:,2:]
        dE, dB, dC = E_-E0, B_-B0, C_-C0
        e1, b1, c1 = map(np.linalg.norm, (dE, dB, dC))
        rec("seed-norm: e1^2+b1^2+c1^2 <= ||s||^2 (completion:seed-norm)", e1**2+b1**2+c1**2 <= np.linalg.norm(s)**2 + 1e-12)
        A0 = np.sqrt(2+abs(1+a)); Astar = min(2, A0+R_s)
        rec("||E'||_op <= A_*", np.linalg.norm(E_, 2) <= Astar + 1e-12)
        PK = KB@KB.conj().T; PB = np.eye(4)-PK; PC = np.eye(4)-KC@KC.conj().T
        rec("P_B = G_B B_0 and P_C = G_C C_0^*", np.allclose(PB, GB@B0) and np.allclose(PC, GC@C0.conj().T))
        rec("||G_B||=1/m_B, ||G_C||=1/m_C", np.isclose(np.linalg.norm(GB,2), 1/mB) and np.isclose(np.linalg.norm(GC,2), 1/mC))
        # Lemma polar
        Y = KC.conj().T @ D_ @ KB
        lhs = 6*np.eye(2) - Y.conj().T@Y
        rhs = KB.conj().T@B_.conj().T@B_@KB + KB.conj().T@D_.conj().T@PC@D_@KB
        rec("deficit identity (completion:deficit-identity)", np.allclose(lhs, rhs))
        q = R_s**2*(1+(Astar**2+6)/mC**2)
        tr = np.trace(lhs).real
        rec("tr(6I-Y*Y) <= b1^2 + (||E'|| b1 + sqrt6 c1)^2/mC^2 <= q", tr <= b1**2 + (np.linalg.norm(E_,2)*b1+np.sqrt(6)*c1)**2/mC**2 + 1e-12 and tr <= q + 1e-12, q)
        rec("q < .32 (rounding.tex:466)", q < .32)
        V = polar(Y); p = q/(np.sqrt(6)+np.sqrt(6-q))
        rec("||Y - sqrt6 V||_F <= p (completion:polar-bound)", np.linalg.norm(Y-np.sqrt(6)*V) <= p + 1e-12)
        # unitary cell at N_u = 64
        al, T, be, ga = u2_params(V)
        rec("U(2) parametrization reproduces V", np.allclose(Umat(al,T,be,ga), V))
        Nu = 64; h = np.pi/Nu
        # cell midpoints: alpha in [0,1/2) cycles with Nu/2 bins; t in [0,1/4] with Nu/4 bins; beta,gamma in [0,1) with Nu bins
        def mid(x_rad, nb, full):   # x in radians; bins of width full/nb (in radians)
            wdt = full/nb; k = min(int(x_rad//wdt), nb-1); return (k+.5)*wdt
        Um = (mid(al, Nu//2, np.pi), mid(T, Nu//4, np.pi/2), mid(be, Nu, TWO_PI), mid(ga, Nu, TWO_PI))
        th = np.array([al,T,be,ga]) - np.array(Um)
        assert np.all(np.abs(th) <= h + 1e-12)
        U = Umat(*Um)
        rec("||V-U||_F <= sqrt6 h", np.linalg.norm(V-U) <= np.sqrt(6)*h + 1e-12)
        rec("||V-U-dU[theta]||_F <= sqrt42/2 h^2", np.linalg.norm(V-U-dU(*Um,th)) <= np.sqrt(42)/2*h**2 + 1e-12)
        # comparison block
        Dspan = -C0@E0.conj().T@GB.conj().T; Dker = np.sqrt(6)*KC@U@KB.conj().T; Dcmp = Dspan + Dker
        Z = D_ - Dcmp
        pieces = [Z@PB, PC@Z@PK, (np.eye(4)-PC)@Z@PK]
        rec("pieces orthogonal and sum to Z", np.allclose(sum(pieces), Z) and all(abs(np.vdot(pieces[i], pieces[j])) < 1e-9 for i in range(3) for j in range(i)))
        rec("piece-B identity", np.allclose(pieces[0], (-dC@E0.conj().T - C_@dE.conj().T - D_@dB.conj().T)@GB.conj().T))
        rec("piece-C identity", np.allclose(pieces[1], GC@(-E_.conj().T@dB - dC.conj().T@D_)@PK))
        rec("piece-K identity", np.allclose(pieces[2], KC@(Y-np.sqrt(6)*U)@KB.conj().T))
        z = np.sqrt((R_s*(np.sqrt(6)+Astar)/m)**2 + (6*h+p)**2)
        rec("||Z||_F <= z (completion:z-bound)", np.linalg.norm(Z) <= z + 1e-12, np.linalg.norm(Z)/z)
        L = max(np.sqrt(6), np.linalg.norm(C0,2)*np.sqrt(1+(A0/mB)**2))
        rec("||[C0 D]||_op <= L and ||D P_K||=sqrt6 and ||D|| <= max(sqrt6, ||C0|| A0/mB)", np.linalg.norm(np.hstack([C0,Dcmp]),2) <= L+1e-12 and np.isclose(np.linalg.norm(Dcmp@PK,2), np.sqrt(6)) and np.linalg.norm(Dcmp,2) <= max(np.sqrt(6), np.linalg.norm(C0,2)*A0/mB)+1e-12)
        rec("L < 6 and ||D|| < 4.6 (rounding.tex:442,461)", L < 6 and np.linalg.norm(Dcmp,2) < 4.6)
        # first variations and linear map J(s,theta)
        E1 = np.zeros((2,2),complex); B1 = np.zeros((2,4),complex); C1 = np.zeros((4,2),complex)
        E1[1,1] = 1j*a*s[0]; B1[1,:] = 1j*b*s[1:5]; C1[:,1] = 1j*c*s[5:9]
        J = (-C1@E0.conj().T - C0@E1.conj().T - Dcmp@B1.conj().T)@GB.conj().T + GC@(-E0.conj().T@B1 - C1.conj().T@Dcmp)@PK + np.sqrt(6)*KC@dU(*Um,th)@KB.conj().T
        e_off = R_s**2*(L+A0+1)/(2*m) + z*R_s/m
        Er = np.sqrt(e_off**2 + (np.sqrt(6*42)/2*h**2 + p)**2)
        rec("||Z - J(s,theta)||_F <= E_r (completion:linear-remainder-bound)", np.linalg.norm(Z-J) <= Er + 1e-12, np.linalg.norm(Z-J)/Er)
        rec("E_r < 1 (rounding.tex:466)", Er < 1)
        # entrywise coefficient formulas (completion.tex:506-529) vs matrix J
        Jc = np.zeros((4,4),complex)
        for i in range(4):
            for j in range(4):
                Jc[i,j] += -c[i]*np.conj(1j*a*GB[j,1])*s[0]
                for k in range(4):
                    Jc[i,j] += (-(GC[i,0]+GC[i,1]*np.conj(a))*1j*b[k]*PK[k,j] - Dcmp[i,k]*np.conj(1j*b[k]*GB[j,1]))*s[1+k]
                    Jc[i,j] += (-(1 if i==k else 0)*1j*c[k]*(np.conj(GB[j,0])+np.conj(a*GB[j,1])) - GC[i,1]*np.conj(1j*c[k])*Dker[k,j])*s[5+k]
        Jc += np.sqrt(6)*KC@dU(*Um,th)@KB.conj().T
        rec("entrywise coefficient formulas == matrix J", np.allclose(Jc, J))
        # norm exclusion test (necessary)
        rec("sum (|D_ij|-1)^2 <= z^2 (completion:norm-test)", ((np.abs(Dcmp)-1)**2).sum() <= z**2 + 1e-12)
        # directional test for random y (necessary inequality)
        cX = (1+a+b.sum())/1j; cY = (1+a+c.sum())/1j
        c_all = np.concatenate([[cX.real, cX.imag, cY.real, cY.imag], ((np.abs(Dcmp)**2-1)/2).ravel()])
        # linear forms: J_X s = a s_a + sum b_j s_bj ; J_Y s = a s_a + sum c_i s_ci ; flat: Re(conj(D_ij) J_ij(s,theta))
        def Jall_apply(sv, thv):
            JX = a*sv[0] + (b*sv[1:5]).sum(); JY = a*sv[0] + (c*sv[5:9]).sum()
            E1 = np.zeros((2,2),complex); B1 = np.zeros((2,4),complex); C1 = np.zeros((4,2),complex)
            E1[1,1] = 1j*a*sv[0]; B1[1,:] = 1j*b*sv[1:5]; C1[:,1] = 1j*c*sv[5:9]
            Jm = (-C1@E0.conj().T - C0@E1.conj().T - Dcmp@B1.conj().T)@GB.conj().T + GC@(-E0.conj().T@B1 - C1.conj().T@Dcmp)@PK + np.sqrt(6)*KC@dU(*Um,thv)@KB.conj().T
            return np.concatenate([[JX.real, JX.imag, JY.real, JY.imag], (np.conj(Dcmp)*Jm).real.ravel()])
        Jall = np.column_stack([Jall_apply(np.eye(9)[k], np.zeros(4)) for k in range(9)] + [Jall_apply(np.zeros(9), np.eye(4)[k]) for k in range(4)])
        # sanity: exact flatness: 0 = c_ij + J_ij + Re(conj D R) + |Z|^2/2
        R_ = Z - J
        flat_check = c_all[4:] + (np.conj(Dcmp)*J).real.ravel() + (np.conj(Dcmp)*R_).real.ravel() + 0.5*(np.abs(Z)**2).ravel()
        rec("flatness expansion identity (completion.tex:607-613)", np.allclose(flat_check, 0))
        rec("seed sum Taylor: |rho_X|,|rho_Y| <= R_s^2/2", abs(cX + a*s[0]+(b*s[1:5]).sum()) <= R_s**2/2+1e-12 and abs(cY + a*s[0]+(c*s[5:9]).sum()) <= R_s**2/2+1e-12)
        okdir = True
        for _ in range(20):
            y = rng.standard_normal(20)
            if rng.random() < .5: y = np.sign(y)
            g = Jall.T@y
            rhs = R_s*np.linalg.norm(g[:9]) + h*np.abs(g[9:]).sum() + Er*np.linalg.norm(Dcmp.ravel()*y[4:]) + .5*z**2*max(0, (-y[4:]).max()) + .5*R_s**2*(np.linalg.norm(y[:2])+np.linalg.norm(y[2:4]))
            okdir &= (y@c_all <= rhs + 1e-10)
        rec("directional test (completion:directional-test) holds for random y", okdir)
        # phase-cover Lemma arc-bounds
        mj = np.abs(Dcmp).ravel(); slack = max(0, z**2 - ((mj-1)**2).sum()); mi = mj.min(); qd = slack/mi
        if qd < 3.8 and mi > .15:
            Dcirc = (Dcmp/np.abs(Dcmp)).ravel(); d = np.angle(D_.ravel()/Dcirc)
            arc = np.sqrt(qd)/np.sqrt(1-qd/4); weighted = np.sqrt(slack)/np.sqrt(1-qd/4)
            rec("arc bounds: ||d|| <= arc, (sum m_j d_j^2)^1/2 <= weighted (phase:arc-weighted)", np.linalg.norm(d) <= arc+1e-12 and np.sqrt((mj*d**2).sum()) <= weighted+1e-12)
            rec("chord identity |D'-D|^2 = (m-1)^2 + m|D'-D°|^2", np.allclose(np.abs(D_.ravel()-Dcmp.ravel())**2, (mj-1)**2 + mj*np.abs(D_.ravel()-Dcirc)**2))
            # raw support representation: u_D = sin d + ... check  Im(conj(D°) Z) = sin d  and |d - sin d| sum <= arc^3/6
            rec("Im(conj(D°)Z_j) = sin d_j and sum|d-sin d| <= arc^3/6", np.allclose((np.conj(Dcirc)*Z.ravel()).imag, np.sin(d)) and np.abs(d-np.sin(d)).sum() <= arc**3/6+1e-12)
        else:
            rec("arc-bound preconditions (qd<3.8, mi>.15) met", False)
        ntests += 1

print(f"exact Hadamard test cases: {ntests}")
allok = True
for k, (okc, tot, vals) in results.items():
    extra = f"  (max observed ratio/value {max(vals):.4g})" if vals else ""
    print(("OK   " if okc == tot else "FAIL ") + f"{k}: {okc}/{tot}{extra}")
    allok &= okc == tot
print("ALL OK" if allok else "SOME FAILURES")
