"""Independent checks: formulas hand-transcribed from the C++ source ("SRC") versus quantities
derived from the TeX ("TEX").  No OpenAI code is compiled or run; SRC blocks are my own Python
transcriptions of what the source lines compute."""
import numpy as np, math, itertools, random
from fractions import Fraction
import mpmath as mp
rng = np.random.default_rng(20261111)
PI = float.fromhex("0x1.921fb54442d18p+1")
ok_all = True
def report(name, val, tol, extra=""):
    global ok_all
    good = bool(val <= tol)
    ok_all &= good
    print(f"[{'OK' if good else 'FAIL'}] {name}: {val:.3e} (tol {tol:.1e}) {extra}")
e = lambda t: np.exp(2j*np.pi*np.asarray(t))

# ---------- 1. expair (SRC seed.cpp:21-31) vs exp(2 pi i x) with mpmath --------------------
def expair(x):
    if not (abs(x) <= 16.): return (math.nan, math.nan)
    k = int(round(4*x))  # nearbyint, ties-to-even == Python round
    y = (x - .25*k)*(2*PI); yy = -y*y
    a = 1.; b = y; ta = 1.; tb = y
    for j in range(1, 10):
        ta = (ta*yy)/((2*j-1.)*(2*j)); tb = (tb*yy)/((2*j)*(2*j+1.))
        a += ta; b += tb
    return [(a, b), (-b, a), (-a, -b), (b, -a)][(k % 4 + 4) % 4]
mp.mp.dps = 40
worst = 0.
xs = list(rng.uniform(-16, 16, 4000)) + [k/8 for k in range(-128, 129)] + [(i+.5)/112 for i in range(112)]
for x in xs:
    c, s = expair(float(x))
    z = mp.e**(2j*mp.pi*mp.mpf(float(x)))
    worst = max(worst, float(abs(mp.mpc(c, s) - z)))
report("1 expair vs e(x), |x|<=16 (tex bound 2e-14, rounding.tex:187-189)", worst, 2e-14)

# ---------- 2. frames (SRC ub.cpp:37-68) vs TEX completion.tex:55-95 ----------------------
def frame(p):
    S = p.sum(); m = math.sqrt(4-abs(S)); top = math.sqrt(4+abs(S))
    inv = 1./(16-abs(S)**2)
    G = np.stack([(4 - p*np.conj(S))*inv, (4*p - S)*inv], axis=1)      # G[j][0], G[j][1]
    basis = [np.full(4, .5+0j), (p - S*.25)/math.sqrt(4.-abs(S)**2/4)]
    Q = np.zeros((4, 2), complex)
    for k in (2, 3):
        best = -1; q = None
        for j in range(4):
            s = 1 - sum(abs(b[j])**2 for b in basis)
            if s < best: continue
            best = s
            q = np.array([(1. if i == j else 0.) - sum(b[i]*np.conj(b[j]) for b in basis) for i in range(4)])
        r = (abs(q)**2).sum()
        if not r > .45*.45: return None
        Q[:, k-2] = q/math.sqrt(r); basis.append(Q[:, k-2].copy())
    P = Q @ Q.conj().T
    return dict(G=G, Q=Q, P=P, m=m, top=top)
def uc(alpha, t, beta, gamma):
    """SRC ub.cpp:159-174 (U, dU) at arbitrary cycle parameters; dU w.r.t. radians."""
    c, s = math.cos(2*np.pi*t), math.sin(2*np.pi*t)
    diagP, diagM, offP, offM = e(alpha+beta), e(alpha-beta), e(alpha+gamma), -e(alpha-gamma)
    U = np.array([diagP*c, offP*s, offM*s, diagM*c])
    dU = [1j*U, np.array([diagP*(-s), offP*c, offM*c, diagM*(-s)]),
          np.array([1j*U[0], 0, 0, -1j*U[3]]), np.array([0, 1j*U[1], -1j*U[2], 0])]
    return U, dU
def Vtex(alpha, t, beta, gamma):
    """TEX completion.tex:211-215."""
    c, s = math.cos(2*np.pi*t), math.sin(2*np.pi*t)
    return e(alpha)*np.array([[e(beta)*c, e(gamma)*s], [-e(-gamma)*s, e(-beta)*c]])
errG = errK = errU = errdU = errD = errJ = 0.; worst1 = worst2 = 0.; mmin = 9; cnt = 0
errPre = 0.
while cnt < 200:
    ph = rng.uniform(0, 1, 9)
    a = e(ph[0]); B = e(ph[1:5]); C = e(ph[5:9])
    FB = frame(np.conj(B)); FC = frame(C)
    if FB is None or FC is None or not (FB['m'] > 1.25 and FC['m'] > 1.25): continue
    cnt += 1
    B0 = np.array([np.ones(4), B]); C0 = np.stack([np.ones(4), C], axis=1); E0 = np.array([[1, 1], [1, a]])
    GB_tex = B0.conj().T @ np.linalg.inv(B0 @ B0.conj().T); GC_tex = C0 @ np.linalg.inv(C0.conj().T @ C0)
    errG = max(errG, abs(FB['G']-GB_tex).max(), abs(FC['G']-GC_tex).max())
    KB, KC = FB['Q'], FC['Q']
    errK = max(errK, abs(B0 @ KB).max(), abs(C0.conj().T @ KC).max(), abs(KB.conj().T@KB-np.eye(2)).max(), abs(KC.conj().T@KC-np.eye(2)).max())
    mmin = min(mmin, abs(FB['m']-np.linalg.svd(B0, compute_uv=False).min()), 9)
    errG = max(errG, abs(FB['m']-np.linalg.svd(B0, compute_uv=False).min()), abs(FC['top']-np.linalg.svd(C0, compute_uv=False).max()))
    par = rng.uniform(0, 1, 4)*np.array([.5, .25, 1, 1])
    U, dU = uc(*par)
    errU = max(errU, abs(U.reshape(2, 2)-Vtex(*par)).max())
    hfd = 1e-6
    for k in range(4):
        d = np.zeros(4); d[k] = hfd/(2*np.pi)          # radians -> cycles
        fd = (Vtex(*(par+d))-Vtex(*(par-d)))/(2*hfd)
        errdU = max(errdU, abs(fd-dU[k].reshape(2, 2)).max())
    # Taylor bounds (TEX completion.tex:251-258) on a random cell of side 1/N
    N = int(rng.choice([16, 32, 64])); h = np.pi/N
    th = rng.uniform(-h, h, 4); th[rng.integers(4)] = rng.choice([-h, h])
    V = Vtex(*(par+th/(2*np.pi)))
    lin = sum(th[k]*dU[k] for k in range(4)).reshape(2, 2)
    worst1 = max(worst1, np.linalg.norm(V-U.reshape(2, 2))/(math.sqrt(6)*h))
    worst2 = max(worst2, np.linalg.norm(V-U.reshape(2, 2)-lin)/(math.sqrt(42)/2*h*h))
    # SRC block::setup (ub.cpp:104-128) + makeD/makeJ (ub.cpp:179-199)
    s6 = math.sqrt(6.)
    base = np.zeros(16, complex); W = np.zeros((16, 4), complex); fixed = np.zeros((9, 16), complex)
    bm = np.zeros((4, 4), complex); cm = np.zeros((4, 4), complex)
    for i in range(4):
        for j in range(4):
            u = 4*i+j
            base[u] = -(1+C[i])*np.conj(FB['G'][j, 0]) - (1+C[i]*np.conj(a))*np.conj(FB['G'][j, 1])
            for l in range(4): W[u, l] = FC['Q'][i, l//2]*np.conj(FB['Q'][j, l % 2])*s6
            fixed[0, u] = -C[i]*np.conj(1j*a*FB['G'][j, 1])
            for k in range(4):
                fixed[5+k, u] = -1j*C[k]*(np.conj(FB['G'][j, 0])+np.conj(a*FB['G'][j, 1])) if i == k else 0
                fixed[1+k, u] = -(FC['G'][i, 0]+FC['G'][i, 1]*np.conj(a))*1j*B[k]*FB['P'][k, j]
    for j in range(4):
        for k in range(4):
            bm[k, j] = -np.conj(1j*B[k]*FB['G'][j, 1]); cm[j, k] = -FC['G'][j, 1]*np.conj(1j*C[k])
    Kc = W @ U; Dc = Kc + base
    jac = np.zeros((13, 16), complex)
    for i in range(4):
        for j in range(4):
            hh = 4*i+j
            jac[0, hh] = fixed[0, hh]
            for k in range(4):
                jac[1+k, hh] = fixed[1+k, hh] + Dc[4*i+k]*bm[k, j]
                jac[5+k, hh] = fixed[5+k, hh] + Kc[4*k+j]*cm[i, k]
                jac[9+k, hh] = (W[hh] * dU[k]).sum()
    # TEX: D = -C0 E0^* GB^* + sqrt6 KC U KB^*   (completion.tex:305-307)
    D_tex = -C0 @ E0.conj().T @ GB_tex.conj().T + s6*KC @ U.reshape(2, 2) @ KB.conj().T
    errD = max(errD, abs(D_tex.reshape(16)-Dc).max())
    # TEX: J(s,theta) (completion.tex:417-425)
    PK = KB @ KB.conj().T
    for trial in range(3):
        s = rng.normal(size=9); th = rng.normal(size=4)
        E1 = np.array([[0, 0], [0, 1j*a*s[0]]]); B1 = np.array([np.zeros(4), 1j*B*s[1:5]]); C1 = np.stack([np.zeros(4), 1j*C*s[5:9]], axis=1)
        dUth = sum(th[k]*dU[k] for k in range(4)).reshape(2, 2)
        Jtex = (-C1@E0.conj().T - C0@E1.conj().T - D_tex@B1.conj().T)@GB_tex.conj().T \
               + GC_tex@(-E0.conj().T@B1 - C1.conj().T@D_tex)@PK + s6*KC@dUth@KB.conj().T
        Jsrc = (np.concatenate([s, th])[:, None]*jac).sum(0)
        errJ = max(errJ, abs(Jtex.reshape(16)-Jsrc).max())
    # pre-equations (ub.cpp:118-127) vs TEX completion.tex:557-562
    SX = 1+a+B.sum(); cX = SX/1j
    errPre = max(errPre, abs(cX.real-SX.imag), abs(cX.imag+SX.real))
report("2a frame G vs B0^*(B0B0^*)^-1, C0(C0^*C0)^-1, m, top", errG, 1e-11)
report("2b K_B in ker B0, K_C in ker C0^*, orthonormal", errK, 1e-11)
report("3a uc::init U vs tex parameterization", errU, 1e-13)
report("3b uc dU vs finite difference (radian derivatives)", errdU, 1e-8)
report("3c ||V-U||_F / (sqrt6 h) max over samples (must be <=1)", worst1, 1.0)
report("3d ||V-U-dU||_F / (sqrt42/2 h^2) max (must be <=1)", worst2, 1.0)
report("4a makeD: base+W.U vs tex D_span+D_ker", errD, 1e-11)
report("4b makeJ: jac columns vs tex linear map J(s,theta)", errJ, 1e-10)
report("4c pre-equation constant c_X=(1+a+sum b)/i", errPre, 1e-15)

# ---------- 5. buildH::gen (SRC patchx.cpp:4-10, 49-72) vs TEX phase-cover.tex:23-38 ------
q = np.zeros((6, 5))
for j in range(5):
    u = 1./math.sqrt((j+1.)*(j+2.)); q[:j+1, j] = u; q[j+1, j] = -(j+1.)*u
Qb = np.einsum('ia,jb->ijab', q, q).reshape(36, 25)
def Ftex(P):
    out = []
    for i in range(6):
        for j in range(i):
            z = e(P[i]-P[j]).sum()/(2j*np.pi); out += [z.real, z.imag]
    return np.array(out)
P = rng.uniform(-1, 1, (6, 6)); unit = e(P)
f = []; J = np.zeros((30, 25)); fullL = np.zeros((30, 36)); at = 0
for i in range(6):
    for j in range(i):
        qq = unit[i]*np.conj(unit[j]); sm = qq.sum()
        f += [sm.imag/(2*PI), -sm.real/(2*PI)]
        cols = qq @ q
        for a_ in range(5):
            u = q[i, a_]-q[j, a_]
            J[at, 5*a_:5*a_+5] = u*cols.real; J[at+1, 5*a_:5*a_+5] = u*cols.imag
        for hh in range(36):
            row, g = divmod(hh, 6)
            fullL[at, hh] = qq[g].real if row == i else (-qq[g].real if row == j else 0)
            fullL[at+1, hh] = qq[g].imag if row == i else (-qq[g].imag if row == j else 0)
        at += 2
report("5a gen residual f vs tex F(P)", abs(np.array(f)-Ftex(P)).max(), 1e-14)
hfd = 1e-6; fdL = np.zeros((30, 36))
for hh in range(36):
    d = np.zeros(36); d[hh] = hfd
    fdL[:, hh] = (Ftex(P+d.reshape(6, 6))-Ftex(P-d.reshape(6, 6)))/(2*hfd)
report("5b fullL vs finite-difference derivative of F", abs(fdL-fullL).max(), 1e-8)
report("5c 25-coordinate J == fullL . Q(Helmert tensor basis)", abs(fullL@Qb-J).max(), 1e-13)
Pi = np.kron(np.eye(6)-1/6, np.eye(6)-1/6)
report("5d Q Q^T == double-centering projector Pi (PI in patch, patchx.cpp:169)", abs(Qb@Qb.T-Pi).max(), 1e-14)
z = rng.normal(size=36); r_ = z.reshape(6, 6)
projsrc = (r_ + r_.sum()/36 - (r_.sum(1)[:, None]+r_.sum(0)[None, :])/6).reshape(36)
report("5e proj() (patchx.cpp:12-16) == Pi", abs(projsrc-Pi@z).max(), 1e-14)
# tensor sign factor (patchx.cpp:145-149) vs (e_i-e_j).(e_i'-e_j')
ei = []; ej = []
for i in range(6):
    for j in range(i): ei += [i, i]; ej += [j, j]
bad = 0
for a_ in range(30):
    for b_ in range(30):
        sign = (ei[a_] == ei[b_])+(ej[a_] == ej[b_])-(ei[a_] == ej[b_])-(ej[a_] == ei[b_])
        va = np.zeros(6); va[ei[a_]] = 1; va[ej[a_]] = -1; vb = np.zeros(6); vb[ei[b_]] = 1; vb[ej[b_]] = -1
        bad += int(sign != va@vb)
report("5f tensor sign factor mismatches", bad, 0)
# remainder constant: sum a_ij = 6 s^2, sum a_ij^2 <= 8 s^4 for horizontal x (phase-cover.tex:205-222)
w8 = 0.
for _ in range(3000):
    x = (Pi@rng.normal(size=36)).reshape(6, 6)
    if rng.random() < .3: x = (Pi@np.outer(rng.normal(size=6), rng.normal(size=6)).reshape(36)).reshape(6, 6)
    a2 = [((x[i]-x[j])**2).sum() for i in range(6) for j in range(i)]
    s2 = (x**2).sum()
    assert abs(sum(a2)-6*s2) < 1e-9*s2
    w8 = max(w8, sum(v*v for v in a2)/s2**2)
report("5g max sum a_ij^2 / s^4 over samples (tex: <= 8)", w8, 8.0)

# ---------- 6. C stage: bin geometry, single-vector and edge linearizations -----------------
a_h = 1.0; best = [0, 0, 0]
for signs in itertools.product([-1, 1], repeat=5):
    d = np.array([0]+list(signs), float); x = d-d.sum()/6
    best = [max(best[0], np.linalg.norm(x)), max(best[1], abs(x).sum()), max(best[2], abs(x).max())]
report("6a X_M/a = sqrt(29/6) (vertex max)", abs(best[0]-math.sqrt(29/6)), 1e-14)
report("6b L_M/a = 5", abs(best[1]-5), 1e-14)
report("6c I_M/a = 3/2", abs(best[2]-1.5), 1e-14)
report("6d source sqrt(5-1./6) == sqrt(29/6)", abs(math.sqrt(5-1./6)-math.sqrt(29/6)), 1e-15)
inv6 = 1/math.sqrt(6.)
def single_src(v, Pm):
    """SRC stageC.cpp:239-257 / 56-80: q[j]=e(v_j - p[j][i])/sqrt6, S_i, c_i, J[i][j] = -Im(conj(S_i) q_ji)."""
    qm = e(v[:, None]-Pm)*inv6            # qm[j,i]
    S = qm.sum(0)
    c = (abs(S)**2-1)/(4*PI)
    Jm = np.array([[-(S[i].real*qm[j, i].imag - S[i].imag*qm[j, i].real) for j in range(6)] for i in range(6)])
    return c, Jm, S, qm
errJ1 = errB = 0.; rowsum = 0.
for _ in range(50):
    Pm = rng.uniform(-1, 1, (6, 6)); v = rng.uniform(0, 1, 6); v[0] = 0
    c, Jm, S, qm = single_src(v, Pm)
    rowsum = max(rowsum, abs(Jm.sum(1)).max())
    hfd = 1e-6
    for j in range(6):
        d = np.zeros(6); d[j] = hfd
        fd = (single_src(v+d, Pm)[0]-single_src(v-d, Pm)[0])/(2*hfd)
        errJ1 = max(errJ1, abs(fd-Jm[:, j]).max())                 # dc_i/dx_j = +J_ij
    # dc_i / d delta_{ji} = -J_ij ; B_ik = sum_j J_ij T[(j,i),k]  (stageC.cpp:95, 254; exclusion.tex:164-166)
    T = rng.normal(size=(36,)); TTm = T.reshape(6, 6)             # TT[k][j][i] = T[6j+i]
    Bsrc = np.array([sum(Jm[i, j]*TTm[j, i] for j in range(6)) for i in range(6)])
    fd = (single_src(v, Pm+hfd*TTm)[0]-single_src(v, Pm-hfd*TTm)[0])/(2*hfd)
    errB = max(errB, abs(fd+Bsrc).max())                           # c(P+eta T) ~ c - B eta
report("6e single J == dc/dx (finite difference)", errJ1, 1e-8)
report("6f rows of J sum to zero", rowsum, 1e-14)
report("6g single B: dc/d(eta) == -B (sign/transposition of T index)", errB, 1e-8)
# edge quantities (SRC stageC.cpp:368-399) vs TEX exclusion.tex:393-441
errDO = errDU = 0.; ratioO = ratioU = 0.; ratioDist = 0.
for _ in range(400):
    u = rng.uniform(0, 1, 6); v = rng.uniform(0, 1, 6)
    av, bv = e(u), e(v)
    Q = np.array([complex(av[k].real*bv[k].real+av[k].imag*bv[k].imag, bv[k].real*av[k].imag-bv[k].imag*av[k].real) for k in range(6)])
    assert abs(Q-e(u-v)).max() < 1e-13
    S = Q.sum(); mag = abs(S)
    A0 = -Q.imag; A1 = Q.real; JU = (S.real*A0+S.imag*A1)/6
    qO = np.array([S.real, S.imag])/(2*PI); qU = (abs(S)**2/6-1)/(4*PI)
    def edge(zv):
        Sp = (Q*e(zv)).sum(); return np.array([Sp.real, Sp.imag])/(2*np.pi), (abs(Sp)**2/6-1)/(4*np.pi), Sp
    hfd = 1e-6
    for k in range(6):
        d = np.zeros(6); d[k] = hfd
        fo = (edge(d)[0]-edge(-d)[0])/(2*hfd); fu = (edge(d)[1]-edge(-d)[1])/(2*hfd)
        errDO = max(errDO, abs(fo[0]-A0[k]), abs(fo[1]-A1[k])); errDU = max(errDU, abs(fu-JU[k]))
    re_ = rng.uniform(.005, .2)
    for _ in range(20):
        zv = rng.normal(size=6); zv -= zv.mean(); zv *= re_/np.linalg.norm(zv)
        o, uu, Sp = edge(zv)
        ratioO = max(ratioO, np.linalg.norm(o+0*qO - (qO+np.array([A0@zv, A1@zv])))/(np.pi*re_**2))
        ratioU = max(ratioU, abs(uu-(qU+JU@zv))/(np.pi*re_**2*(1+mag/6)))
        ratioDist = max(ratioDist, abs(Sp-S)/(2*np.pi*math.sqrt(6)*re_))
report("6h edge D_O == d(S'/2pi)/dz", errDO, 1e-8)
report("6i edge D_U == d((|S'|^2/6-1)/4pi)/dz", errDU, 1e-8)
report("6j O Taylor remainder / (pi r^2) max (must be <=1)", ratioO, 1.0)
report("6k U Taylor remainder / (pi r^2 (1+|S|/6)) max (must be <=1)", ratioU, 1.0)
report("6l |S'-S| / (2 pi sqrt6 r_e) max (must be <=1)", ratioDist, 1.0)
lhs = 6-2*math.pi*math.sqrt(6)*.2
report("6m Lem 6.5: 6-2pi sqrt6 *.2 > sqrt6 (print margin as negative slack)", math.sqrt(6)-lhs, 0.0, f"lhs={lhs:.4f}")

# ---------- 7. covout scalars (SRC covout.cpp:100-103) and seed constants -------------------
R = .032; n = 112
Rn = R*n; RN2 = Rn*Rn - .0001; T = (math.sqrt(RN2)+1.5)**2
exact = (Fraction(32, 1000)*112)**2
print(f"     RN2={RN2!r}  exact (R n)^2={float(exact)!r}  exact-RN2={float(exact-Fraction(RN2)):.3e}  T={T:.6f} int(T)={int(T)}  K={max(1,int(R*n+.51))}")
report("7a RN2 < exact (.032*112)^2 (print exact-RN2 must be >0)", float(Fraction(RN2)-exact), 0.0)
print(f"     math.pow(R*n,2) - (R*n)*(R*n) = {math.pow(Rn,2)-Rn*Rn:.3e} (info only; this platform's libm)")
report("7b 3/(2*112) < .032", 3/224-.032, 0.0)
report("7c 1.251 pi/N >= 5 pi/(4N)", 1.25-1.251, 0.0)
# J_* layout of supdual rows (covout.cpp:43-46) vs TEX seed.tex:372-378: identity 2k.h = 2(k-J*^T v).h + 2 v.(J* h)
errS = 0.
for _ in range(200):
    s = int(rng.integers(5)); JX = e(rng.uniform(0, 1, 5)); JY = e(rng.uniform(0, 1, 5))
    rows = np.zeros((9, 4))
    for j in range(5): rows[j, 0] = JX[j].real; rows[j, 1] = JX[j].imag
    for j in range(5):
        hh = 4+j if j else s; rows[hh, 2] = JY[j].real; rows[hh, 3] = JY[j].imag
    h = rng.uniform(-.5, .5, 9)
    zX = (JX*h[:5]).sum(); zY = (JY*np.concatenate([[h[s]], h[5:]])).sum()
    Jh_tex = np.array([zX.real, zX.imag, zY.real, zY.imag])
    errS = max(errS, abs(rows.T@h-Jh_tex).max())
report("7d supdual rows[] == J_* of tex (J_* h)", errS, 1e-14)
# seed residual c (seed.cpp:36-39) vs TEX (1+sum e(u))/(2 pi i)
uu = rng.uniform(0, 1, 5); z = (1+e(uu).sum())/(2j*np.pi)
c0 = sum(math.sin(2*np.pi*t) for t in uu)/(2*PI); c1 = -1/(2*PI)-sum(math.cos(2*np.pi*t) for t in uu)/(2*PI)
report("7e sumok c vs tex c", max(abs(c0-z.real), abs(c1-z.imag)), 1e-14)

# ---------- 8. scale claims used by tex margins ----------------------------------------------
Rs = 2*PI*.0320001+1e-9; mlow = 1.25-1e-8
qmax = Rs*Rs*(1+(4+6)/mlow**2)+1e-7
print(f"     Rs={Rs:.6f}  q_max(m=1.25,Amax=2)={qmax:.5f} (tex <.32)  pol={qmax/(math.sqrt(6)+math.sqrt(6-qmax)):.5f} (tex <.068)")
report("8a R_s < .203", Rs-.203, 0.0)
report("8b q < .32 at worst accepted frame", qmax-.32, 0.0)
arc = math.sqrt(3.8)*(1/math.sqrt(1-3.8/4)+1e-7)+1e-8
report("8c hor*sqrt(Rs^2+arc^2) < 8.81e-5 with hor<1e-5, q<3.8", 1e-5*math.hypot(Rs, arc)-8.81e-5, 0.0, f"arc<{arc:.4f}")
report("8d K=1.00001 > sqrt(1+2pi sqrt2/6 *1e-6)", math.sqrt(1+2*math.pi*math.sqrt(2)/6*1e-6)-1.00001, 0.0)
# ellipsoid lemma constants: (1+b) and (1+1/b) with b=.45, .8 (patchx.cpp:336,343)
report("8e 1.45 == 1+.45 and 1.8 == 1+.8 (as doubles, within 1ulp)", max(abs(1.45-(1+.45)), abs(1.8-(1+.8))), 3e-16)
# B grid: L=ceil(E/w+.500001) => (L+1/2) w >= E
wbad = 0
for _ in range(20000):
    E = rng.uniform(0, .127); w = rng.choice([.013, .005, .0038])
    L = int(math.ceil(E/w+.500001)); wbad += int((L+.5)*w < E)
report("8f B grid covers [-E,E] (violations)", wbad, 0)
print("\nALL OK" if ok_all else "\nSOME CHECK FAILED")
