"""Audit 268, companion paper "A boundary-field gap for the spin-one Heisenberg chain".
(a) ED of the open chains H_{n,g} (n = 2..6, g = 1,-1,0) at b_0 = 49/4: compare the
    15 thermal entries k_{n,g} (thermal.tex, eq. appth:integer-table), the I(n) centers,
    Y_6 < .027 and the polynomial bound c_0 P_*(lambda_0)^2 < .003745.
(b) Exact-integer open-chain trial contraction (trials.tex, Lemma apptr:contractions)
    for m = 34, 78, 142, 240 with the 26x26 periodic matrices and endpoint vectors.
(c) Exact rational check of the scalar chain in periodic.tex / initialization.tex /
    conclusion.tex.
"""
import numpy as np, functools, itertools
from fractions import Fraction as Fr
from math import isqrt, factorial, exp, log

fails = []
def check(label, ok, extra=""):
    print(f"[{'OK ' if ok else 'FAIL'}] {label} {extra}")
    if not ok: fails.append(label)
def fl(x): return f"{float(x):.10g}"

a = Fr(700741, 500000); b0 = Fr(49, 4); h = Fr(3, 5); af = float(a); b0f = float(b0)

# ---------- (a) open-chain thermal traces ----------
print("== (a) companion thermal appendix: open chains n=2..6 ==")
Sz = np.diag([1.0, 0.0, -1.0]); Sp = np.zeros((3,3)); Sp[0,1] = Sp[1,2] = np.sqrt(2); Sm = Sp.T; Sx = (Sp+Sm)/2
def op(o, j, n):
    mats = [np.eye(3)]*n; mats[j] = o
    return functools.reduce(np.kron, mats)
def Hopen(n, g):
    H = np.zeros((3**n, 3**n))
    for j in range(n-1):
        H += op(Sz,j,n)@op(Sz,j+1,n) + 0.5*(op(Sp,j,n)@op(Sm,j+1,n) + op(Sm,j,n)@op(Sp,j+1,n))
    H += -0.6*op(Sz,0,n)
    R = {1: op(Sz,n-1,n), -1: -op(Sz,n-1,n), 0: op(Sx,n-1,n)}[g]
    H += 0.6*R
    return H
ktab = {2: (10185296075298, 4160017987, 583578389575), 3: (259254948723, 571111660440989, 47796720049308),
        4: (318278534893, 661347985, 43971049724), 5: (260908126900, 29917222640283, 5659378873682),
        6: (16467340956, 328219854, 4086161507)}
Zc = {}
Icenters = {}
for n in range(2, 7):
    Cn = 3*((n-1)//2) + 2*((n-1) % 2) + Fr(6,5)
    for gi, g in enumerate((1, -1, 0)):
        ev = np.linalg.eigvalsh(Hopen(n, g))
        assert ev[0] + float(Cn) >= -1e-9 and ev[-1] + float(Cn) <= 16 + 1e-9, "contraction bound A"
        Z = float(np.sum(np.exp(-b0f*(ev + af*n))))
        Zc[(n,g)] = Z
        # paper: Z = e^{b0 (C_n - a n)} W,  W in [k/1e16, (k+1)/1e16)
        W = Z*exp(-b0f*(float(Cn) - af*n))
        k = ktab[n][gi]
        ok = k/1e16 - 1e-9 <= W <= (k+1)/1e16 + 1e-9   # allow float tolerance ~1e-9 on W ~ 1e-3..1e-1
        check(f"n={n} g={g:2d}: Z={Z:.10f}  W=Z e^-b0(Cn-an)={W:.14f}  table k/1e16={k/1e16:.14f}", ok, f"diff {W-k/1e16:+.2e}")
    I_n = (Zc[(n,1)] + Zc[(n,-1)] + 4*Zc[(n,0)])/6
    Icenters[n] = I_n
cen = {2: .02703074, 3: .01203562, 4: .01193511, 5: .00931560, 6: .00893272}
for n in range(2,7):
    check(f"I({n}) = {Icenters[n]:.8f} vs center {cen[n]} (|diff|<1e-7)", abs(Icenters[n]-cen[n]) < 1e-7, f"diff {Icenters[n]-cen[n]:+.2e}")
check("Y_6(b0) = Z_{6,1} < .027", Zc[(6,1)] < .027, f"{Zc[(6,1)]:.10f} (paper interval .0266319440..0266319607)")
poly = Icenters[6] - .28*Icenters[5] - .3404*Icenters[4] + .0504*Icenters[3] + .0324*Icenters[2]
check("I(6)-.28I(5)-.3404I(4)+.0504I(3)+.0324I(2) < .003745", poly < .003745, f"{poly:.12f} (paper .003744031780 from centers)")
# P_*(x)^2 expansion check: P_* = x^3 - .14x^2 - .18x -> P^2 = x^6 - .28 x^5 + (.0196-.36) x^4 + 2(.14)(.18) x^3 + .0324 x^2
P2 = np.polymul([1,-.14,-.18,0],[1,-.14,-.18,0])
check("P_*(x)^2 coefficients = (1, -.28, -.3404, .0504, .0324, 0, 0)", np.allclose(P2, [1,-.28,-.3404,.0504,.0324,0,0]), str(P2))

# ---------- (b) open-chain trial contraction ----------
print("\n== (b) companion trial appendix: exact integer contraction ==")
ell = [1,1,1,1,3,3,3,5]
Dpar = {1:-289,2:-84,3:-9,4:33,5:-22,6:27,7:71,8:-11}
Ppar = {(1,5):-11,(1,6):-35,(1,7):-28,(2,5):-23,(2,6):16,(2,7):-49,(3,5):-6,(3,6):-25,(3,7):-20,(4,5):-34,(4,6):9,(4,7):-32,(5,8):1,(6,8):2,(7,8):-16}
beta = {1:10000,2:5620,3:2456,4:4868,5:82,6:-809,7:12475,8:-9645}
kappa = {1:12,2:12,3:12,4:12,5:6,6:6,7:6,8:1}
labels = [(i,d) for i in range(1,9) for d in range(-ell[i-1], ell[i-1]+1, 2)]
idx = {lab:k for k,lab in enumerate(labels)}
def entry(al,d,be,e,s):
    l,k = ell[al-1], ell[be-1]
    if e != d+2*s or abs(k-l) > 2 or (k == l and al != be): return 0
    r = (d+l)//2
    if k == l+2: p = Ppar[(al,be)]; fac = {-1:2*(l+1-r)*(l+2-r), 0:2*(r+1)*(l-r+1), 1:(r+1)*(r+2)}[s]
    elif k == l: p = Dpar[al]; fac = {-1:2*(l-r+1), 0:d, 1:-(r+1)}[s]
    else: p = -Ppar[(be,al)]; fac = {-1:2, 0:-2, 1:1}[s]
    return p*fac
A = {s: np.array([[entry(al,d,be,e,s) for (be,e) in labels] for (al,d) in labels], dtype=object) for s in (-1,0,1)}
bvec = np.array([beta[i] if d == -1 else 0 for (i,d) in labels], dtype=object)
cvec = np.array([kappa[i]*beta[i] if d == -1 else 0 for (i,d) in labels], dtype=object)
T = sum(np.kron(A[s],A[s]) for s in (-1,0,1)); G = sum(s*np.kron(A[s],A[s]) for s in (-1,0,1))
B = {(s,t): A[s].dot(A[t]) for s in (-1,0,1) for t in (-1,0,1)}
O = sum(s*t*np.kron(B[(s,t)],B[(s,t)]) for s in (-1,0,1) for t in (-1,0,1))
for s in (-1,0):
    for t in (0,1):
        O = O + np.kron(B[(s,t)],B[(s+1,t-1)]) + np.kron(B[(s+1,t-1)],B[(s,t)])
pairs = [(i,j) for i in range(26) for j in range(26)]
K0 = [k for k,(i,j) in enumerate(pairs) if labels[i][1] == labels[j][1]]
check("dim K_0 = 162", len(K0) == 162)
Tr, Gr, Or = T[np.ix_(K0,K0)], G[np.ix_(K0,K0)], O[np.ix_(K0,K0)]
L0 = np.kron(bvec,bvec)[K0]; u = np.kron(cvec,cvec)[K0]
x = L0.copy(); y = L0.dot(Tr); w = -300000*L0.dot(Gr)
tab = {34:(202,715867,4031974), 78:(440,747265,4041175), 142:(787,121288,4041110), 240:(1317,165352,4041006)}
# brute-force check of the contraction formula at small m against explicit psi and H^opp
Af = {s: A[s].astype(float) for s in A}
for m in (3, 4):
    psi = np.zeros(3**m)
    for conf in itertools.product((1,0,-1), repeat=m):
        M = bvec.astype(float)
        for s in conf: M = M @ Af[s]
        amp = M @ cvec.astype(float)
        ind = 0
        for s in conf: ind = 3*ind + (1-s)
        psi[ind] = amp
    H = Hopen(m, 1)
    # recurrence to m
    xx = L0.copy(); yy = L0.dot(Tr); ww = -300000*L0.dot(Gr)
    for mm in range(2, m+1):
        ww, xx, yy = ww.dot(Tr) + 500000*xx.dot(Or) + 700741*xx.dot(Tr).dot(Tr), yy, yy.dot(Tr)
    V = int(yy.dot(u)); Hs = int(ww.dot(u) + 300000*xx.dot(Gr.dot(u)) + 700741*V)
    lhs = 500000*(psi @ H @ psi + af*m*(psi@psi))
    check(f"m={m}: brute force <psi,(H+am)psi>*5e5 = {lhs:.6e} vs H_m^* = {Hs} ; |psi|^2 = {psi@psi:.6e} vs V_m = {V}", abs(lhs-Hs) < 1e-6*abs(Hs) and abs(psi@psi - V) < 1e-9*V)
for m in range(2, 241):
    w, x, y = w.dot(Tr) + 500000*x.dot(Or) + 700741*x.dot(Tr).dot(Tr), y, y.dot(Tr)
    if m in tab:
        V = int(y.dot(u)); Hs = int(w.dot(u) + 300000*x.dot(Gr.dot(u)) + 700741*V)
        p, q, k = tab[m]
        ok = (10**(p-1) <= V < 10**p) and (q*10**(p-6) <= V < (q+1)*10**(p-6)) and (k*V <= 20*Hs < (k+1)*V) and (202100*V - Hs > 0)
        check(f"m={m}: p,q,k certificate and 202100V-H*>0", ok, f"H*/(5e5 V) = {Hs/(500000*V):.9f} < alpha=.4042")

# ---------- (c) scalar chain ----------
print("\n== (c) scalar chain (periodic.tex, initialization.tex, conclusion.tex) ==")
l, U, V, mJ, mN = Fr('.999'), Fr('1.00139'), Fr('.872'), Fr('1.0657205'), Fr('.2783155')
def f(u): return u*u/(2*(1-u)**2)
def d(u): return u/(2-3*u)
check("m_J - l^12 < (.874 l)^12 and V < .874 l", mJ-l**12 < (Fr('.874')*l)**12 and V < Fr('.874')*l)
def D(k):  # k multiple of 12 only for exact rationals
    return ((mJ-l**12)**Fr(k,12) if k % 12 == 0 else None)
# D(32) via cube-root bounds as in the paper
c1, c2 = Fr('.426635'), Fr('.439736')
check("(m_J-l^12)^(1/3) < .426635 (cube)", mJ-l**12 < c1**3)
check("(m_N-V^12)^(1/3) < .439736 (cube)", mN-V**12 < c2**3)
D32 = (c1**8 + 3*(V**32 + c2**8))/l**32
check("D(32) < .0442", D32 < Fr('.0442'), fl(D32))
D120 = ((mJ-l**12)**10 + 3*(V**120 + (mN-V**12)**10))/l**120
check("D(120) < .00000025", D120 < Fr('.00000025'), fl(D120))
check("(.9856)^9 < (.968)^4  [i.e. .968^(-32/72) < 1/.9856]", Fr('.9856')**9 < Fr('.968')**4)
t32 = U**64*(1+Fr('.0442'))**2/Fr('.9856') - 1
check("U^64(1.0442)^2/.9856 - 1 < .2092 < .9525^32", t32 < Fr('.2092') < Fr('.9525')**32, f"{fl(t32)}, .9525^32={fl(Fr('.9525')**32)}")
check(".0604996 < 2(.032)(.968)", Fr('.0604996') < 2*Fr('.032')*Fr('.968'))
s_ = (1+Fr('.9525')**60)**-2; R_ = U**120*(1+Fr('.00000025'))-1; t_ = (1+R_**2)**-2
check("s > .9002", s_ > Fr('.9002'), fl(s_)); check("t > .9372", t_ > Fr('.9372'), fl(t_))
check("f(1-.9002^2 .9372) < .0502", f(1-Fr('.9002')**2*Fr('.9372')) < Fr('.0502'), fl(f(1-Fr('.9002')**2*Fr('.9372'))))
check("f(1-.9372^2(1-.0502)) < .0198", f(1-Fr('.9372')**2*(1-Fr('.0502'))) < Fr('.0198'))
def Gf(u): return Fr(1,2)*((1-u)**-1+(1-u)**-2+(1-u)**-3)**2
check("G(.052) < 5.583 < 6 and 6u^2 <= u on [0,.052]", Gf(Fr('.052')) < Fr('5.583') and 6*Fr('.052')**2 <= Fr('.052'), fl(Gf(Fr('.052'))))
check("d(.052) < .971^120", d(Fr('.052')) < Fr('.971')**120, f"{fl(d(Fr('.052')))} vs {fl(Fr('.971')**120)}")
# initialization
alpha = Fr('.4042'); eb0a = exp(float(b0*alpha))
def eps(m): return Fr('.027')*eb0a*float(U**m)*float(Fr('.874')**(m-6))/float(l**6)
def Rm(m): return exp(.132)*float(U**m) + eps(m) - 1
Pst = lambda x: x**3 - Fr('.14')*x*x - Fr('.18')*x
check("P_*(x) = x(x-1/2)(x+9/25)", all(Pst(Fr(k,7)) == Fr(k,7)*(Fr(k,7)-Fr(1,2))*(Fr(k,7)+Fr(9,25)) for k in range(1,5)))
lc = float(Fr('.003745')*Fr(1))*eb0a/float(Pst(U)**2)
check(".003745 e^{b0 alpha}/P_*(U)^2 < 1.133226 and e^.132 > 1.141108", lc < 1.133226 and exp(.132) > 1.141108, f"{lc:.6f}")
check("eps(34) < .0928", eps(34) < .0928, f"{eps(34):.5f}")
def F(q, e, t): return q/(1-q)*(1/((1-e)**2*t) - 1)
e78 = F(.953**44, .093, (1+Rm(78))**-2)
e142 = F(.971**64, .133, (1+Rm(142)**2)**-2)
A240 = (1+Rm(240)**4)**-2
check("(1+R(78))^-2 > .6180", (1+Rm(78))**-2 > .6180, f"{(1+Rm(78))**-2:.5f}")
check("(1+R(142)^2)^-2 > .7533", (1+Rm(142)**2)**-2 > .7533, f"{(1+Rm(142)**2)**-2:.5f}")
check("(1+R(240)^4)^-2 > .7923", A240 > .7923, f"{A240:.5f}")
check("e_78(1) bound < .132129 < .133", e78 < .132129, f"{e78:.6f}")
check("e_142(2) bound < .137301 < .14", e142 < .137301, f"{e142:.6f}")
check("1-A(240,2) < .207518 < .21", 1-A240 < .207518, f"{1-A240:.6f}")
# coupled boundary updates table (in:table): u_2=.052, u_{j+1}=6u_j^2 ; M_j = 240*2^(j-2), k_2=142, k_{j+1}=M_j
u = {2: Fr('.052')}
for j in range(2, 6): u[j+1] = 6*u[j]**2
E = {2: Fr('.14')}; eta = {2: Fr('.21')}
Mj = lambda j: 240*2**(j-2)
kj = {2: 142}
for j in range(2, 5):
    kj[j+1] = Mj(j)
    expo = Fr(Mj(j)-kj[j], Mj(j))
    q = float(d(u[j+1]))**float(expo)
    qbound = {2: .142, 3: .0282, 4: .0028}[j]
    Enew = F(q, float(E[j]), float(1-eta[j]))
    Etab = {2: Fr('.12'), 3: Fr('.016'), 4: Fr('.001')}[j]
    etanew = float(f(1-(1-eta[j])*(1-u[j])*(1-u[j+1])*(1-Etab)))
    etatab = {2: Fr('.155'), 3: Fr('.026'), 4: Fr('.001')}[j]
    check(f"step {j}->{j+1}: q={q:.5f} <= {qbound}; E_out={Enew:.6f} <= {float(Etab)}; eta_out={etanew:.6f} <= {float(etatab)}; exponent {expo}",
          q <= qbound and Enew <= float(Etab) and etanew <= float(etatab))
    E[j+1] = Etab; eta[j+1] = etatab
check("u_5 < .000015", u[5] < Fr('.000015'), fl(u[5]))
check("d(u_3)^49 < .142^120 (first q bound, exact)", d(u[3])**49 < Fr('.142')**120)
check("42430000/6792093 < 10 and 3645000/1495729 < 10", Fr(42430000,6792093) < 10 and Fr(3645000,1495729) < 10)
xx = Fr('.01')
check("F(2x,x,1-x)/x^2 endpoint value at x=.01 equals 42430000/6792093",
      Fr(2,(1-2*xx))*((1-xx)**-1+(1-xx)**-2+(1-xx)**-3) == Fr(42430000,6792093) and F(2*xx,xx,1-xx)/xx**2 == Fr(2,(1-2*xx))*((1-xx)**-1+(1-xx)**-2+(1-xx)**-3))
check("f(2x+16x^2)/x^2 at x=.01 equals 3645000/1495729", f(2*xx+16*xx*xx)/xx**2 == Fr(3645000,1495729))
check("sqrt(6x^2/(2-18x^2)) <= 2x at x=.01", Fr(6)*xx*xx/(2-18*xx*xx) <= 4*xx*xx)
check("b_5 = 392, M_5 = 1920, 4/(10(1-.04)) < 1", b0*32 == 392 and Mj(5) == 1920 and Fr(4,10*(1-Fr('.04'))) < 1)
print(f"   gap bound log(10)/392 = {log(10)/392:.6e}")
print("\nSUMMARY:", "ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURES: {fails}")
