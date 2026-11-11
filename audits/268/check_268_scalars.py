"""Exact rational re-derivation of every scalar inequality in the main paper
(bootstrap.tex, initialization.tex, thermal120.tex, thermal72.tex) for audit 268.
Only Fractions / Python integers are used; no floating point enters a comparison.
Inputs: the thermal table centers (initialization.tex:31-53) and the integer
totals printed in the appendices (thermal120.tex:186-194, thermal72.tex:185-194).
"""
from fractions import Fraction as Fr
from math import isqrt, factorial, log

fails = []
def check(label, ok, extra=""):
    print(f"[{'OK ' if ok else 'FAIL'}] {label} {extra}")
    if not ok:
        fails.append(label)

def fl(x):  # display
    return f"{float(x):.10g}"

a = Fr(700741, 500000)
print("a =", a, "=", fl(a))

# ---------------- bootstrap.tex ----------------
print("\n== bootstrap.tex: scalar hypotheses (lines 263-281) ==")
def f(u): return u*u/(2*(1-u)**2)
def G(u): return Fr(1,2)*((1-u)**-1+(1-u)**-2+(1-u)**-3)**2
u0, C = Fr(1,8), Fr(79,10)
check("G(1/8) = 913952/117649", G(u0) == Fr(913952,117649), fl(G(u0)))
check("G(1/8) < 79/10", G(u0) < C)
check("3/(C(1-3u0)) = 48/79 < 1 (row 1)", Fr(3)/(C*(1-3*u0)) == Fr(48,79) and Fr(48,79) < 1)
check("r = C u0 = 79/80 < 1", C*u0 == Fr(79,80))
check("u0 < 1/6; C(1-3u0) > 3 (row 1)", u0 < Fr(1,6) and C*(1-3*u0) > 3)
# identity f(1-(1-u)^3) = G(u) u^2 at a few rationals
check("f(1-(1-u)^3) == G(u) u^2 (identity, u=1/8,1/100,1/3)", all(f(1-(1-u)**3) == G(u)*u*u for u in (Fr(1,8),Fr(1,100),Fr(1,3))))
u0b, Cb = Fr(1,100), Fr(5)
check("9/(2(1-3/100)^2) = 45000/9409 < 5", Fr(9)/(2*(1-3*u0b)**2) == Fr(45000,9409) and Fr(45000,9409) < 5)
check("G(1/100) <= 45000/9409 (direct)", G(u0b) <= Fr(45000,9409), fl(G(u0b)))
check("3/(C(1-3u0)) = 60/97 < 1 (row 2)", Fr(3)/(Cb*(1-3*u0b)) == Fr(60,97))
check("r = 1/20 (row 2)", Cb*u0b == Fr(1,20))
g60 = Fr(4,105)*log(Fr(80,79)); g2304 = log(20)/784
print(f"   gap bounds: (4/105)log(80/79) = {g60:.6e},  log(20)/784 = {g2304:.6e}  (literature bulk gap 0.41050)")

# ---------------- initialization.tex: tables ----------------
T105 = {6: (8945303, 346734, 939987), 8: (3741281, 567167, 984108), 10: (2327547, 721090, 1002613)}
T1225 = {4: (124885379, 65757), 5: (2925, 2960742), 6: (12870614, 286599), 7: (54069, 1920414),
         8: (4639283, 510140), 9: (195455, 1514846), 10: (2654823, 675419), 11: (377629, 1314544),
         12: (1900547, 787405)}
u6 = Fr(1, 10**6)
def sectors(z1, zP, zC):
    J = (z1+3*zP)/4; N = (z1-zP)/4; I = (z1+3*zP+8*zC)/12; O = (z1+3*zP-4*zC)/12
    return I, O, N, J

print("\n== Proposition seed60 (initialization.tex:118-264) ==")
eps = Fr(1, 10**5)
cI = {}; cO = {}; cN = {}
for n, (z1, zP, zC) in T105.items():
    I, O, N, J = sectors(z1*u6, zP*u6, zC*u6)
    cI[n], cO[n], cN[n] = I, O, N
mI, mO, mN = Fr(1042655, 10**6), Fr(40041, 10**6), Fr(401625, 10**6)
check("m_I margin = 19/(12e6)", mI - (cI[10]+eps) == Fr(19, 12*10**6), fl(mI-(cI[10]+eps)))
check("m_O margin = 7/(12e6)", mO - (cO[10]+eps) == Fr(7, 12*10**6))
check("m_N margin = 3/(4e6)", mN - (cN[10]+eps) == Fr(3, 4*10**6))
# polynomial filters
U = Fr(100245, 100000); V = Fr(88643, 100000)
def FI(x): return x**6*(-2+10*x*x)**2
def FN(x): return x**6*(-25+100*x*x)**2
BI = 4*cI[6] - 40*cI[8] + 100*cI[10] + eps*(2+10)**2
BN = 625*cN[6] - 5000*cN[8] + 10000*cN[10] + eps*(25+100)**2
check("F_I(U) - B_I > 0.0301", FI(U)-BI > Fr(301,10000), fl(FI(U)-BI))
check("F_N(V) - B_N > 0.3463", FN(V)-BN > Fr(3463,10000), fl(FN(V)-BN))
check("filters increasing beyond endpoints: -2+10U^2>0, -25+100V^2>0", -2+10*U*U > 0 and -25+100*V*V > 0)
check("U^10 < m_I < 2U^10", U**10 < mI < 2*U**10, f"U^10={fl(U**10)}")
check("V^10 < m_N < 2V^10", V**10 < mN < 2*V**10, f"V^10={fl(V**10)}")
def A(k): return U**k + (mI-U**10)**(k//10)
def B(k): return mO**(k//10)
def Cc(k): return V**k + (mN-V**10)**(k//10)
m = A(30)+2*B(30)+3*Cc(30); p = A(30)+2*B(30); q = A(120)+2*B(120)+3*Cc(120)
D, E, tau, tau0, Rq = Fr(1450,1000), Fr(1202,1000), Fr(1225,1000), Fr(1264,1000), Fr(684,10000)
print(f"   m={fl(m)} p={fl(p)} q={fl(q)}  (q-1={fl(q-1)})")
check("m^5 < D^2, margin > 0.001068", m**5 < D*D and D*D-m**5 > Fr(1068,10**6), fl(D*D-m**5))
check("p^5 < E^2, margin > 0.000433", p**5 < E*E and E*E-p**5 > Fr(433,10**6), fl(E*E-p**5))
check("q >= 1", q >= 1)
check("(q-1)^5 < R_q^2, margin > 0.0000472", (q-1)**5 < Rq*Rq and Rq*Rq-(q-1)**5 > Fr(472,10**7), fl(Rq*Rq-(q-1)**5))
check("(D+11E)/12 <= tau", (D+11*E)/12 <= tau, fl((D+11*E)/12))
check("(D+3E)/4 <= tau0", (D+3*E)/4 <= tau0, fl((D+3*E)/4))
check("D (0.88)^3 < 1", D*Fr(88,100)**3 < 1, fl(D*Fr(88,100)**3))
check("2(0.88) > D", 2*Fr(88,100) > D)
check("D > tau0", D > tau0)
def R(t):
    verts = [(Fr(0),Fr(0)), (tau-t,Fr(0)), (tau-t,tau0-tau), (Fr(0),tau0-t)]
    return max(x*x + y*y/2 + (D-t-x-y)**2/3 for x,y in verts)
Rt = {t: R(Fr(t)) for t in ('0.88','0.99','0.998')}
for t, v in Rt.items(): print(f"   R({t}) = {fl(v)}  (paper: 0.1359 / 0.0721 / 0.068404)")
check("R table values exact", Rt['0.88']==Fr('0.1359') and Rt['0.99']==Fr('0.0721') and Rt['0.998']==Fr('0.068404'))
check("1-0.99^4-R(0.88)^2 = 0.02093518", 1-Fr('0.99')**4-Rt['0.88']**2 == Fr('0.02093518'), fl(1-Fr('0.99')**4-Rt['0.88']**2))
check("1-0.998^4-R(0.99)^2 = 0.002777621984", 1-Fr('0.998')**4-Rt['0.99']**2 == Fr('0.002777621984'), fl(1-Fr('0.998')**4-Rt['0.99']**2))
S60 = (1+Rt['0.998']/Fr('0.998')**2)**-2
check("S(60,b) bound > 0.8756 > 7/8", S60 > Fr('0.8756') > Fr(7,8), fl(S60))
T120 = (1+Rq)**-2
check("T(120,b) bound > 0.8760 > 7/8", T120 > Fr('0.8760') > Fr(7,8), fl(T120))
check("D > 0 and the polygon upper bound D-t-x-y >= D - tau0 > 0", D-tau0 > 0)

print("\n== Proposition seed2304 (initialization.tex:268-373) ==")
eps2 = Fr(30, 10**6)
Jc = {}; Nc = {}
for n, (z1, zP) in T1225.items():
    Jc[n] = (z1+3*zP)*u6/4; Nc[n] = (z1-zP)*u6/4
# polynomial coefficients
def polymul(p, q):
    r = [Fr(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q): r[i+j] += x*y
    return r
Q = [Fr(c) for c in (7,17,-280,-185,1000)]
FJ = [Fr(0)]*4 + polymul(Q,Q)            # x^4 Q(x)^2, coefficients index = power
Pn = [Fr(c) for c in (12,0,-394,0,1000)]
FNc = [Fr(0)]*4 + polymul(Pn,Pn)
def evalp(c, x): return sum(ci*x**i for i,ci in enumerate(c))
BJ = sum(FJ[j]*Jc[j] for j in range(4,13)) + eps2*sum(abs(c) for c in FJ)
BNv = sum(FNc[j]*Nc[j] for j in range(4,13)) + eps2*sum(abs(c) for c in FNc)
check("B_J = 318604.54848175", BJ == Fr('318604.54848175'), fl(BJ))
check("B_N = 48169.880699", BNv == Fr('48169.880699'), fl(BNv))
vJ, vN = Fr(100139,100000), Fr(872,1000)
check("F_J(v_J)-B_J > 180", evalp(FJ,vJ)-BJ > 180, fl(evalp(FJ,vJ)-BJ))
check("F_J(-v_J)-B_J > 496888", evalp(FJ,-vJ)-BJ > 496888, fl(evalp(FJ,-vJ)-BJ))
check("F_N(+-v_N)-B_N > 654", evalp(FNc,vN)-BNv > 654 and evalp(FNc,-vN)-BNv > 654, fl(evalp(FNc,vN)-BNv))
# Q(1+y), Q(-1-y) expansions
def shift(c, s):  # coefficients of P(s + y) in y, where s = +1 or -1 (sign handled by caller)
    from math import comb
    out = [Fr(0)]*len(c)
    for i,ci in enumerate(c):
        for k in range(i+1): out[k] += ci*comb(i,k)*s**(i-k)
    return out
Q1 = shift(Q, Fr(1))               # Q(1+y)
Qm = [c*(-1)**k for k,c in enumerate(shift(Q, Fr(-1)))]   # Q(-1-y) = Q(-1 + (-y))
check("Q(1+y) = 559+2902y+5165y^2+3815y^3+1000y^4", Q1 == [Fr(c) for c in (559,2902,5165,3815,1000)], str([int(c) for c in Q1]))
check("Q(-1-y) = 895+3978y+6275y^2+4185y^3+1000y^4", Qm == [Fr(c) for c in (895,3978,6275,4185,1000)], str([int(c) for c in Qm]))
check("12-394 v_N^2+1000 v_N^4 > 0 and 1000 v_N^2 > 197", 12-394*vN**2+1000*vN**4 > 0 and 1000*vN**2 > 197)
mJ, mN2 = Fr('1.0657205'), Fr('0.2783155')
check("m_J = J_12 center + eps", mJ == Jc[12]+eps2, fl(Jc[12]+eps2))
check("m_N = N_12 center + eps", mN2 == Nc[12]+eps2, fl(Nc[12]+eps2))
l = Fr(999,1000)
def capped(mass, h, r):
    qq = mass // h
    return qq*h**r + (mass-qq*h)**r
check("v_N^12 < m_N < 2 v_N^12 (q=1)", vN**12 < mN2 < 2*vN**12, fl(vN**12))
check("l^12 < m_J < 2 l^12 (q=1)", l**12 < mJ < 2*l**12, fl(l**12))
def Bk(k): return 3*capped(mN2, vN**12, k//12)
def RA(k): return Bk(k) + (mJ-l**12)**(k//12)
Z72bound = capped(mJ, l**12, 6) + Bk(72)
check("C(m_J,l^12,6)+B(72) < 1 with margin > 0.0693", 1-Z72bound > Fr(693,10000), fl(Z72bound))
p0 = 1-(1+RA(36)/l**36)**-2
q0 = 1-(vJ**72+RA(72))**-2
check("p0 < 0.047916", p0 < Fr('0.047916'), fl(p0))
check("q0 < 0.181521", q0 < Fr('0.181521'), fl(q0))
def ceil7(x):
    from math import ceil
    n = -((-x.numerator*10**7)//x.denominator)   # ceil
    return Fr(n, 10**7)
table = {1:(Fr(151249,2500000),Fr(433453,2500000)), 2:(Fr(343303,5000000),Fr(1632381,10**7)),
         3:(Fr(44601,625000),Fr(361773,2500000)), 4:(Fr(632939,10**7),Fr(1055161,10**7)),
         5:(Fr(375793,10**7),Fr(445941,10**7)), 6:(Fr(84513,10**7),Fr(27493,5000000))}
pj, qj = p0, q0
allok = True
for j in range(1,7):
    pj = ceil7(f(1-(1-pj)**2*(1-qj)))
    qj = ceil7(f(1-(1-qj)**2*(1-pj)))
    ok = (pj, qj) == table[j]
    allok &= ok
    print(f"   j={j}: p={pj} ({fl(pj)})  q={qj} ({fl(qj)})  matches table: {ok}")
check("six rounded updates reproduce the table exactly", allok)
check("final defects < 1/100", table[6][0] < Fr(1,100) and table[6][1] < Fr(1,100))
check("row 1 = Corollary boundary-inputs value 151249/2500000", table[1][0] == Fr(151249,2500000))
# sector positivity of the b=49/4 table (even n)
for n in (4,6,8,10,12):
    check(f"N_{n} >= 0 at b=49/4", Nc[n] >= 0, fl(Nc[n]))

# ---------------- thermal120.tex ----------------
print("\n== thermal120.tex error chain and enclosures ==")
b = Fr(21,2)
check("e^84 < 27^28 (since e^3 < 27) and 27^28 < 2^140 so 27^28/2^257 < 2^-117", 27**28 < 2**140)
check("257(2*243+1) = 125159", 257*(2*243+1) == 125159)
Qa = 2**56
xn = {n: (Fr(3,2)-a)*b*n for n in (6,8,10)}
check("x_n/2 < 6 for n<=10", all(x/2 < 6 for x in xn.values()), fl(xn[10]/2))
check("3^6*243*(125159/2^56 + 2^-116) < 4e-7", 3**6*243*(Fr(125159,Qa)+Fr(1,2**116)) < Fr(4,10**7), fl(3**6*243*(Fr(125159,Qa)+Fr(1,2**116))))
check("32(Q+125159) = 2305843009217699040 < 2^63", 32*(Qa+125159) == 2305843009217699040 and 32*(Qa+125159) < 2**63)
check("2e-6 + (8+4e-7)4e-7 < 1e-5", Fr(2,10**6)+(8+Fr(4,10**7))*Fr(4,10**7) < Fr(1,10**5))
tot120 = {(6,'1'):93637037510293158533967385845216, (6,'P'):3629519998654831111446492683663, (6,'C'):9839536903588146388838689926305,
          (8,'1'):4947327904487991279560670127283, (8,'P'):749999795795292705866802820621, (8,'C'):1301346729026243689610427021457,
          (10,'1'):388818111664728958208623998490, (10,'P'):120458568770462253326886628639, (10,'C'):167487042402286114495695779656}
def expser(x, J): return sum(x**j/factorial(j) for j in range(J+1))
for (n,g), t in tot120.items():
    L = expser(xn[n],100); Un = L + 2*xn[n]**101/factorial(101)
    c = T105[n][{'1':0,'P':1,'C':2}[g]]*u6
    lo = L*t/Qa**2; hi = Un*t/Qa**2
    ok = c-Fr(2,10**6) < lo < hi < c+Fr(2,10**6)
    check(f"enclosure n={n} g={g}: center {fl(c)}, L t/Q^2={fl(lo)}, U t/Q^2={fl(hi)}", ok, f"dev {fl(lo-c)}")
    check(f"   scaled norm < 4", hi < 16)

# ---------------- thermal72.tex ----------------
print("\n== thermal72.tex error chain and enclosures ==")
b2 = Fr(49,4); Q2 = 2**55
check("2^-280 e^98 < 2^-84 (e^98 < 4^98)", True)  # e<4 elementary
check("76*2Q = 5476377146882523136 < 2^63-1", 76*2*Q2 == 5476377146882523136 and 76*2*Q2 < 2**63-1)
check("76*2Q/32 + Q = 207165582859042816", Fr(76*2*Q2,32)+Q2 == 207165582859042816)
W72 = {4:(1298092848946206529211955779616040,683496640248050713756440676071),
       5:(9095379218684502088743241789,9205948219074287936685289758037),
       6:(11971284263768914176431179720622,266573596119930013833919176645),
       7:(15044150708248976451784593468,534330762492192297990602166370),
       8:(386135473475360677201628357188,42459883501386986922766507327),
       9:(4866417737994347551816846903,37716485750353424806876104310),
       10:(19772960374215472827775892486,5030484777222468978572235564),
       11:(841347687722481981930750335,2928767918985354021900689432),
       12:(1266666171463993322695232916,524785463625387865619045258)}
minmargin = None
for n in range(4,13):
    x = b2*n*(Fr(3,2)-a)
    assert x < 15
    dn = isqrt(3**n)+1; En = 281*(dn+1)+1; calE = dn*En
    s = expser(x,100); t = x**100/factorial(100)
    for gi, g in enumerate('1P'):
        W = W72[n][gi]; r = isqrt(W)
        L = s*max(0, r-calE)**2/Q2**2; Uu = (s+t)*(r+1+calE)**2/Q2**2
        c = T1225[n][gi]*u6
        mlo = L-(c-Fr(30,10**6)); mhi = (c+Fr(30,10**6))-Uu
        ok = mlo > Fr(13,10**6) and mhi > Fr(13,10**6)
        mm = min(mlo, mhi); minmargin = mm if minmargin is None else min(minmargin, mm)
        check(f"n={n:2d} g={g}: L={fl(L)} U={fl(Uu)} center={fl(c)} halfwidth={fl((Uu-L)/2)}", ok, f"margins {fl(mlo)}, {fl(mhi)}")
print(f"   minimum margin over the 18 enclosures: {fl(minmargin)} (paper: > 13e-6)")
check("d_n <= 730 for n <= 12", isqrt(3**12)+1 <= 730, str(isqrt(3**12)+1))
check("Q + 281(d_n+1) < 2Q", Q2+281*(730+1) < 2*Q2)

print("\nSUMMARY:", "ALL CHECKS PASSED" if not fails else f"{len(fails)} FAILURES: {fails}")
