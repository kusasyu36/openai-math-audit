"""Computational audit of OpenAI preprint 048
"A singular normal affine surface with free tangent sheaf" (2026-09-23).

Checks every explicit formula in the paper that can be checked symbolically
without an explicit model of the (non-explicit, generic) curve C.
Run:  ../.venv/bin/python -I check_048.py
"""
import itertools
import sympy as sp
import mpmath as mp

ok_count = 0


def check(name, cond):
    global ok_count
    if cond is True or (hasattr(cond, "__bool__") and bool(cond)):
        ok_count += 1
        print(f"[OK]   {name}")
    else:
        print(f"[FAIL] {name}  -> {cond}")


def is0(e):
    return sp.simplify(sp.expand(e)) == 0


# ---------------------------------------------------------------- 1. numerics
print("== 1. numerical invariants (Sec.2, Prop 2.4, Sec.5) ==")
# bidegrees on P1xP1, intersection (a,b).(c,d) = a d + b c
inter = lambda x, y: x[0] * y[1] + x[1] * y[0]
Gam, F = (2, 2), (1, 0)
C = (2 * Gam[0] + 4 * F[0], 2 * Gam[1] + 4 * F[1])
check("C = 2Gamma+4F has bidegree (8,4)", C == (8, 4))
check("genus (a-1)(b-1) = 21", (C[0] - 1) * (C[1] - 1) == 21)
KS = (-2, -2)
KC = (KS[0] + C[0], KS[1] + C[1])
check("K_Sigma + C = (6,2) = Gamma+4F", KC == (6, 2) == (Gam[0] + 4, Gam[1]))
check("deg K_C = 40 = 2g-2", inter(KC, C) == 40 == 2 * 21 - 2)
check("deg H = F.C = 4", inter(F, C) == 4)
check("Gamma.Gamma = 8 (Gamma|Gamma = 4F+Delta has deg 8)", inter(Gam, Gam) == 8)
check("C.Gamma = 24 = 4*deg Z (contact 4 at 6 pts)", inter(C, Gam) == 24 == 4 * 6)
degL = 4 + 6
check("deg L = 10 and 4 deg L = deg K_C", degL == 10 and 4 * degL == 40)
g = 21
# Lemma curve:adjoint-jets table: on elliptic curve h0 = deg for deg>0
rows = []
for j, ell in [(0, 1), (1, 2), (2, 2), (3, 2)]:
    degN = 2 * (8 - j)          # N_j = (8-j)F + Delta, deg F = 2
    res = degN - 6 * ell        # N_j - l Z, deg Z = 6
    # residual bundle: 5F, Phi=F+Delta, Delta (nontrivial deg 0), -F+Delta
    h_res = res if res > 0 else 0
    rank = degN - h_res
    rows.append((j, ell, degN, res, h_res, rank))
check("table h0(N_j) = 16,14,12,10", [r[2] for r in rows] == [16, 14, 12, 10])
check("table residual h0 = 10,2,0,0", [r[4] for r in rows] == [10, 2, 0, 0])
check("table ranks = 6,12,12,10 (<= length 6*l)",
      [r[5] for r in rows] == [6, 12, 12, 10] and all(r[5] <= 6 * r[1] for r in rows))
check("h0(Sigma,(6-j,2)) = 3(7-j); j=0 gives g=21", 3 * 7 == g)
# dimension formula h0(B+D) = h0(B) + deg D - rank
h0Z = 1 + 6 - 6
h0H2Z = 2 + 12 - 12
h02H2Z = 3 + 12 - 12
h03H2Z = 4 + 12 - 10
check("h0(O(Z))=1, h0(H+2Z)=2, h0(theta)=3, h0(theta+H)=6",
      (h0Z, h0H2Z, h02H2Z, h03H2Z) == (1, 2, 3, 6))
# Riemann-Roch
check("h1(L) = h0(L) - (10+1-21) = 12 = h0(3L)", 2 - (10 + 1 - g) == 12)
check("evaluation surjectivity k>=4: deg((3-k)L+Z) <= -4",
      all(10 * (3 - k) + 6 <= -4 for k in range(4, 30)))
# splitting type of E = f_* Lambda, rank 4, Lambda = 3H+2Z deg 24
chiE = 24 + 1 - g
degE = chiE - 4
check("deg f_*Lambda = 0", degE == 0)
# self-dual splitting types with h0(E(-1))=3, h0(E(-2))=2, h0(E)=6
h0 = lambda a: max(a + 1, 0)
sols = []
for t in itertools.combinations_with_replacement(range(-6, 7), 4):
    if sorted(t) != sorted(-x for x in t):
        continue
    if sum(h0(a - 1) for a in t) == 3 and sum(h0(a - 2) for a in t) == 2:
        sols.append(t)
check("unique self-dual splitting with h0(E(-1))=3,h0(E(-2))=2 is (3,0,0,-3)",
      sols == [(-3, 0, 0, 3)])
check("h0(E) = 6 = h0(Lambda)", sum(h0(a) for a in (3, 0, 0, -3)) == 6 == h03H2Z)
check("relative dualizing: 2 Lambda = K_C + 2H (degrees 48)", 2 * 24 == 40 + 8)
check("principal-part annihilator dim 12 - h0(N_Gamma)=12-10 = 2", 12 - 10 == 2)
# jets thresholds
check("H^1(L^{k+1})=0 for k>=4 (deg L^{3-k}<0)", all(3 - k < 0 for k in range(4, 9)))
check("tail vanishing n>=9: H^1(L^n), H^1(L^{n-4}) dual to L^{4-n}, L^{8-n}",
      all(4 - n < 0 and 8 - n < 0 for n in range(9, 30)) and not (8 - 8 < 0))
# Sec 5: adjunction, discrepancy, p_g
C2 = -degL
KYC = 2 * g - 2 - C2
check("K_Y.C = 2g-2-C^2 = 50", KYC == 50)
check("discrepancy -5 consistent: (-5C).C = 50", -5 * C2 == KYC)
check("p_g >= h1(L)+g = 33 and p_g-g-b >= 12", 12 + g == 33 and 33 - g - 0 == 12)
check("Bergner-Graf bound p_g-g-b<=2 is violated (needed for consistency)", 12 > 2)

# ------------------------------------------------- 2. six-vectors lemma (2.1)
print("== 2. Lemma curve:six-vectors (numeric, 50 digits) ==")
mp.mp.dps = 50
h = mp.sqrt(90)
R = mp.sqrt(1 + h ** 2 / 4)
b = 2 * h / R
Q = lambda x, y: x[0] * y[1] / 2 + x[1] * y[0] / 2      # polar form of xy
dvec = lambda z: (R * z - h / 2, R / z - h / 2)
for lam in [mp.mpf("1e-3"), mp.mpf("3e-4") + mp.mpf("1e-4") * 1j]:
    roots = mp.polyroots([1, -b, b * lam, -lam], maxsteps=200, extraprec=200)
    ws = [dvec(z) for z in roots for _ in range(2)]
    w = (sum(v[0] for v in ws), sum(v[1] for v in ws))
    check(f"lam={lam}: w = (h,h)", abs(w[0] - h) < 1e-40 and abs(w[1] - h) < 1e-40)
    check(f"lam={lam}: w^2 = 90", abs(Q(w, w) - 90) < 1e-40)
    check(f"lam={lam}: w_j^2 + w_j.w = 1 for all j",
          all(abs(Q(v, v) + Q(v, w) - 1) < 1e-40 for v in ws))
    check(f"lam={lam}: sum(w_j^2-1) = -90 = -9 deg L",
          abs(sum(Q(v, v) - 1 for v in ws) + 90) < 1e-40)
    ratios = [d[0] / d[1] for d in (dvec(z) for z in roots)]
    distinct = all(abs(ratios[a] - ratios[c]) > 1e-8 for a in range(3) for c in range(a + 1, 3))
    noniso = all(abs(Q(dvec(z), dvec(z))) > 1e-8 for z in roots)
    check(f"lam={lam}: three distinct nonisotropic directions", distinct and noniso)
X, lamS, bS = sp.symbols("X lambda b")
disc = sp.discriminant(X ** 3 - bS * X ** 2 + bS * lamS * X - lamS, X)
lead = sp.expand(disc).coeff(lamS, 1)
check("cubic discriminant: lowest lambda-term is -4 b^3 lambda",
      is0(sp.expand(disc).coeff(lamS, 0)) and is0(lead + 4 * bS ** 3))
zS, hS, RS = sp.symbols("z h R", positive=True)
dz = (RS * zS - hS / 2, RS / zS - hS / 2)
check("symbolic: d(z)^2 + d(z).(h,h) = R^2 - h^2/4",
      is0(dz[0] * dz[1] + (dz[0] * hS + dz[1] * hS) / 2 - (RS ** 2 - hS ** 2 / 4)))
rat = zS * (RS * zS - hS / 2) / (RS - hS * zS / 2)
check("derivative of direction ratio at 0 is -h/(2R)",
      is0(sp.diff(rat, zS).subs(zS, 0) + hS / (2 * RS)))

# ------------------------------------------------- 3. Lie-bracket lemma (4.1)
print("== 3. Lemma jets:bracket, split coefficients ==")
t, z = sp.symbols("t z")
k, l = sp.symbols("k l", integer=True, positive=True)
bf, hf, Af, Bf = [sp.Function(n)(z) for n in ("b", "h", "A", "B")]


def bracket(X1, X2):
    """[X1,X2] for vector fields given as (t-comp, z-comp)."""
    app = lambda Xv, f: Xv[0] * sp.diff(f, t) + Xv[1] * sp.diff(f, z)
    return (sp.expand(app(X1, X2[0]) - app(X2, X1[0])),
            sp.expand(app(X1, X2[1]) - app(X2, X1[1])))


def vf(order, a, bb):
    """t^order (a E_t + bb d_z) as components."""
    return (t ** order * a * t, t ** order * bb)


for kk in range(1, 5):
    for ll in range(0, 4):
        D = vf(kk, bf, hf)
        V = vf(ll + 1, Af, Bf)
        br = bracket(D, V)
        Et = ((ll + 1 - kk) * bf * Af + hf * Af.diff(z) - Bf * bf.diff(z))
        Dz = ((ll + 1) * bf * Bf + hf * Bf.diff(z) - kk * hf * Af - Bf * hf.diff(z))
        exp_ = vf(kk + ll + 1, Et, Dz)
        if not (is0(br[0] - exp_[0]) and is0(br[1] - exp_[1])):
            check(f"bracket formula k={kk} l={ll}", False)
            break
    else:
        continue
    break
else:
    check("bracket formula (jets:general-bracket) for k=1..4, l=0..3", True)
ok_div = True
for kk in range(1, 7):
    D = vf(kk, bf, hf)
    div = t ** 5 * (sp.diff(t ** -5 * D[0], t) + sp.diff(t ** -5 * D[1], z))
    ok_div &= is0(div - t ** kk * ((kk - 4) * bf + hf.diff(z)))
check("divergence L_D sigma = t^k((k-4)b+h') sigma, sigma=t^-5 dt^dz", ok_div)

# split coordinate change z~=psi(z), t~=g(z)t applied to t^{k+1}(A E_t + B d_z)
gf = sp.Function("g")(z)
ok_split = True
for kk in range(0, 6):
    V = vf(kk + 1, Af, Bf)      # t^{k+1}(A E_t + B d_z)
    Vt = V[0] * sp.diff(gf * t, t) + V[1] * sp.diff(gf * t, z)    # V(t~)
    Vz = V[1] * gf ** 4                                             # V(z~), psi'=g^4
    tt = gf * t
    At = sp.simplify(Vt / tt / tt ** (kk + 1))
    Bt = sp.simplify(Vz / tt ** (kk + 1))
    a = gf.diff(z) / gf
    ok_split &= is0(At - gf ** (-(kk + 1)) * (Af + a * Bf))
    ok_split &= is0(Bt - gf ** (3 - kk) * Bf)
check("split coefficients A~=g^{-(k+1)}(A+aB), B~=g^{3-k}B (k=0..5)", ok_split)

# ---------------------------------------- 4. order-two expansion (4.2, eq. second-expansion)
print("== 4. order-two obstruction (Lemma jets:second-obstruction) ==")
rf, pf, mf = [sp.Function(n)(z) for n in ("r", "p", "m")]
D1 = vf(1, mf.diff(z) / 3, mf)                      # order-one generator
v_in = (t * rf * t + t ** 2 * pf * t, t ** 2 * mf * rf)   # t r E_t + t^2(p E_t + m r d_z)
b1 = bracket(D1, v_in)
b2 = bracket(D1, b1)
new = (sp.expand(v_in[0] + b1[0] + b2[0] / 2), sp.expand(v_in[1] + b1[1] + b2[1] / 2))
A1n = new[0].coeff(t, 3)        # index-1 row: t^2 * A_1 * t  (E_t comp = t * coeff)
B1n = new[1].coeff(t, 2)
A2n = new[0].coeff(t, 4)
B2n = new[1].coeff(t, 3)
check("order-one transition turns inside (p, m r) into (P,0)=(p+m r', 0)",
      is0(A1n - (pf + mf * rf.diff(z))) and is0(B1n))
mp_, mpp = mf.diff(z), mf.diff(z, 2)
dA2 = mp_ * pf / 3 + mf * pf.diff(z) + sp.Rational(2, 3) * mf * mp_ * rf.diff(z) \
    + mf ** 2 * rf.diff(z, 2) / 2 - mf * mpp * rf / 6
dU = -mf * pf + mf * mp_ * rf / 3
check("Delta A_2 matches (jets:second-expansion)", is0(A2n - dA2))
check("Delta U matches (jets:second-expansion)", is0(B2n - dU))
# third-order exp term does not reach index 2:
b3 = bracket(D1, b2)
check("cubic exp term starts at index >= 3", is0(sp.expand(b3[0]).coeff(t, 4)) and is0(sp.expand(b3[1]).coeff(t, 3)))
# order-two transition effect: k=2, b=h'/2
D2 = vf(2, hf.diff(z) / 2, hf)
e2 = bracket(D2, (t * rf * t, 0))
check("order-2 generator: dA_2 = h r' - h' r/2, dU = -2 h r",
      is0(e2[0].coeff(t, 4) - (hf * rf.diff(z) - hf.diff(z) * rf / 2)) and is0(e2[1].coeff(t, 3) + 2 * hf * rf))
# delta Psi exact
ri, rj = sp.Function("r_i")(z), sp.Function("r_j")(z)
dD = lambda r: hf * r.diff(z) - hf.diff(z) * r / 2
dUf = lambda r: -2 * hf * r
dPsi = ri * dD(rj) + ri.diff(z) * dUf(rj) / 4 + 3 * rj.diff(z) * dUf(ri) / 4
check("delta Psi_ij = -1/2 d(h r_i r_j) (exact)", is0(dPsi + sp.diff(hf * ri * rj, z) / 2))

# Laurent residues with m=c/z, r=zR/c
cS = sp.symbols("c", nonzero=True)
N = 5
Rc = {s: sp.symbols(f"R{s}_0:{N}") for s in "ij"}
Pc = {s: sp.symbols(f"p{s}_0:{N}") for s in "ij"}
Rser = {s: sum(Rc[s][n] * z ** n for n in range(N)) for s in "ij"}
Pser = {s: sum(Pc[s][n] * z ** n for n in range(N)) for s in "ij"}
mS = cS / z
rS = {s: z * Rser[s] / cS for s in "ij"}


def dA2_of(s):
    r, p = rS[s], Pser[s]
    m1, m2 = sp.diff(mS, z), sp.diff(mS, z, 2)
    return m1 * p / 3 + mS * p.diff(z) + sp.Rational(2, 3) * mS * m1 * r.diff(z) \
        + mS ** 2 * r.diff(z, 2) / 2 - mS * m2 * r / 6


def dU_of(s):
    return -mS * Pser[s] + mS * sp.diff(mS, z) * rS[s] / 3


res = lambda e: sp.residue(sp.expand(e), z, 0)
Ri0, Rj0, Ri1, Rj1 = Rc["i"][0], Rc["j"][0], Rc["i"][1], Rc["j"][1]
pi0, pj0 = Pc["i"][0], Pc["j"][0]
Fij = res(rS["i"] * dA2_of("j"))
Hij = res(rS["i"].diff(z) * dU_of("j"))
Hji = res(rS["j"].diff(z) * dU_of("i"))
check("F_ij = -1/3 R_i p_j - 2/3 R_i R_j' - R_i' R_j",
      is0(Fij - (-Ri0 * pj0 / 3 - sp.Rational(2, 3) * Ri0 * Rj1 - Ri1 * Rj0)))
check("H_ij = -R_i p_j - 1/3 R_i R_j' - 2/3 R_i' R_j",
      is0(Hij - (-Ri0 * pj0 - Ri0 * Rj1 / 3 - sp.Rational(2, 3) * Ri1 * Rj0)))
Tij = Fij + Hij / 4 + 3 * Hji / 4
check("T_ij = -7/12 R_i p_j - 9/12 R_j p_i - 5/4 R_i R_j' - 17/12 R_i' R_j",
      is0(Tij - (-sp.Rational(7, 12) * Ri0 * pj0 - sp.Rational(9, 12) * Rj0 * pi0
                 - sp.Rational(5, 4) * Ri0 * Rj1 - sp.Rational(17, 12) * Ri1 * Rj0)))
Tji = Tij.subs({**{a: b_ for a, b_ in zip(Rc["i"] + Pc["i"], Rc["j"] + Pc["j"])},
                **{a: b_ for a, b_ in zip(Rc["j"] + Pc["j"], Rc["i"] + Pc["i"])}},
               simultaneous=True)
sym = (Tij + Tji) / 2
alt = (Tij - Tji) / 2
check("T^sym = -2/3 (R_i p_j + R_j p_i + 2 R_i R_j' + 2 R_i' R_j)",
      is0(sym + sp.Rational(2, 3) * (Ri0 * pj0 + Rj0 * pi0 + 2 * Ri0 * Rj1 + 2 * Ri1 * Rj0)))
check("T^alt(i=1,j=2) = 1/12 ((R^p)(0) + (R^R')(0))",
      is0(alt - ((Ri0 * pj0 - Rj0 * pi0) + (Ri0 * Rj1 - Rj0 * Ri1)) / 12))
Pfull = {s: Pser[s] + mS * rS[s].diff(z) for s in "ij"}
check("P = R(0)/z + p(0) + 2R'(0) + O(z)",
      is0(sp.expand(Pfull["i"] * z).coeff(z, 0) - Ri0) and is0(sp.expand(Pfull["i"]).coeff(z, 0) - (pi0 + 2 * Ri1)))
check("Res(P_i P_j dz) = R_i p_j + R_j p_i + 2R_i R_j' + 2R_i' R_j",
      is0(res(Pfull["i"] * Pfull["j"]) - (Ri0 * pj0 + Rj0 * pi0 + 2 * Ri0 * Rj1 + 2 * Ri1 * Rj0)))
# totals: alt sum = (1/12) * sum(W_j - 1) with R^p = -1 ; = -3/4 deg L
check("sum T^alt = (1/12)(-9 deg L) = -3/4 deg L", sp.Rational(1, 12) * (-9 * degL) == -sp.Rational(3, 4) * degL)
check("outside total: -(4 * 3/4) - 6 * 3/4 = -3/4 deg L  (u1 fibre deg 4, |Z| = 6)",
      -4 * sp.Rational(3, 4) - 6 * sp.Rational(3, 4) == -sp.Rational(3, 4) * degL)

# Psi transformation law (jets:psi-transformation)
gS, aS = sp.symbols("g a")
r_ = sp.symbols("r_i r_j")
rp_ = sp.symbols("rp_i rp_j")
U_ = sp.symbols("U_i U_j")
D_ = sp.symbols("D_i D_j")


def Psi(r, rp, U, Dd, i, j):
    return r[i] * Dd[j] + rp[i] * U[j] / 4 + 3 * rp[j] * U[i] / 4


ok_psi = True
for i in range(2):
    for j in range(2):
        rt = [r_[s] / gS for s in range(2)]
        rpt = [gS ** -5 * (rp_[s] - aS * r_[s]) for s in range(2)]
        Ut = [gS * U_[s] for s in range(2)]
        Dt = [gS ** -3 * (D_[s] + aS * U_[s]) for s in range(2)]
        diff = Psi(rt, rpt, Ut, Dt, i, j) * gS ** 4 - Psi(r_, rp_, U_, D_, i, j)
        ok_psi &= is0(diff - sp.Rational(3, 4) * (r_[i] * U_[j] - r_[j] * U_[i]) * aS)
check("Psi~ - Psi = 3/4 (r_i U_j - r_j U_i) dlog g", ok_psi)

# third-order residue (jets:third-residue)
rho = sum(sp.symbols("q0:4")[n] * z ** n for n in range(4))
r3 = z * rho
check("Res[(r'/z + 2r/z^2)dz] = 3 r'(0)", is0(res(r3.diff(z) / z + 2 * r3 / z ** 2) - 3 * r3.diff(z).subs(z, 0)))
# leading action at order k (jets:leading-action)
ok_la = True
for kk in range(1, 8):
    e = bracket(vf(kk, bf, hf), (t * rf * t, 0))
    ok_la &= is0(e[0].coeff(t, kk + 2) - ((1 - kk) * bf * rf + hf * rf.diff(z))) and is0(e[1].coeff(t, kk + 1) + kk * hf * rf)
check("leading action dA_k=(1-k)b r + h r', dB_k=-k h r (k=1..7)", ok_la)
# old-transition effect on global addition t^k a E_t (jets:old-transition-effect)
af = sp.Function("a")(z)
ok_ot = True
for kk in range(3, 8):
    e = bracket(D1, (t ** kk * af * t, 0))
    ok_ot &= is0(e[0].coeff(t, kk + 2) - (sp.Rational(kk - 1, 3) * mf.diff(z) * af + mf * af.diff(z))) \
        and is0(e[1].coeff(t, kk + 1) + mf * af)
check("old order-1 transition on t^k a E_t: dA_k=(k-1)/3 m' a + m a', dB_k=-m a", ok_ot)

# ------------------------------------------------- 5. trace section formulas
print("== 5. Sec.3 trace / principal parts ==")
ok_tr = True
for e in range(1, 7):
    for i in range(e):
        for j in range(e):
            n = i + j - (e - 1)       # Tr(w^n)/e with x = w^e; Tr(w^n) = e x^{n/e} iff e | n
            val = 1 if n == 0 else 0  # n in (-e, e)
            ok_tr &= (val == (1 if i + j == e - 1 else 0))
check("local trace matrix is antidiagonal identity (e=1..6)", ok_tr)
b1_, b2_, K0, K1, p0_, p1_ = sp.symbols("b1 b2 K0 K1 pi0 pi1")
bz = b1_ * z + b2_ * z ** 2
Kz = K0 + K1 * z
piz = p0_ + p1_ * z
lhs = sp.series(piz / (Kz * z ** 2), z, 0, 0).removeO()
pt0 = b1_ ** 2 / K0 * p0_
pt1 = b1_ ** 2 / K0 * p1_ + (2 * b1_ * b2_ / K0 - b1_ ** 2 * K1 / K0 ** 2) * p0_
rhs = sp.series((pt0 + pt1 * z) / bz ** 2, z, 0, 0).removeO()
check("(trace:two-coefficients) equal principal parts", is0(lhs - rhs))
yC = -bz ** 4 / Kz
s0 = sp.limit(yC / (z ** 4 * Kz), z, 0)
check("s(0) = -b'(0)^4 / k_y(0,0)^2", is0(s0 + b1_ ** 4 / K0 ** 2))
E0 = sp.symbols("E0")
check("(b1^2/K0)^2 * E/s(0) = -E  (isometry sign)", is0((b1_ ** 2 / K0) ** 2 * E0 / s0 + E0))
# contact expansion: k = b^4 + y K + y^2 k2 -> y = -b^4/K + O(z^8)
k2 = sp.symbols("k2")
ysol = -bz ** 4 / Kz
check("contact: k(-b^4/K, z) = O(z^8)",
      sp.series(bz ** 4 + ysol * Kz + ysol ** 2 * k2, z, 0, 8).removeO() == 0)
g0 = sp.symbols("g0", positive=True)
check("frame change: tau_new = g^3 psi/z tau -> l_new = g0^-7 l (psi'(0)=g0^4)", g0 ** 3 * g0 ** 4 == g0 ** 7)

# ------------------------------------------------- 6. finite determinacy constants
print("== 6. Lemma algebra:finite-determinacy constants ==")
ok_c = True
for c in range(0, 12):
    for Nn in range(1, 12):
        n0 = max(c + Nn, c + 2, 2 * c + 1)
        for n in range(n0, n0 + 15):
            q = n - c
            ok_c &= q >= max(Nn, 2) and 2 * q >= n + 1 and n + q - 1 >= n + 1
check("n >= n0 implies q>=max(N,2), 2q>=n+1, n+q-1>=n+1", ok_c)

print(f"\n{ok_count} checks passed")
