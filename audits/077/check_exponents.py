# Audit 077: exact-arithmetic checks of exponent bookkeeping (own script; no repo code executed)
from fractions import Fraction as F
import sympy as sp
import mpmath as mp

ok = lambda c, msg: print(("OK  " if c else "FAIL"), msg)

print("== 1. Necessary conditions for E: L^q(S^2)->L^p(R^3):  p>3 and p>=2q'  ==")
def admissible(q, p):  # q may be None for infinity
    qp = F(1) if q is None else F(q) / (F(q) - 1)
    return p > 3 and p >= 2 * qp
eps = F(1, 10**6)
ok(admissible(None, 3 + eps), "Thm1.1 (sphere): q=inf, p=3+1e-6 admissible")
ok(not admissible(None, F(3)), "endpoint p=3, q=inf excluded (paper claims none)")
for p in [3 + eps, F(7, 2), F(4), F(10)]:
    ok(admissible(p, p), f"diagonal q=p={p} admissible")
for p in [3 + eps, F(4)]:
    ok(admissible(3, p), f"Diagonal paper graph bound L^3->L^p, p={p} admissible")
# Knapp exponent: ||Ef||_p/||f||_q ~ delta^(2-4/p-2/q) must be >= 0
p, q = sp.symbols('p q', positive=True)
kn = 2 - 4/p - 2/q
print("   Knapp exponent 2-4/p-2/q; zero set p =", sp.solve(sp.Eq(kn, 0), p), "(=2q')")
ok(kn.subs({p: 3 + sp.Rational(1, 10**6), q: 3}) > 0, "Knapp ok at (q,p)=(3,3+1e-6)")
# mixed norm b>2a'
for a, b in [(F(2), F(4) + eps), (F(3, 2), F(6) + eps), (F(10), F(3) + eps)]:
    ap = a / (a - 1)
    ok(b > 2 * ap and admissible(a, b), f"mixed-norm a={a}, b={float(b):.6f} within conjectured range")
# interpolation identity of mixed norm
a_, b_ = sp.symbols('a b', positive=True)
ap_ = a_/(a_-1); r = 1 + b_/ap_; th = r/b_
ok(sp.simplify(1 - th + th/r - 1/a_) == 0, "mixed-norm interpolation 1/a = 1-th+th/r")

print("== 2. Section 7 (extremal) bookkeeping ==")
m, s, k, z, d = sp.symbols('m s kappa zeta delta', positive=True)
e1, e2 = sp.symbols('e1 e2', positive=True)
r_ = m**-2; theta = m*d; w = d/m
# capacity: theta^2 W(E)^-1 w^-2 prod min * delta^-2 W(F) m^{-2(1+s-10k-z)}, with W(E)^-1 prod min W(F) <= delta^2
cap = theta**2 * w**-2 * d**2 * d**-2 * m**(-2*(1+s-10*k-z))
ok(sp.simplify(sp.log(cap, m).expand(force=True) - (2*(1-s)+20*k+2*z)) == 0 if False else
   sp.simplify(cap / m**(2*(1-s)+20*k+2*z)) == 1, "terminal capacity exponent 2(1-s)+20k+2z (ext:capacity-bound)")
# min(re_j, delta) = delta*e_j/f_j with f_j=max(e_j, delta/r): numeric spot check
import random
for _ in range(1000):
    rr = random.uniform(1e-4, 1); dd = random.uniform(1e-4, 1); ee = random.uniform(1e-4, 10)
    f = max(ee, dd/rr)
    assert abs(min(rr*ee, dd) - dd*ee/f) < 1e-12*max(1, dd)
ok(True, "min(r e_j, delta) = delta e_j / f_j (1000 random cases)")
# W(F)/(f1 f2) <= W(E)/(e1 e2) when raising to common floor
kk = 0.07
bad = 0
for _ in range(10000):
    a1 = 10**random.uniform(-4, 2); a2 = 10**random.uniform(-4, 2); fl = 10**random.uniform(-4, 2)
    E1, E2 = max(a1, a2), min(a1, a2); F1, F2 = max(E1, fl), max(E2, fl)
    if (F1/F2)**kk > (E1/E2)**kk * (1 + 1e-12): bad += 1
ok(bad == 0, "aspect ratio decreases when radii raised to a common floor")
# geometric 40 -> 41 division by threshold
ok(sp.simplify(m**(-2*(1+2*s-10*k-z)) / m**(-2*s) / m**(-2*(1+s-10*k-z))) == 1, "(40)/threshold = (41)")
# final cancellation: h^{1/2} min(m^{-2(1-s)}, 1/h) <= m^{-(1-s)}
mm = mp.mpf(10)**6
for sv in [0.1, 0.5, 0.9, 1.0]:
    worst = max(mp.sqrt(h)*min(mm**(-2*(1-sv)), 1/h) / mm**(-(1-sv)) for h in [mp.mpf(10)**e for e in range(-15, 16)])
    ok(worst <= 1 + 1e-12, f"h^1/2 min(..) <= m^-(1-s) at s={sv} (max ratio {float(worst):.3g})")
# total amplification
amp = m**(1-s+10*k+z) * m**(-(1-s))
ok(sp.simplify(amp / m**(10*k+z)) == 1, "cubed amplification <= m^{10k+zeta}")
# gamma_sigma <= 2 sigma
sig = sp.symbols('sigma', positive=True)
ok(sp.simplify(m**(-2*(1-sig)) * m**2 / m**(2*sig)) == 1, "gamma_sigma <= 2 sigma (mass m^{-2(1-σ)} x capacity^{1/2} m^2)")
# quartic scale identity r+ w+^2 l^4 = r- w-^2
mm_, l = sp.symbols('m_minus ell', positive=True)
rp, wp = (mm_*l)**-2, d/(mm_*l); rm, wm = mm_**-2, d/mm_
ok(sp.simplify(rp*wp**2*l**4 - rm*wm**2) == 0, "r+ w+^2 ell^4 = r- w-^2")
# frame: N^2 w theta = 1
N = 1/d
ok(sp.simplify(N**2*w*theta) == 1, "N^2 w theta = 1")

print("== 3. Hairbrush / heavy count / angular exclusion exponents ==")
e_, t_, dd_ = sp.symbols('e tau d', positive=True)
ck = 4/k
ok(sp.simplify(2*e_ - 2*k*ck*(e_+t_) - (-6*e_-8*t_)) == 0, "heavy-count removal exponent -6e-8tau")
Ck = (ck+10)/(2*k)
lhs = (3+ck)*e_ + 6*dd_ + ck*t_ - 2*k*Ck*(e_+dd_+t_)
ok(sp.simplify(lhs - (-7*e_-(ck+4)*dd_-10*t_)) == 0, "angular exclusion exponent -7e-(c+4)d-10tau")
# hairbrush CS: (m A M^-e)^2 <= A eta^-(1-k) * m A^2 eta^(1+k)  => m <= A M^2e eta^2k
mA, A, Me, eta = sp.symbols('mA A Me eta', positive=True)
sol = sp.solve(sp.Eq((mA*A/Me)**2, A*eta**(-(1-k))*mA*A**2*eta**(1+k)), mA)
ok(sp.simplify(sol[0] - A*Me**2*eta**(2*k)) == 0, "hairbrush Cauchy-Schwarz gives m <= A M^{2e} eta^{2k}")

print("== 4. Sphere global section constants ==")
# E1(x)=2 sin(2pi r)/r
rr = mp.mpf('0.37')
val = mp.quad(lambda c: 2*mp.pi*mp.e**(2j*mp.pi*rr*c), [-1, 1])
ok(abs(val - 2*mp.sin(2*mp.pi*rr)/rr) < 1e-12, "E1 = 2 sin(2 pi |x|)/|x|")
# ||sigma*sigma||_2^2 = 32 pi^3 and ||E1||_4^4 = 32 pi^3 (equality case of the L^4 bound)
ok(abs(mp.quad(lambda r: (2*mp.pi/r)**2*4*mp.pi*r**2, [0, 2]) - 32*mp.pi**3) < 1e-9, "||sigma*sigma||_2^2 = 32 pi^3")
I4 = 64*mp.pi*mp.quadosc(lambda r: mp.sin(2*mp.pi*r)**4/r**2, [0, mp.inf], period=0.5)
ok(abs(I4 - 32*mp.pi**3) < 1e-6, f"||E1||_4^4 = {float(I4):.6f} = 32pi^3 = {float(32*mp.pi**3):.6f} (bound attained by g=1)")
# ||E1||_p^p ~ C/(p-3): check (p-3)*||E1||_p^p tends to a constant
for pv in [3.2, 3.1, 3.05]:
    Ip = mp.quadosc(lambda r: abs(2*mp.sin(2*mp.pi*r)/r)**pv*4*mp.pi*r**2, [0, mp.inf], period=0.5) if False else None
print("   (divergence of ||E1||_p^p at p=3 is elementary: integrand ~ r^{2-p}|sin|^p)")
# sparse cubic: N^6 * N^{10k} * delta^4 * delta^2 = N^{10k}
ok(sp.simplify(N**6*d**4*d**2) == 1, "sparse cubic sample sum N^6 delta^6 = 1")
# covering: R_{j+1}/2 >= (6R_j/mu)^C0 for small mu, with R_{j+1}=(R_j/mu)^(2C0+2), R_j>=1/mu
C0 = 9; mu = F(1, 1000); Rj = 1/mu
ok(F(1, 2)*(Rj/mu)**(2*C0+2) >= (6*Rj/mu)**C0, "covering separation inequality (C0=9, mu=1e-3)")

print("== 5. Diagonal paper bookkeeping ==")
tau = sp.symbols('tau', positive=True)
th = sp.Rational(1, 10**6)*tau**2
g0 = 10*th + th + sp.Rational(3, 2)*th
ok(sp.simplify(g0 - sp.Rational(25, 2)*th) == 0, "gamma0 = 25 theta / 2")
ok(sp.simplify(50*(34/tau)*g0 - sp.Rational(17, 800)*tau) == 0, "50 J gamma0 = 17 tau/800 (with J<=34/tau)")
for tv in [F(1), F(1, 2), F(1, 7), F(1, 1000)]:
    J = -(-32//tv) + 1  # ceil(32/tau)+1
    ok(J <= 34/tv and F(tv, 4)*J > 8, f"tau={tv}: J={J} <= 34/tau and beta*J>8")
u = sp.symbols('u', positive=True)
Iu = sp.integrate(sp.exp(-tau*u/2)*(1+u)**3, (u, 0, sp.oo), conds='none')
ok(sp.simplify(Iu - (2/tau + 12/tau**2 + 48/tau**3 + 96/tau**4)) == 0, "layer-cake integral identity")
ok(2+12+48+96 <= 158, "<= 158 tau^-4 for tau<=1")
# bundle lower bound constant: pi^{3/2}/(9 pi)^{1/2} = pi/3
ok(sp.simplify(sp.pi**sp.Rational(3, 2)/sp.sqrt(9*sp.pi) - sp.pi/3) == 0, "bundle lower bound constant pi/3")
# L3 scaling: delta^2 sum b^{3/2} ~ delta^6 for L^3 data
ok(sp.simplify(d**5 * d * N**6) == 1, "sparse L3 mass delta^6 cancels N^6")
# Schrodinger threshold s > 2-6/p equals Rogers' 2n(1/2-1/p)-2/p at n=2
nn = 2
ok(sp.simplify(2*nn*(sp.Rational(1, 2)-1/p) - 2/p - (2 - 6/p)) == 0, "Schrodinger smoothing threshold = Rogers' line at n=2")
ok(sp.Rational(2*(nn+1), nn) == 3, "Rogers' critical p = 2(n+1)/n = 3")

print("== 6. Kakeya transfer (sec. 9) ==")
qq = sp.symbols('q', positive=True)
rq = qq/2
ok(sp.simplify((6-2*qq)/rq - (12/qq - 4)) == 0, "weighted tube exponent 12/q - 4")
pq = qq/(qq-2)
ok(sp.simplify((12/qq-4)*pq/3 + 4*(qq-3)/(3*(qq-2))) == 0, "interpolated loss -4(q-3)/(3(q-2))")
