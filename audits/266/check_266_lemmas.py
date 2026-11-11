"""Exact / high-precision checks of small lemmas and constants in OpenAI preprint 266
(main paper) and in the companion paper. Independent code; uses sympy/mpmath only."""
import sympy as sp, mpmath as mp, itertools, random
from fractions import Fraction as Fr
mp.mp.dps = 40
ok_all = True
def report(name, ok, detail=""):
    global ok_all
    ok_all &= bool(ok)
    print(("OK   " if ok else "FAIL ") + name + (("  " + detail) if detail else ""))

# ---- introduction.tex:166-210  three MUBs in C^6 (tensor construction) ----
I2 = sp.eye(2); A1 = sp.Matrix([[1,1],[1,-1]])/sp.sqrt(2); A2 = sp.Matrix([[1,1],[sp.I,-sp.I]])/sp.sqrt(2)
w = (-1 + sp.sqrt(3)*sp.I)/2
F3 = sp.Matrix(3,3,lambda j,k: w**(j*k))/sp.sqrt(3); D3 = sp.diag(1,w,w)
Us = [sp.kronecker_product(a,b) for a,b in [(I2,sp.eye(3)),(A1,F3),(A2,D3*F3)]]
good = True
for U in Us:
    good &= sp.expand(U.H*U - sp.eye(6)) == sp.zeros(6)
for r in range(3):
    for s in range(r+1,3):
        M = Us[r].H*Us[s]
        good &= all(sp.expand(M[i,j]*sp.conjugate(M[i,j])) == sp.Rational(1,6) for i in range(6) for j in range(6))
report("Lemma three-mub-six: U_r unitary and |(U_r^*U_s)_ij|^2 = 1/6", good)
# the F^*DF entry claim: (1+2w) for t=0, (2+w^2) for t=1,2, squared modulus 3
report("F^*DF sums: |1+2w|^2=3, |2+w^2|^2=3", sp.expand((1+2*w)*sp.conjugate(1+2*w))==3 and sp.expand((2+w**2)*sp.conjugate(2+w**2))==3)

# ---- rounding.tex:151-179  pi literal and Machin certificate ----
pihat = sp.Rational(884279719003555, 281474976710656)
report("0x1.921fb54442d18p+1 == 884279719003555/281474976710656", sp.Rational(0x1921fb54442d18, 2**51) == pihat)
def Aser(x, n): return sum(sp.Rational((-1)**j) * x**(2*j+1) / (2*j+1) for j in range(n))
s = 16*Aser(sp.Rational(1,5),12) - 4*Aser(sp.Rational(1,239),4)
diff = pihat - s
report("pihat - s rational (rounding.tex:171-173)", diff == sp.Rational(-1204334632739641172855772829555432415587, 10003273510806667494906565283020800000000000000000000000), str(diff)[:60])
rem = sp.Rational(16, 25*5**25) + sp.Rational(4, 9*239**9)
bound = abs(diff) + rem
report("|pihat - pi| < 1.225418e-16 (rounding.tex:178)", bound < sp.Rational(1225418, 10**22), f"certified bound {float(bound):.7e}; true {float(abs(pihat - sp.pi.evalf(50))):.7e}")
# series tails (rounding.tex:231-233)
y = mp.mpf('0.786')
report("cos tail y^20/20! < 3.330e-21, sin tail y^21/21! < 1.247e-22 at y=.786", y**20/mp.factorial(20) < mp.mpf('3.330e-21') and y**21/mp.factorial(21) < mp.mpf('1.247e-22'), f"{mp.nstr(y**20/mp.factorial(20),5)} {mp.nstr(y**21/mp.factorial(21),5)}")
report("(1+u)^27-1 < 5.5e-15, (1+u)^9-1 < 1.81e-15 with u=2e-16", (1+mp.mpf('2e-16'))**27-1 < mp.mpf('5.5e-15') and (1+mp.mpf('2e-16'))**9-1 < mp.mpf('1.81e-15'))
report("gamma_m: (1+u)^(m+1)-1 < 2.1e-16(m+1) for m<=2000", all((1+mp.mpf('2e-16'))**(m+1)-1 < mp.mpf('2.1e-16')*(m+1) for m in (1,10,100,1000,2000)))
report("sqrt2*(2u+u^2) < 1e-15", mp.sqrt(2)*(2*mp.mpf('2e-16')+mp.mpf('4e-32')) < mp.mpf('1e-15'))

# ---- seed.tex: finite counts (execution.tex:191-199) ----
report("C(116,5)=160,389,488 and C(115,4)=6,913,340", sp.binomial(116,5)==160389488 and sp.binomial(115,4)==6913340)
report("3/224 < 0.032 (seed.tex:477)", sp.Rational(3,224) < sp.Rational(32,1000))
# Lemma sum-test constant: 5 remainders pi*Delta^2 with |Delta|<=w/2 -> 5 pi w^2/4
report("sum-test remainder 5*pi*(w/2)^2 = 5 pi w^2/4", sp.simplify(5*sp.pi*(sp.Symbol('w')/2)**2 - 5*sp.pi*sp.Symbol('w')**2/4)==0)

# ---- completion.tex: A_0 = ||E_0||_op = sqrt(2+|1+a|) ----
import numpy as np
okE = True
for thv in np.linspace(0, 2*np.pi, 37):
    a_ = np.exp(1j*thv); E0n = np.array([[1,1],[1,a_]]); evn = np.sort(np.linalg.eigvalsh(E0n@E0n.conj().T))
    okE &= np.allclose(evn, [2-abs(1+a_), 2+abs(1+a_)])
report("E_0E_0^* eigenvalues are 2 +- |1+a| (so A_0 = sqrt(2+|1+a|))", okE)
# Lemma frames: M=(1 p), M^*M = [[4,S],[conj S,4]], and M(M^*M)^{-1} formula
random.seed(1)
p = sp.Matrix([sp.exp(sp.I*sp.Rational(random.randint(1,999),1000)*2*sp.pi) for _ in range(4)])
one = sp.ones(4,1); M = one.row_join(p); S = sum(p)
MM = sp.simplify(M.H*M)
report("M^*M = [[4,S],[S*,4]]", sp.simplify(MM - sp.Matrix([[4,S],[sp.conjugate(S),4]])) == sp.zeros(2))
lhs = sp.simplify(M*MM.inv()); rhs = sp.simplify(((4*one - p*sp.conjugate(S)).row_join(4*p - one*S))/(16-abs(S)**2))
report("M(M^*M)^{-1} = (4*1 - p S*, 4p - 1 S)/(16-|S|^2)", sp.simplify(lhs-rhs).norm() < 1e-20 or sp.simplify((lhs-rhs).evalf(30)).norm() < 1e-25)
# Lemma polar: sqrt6 - sqrt(6-l) = l/(sqrt6+sqrt(6-l))
l = sp.Symbol('l', positive=True)
report("sqrt6 - sqrt(6-l) = l/(sqrt6+sqrt(6-l))", sp.simplify((sp.sqrt(6)-sp.sqrt(6-l)) - l/(sp.sqrt(6)+sp.sqrt(6-l))) == 0)
# rounding.tex:469-471 scale examples
report("q<.32: .203^2(1+10/1.239^2)+1e-7 < .32", mp.mpf('.203')**2*(1+10/mp.mpf('1.239')**2)+1e-7 < 0.32, mp.nstr(mp.mpf('.203')**2*(1+10/mp.mpf('1.239')**2),6))
report("off<.74: .203(sqrt6+2)/1.239+1e-7 < .74", mp.mpf('.203')*(mp.sqrt(6)+2)/mp.mpf('1.239')+1e-7 < .74)
report("pol<.068: .32/(sqrt6+sqrt(6-.32))+1e-7 < .068", mp.mpf('.32')/(mp.sqrt(6)+mp.sqrt(6-mp.mpf('.32')))+1e-7 < .068)
report("R_s = 2*pi*.0320001 < .203 and h=pi/16 < .197", 2*mp.pi*mp.mpf('.0320001') < .203 and mp.pi/16 < .197, mp.nstr(2*mp.pi*mp.mpf('.0320001'),6))
# at an exact Hadamard, |sum b_j| = |1+a| <= 2 so sigma_min(B_0)^2 = 4-|S| >= 2 > 1.24^2
report("exact seed: sigma_min(B_0) = sqrt(4-|1+a|) >= sqrt2 > 1.24", mp.sqrt(2) > 1.24)

# ---- completion.tex:212-251  U(2) cell Taylor bounds (numerical sup over random cells) ----
import numpy as np
rng = np.random.default_rng(0)
def Umat(al, T, be, ga):
    return np.exp(1j*al)*np.array([[np.exp(1j*be)*np.cos(T), np.exp(1j*ga)*np.sin(T)],[-np.exp(-1j*ga)*np.sin(T), np.exp(-1j*be)*np.cos(T)]])
def dU(al,T,be,ga,th):
    x,y,z,wv = th
    U = Umat(al,T,be,ga)
    dA = 1j*U
    dT = np.exp(1j*al)*np.array([[-np.exp(1j*be)*np.sin(T), np.exp(1j*ga)*np.cos(T)],[-np.exp(-1j*ga)*np.cos(T), -np.exp(-1j*be)*np.sin(T)]])
    dB = np.array([[1j*U[0,0],0],[0,-1j*U[1,1]]])
    dG = np.array([[0,1j*U[0,1]],[-1j*U[1,0],0]])
    return x*dA + y*dT + z*dB + wv*dG
worst1 = worst2 = 0
for _ in range(20000):
    h = rng.uniform(0.001, np.pi/16)
    base = rng.uniform(0, 2*np.pi, 4); th = rng.uniform(-h, h, 4)
    if rng.random() < 0.3: th = h*rng.choice([-1,1],4)
    U0 = Umat(*base); V = Umat(*(base+th))
    worst1 = max(worst1, np.linalg.norm(V-U0)/(np.sqrt(6)*h))
    worst2 = max(worst2, np.linalg.norm(V-U0-dU(*base,th))/(np.sqrt(42)/2*h**2))
report("||V-U||_F <= sqrt6 h and ||V-U-dU||_F <= sqrt42/2 h^2 (20000 random cells)", worst1 <= 1 and worst2 <= 1, f"max ratios {worst1:.4f}, {worst2:.4f}")
# second-derivative constant 42: sup ||U''||^2/h^4 along straight segments with |theta_j|<=h
worst = 0
for _ in range(20000):
    base = rng.uniform(0,2*np.pi,4); th = rng.choice([-1.0,1.0],4)*rng.uniform(0.5,1,4)
    eps = 1e-4
    Upp = (Umat(*(base+eps*th)) - 2*Umat(*base) + Umat(*(base-eps*th)))/eps**2
    worst = max(worst, np.linalg.norm(Upp)**2/np.max(np.abs(th))**4)
report("sup ||U''||_F^2 / h^4 <= 42 (finite differences)", worst <= 42*(1+1e-6), f"observed sup {worst:.4f}")

# ---- phase-cover.tex:206-222  horizontal-row identities ----
worst6 = 0; worst8 = 0
for _ in range(5000):
    x = rng.standard_normal((6,6)); x -= x.mean(axis=0, keepdims=True); x -= x.mean(axis=1, keepdims=True)
    s2 = (x**2).sum(); a = [((x[i]-x[j])**2).sum() for i in range(6) for j in range(i)]
    worst6 = max(worst6, abs(sum(a) - 6*s2)/s2**1); worst8 = max(worst8, sum(v*v for v in a)/s2**2)
report("sum a_ij = 6 s^2 and sum a_ij^2 <= 8 s^4 for horizontal x", worst6 < 1e-9 and worst8 <= 8, f"max sum a^2/s^4 = {worst8:.4f}")
# residual identities ||AA^*-6I||_F = 2 pi sqrt2 ||F(P)|| and ||F||^2 = (||AA^*||_F^2-216)/(8 pi^2)
P = rng.uniform(0,1,(6,6)); A = np.exp(2j*np.pi*P)
Fv = np.array([(A[i]@A[j].conj())/(2j*np.pi) for i in range(6) for j in range(i)])
nF = np.sqrt((abs(Fv)**2).sum())
report("||AA*-6I||_F = 2 pi sqrt2 ||F||  and  ||F||^2 = (||AA*||_F^2-216)/(8pi^2)", abs(np.linalg.norm(A@A.conj().T-6*np.eye(6)) - 2*np.pi*np.sqrt(2)*nF) < 1e-10 and abs(nF**2 - (np.linalg.norm(A@A.conj().T)**2-216)/(8*np.pi**2)) < 1e-10)
report("transpose invariance ||F(P^T)|| = ||F(P)||", abs(np.linalg.norm(A.T@A.conj()-6*np.eye(6)) - np.linalg.norm(A@A.conj().T-6*np.eye(6))) < 1e-10)

# ---- exclusion.tex: bin geometry constants ----
best = {"X":0,"L":0,"I":0}
aa = 1.0
for signs in itertools.product([-1,1], repeat=5):
    d = np.array((0,)+signs, float)*aa
    x = d - d.sum()/6
    best["X"] = max(best["X"], np.linalg.norm(x)); best["L"] = max(best["L"], abs(x).sum()); best["I"] = max(best["I"], abs(x).max())
report("X_M=sqrt(29/6)a, L_M=5a, I_M=1.5a attained at vertices", abs(best["X"]-np.sqrt(29/6))<1e-12 and abs(best["L"]-5)<1e-12 and abs(best["I"]-1.5)<1e-12, str(best))
report("distinct labels: 6-2pi sqrt6 *0.2 > 20/7 > sqrt6", 6-2*mp.pi*mp.sqrt(6)*mp.mpf('0.2') > mp.mpf(20)/7 > mp.sqrt(6), mp.nstr(6-2*mp.pi*mp.sqrt(6)*mp.mpf('0.2'),6))
report("K: sqrt(1+(2 pi sqrt2/6)1e-6) < 1.00000075 < 1.00001", mp.sqrt(1+(2*mp.pi*mp.sqrt(2)/6)*mp.mpf('1e-6')) < mp.mpf('1.00000075'))
# Edge Taylor constants: R <= 2 pi^2 r^2 ; pairing with S: 2|S||R|/(24 pi) <= pi |S| r^2/6 ; increment^2/(24 pi) <= pi r^2
r_, S_ = sp.symbols('r S', positive=True)
report("E_U pieces: 2S(2pi^2r^2)/(24pi) = pi S r^2/6 and (2pi sqrt6 r)^2/(24pi) = pi r^2", sp.simplify(2*S_*2*sp.pi**2*r_**2/(24*sp.pi) - sp.pi*S_*r_**2/6)==0 and sp.simplify((2*sp.pi*sp.sqrt(6)*r_)**2/(24*sp.pi) - sp.pi*r_**2)==0)
report("single-vector Taylor: |e(t)-1-2 pi i t| <= 2 pi^2 t^2 (numeric)", all(abs(mp.expj(2*mp.pi*t)-1-2j*mp.pi*t) <= 2*mp.pi**2*t**2 for t in [mp.mpf(k)/97 for k in range(-40,41)]))
report("arc: |d| = 2 arcsin(c/2) <= c (1-q/4)^(-1/2) for c^2<=q<3.8", all(2*mp.asin(c/2) <= c/mp.sqrt(1-q/4) for q in [0.1,1,2,3,3.79] for c in [mp.sqrt(q)*k/50 for k in range(1,51)]))

# ---- companion paper: final-line arithmetic and misc ----
f = Fr(-613302797399911, 6480); v = Fr(179721388988719, 6480)
report("(-613302797399911+179721388988719)/6480 = -2007321335237/30 < 0", f+v == Fr(-2007321335237,30) and f+v < 0, str(f+v))
report("|O| = 30*C(4,2) = 180 (03-complete-family.tex:320)", 30*sp.binomial(4,2) == 180)
report("36*2*10 = 720 symmetries (Lemma cubic-symmetry)", 36*2*10 == 720)
print("ALL OK" if ok_all else "SOME CHECKS FAILED")
