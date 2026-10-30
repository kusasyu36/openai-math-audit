# Audit 029: independent checks (own code; no repository code executed)
from fractions import Fraction as Fr
import cmath, math, itertools
ok = 0; bad = 0
def chk(name, cond):
    global ok, bad
    print(("OK  " if cond else "FAIL"), name)
    if cond: ok += 1
    else: bad += 1

# ---------- 1. exponent bookkeeping (Sections 4, 7, 8) ----------
a, b = Fr(9,10), Fr(133,1000); M0 = a+b; h0 = 1-a
chk("M0=1033/1000", M0 == Fr(1033,1000))
refl = (2*M0-1-b)/2
chk("reflection exponent 933/2000", refl == Fr(933,2000))
Cstar1 = a/2 + 1 - 1 + h0/6
chk("C*(1)=a/2+h0/6=7/15", Cstar1 == Fr(7,15))
chk("C*(s)=s-8/15", a/2 - 1 + h0/6 == Fr(-8,15))
chk("margin 7/15-933/2000 = 1/6000", Cstar1 - refl == Fr(1,6000))
chk("condition b < 4(1-a)/3", b < Fr(4,3)*(1-a))
chk("1/6000 - 1e-5 > 1/20000", Fr(1,6000)-Fr(1,100000) > Fr(1,20000))
chk("1/20000 > 1e-6 (Thm 1.2 width weaker than proved)", Fr(1,20000) > Fr(1,10**6))
# X^(1/2-z) Z^(xi+z-1) Y^(w-1)
def E(sx, sw, sz): return a*(Fr(1,2)-sz) + sx + sz - 1 + b*(sw-1)
chk("E* formula equals a/2+sx-1+h0 sz+b(sw-1)", all(E(x,y,z)==a/2+x-1+h0*z+b*(y-1) for x,y,z in [(Fr(1),Fr(1),Fr(1,6)),(Fr(3),Fr(2),Fr(5))]))
d1 = Fr(1,10**7); s1 = 1+d1; z1 = Fr(1,6)+d1
Rw = E(s1, Fr(95,100), z1) - Cstar1
Rz = E(s1, Fr(1), Fr(16,100)) - Cstar1
print("   Rw rel =", float(Rw), " Rz rel =", float(Rz))
chk("Rw rel = -0.00664989", Rw == Fr(-664989,10**8))
chk("Rz rel ~ -0.000666566666", abs(float(Rz) + 0.000666566666) < 1e-11)
chk("Rw,Rz < -0.0005", Rw < Fr(-5,10000) and Rz < Fr(-5,10000))
base = E(s1, Fr(3,1000), Fr(4,10)) - Cstar1
print("   base rel =", float(base))
chk("base rel = 1e-7 + 7/300 - .132601", base == d1 + Fr(7,300) - Fr(132601,10**6))
r1 = base + Fr(9,100)*Fr(11,10); r2 = base + Fr(10001,100000)*Fr(6,10)
r3 = base - (s1-Fr(998,1000)) + Fr(10001,100000)*Fr(11,10)
print("   rows: small", float(r1), " exceptional", float(r2), " good", float(r3))
chk("small rows -0.010267566..", abs(float(r1)+0.010267566666)<1e-11)
chk("exceptional rows -0.049261566..", abs(float(r2)+0.049261566666)<1e-11)
chk("good rows -377/300000", r3 == Fr(-377,300000))
R0 = 350000
big = E(Fr(4),Fr(4),Fr(R0)) + Fr(10001,100000)*(1-R0)
chk("E*(4,4,R0)=3.849+.1R0", E(Fr(4),Fr(4),Fr(R0)) == Fr(3849,1000)+Fr(R0,10))
chk("large rows exponent .44901 < C*(1)-.01", big == Fr(44901,100000) and big < Cstar1-Fr(1,100))
# Euler-product margins (7.20),(7.21)
xi, w1, z1r = Fr(998,1000), Fr(3,1000), Fr(4,10)
w2, z2 = Fr(95,100), Fr(16,100)
chk("|D1W| exponent region1 = -1.001", -xi-w1 == Fr(-1001,1000))
chk("1-xi-w-6z max = -1.908", max(1-xi-w1-6*z1r, 1-xi-w2-6*z2) == Fr(-1908,1000))
chk("4-6xi-6z max = -2.948", max(4-6*xi-6*z1r, 4-6*xi-6*z2) == Fr(-2948,1000))
e = Fr(1,10**5)
chk("enlarged margins .00098,.90792,1.94788",
    (-(xi-e)-(w1-e)+1 == Fr(-98,100000)) and (1-(xi-e)-(w2-e)-6*(z2-e)+1 == Fr(-90792,100000)) and (4-6*(xi-e)-6*(z2-e)+1 == Fr(-194788,100000)))
chk("(7.21) l=1,e0=0 region1 = -0.497", (2-3*xi) + (Fr(1,2)-w1) == Fr(-497,1000))
chk("t=1,k=1 exponent 1-xi-w = -.001/-.948", 1-xi-w1 == Fr(-1,1000) and 1-xi-w2 == Fr(-948,1000))
# Prop 8.2 numerics
chk("(8.7) 10(1/2-.996)+1/2+1 = -3.46", 10*(Fr(1,2)-Fr(996,1000))+Fr(3,2) == Fr(-346,100))
chk("(8.9) 11*(2-2*.996) = .088", 11*(2-2*Fr(996,1000)) == Fr(88,1000))
th = math.log(1.002/.999)/math.log(1.003/.999)
chk("three-circles theta<1 (%.4f); disk edges .9965,.998,1.001" % th, th<1 and abs(2-1.0035-.9965)<1e-12 and abs(2-1.002-.998)<1e-12 and abs(2-.999-1.001)<1e-12)
# cross-prime cancellation in Prop 7.3 and reflection exponent
chk("tk'+t'k+(t-k)t'+(t'-k')t = 2tt'", all(t*k2+t2*k+(t-k)*t2+(t2-k2)*t == 2*t*t2 for t,k,t2,k2 in itertools.product(range(5),repeat=4)))
chk("2(e0+3l)(e0'+3l')+4e0e0' = 0 mod 6", all((2*(e0+3*l)*(f0+3*m)+4*e0*f0)%6==0 for e0,f0 in itertools.product((0,1),repeat=2) for l in range(4) for m in range(4)))
chk("row exponent 8-j at j=1 is 1; 1-4 = 3 mod 6 (quadratic)", (8-1)%6==1 and (1-4)%6==3)
# sieve depth in 12.3 (TeX: H = 2 floor(1/(40b)))
for bb in [Fr(1,200), Fr(1,300), Fr(1,1000), Fr(3,1000)]:
    H = 2*math.floor(1/(40*bb))
    chk("b=%s: b(4H+2)<=0.2+2b and e^-H<=e^2 e^{-1/(20b)}" % bb, bb*(4*H+2) <= Fr(2,10)+2*bb and -H <= 2 - 1/(20*float(bb)))
# 12.6 final margin
chk("0.9-0.2-1/4 > 1/4", Fr(9,10)-Fr(2,10)-Fr(1,4) > Fr(1,4))

# ---------- 2. Lemma 2.4 (progression of quadratic non-residues) ----------
def legendre(a_, p):
    a_ %= p
    if a_ == 0: return 0
    return 1 if pow(a_, (p-1)//2, p) == 1 else -1
def primes_upto(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
PR = primes_upto(200000)
table = {-1:3, 2:3, -2:7, 3:5, -3:5, 6:3, -6:7}
for d0,u8 in table.items():
    cand = [u for u in range(24) if u%8==u8 and u%3==2]
    u = cand[0]
    good = all(legendre(d0,p)==-1 for p in PR if p>3 and p%24==u)
    chk("table d0=%d u mod 8=%d (mod 24 class %d)" % (d0,u8,u), good)
def lemma24(aint):
    M = 8
    for l in PR:
        if l>2 and aint % l == 0: M *= l
    found = None
    for u in range(1, M):
        if math.gcd(u,M)!=1 or u%8 not in (3,5,7): continue
        if any(u%l==1 for l in PR if l>2 and M%l==0): continue
        c = 2 if u%8 in (3,7) else 4
        if (u-1)%c or math.gcd((u-1)//c, M//c)!=1: continue
        ps = [p for p in PR if p%M==u][:300]
        if ps and all(aint%p!=0 and legendre(aint,p)==-1 for p in ps):
            found=(M,c,u); break
    return found
allok = True
for aint in list(range(-60,61)):
    if aint in (-1,0,1) or (aint>0 and int(math.isqrt(aint))**2==aint): continue
    if lemma24(aint) is None: allok=False; print("   no (M,c,u) for a=",aint)
chk("Lemma 2.4 brute force for all admissible |a|<=60", allok)

# ---------- 3. Finite-field Gauss/Jacobi identities (Lemma 3.4) ----------
def find_gen(p):
    for g in range(2,p):
        if all(pow(g,(p-1)//q,p)!=1 for q in set(f for f in PR if (p-1)%f==0 and f<p)):
            return g
w = cmath.exp(2j*math.pi/3)
def ff_checks(p):
    g = find_gen(p); ind = {}
    x=1
    for k in range(p-1): ind[x]=k; x=x*g%p
    z6 = cmath.exp(2j*math.pi/6)
    chi = lambda x: 0 if x%p==0 else z6**(ind[x%p]%6)
    e = lambda x: cmath.exp(2j*math.pi*(x%p)/p)
    G = {j: sum((chi(x)**j if x%p else 0)*e(x) for x in range(p)) for j in range(1,6)}
    Q = p
    r1 = abs(G[2]*G[4]-Q)<1e-6
    rho = lambda x: chi(x)**2 if x%p else 0
    J = sum(rho(x)*rho(1-x) for x in range(p))
    r2 = abs(G[2]**2 - J*G[4])<1e-6
    r3 = abs((G[2]/math.sqrt(Q))**3 - J/math.sqrt(Q))<1e-6
    # J in Z[omega]: J = A + B w
    B = J.imag/ w.imag; A = J.real - B*w.real
    A, B = round(A), round(B)
    r4 = abs(A + B*w - J)<1e-6 and (A - (-1)) % 3 == 0 and B % 3 == 0 and A*A - A*B + B*B == Q
    # primary pi_+ = -J ; check -J ≡ 1 mod 3 and rho(x) ≡ x^{(p-1)/3} mod (pi) where pi=-J:
    # reduce: find integer m with w ≡ m mod pi (m^2+m+1 ≡ 0 mod p) and pi | (A + B m)
    ms = [m for m in range(p) if (m*m+m+1)%p==0]
    mm = [m for m in ms if (A+B*m)%p==0]
    r5 = False
    if mm:
        m = mm[0]
        # rho(x) = w^{k}; check w^k ≡ x^{(p-1)/3} mod p via w->m
        r5 = all(pow(m, (round((cmath.phase(rho(x))/(2*math.pi/3)))%3), p) == pow(x,(p-1)//3,p) for x in range(1,p))
    c4 = chi(4)
    r6 = abs(G[1]*G[4] - c4.conjugate()*G[2]*G[3])<1e-6
    gam = {j: G[j]/math.sqrt(Q) for j in G}
    r7 = abs(gam[1]*gam[2] - c4.conjugate()*gam[3]*gam[2]**3)<1e-6
    r8 = abs(gam[3]**2-1)<1e-9
    # digit identity (6.4): sum_v alpha_{j}(v) e(t v) = chi(t)^j  (t mod p)
    r9 = True
    for j in range(0,6):
        tau = sum((chi(v)**(-j) if True else 0)*e(v) for v in range(1,p)) if j else None
        def alpha(v):
            if j==0: return (1-1/Q) if v==0 else -1/Q
            return 0 if v==0 else chi(v)**(-j)/tau
        for t in range(p):
            lhs = sum(alpha(v)*e(t*v) for v in range(p))
            rhs = (chi(t)**j if t%p else 0) if j else (1 if t%p else 0)
            if abs(lhs-rhs)>1e-6: r9=False
    # Lemma 6.2 e=0, j=1: sum_v alpha_1(v) chi(v)^4 e(t/v) = chi(t)^{8-1-4}*omega ... check unit dependence chi(t)^{4-j}
    tau1 = sum(chi(v)**(-1)*e(v) for v in range(1,p))
    S = lambda t: sum(chi(v)**(-1)/tau1*chi(v)**4*e(t*pow(v,p-2,p)) for v in range(1,p))
    r10 = all(abs(S(t) - chi(t)**3*S(1))<1e-6 for t in range(1,p))  # chi(t)^{4-j}=chi^3, so total 4+3=7≡1
    return r1,r2,r3,r4,r5,r6,r7,r8,r9,r10
for p in [13,37,61,73,97,109]:
    res = ff_checks(p)
    chk("F_%d: G2G4=Q, G2^2=J G4, (G2/sqrtQ)^3=J/sqrtQ, J≡-1 mod 3 & N=Q, J=-pi+ w/ rho≡x^((p-1)/3) mod pi, G1G4=conj(chi4)G2G3, g1g2=conj(chi4)g3 g2^3, g3=±1, digit id (6.4), Lemma6.2 unit power" % p, all(res))
    if not all(res): print("   ", res)

# ---------- 4. Buchstab identity (11.6) via omega ----------
# D_g(w) = omega(w/g)/g ; check int_b^{1/2} D_t(1-t) dt/t = D_b(1)-1 for b=0.15
import functools
def omega_table(U, h=1e-4):
    n = int(U/h)+2; om=[0.0]*n; us=[1+i*h for i in range(n)]
    for i,u in enumerate(us):
        if u<=2: om[i]=1/u
        else:
            # u om(u) = 1 + int_1^{u-1} om
            k = int(round((u-1-1)/h))
            s = sum(om[1:k]) * h + (om[0]+om[k])*h/2
            om[i] = (1+s)/u
    return us, om
us, om = omega_table(8.0, 2e-3)
def omega(u):
    i = (u-1)/2e-3; i0=int(i); fr=i-i0
    return om[i0]*(1-fr)+om[i0+1]*fr
bq = 0.15
N=20000; tot=0
for k in range(N):
    t = bq + (0.5-bq)*(k+0.5)/N
    tot += omega((1-t)/t)/t * (1/t) * (0.5-bq)/N
lhs, rhs = tot, omega(1/bq)/bq - 1
print("   Buchstab: lhs=%.5f rhs=%.5f" % (lhs, rhs))
chk("(11.6) Buchstab identity numerically (b=0.15)", abs(lhs-rhs)<5e-3)
chk("omega(u)->e^-gamma (u=7.9)", abs(omega(7.9)-math.exp(-0.5772156649))<1e-4)
print("TOTAL OK=%d FAIL=%d" % (ok,bad))
