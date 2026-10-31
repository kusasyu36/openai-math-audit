"""
Audit 138: independent recomputation of Proposition "Filtering moments"
(filtering.tex:90-98) and of the certificate tables in numerical.tex.

We recompute everything from the definitions in filtering.tex:18-40
(not from the appendix identities), using mpmath at 40 digits and
periodic trapezoid rules with 64 / 512 / 4096 nodes.
"""
import mpmath as mp
mp.mp.dps = 40

# --- constants, filtering.tex:26, 61-66, 76-80 -------------------------
kappa = mp.mpf('1.0315')
eta, tau = mp.mpf('0.51772'), mp.mpf('0.48996')
Delta, Gamma = eta - tau, 2*tau - eta
a = Delta/Gamma
r_c = [mp.mpf(x) for x in ('0.851','0.466','0.325','0.253','0.224','-0.314')]
s_c = [mp.mpf(x) for x in ('1.281','-0.938','0.731','-0.559','0.391','0.110')]
LOG2 = mp.log(2)

def r(x): return sum(c*mp.cos((j+1)*x) for j,c in enumerate(r_c))
def s(x): return sum(c*mp.cos((j+1)*x) for j,c in enumerate(s_c))
def F(x): return mp.power(2, r(x))
def G(x): return mp.power(2, s(x))
def M2(x): return 2+2*mp.cos(x)
def log2(y): return mp.log(y)/LOG2

def nodes(N): return [2*mp.pi*j/N for j in range(N)]
def avg(vals): return mp.fsum(vals)/len(vals)

def scores(N):
    """finite-circle (N-point) version of filtering.tex:18-40, exactly as
    Lemma finite-filtering says: every integral -> N-point average."""
    xs = nodes(N)
    Fv = [F(x) for x in xs]; Gv = [G(x) for x in xs]
    Zf, Zg = avg(Fv), avg(Gv)
    f = [v/Zf for v in Fv]; g = [v/Zg for v in Gv]
    fk = [mp.power(v,kappa) for v in f]; gk = [mp.power(v,kappa) for v in g]
    Zfk, Zgk = avg(fk), avg(gk)
    fs = [v/Zfk for v in fk]; gs = [v/Zgk for v in gk]
    lfs = [log2(v) for v in fs]; lgs = [log2(v) for v in gs]
    out = []
    for jx in range(N):
        prod = [f[(jx+i)%N]*g[i] for i in range(N)]
        Q = avg(prod)
        K = [p/Q for p in prod]
        D = avg([K[i]*log2(K[i]) for i in range(N)])
        A = avg([K[i]*lfs[(jx+i)%N] for i in range(N)])
        E = avg([K[i]*lgs[i] for i in range(N)])
        u = A+E-D; v = 2*D-A-E
        out.append((xs[jx], Q, D, A, E, u, v))
    return out, (Zf, Zg, Zfk*Zf**kappa, Zgk*Zg**kappa)

def moments(sc):
    X1 = [M2(x)*mp.power(2, mp.mpf('8.5')*(Delta-u)) for (x,Q,D,A,E,u,v) in sc]
    X2 = [M2(x)*mp.power(2, mp.mpf('8.3')*(a*v-u)) for (x,Q,D,A,E,u,v) in sc]
    return avg(X1), avg(X2), X1, X2

print("Delta =", Delta, " Gamma =", Gamma, " a =", a)
print("target 2^0.5168 =", mp.power(2, mp.mpf('0.5168')), " 2^0.5169 =", mp.power(2, mp.mpf('0.5169')))
print("bounds |r|<=2.433 ?", max(abs(r(x)) for x in nodes(4096)) , " |s|<=4.010 ?", max(abs(s(x)) for x in nodes(4096)))

res = {}
for N in (64, 512, 4096):
    sc, Z = scores(N)
    m1, m2, X1, X2 = moments(sc)
    res[N] = (m1, m2, sc, Z, X1, X2)
    print(f"N={N:5d}: int M^2 2^(8.5(Delta-u)) = {mp.nstr(m1,12)}   int M^2 2^(8.3(a v-u)) = {mp.nstr(m2,12)}")
    print(f"          Z_f={mp.nstr(Z[0],9)} Z_g={mp.nstr(Z[1],9)} Z_f,kappa(unnormalised F^kappa)={mp.nstr(Z[2],9)} Z_g,kappa={mp.nstr(Z[3],9)}")

# entropy order (filtering.tex:44-46) and the Jensen bound 2^u <= Q_*(x)
sc = res[4096][2]
print("min D, min D-A, min D-E over 4096 nodes:",
      mp.nstr(min(t[2] for t in sc),6), mp.nstr(min(t[2]-t[3] for t in sc),6), mp.nstr(min(t[2]-t[4] for t in sc),6))
print("mean of 2^u (must be <= 1):", mp.nstr(avg([mp.power(2,t[5]) for t in sc]),8))
print("max u, min u, max v:", mp.nstr(max(t[5] for t in sc),6), mp.nstr(min(t[5] for t in sc),6), mp.nstr(max(t[6] for t in sc),6))

# compare with Table 1 (numerical.tex:96-128): hat X_j, hat X'_j at 64 nodes
tab = """0 1.71479 1.72487
1 1.70072 1.70914
2 1.66451 1.66858
3 1.62068 1.61932
4 1.58491 1.57884
5 1.56813 1.55943
6 1.57405 1.56556
7 1.59915 1.59372
8 1.63351 1.63323
9 1.66273 1.66812
10 1.67219 1.68188
11 1.65371 1.66486
12 1.61070 1.62013
13 1.55729 1.56268
14 1.51133 1.51181
15 1.48634 1.48247
16 1.48675 1.48020
17 1.50690 1.50018
18 1.53309 1.52883
19 1.54784 1.54794
20 1.53630 1.54099
21 1.49268 1.50055
22 1.42405 1.43289
23 1.34881 1.35678
24 1.29024 1.29629
25 1.26827 1.27187
26 1.29131 1.29177
27 1.34397 1.33998
28 1.36465 1.35461
29 1.23032 1.21474
30 0.82842 0.81341
31 0.27604 0.26987
32 0 0"""
X1_64, X2_64 = res[64][4], res[64][5]
worst = mp.mpf(0)
for line in tab.splitlines():
    j, x1, x2 = line.split(); j = int(j)
    d1 = abs(X1_64[j]-mp.mpf(x1)); d2 = abs(X2_64[j]-mp.mpf(x2))
    worst = max(worst, d1, d2)
print("max |my 64-node hat X_j - table| (table tolerance 0.0004):", mp.nstr(worst, 5))
print("64-node averages (paper eq hatted-average-bound: <1.4295):", mp.nstr(res[64][0],8), mp.nstr(res[64][1],8))

# Check the 1.4295 average via the symmetric weighting the paper describes
w = lambda j: 1 if j in (0,32) else 2
avg1 = sum(w(j)*mp.mpf(l.split()[1]) for j,l in enumerate(tab.splitlines()))/64
avg2 = sum(w(j)*mp.mpf(l.split()[2]) for j,l in enumerate(tab.splitlines()))/64
print("table-based 64-averages:", mp.nstr(avg1,7), mp.nstr(avg2,7), " +0.0004 tol:", mp.nstr(avg1+mp.mpf('0.0004'),7), mp.nstr(avg2+mp.mpf('0.0004'),7))

# Final slack: the exponent used in checksum.tex:488-495
eps = mp.mpf('1e-8')
print("exponent -tau-Delta+0.5169+48eps =", mp.nstr(-tau-Delta+mp.mpf('0.5169')+48*eps, 10), " (paper: -0.00081952)")
print("log2 of computed moment1, moment2 (4096 nodes):", mp.nstr(log2(res[4096][0]),8), mp.nstr(log2(res[4096][1]),8), " need < eta =", eta)

# "Ideal" lower bound: min over u with E 2^u <= 1 of E[M^2 2^{lam(Delta-u)}] = 2^{lam Delta} (E (M^2)^{1/(lam+1)})^{lam+1}
for lam in (mp.mpf('8.5'), mp.mpf(4), mp.mpf(20)):
    p = 1/(lam+1)
    I = mp.quad(lambda x: mp.power(M2(x), p), [0, mp.pi])/mp.pi
    print(f"ideal bound at lambda={lam}: 2^(lam*Delta)*(E (M^2)^(1/(lam+1)))^(lam+1) = {mp.nstr(mp.power(2,lam*Delta)*mp.power(I, lam+1), 8)}")

# Analytic-strip sanity (numerical.tex:290-297): |Q(z)-1| on |Im z|<=0.38 (continuous, 512 nodes)
sc512, Z512 = res[512][2], res[512][3]
Zf, Zg = Z512[0], Z512[1]
def Qz(z, N=512):
    xs = nodes(N)
    return avg([mp.power(2, r(z+w))/Zf * mp.power(2, s(w))/Zg for w in xs])
worstQ = max(abs(Qz(mp.mpc(x, y))-1) for x in [2*mp.pi*j/32 for j in range(32)] for y in (mp.mpf('0.38'), mp.mpf('-0.38'), mp.mpf('0.19')))
print("max |Q(z)-1| on sampled strip points |Im z|=0.38,0.19 (paper: <0.54):", mp.nstr(worstQ, 6))
