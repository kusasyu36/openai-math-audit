"""Audit 138: sensitivity of the two filtering moments (filtering.tex:94-95) to the
coefficients r_j, s_j (filtering.tex:78-79), kappa and lambda.  64-node periodic
trapezoid rule (already 12-digit converged, see check_138_filtering.out)."""
import mpmath as mp
mp.mp.dps = 25
eta, tau = mp.mpf('0.51772'), mp.mpf('0.48996'); Delta = eta-tau; a = Delta/(2*tau-eta)
R0 = [mp.mpf(x) for x in ('0.851','0.466','0.325','0.253','0.224','-0.314')]
S0 = [mp.mpf(x) for x in ('1.281','-0.938','0.731','-0.559','0.391','0.110')]
TARGET = mp.power(2, mp.mpf('0.5168'))
LOG2 = mp.log(2)
def moments(rc, sc, kappa=mp.mpf('1.0315'), l1=mp.mpf('8.5'), l2=mp.mpf('8.3'), N=64):
    xs = [2*mp.pi*j/N for j in range(N)]
    r = lambda x: sum(c*mp.cos((j+1)*x) for j,c in enumerate(rc))
    s = lambda x: sum(c*mp.cos((j+1)*x) for j,c in enumerate(sc))
    Fv = [mp.power(2, r(x)) for x in xs]; Gv = [mp.power(2, s(x)) for x in xs]
    f = [v/(mp.fsum(Fv)/N) for v in Fv]; g = [v/(mp.fsum(Gv)/N) for v in Gv]
    fk = [v**kappa for v in f]; gk = [v**kappa for v in g]
    lfs = [mp.log(v/(mp.fsum(fk)/N))/LOG2 for v in fk]; lgs = [mp.log(v/(mp.fsum(gk)/N))/LOG2 for v in gk]
    m1 = m2 = mp.mpf(0)
    for jx in range(N):
        prod = [f[(jx+i)%N]*g[i] for i in range(N)]; Q = mp.fsum(prod)/N; K = [p/Q for p in prod]
        D = mp.fsum(K[i]*mp.log(K[i]) for i in range(N))/N/LOG2
        A = mp.fsum(K[i]*lfs[(jx+i)%N] for i in range(N))/N
        E = mp.fsum(K[i]*lgs[i] for i in range(N))/N
        u = A+E-D; v = 2*D-A-E; M2 = 2+2*mp.cos(xs[jx])
        m1 += M2*mp.power(2, l1*(Delta-u)); m2 += M2*mp.power(2, l2*(a*v-u))
    return m1/N, m2/N
def show(tag, m):
    print(f"{tag:42s} I1={mp.nstr(m[0],8)} I2={mp.nstr(m[1],8)}  {'both < 2^0.5168' if max(m) < TARGET else 'VIOLATED'}")
show("paper values", moments(R0, S0))
for d in (mp.mpf('0.001'), mp.mpf('-0.001'), mp.mpf('0.01'), mp.mpf('-0.01')):
    for j in (0, 5):
        rc = list(R0); rc[j] += d; show(f"r_{j+1} {'+' if d>0 else ''}{d}", moments(rc, S0))
        sc = list(S0); sc[j] += d; show(f"s_{j+1} {'+' if d>0 else ''}{d}", moments(R0, sc))
for k in ('1.0','1.02','1.0315','1.05','1.1'):
    show(f"kappa={k}", moments(R0, S0, kappa=mp.mpf(k)))
for l in ('6','8','8.5','9','10','12'):
    show(f"lambda1=lambda2={l}", moments(R0, S0, l1=mp.mpf(l), l2=mp.mpf(l)))
# the actual requirement is I < 2^eta; how much room would 2^eta allow?
print("2^0.5168 =", mp.nstr(TARGET,8), " 2^eta =", mp.nstr(mp.power(2,eta),8))
