"""Random pointwise tests of Lemma 3.4 (3.21)/(3.22), block Cauchy-Schwarz (3.23), weighted square (3.25)."""
import mpmath as mp, random
mp.mp.dps=30
random.seed(7)
al=mp.mpf(3)/40; s=1/(1-al); eps=mp.mpf(3)/20; kap=mp.mpf(49)/100
def rnd(m,sym=False,scale=1):
    A=mp.matrix(m,m)
    for i in range(m):
        for j in range(m): A[i,j]=mp.mpf(random.gauss(0,scale))
    return (A+A.T)/2 if sym else A
tr=lambda A: sum(A[i,i] for i in range(A.rows))
def fro2(A): return sum(A[i,j]**2 for i in range(A.rows) for j in range(A.cols))
def mfun(A,f):
    E,Q=mp.eigsy(A); return Q*mp.diag([f(e) for e in E])*Q.T
I=lambda m: mp.eye(m)
w34=w323=w325=-1; wpt=0
for trial in range(1500):
    m=random.choice([2,3,4]); sc=random.choice([.2,1,3])
    X=[rnd(m,True,sc) for _ in range(m)]
    B=sum((Xr*Xr for Xr in X), mp.zeros(m,m))
    S=mp.inverse(I(m)/s+al*B); T=S-I(m)
    # Lemma 3.4
    w=[rnd(m) for _ in range(m)]
    b0=rnd(m,True)
    Fm=sum((X[r]*w[r] for r in range(m)), mp.zeros(m,m))
    PT=lambda b: tr(T*b*B*b - T*B*b*b)
    comm=lambda b: T*b-b*T
    # check (3.22)
    S12=mfun(S,lambda x: 1/mp.sqrt(x)); C=S12*comm(b0)*S12
    wpt=max(wpt, abs(PT(b0)-fro2(C)/(2*al)))
    P=PT(b0); Lin=tr(comm(b0)*Fm)
    c0=-tr(T*sum((w[r].T*w[r] for r in range(m)), mp.zeros(m,m)))
    best=c0+(Lin**2/P if P>0 else 0)   # max over t of -t^2 P - 2 t Lin + c0
    rhs=(1+s/8)*sum(fro2(wr) for wr in w)
    w34=max(w34,best/rhs)
    # (3.23)
    a=rnd(m); L=[rnd(m) for _ in range(m)]
    G=a+sum((X[r]*L[r] for r in range(m)), mp.zeros(m,m))
    lhs=tr(G.T*S*G); rhs=s*fro2(a)+sum(fro2(Lr) for Lr in L)/al
    w323=max(w323,lhs/rhs)
    # (3.25): p,h0 functions of a random symmetric A
    A=rnd(m,True,2); E,Q=mp.eigsy(A)
    pv=[mp.mpf(random.uniform(.01,.5)) for _ in range(m)]; hv=[mp.mpf(random.gauss(0,.5)) for _ in range(m)]
    p=Q*mp.diag(pv)*Q.T; h0=Q*mp.diag(hv)*Q.T; pinvh=Q*mp.diag([1/mp.sqrt(x) for x in pv])*Q.T
    Y=rnd(m,True); M=[rnd(m) for _ in range(m)]
    H=al*S*(B-I(m))*(B-I(m)); rho=eps*I(m)+kap*S/s
    Fv=Y*h0+sum((X[r]*M[r] for r in range(m)), mp.zeros(m,m))
    # optimise over scaling of (Y,M): LHS linear t*tr(TF), RHS = tr(H(I-rho)p) + t^2 * quad -> worst t
    lin=tr(T*Fv); quad=sum(fro2(Mr*pinvh) for Mr in M)/(4*(1-eps))+al*s*fro2(Y*h0*pinvh)/(4*(1-eps-kap))
    c=tr(H*(I(m)-rho)*p)
    # need t*lin <= c + t^2 quad for all t  <=> lin^2 <= 4 c quad
    ratio = lin**2/(4*c*quad) if c>0 else (0 if abs(lin)<1e-20 else mp.inf)
    w325=max(w325,ratio)
print('(3.22) identity max abs error:',mp.nstr(wpt,5))
print('Lemma 3.4: max over trials of [max_b LHS]/RHS (<=1 required):',mp.nstr(w34,10))
print('(3.23): max LHS/RHS (<=1 required):',mp.nstr(w323,10))
print('(3.25): max lin^2/(4 c quad) (<=1 required):',mp.nstr(w325,10))
