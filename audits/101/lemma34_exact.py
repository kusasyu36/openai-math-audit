"""Lemma 3.4: exact maximisation over symmetric b (closed form in S-eigenbasis) and over w (eigenvalue)."""
import mpmath as mp, random
mp.mp.dps=25; random.seed(3)
al=mp.mpf(3)/40; s=1/(1-al)
def rnd(m,sc): 
    A=mp.matrix(m,m)
    for i in range(m):
        for j in range(m): A[i,j]=mp.mpf(random.gauss(0,sc))
    return (A+A.T)/2
worst=0
for trial in range(300):
    m=random.choice([2,3]); sc=random.choice([.1,.5,1,2,5])
    X=[rnd(m,sc) for _ in range(m)]
    B=sum((x*x for x in X), mp.zeros(m,m))
    Sv,U=mp.eigsy(mp.inverse(mp.eye(m)/s+al*B))
    Xs=[U.T*x*U for x in X]   # work in S-eigenbasis
    n=m*m*m
    def form(wv):
        w=[mp.matrix(m,m) for _ in range(m)]
        for r in range(m):
            for i in range(m):
                for j in range(m): w[r][i,j]=wv[r*m*m+i*m+j]
        F=sum((Xs[r]*w[r] for r in range(m)), mp.zeros(m,m))
        val=sum(al*Sv[i]*Sv[j]*(F[i,j]-F[j,i])**2 for i in range(m) for j in range(i+1,m))
        val-=sum((Sv[j]-1)*w[r][i,j]**2 for r in range(m) for i in range(m) for j in range(m))  # tr(T w^t w): T acts on column index
        return val
    Qm=mp.matrix(n,n)
    e=lambda k:[mp.mpf(1) if t==k else mp.mpf(0) for t in range(n)]
    diag=[form(e(k)) for k in range(n)]
    for k in range(n):
        Qm[k,k]=diag[k]
        for l in range(k+1,n):
            v=[a+b for a,b in zip(e(k),e(l))]
            Qm[k,l]=Qm[l,k]=(form(v)-diag[k]-diag[l])/2
    ev=max(mp.eigsy(Qm,eigvals_only=True))
    worst=max(worst,ev/(1+s/8))
print('Lemma 3.4 exact sup ratio (<=1 required):',mp.nstr(worst,10))
