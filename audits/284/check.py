import math, random
from fractions import Fraction as Fr
# 1. Grover average bound: avg >= 1/2 - 1/(2M sin2θ) >= 9/20 for gamma<=p<=1/2, M=ceil(10/sqrt(gamma))
worst=1
for gamma in [1/(2*m) for m in [2,4,16,64,1024,4096]]:
    M=math.ceil(10/math.sqrt(gamma))
    for p in [gamma*(1+i/50) for i in range(50)]+[0.5,0.25,0.1]:
        if p<gamma or p>0.5: continue
        th=math.asin(math.sqrt(p))
        avg=sum(math.sin((2*l+1)*th)**2 for l in range(M))/M
        worst=min(worst,avg)
        assert avg>=0.5-1/(2*M*math.sin(2*th))-1e-12
print("grover avg min",worst, ">=0.45?",worst>=0.45)
# 2. drift identity with random complementary a
for _ in range(1000):
    n=random.randint(2,7); w=[random.random() for _ in range(n)]
    a=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            x=random.random(); a[i][j]=x; a[j][i]=1-x
    W=sum(w); E=sum(w[p]/W*sum(w[i]*a[i][p] for i in range(n)) for p in range(n))
    assert abs(E-(W*W-sum(x*x for x in w))/(2*W))<1e-12
print("drift identity ok")
# 3. quantum error budget
for m in [4**j for j in range(1,30)]:
    T=math.ceil(math.log2(1000*m)); eps=Fr(1,1000*T)
    tot=Fr(2*m,2**T)+T*eps+eps
    assert tot<Fr(1,3), m
    assert (1-eps)**T > Fr(1,2) if T<40 else True
print("quantum budget ok; e.g. m=4^10:", float(Fr(2*4**10,2**math.ceil(math.log2(1000*4**10)))+Fr(1,1000)+Fr(1,1000*math.ceil(math.log2(1000*4**10)))))
# 4. classical exponents
print("20*log2(3/2)=",20*math.log2(1.5), " => G(2/3)^r <= m^", 2-20*math.log2(1.5)/1)
# exponent of separation: R ~ m^{2-1/k}, Q ~ m^{1/2} => R ~ Q^{4-2/k}
for k in range(2,8): print(k, Fr(2*k-1,k)/Fr(1,2))
# 5. circuit comparisons: Phi = N/M, M=m*sqrt(m)^(k-1)
for k in [2,3,4]:
  for m in [4,16,64]:
    Mk=m*int(math.isqrt(m))**(k-1); delta=Fr(1,2**(5*k))
    for N in range(-Mk,Mk+1,max(1,Mk//997)):
        Phi=Fr(N,Mk)
        assert (Phi>=delta)==(2**(5*k)*N>=Mk)
        assert (abs(Phi)<=delta/2)==(2**(5*k+1)*abs(N)<=Mk)
print("threshold comparisons ok")
# 6. Forrelation normalization: u^T D1 H D2 ... H Dk u vs BS (1/N) z1^T H Z2 ... H zk
def had(n):
    H=[[1]]
    while len(H)<n: H=[r+r for r in H]+[r+[-x for x in r] for r in H]
    return H
import itertools
n=8;H=had(n);s=1/math.sqrt(n)
for _ in range(20):
    k=3;z=[[random.choice([1,-1]) for _ in range(n)] for _ in range(k)]
    v=[s*z[k-1][i] for i in range(n)]
    for h in range(k-2,-1,-1):
        v=[s*sum(H[i][j]*v[j] for j in range(n)) for i in range(n)]
        v=[z[h][i]*v[i] for i in range(n)]
    phi=sum(s*x for x in v)
    # BS: (1/N^{(k+1)/2}) sum over ... same as (1/N) z1^T Hn Z2 Hn z3 with Hn normalized
    w=list(z[2]); w=[s*sum(H[i][j]*w[j] for j in range(n)) for i in range(n)]; w=[z[1][i]*w[i] for i in range(n)]
    w=[s*sum(H[i][j]*w[j] for j in range(n)) for i in range(n)]
    bs=sum(z[0][i]*w[i] for i in range(n))/n
    assert abs(phi-bs)<1e-12
print("normalization matches BS21")
# 7. BS bound vs L: (m/log(km))^{1-1/k} / (log2 m)^{-2}... ratio grows
for k in [2,5]:
  for e in [10,40,100]:
    m=4**e; bs=(m/math.log2(k*m))**(1-1/k); L=m**(1-1/k)/math.log2(m)**2
    print(k,e,"BS/L ratio",bs/L)
