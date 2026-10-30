# Sec 13: dyadic prefix allocation (prop:patch-extraction-allocation) brute-force check
from fractions import Fraction as F
import random, itertools, math
random.seed(1)
def alloc(M):
    S=sum(M)
    r=[]
    for m in M:
        p=F(m,S); l=0
        while F(1,2**(l+1))>=p: l+=1   # largest dyadic 2^-l >= p
        assert p<=F(1,2**l)<2*p
        r.append((F(1,2**l),l,m))
    r.sort(key=lambda t:-t[0])
    J=[];s=F(0)
    for t in r:
        J.append(t); s+=t[0]
        if s>=1: break
    assert s==1, s
    assert sum(t[2] for t in J)*2>S
    return J,S
def prefix_code(J):
    # canonical complete prefix code from lengths sorted increasing
    code=[];c=0;prev=0
    for (_,l,m) in J:
        c<<= (l-prev); code.append((format(c,'0%db'%l) if l>0 else '',m)); c+=1; prev=l
    return code
bad=0
for trial in range(3000):
    k=random.randint(1,12)
    M=[random.randint(1,10**random.randint(1,6)) for _ in range(k)]
    J,S=alloc(M)
    code=prefix_code(J)
    n=max(len(w) for w,_ in code)+2
    words=[w for w,_ in code]
    # prefix-free and complete
    for a,b in itertools.permutations(words,2):
        assert not b.startswith(a)
    assert sum(F(1,2**len(w)) for w in words)==1
    # crossing check: flipping coord j<=l of vertex in leaf w lands in distinct other leaves
    for w in words:
        hit=[]
        for j in range(len(w)):
            x=w[:j]+('1' if w[j]=='0' else '0')+w[j+1:]+'0'*n
            leaf=[u for u in words if x.startswith(u)]
            assert len(leaf)==1 and leaf[0]!=w
            hit.append(leaf[0])
        assert len(set(hit))==len(hit)
print("allocation: all trials OK (sum r=1 exactly, sum_J M > S/2, complete prefix code, l_i crossings to distinct leaves)")
# column-sum factor |E_i|/M_i <= 2^n/S
for trial in range(1000):
    k=random.randint(1,8); M=[random.randint(1,1000) for _ in range(k)]
    J,S=alloc(M)
    for (r,l,m) in J:
        assert F(1,2**(l+1))/m < F(1,S)   # 2^{n-l-1}/M_i < 2^n/S  (divide by 2^n)
print("|E_i|/M_i < 2^n/S verified")
