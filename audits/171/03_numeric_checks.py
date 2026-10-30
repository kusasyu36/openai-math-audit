# Section 3 audit: brute-force check of Lemma scattered-moments on small instances,
# and Monte Carlo of the sequential sampler of Lemma near-product-injection (Steps 1-2).
from fractions import Fraction as F
import itertools, random, math

print("== scattered moments: exact brute force, independent Z_v (K=1), near set = {v} or cyclic windows")
random.seed(3)
worst=0
for trial in range(300):
    U=random.randint(3,5); n=random.randint(2,4); w=random.randint(0,1)  # near window radius
    near=[set((v+o)%U for o in range(-w,w+1)) for v in range(U)]
    f=F(max(len(s) for s in near),U)
    # each Z_v independent, values in {0,1,..,L}
    L=random.randint(1,3)
    dists=[]
    for v in range(U):
        wts=[random.randint(0,4) for _ in range(L+1)]
        if sum(wts)==0: wts[0]=1
        S=sum(wts); dists.append([F(x,S) for x in wts])
    dmean=[sum(k*pk for k,pk in enumerate(dv)) for dv in dists]
    # hypothesis: separated sequences -> distinct indices? (separated => outside near set => distinct). Independent => E prod = prod d. K=1
    # LHS exact
    lhs=F(0)
    for vals in itertools.product(*[range(L+1)]*U):
        pr=F(1)
        for v,x in enumerate(vals): pr*=dists[v][x]
        if pr==0: continue
        lhs+=pr*(F(sum(vals),U))**n
    rhs=(sum(dmean)/U + n*f*L)**n
    assert lhs<=rhs, (lhs,rhs)
    worst=max(worst,float(lhs/rhs))
print("all 300 instances satisfy the bound; max lhs/rhs =",worst)

print("== near-product Steps 1-2: sequential sampler tracking simulation")
def run(d,m,atom_exp,seed):
    rng=random.Random(seed)
    # real rows: each uniform on a random support of size s ~ d^{atom_exp}; column sums kept <= .4 roughly
    s=int(math.ceil(d**atom_exp))
    rows=[]; col=[0.0]*d
    for i in range(m):
        sup=rng.sample(range(d),s)
        rows.append({y:1.0/s for y in sup})
        for y in sup: col[y]+=1.0/s
    cmax=max(col)
    t=math.ceil(2*d/3)
    dummy={y:(t/d-col[y])/(t-m) for y in range(d)}
    assert min(dummy.values())>=0
    allrows=rows+[dummy]*(t-m)
    order=list(range(t)); rng.shuffle(order)
    allrows=[allrows[i] for i in order]
    used=set(); D={}  # track D for real rows only (sufficient for illustration)
    realidx=[i for i,o in enumerate(order) if o<m]
    maxerr=0.0
    Dvals={a:0.0 for a in realidx}
    for jj,row in enumerate(allrows):
        sj=jj/d
        for a in realidx: maxerr=max(maxerr,abs(Dvals[a]-sj))
        free=[(y,pv) for y,pv in row.items() if y not in used]
        tot=sum(pv for _,pv in free)
        r=rng.random()*tot; acc=0
        for y,pv in free:
            acc+=pv
            if acc>=r: break
        used.add(y)
        for a in realidx: Dvals[a]+=allrows[a].get(y,0.0)
    return cmax,maxerr
for d,m in [(2000,600),(4000,1200)]:
    for seed in range(2):
        cmax,me=run(d,m,0.95,seed)
        print(f"d={d} m={m} colmax={cmax:.3f} max_{{a,j}}|D_j^a-s_j|={me:.4f}  (d^-1/8={d**-0.125:.3f}, stop threshold 1/20)")
