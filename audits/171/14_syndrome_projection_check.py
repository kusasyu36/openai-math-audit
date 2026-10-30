# Sec 14: psi(z)=z xor e_{s(z)}, s(z)=sum l z_l in (Z/2)^m, h=2^m coordinates indexed by group elements
import itertools, collections
for m in (2,3,4):
    h=2**m
    def s(z):
        acc=0
        for l in range(h):
            if (z>>l)&1: acc^=l
        return acc
    psi=lambda z: z ^ (1<<s(z))
    fib=collections.Counter()
    for z in range(2**h):
        w=psi(z); assert s(w)==0; assert bin(z^w).count('1')==1; fib[w]+=1
    assert set(fib.values())=={h}, set(fib.values())
    # star multiplicities
    minj=10**9
    for v in range(2**h):
        if bin(v).count('1')%2: continue
        c=collections.Counter(psi(v^(1<<l)) for l in range(h))
        assert sum(c.values())==h and min(c.values())>=2
        minj=min(minj,min(c.values()))
    # site neighborhood distances and incidence count
    odd=[z for z in range(2**h) if bin(z).count('1')%2]
    groups=collections.defaultdict(list)
    for z in odd: groups[psi(z)].append(z)
    inc=collections.Counter(); maxd=0
    for w,G in groups.items():
        sites=set()
        for z in G:
            for l in range(h):
                u=z^(1<<l); sites.add(psi(u))
        for a in sites:
            for b in sites: maxd=max(maxd,bin(a^b).count('1'))
        for a in sites: inc[a]+=1
    print(f"h={h}: fibers size h OK, star multiplicities >=2 (min {minj}), max site distance {maxd}, max neighborhoods per site {max(inc.values())} (h^3={h**3})")
