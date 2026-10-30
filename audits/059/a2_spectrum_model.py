# Independent audit script. End-to-end test of the spectral comparison in Prop. 3.3
# on a small model (d = 3^k small) with ARBITRARY integer e_b (no product relation),
# using actual chain polynomials built by Lemma 3.2, their exact spectra
# (Steenbrink formula prod (t^w - t)/(1 - t^w), divided exactly), reduction mod 2,
# and exact character evaluation in F_p.  Prediction from the paper's Cases 1-3:
#   odd k  : chi_k(P2) == chi_k(P1)
#   k = 2m : chi_k(P2) == chi_k(P1) * [prod A_b^{2e_b} if d|m] * [prod (R_b/A_b)^{2e_b} if 3d|m]
# Then the real data: per-row cotangent identity in F_p for d = 3^11.
import random, sys
from fractions import Fraction

def is_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    dd, s = n - 1, 0
    while dd % 2 == 0: dd //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, dd, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def find_prime(mod, start):
    p = (start // mod) * mod + 1
    while not is_prime(p): p += mod
    return p

def prim_root_of_order(p, n):
    # element of exact order n in F_p^*, n | p-1
    fac = []
    x = n; q = 2
    while q*q <= x:
        if x % q == 0:
            fac.append(q)
            while x % q == 0: x //= q
        q += 1
    if x > 1: fac.append(x)
    for g in range(2, 10**6):
        z = pow(g, (p-1)//n, p)
        if all(pow(z, n//q, p) != 1 for q in fac):
            return z
    raise RuntimeError

# ---- polynomials in s = t^{1/(2D)} as dict-free lists ----
def pmul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: r[i+j] += x*y
    return r

def pdiv_exact(a, b):
    a = a[:]; q = [0]*(len(a)-len(b)+1)
    lb = b[-1]
    for i in range(len(a)-len(b), -1, -1):
        c = a[i+len(b)-1]
        if c:
            assert c % lb == 0
            c //= lb; q[i] = c
            for j, y in enumerate(b): a[i+j] -= c*y
    assert not any(a), "non-exact division"
    return q

def chain_weights(D, c):
    ws = []
    while c:
        a = D // c; r = D - a*c
        ws.append(c); c = r
    return ws  # numerators over D

def spectrum_poly(D, nums):
    # nums: weight numerators over D  -> polynomial in s with s^{2D}=t
    # factor (t^w - t)/(1 - t^w) = (s^{2c} - s^{2D})/(1 - s^{2c})
    num = [1]; den = [1]
    for c in nums:
        f = [0]*(2*D+1); f[2*c] += 1; f[2*D] -= 1
        g = [0]*(2*c+1); g[0] = 1; g[2*c] = -1
        num = pmul(num, f); den = pmul(den, g)
    q = pdiv_exact(num, den)
    assert all(x >= 0 for x in q)
    return q

def reduce_mod(q, n):
    r = [0]*n
    for i, x in enumerate(q):
        if x: r[i % n] += x
    return r

def conv_cyc(a, b, n):
    r = [0]*n
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: r[(i+j) % n] += x*y
    return r

def model(k3, seed):
    d = 3**k3; D = 3*d; n = 4*D
    rnd = random.Random(seed)
    B = [b for b in range(1, (d+1)//2) if b % 3]
    es = {b: rnd.choice([-2,-1,1,2]) for b in B}
    E = sum(es.values())
    if E <= 0:
        es[B[0]] += 1 - E; E = sum(es.values())
    sides = [[1] + [0]*(n-1), [1] + [0]*(n-1)]  # group ring elements (index = 2D * alpha mod 4D)
    nvars = [0, 0]
    cache = {}
    def add(side, nums, times):
        key = tuple(nums)
        if key not in cache:
            cache[key] = reduce_mod(spectrum_poly(D, list(nums)), n)
        for _ in range(times):
            sides[side] = conv_cyc(sides[side], cache[key], n)
        nvars[side] += times*len(nums)
    for b, e in es.items():
        for c, cnt in ((3*b, 2*e), (b, -2*e), (d-b, 2*e), (d+b, 2*e)):
            ws = chain_weights(D, c)
            own = 1 if cnt > 0 else 0
            add(own, ws, abs(cnt))
            if len(ws) > 1:
                add(1-own, ws[1:], abs(cnt))
    # 4E squares on side 1 (index 0): spectrum t^{1/2} = s^{D}
    sq = [0]*n; sq[D] = 1
    for _ in range(4*E):
        sides[0] = conv_cyc(sides[0], sq, n); nvars[0] += 1
    assert nvars[0] == nvars[1], nvars
    # two cubes each side: t^{1/3}+t^{2/3}
    cube = [0]*n; cube[2*D//3] += 1; cube[4*D//3] += 1
    for j in (0,1):
        for _ in range(2):
            sides[j] = conv_cyc(sides[j], cube, n)
    # character check over two primes
    A = Fraction(1); RA = Fraction(1)
    for b, e in es.items():
        A *= Fraction(d-b, b)**(2*e)
        RA *= (Fraction((2*d+b)*(2*d-b), (3*d-b)*(d+b)) / Fraction(d-b, b))**(2*e)
    for p in (find_prime(n, 10**12), find_prime(n, 3*10**12)):
        z = prim_root_of_order(p, n)
        pw = [pow(z, i, p) for i in range(n)]
        def chi(vec, k):
            return sum(x*pw[(i*k) % n] for i, x in enumerate(vec) if x) % p
        def fr(q): return q.numerator % p * pow(q.denominator % p, p-2, p) % p
        bad = 0
        for k in range(n):
            v1, v2 = chi(sides[0], k), chi(sides[1], k)
            if k % 2:
                pred = 1
            else:
                m = k//2
                pred = 1
                if m % d == 0: pred = pred * fr(A) % p
                if m % (3*d) == 0: pred = pred * fr(RA) % p
            if v2 != pred * v1 % p:
                bad += 1
        print(f"  d=3^{k3} seed={seed} E={E} vars/side={nvars[0]+2} p={p}: mismatching characters = {bad} / {n}")
        assert bad == 0
    # sanity: the comparison is non-trivial (sides differ in Z[C_n] when A,RA != 1)
    print("    sides equal as mod-2 multisets?", sides[0] == sides[1], " prodA=1?", A == 1, " prod(R/A)=1?", RA == 1)
    print("    Milnor numbers mu1, mu2 =", sum(sides[0]), sum(sides[1]), " ratio", Fraction(sum(sides[1]), sum(sides[0])), "pred", A*RA)

print("[small-model end-to-end test of Cases 1-3]")
for k3, seed in ((2, 1), (2, 2), (3, 3)):
    model(k3, seed)

# ---- real data: per-row cot identity in F_p, d = 3^11 ----
import json
data = json.load(open(sys.argv[1]))
d = int(data["d"]); D = 3*d; n = 4*D
bs = [int(r["b"]) for r in data["relation"]]
p = find_prime(n, 10**13)
z = prim_root_of_order(p, n)
def V(c, k):  # value of (t^w - t)/(1-t^w) at odd character k, w=c/D; t^w = z^{2kc}, t=-1
    x = pow(z, (2*k*c) % n, p)
    assert x != 1
    return (x + 1) * pow((1 - x) % p, p-2, p) % p
rnd = random.Random(7)
cnt = 0
for b in bs:
    ks = [rnd.randrange(n) | 1 for _ in range(300)]
    for k in ks:
        if k % 3 == 0: continue
        r = V(3*b, k) * V(d-b, k) % p * V(d+b, k) % p * pow(V(b, k), p-2, p) % p
        assert r == p - 1, (b, k)   # = -1, so its 2e_b-th power is 1
        cnt += 1
print(f"[real data] per-row identity V(3b)V(d-b)V(d+b)/V(b) = -1 checked at {cnt} (b, odd k, 3∤k) pairs, all 849 b")
# exhaustive over all odd k for a few b
for b in bs[:3]:
    for k in range(1, n, 2):
        if k % 3 == 0: continue
        r = V(3*b, k) * V(d-b, k) % p * V(d+b, k) % p * pow(V(b, k), p-2, p) % p
        assert r == p - 1
print("[real data] exhaustive over all odd k (3∤k) for b =", bs[:3], ": OK")
# cube vanishing at odd 3|k
for k in range(3, n, 6):
    pass
zz = [pow(z, (k*2*D//3) % n, p) + pow(z, (k*4*D//3) % n, p) for k in range(3, n, 6)]
assert all(x % p == 0 for x in zz)
print("[real data] cube spectrum vanishes at all odd k with 3|k: OK")
