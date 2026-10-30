"""Numerical spot checks of explicit claims in the pointwise 2-converse paper (audit 004b)."""
import random, itertools, math
from fractions import Fraction

def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def count_points(a1, a2, a3, a4, a6, p):
    # y^2 + a1xy + a3y = x^3 + a2x^2 + a4x + a6 over F_p, p odd
    n = 1
    for x in range(p):
        # solve y^2 + (a1x+a3)y - rhs = 0 ; disc = (a1x+a3)^2 + 4 rhs
        rhs = (x**3 + a2*x*x + a4*x + a6) % p
        d = ((a1*x + a3)**2 + 4*rhs) % p
        if d == 0: n += 1
        elif pow(d, (p-1)//2, p) == 1: n += 2
    return n

def leg(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

PR = [p for p in primes(1500) if p > 3]

# (a) pw:good-simple : E: y^2 = x(x^2+ax+b), F=Q(sqrt(a^2-4b)), F'=Q(sqrt(b)) [E': y^2=x(x^2-2ax+a^2-4b)]
def check_a():
    bad = 0; tested = 0
    for (a, b) in [(1, -1), (3, 1), (2, 3), (-1, 5), (5, -3), (7, 11), (1, 2)]:
        D = a*a - 4*b
        for p in PR:
            if (b * D) % p == 0: continue
            if leg(D, p) == -1 and leg(b, p) == -1:
                n = count_points(0, a, 0, b, 0, p)
                tested += 1
                if n % 4 != 2: bad += 1
            # also: F' generator check: E' 2-torsion field
    print(f"(a) #E(F_p)=2 mod 4 at p inert in F and F': tested {tested}, failures {bad}")

# (b) co:cp recurrence; split p (full 2-torsion over F_p) => a_p = p+1 mod 4, Delta_s/2 constant mod 2
def check_b():
    bad = 0; tested = 0
    for es in [(0, 1, -1), (0, 5, -5), (0, 3, -6), (1, 4, -9), (0, 7, 16)]:
        for p in PR[:150]:
            e1, e2, e3 = es
            if ((e1-e2)*(e1-e3)*(e2-e3)) % p == 0: continue
            # y^2=(x-e1)(x-e2)(x-e3)
            a2 = -(e1+e2+e3); a4 = e1*e2+e1*e3+e2*e3; a6 = -e1*e2*e3
            ap = p + 1 - count_points(0, a2, 0, a4, a6, p)
            assert (ap - (p+1)) % 4 == 0
            for kappa in (1, -1):
                b = [1, ap - kappa]
                for s in range(1, 8): b.append(ap*b[-1] - p*b[-2])
                Ds = [b[i+1]-b[i] for i in range(8)]
                target = ((ap - kappa - 1)//2) % 2
                tested += 1
                if any(d % 2 for d in Ds) or any((d//2) % 2 != target for d in Ds): bad += 1
    print(f"(b) co:cp residues Delta_s/2 = (a_p-kappa-1)/2 mod 2 for all s<=7: tested {tested}, failures {bad}")

# (c) co:symbol-interface final: X = a + kappa(q-eps) - q eps^i (1+eps); X/2 = #E(F_q)/2 mod 2
def check_c():
    bad = 0; tested = 0
    for q in PR[:200]:
        eps = 1 if q % 4 == 1 else -1
        for a in range(-2*int(q**.5)-2, 2*int(q**.5)+3):
            if a % 2: continue
            for kappa in (1, -1):
                for i in (0, 1):
                    X = a + kappa*(q-eps) - q*(eps**i)*(1+eps)
                    tested += 1
                    if X % 2 or ((X//2) - (q+1-a)//2) % 2: bad += 1
    print(f"(c) squared-trace residue = #E(F_q)/2 mod 2: tested {tested}, failures {bad}")

# (d) cy:binary : digits of z(x) are Boolean polys of degree <= 2^j ; and Chevalley-Warning existence
def anf_degree(vals, b):
    f = list(vals)
    for i in range(b):
        for m in range(1 << b):
            if m >> i & 1: f[m] ^= f[m ^ (1 << i)]
    return max((bin(m).count('1') for m in range(1 << b) if f[m]), default=-1)

def check_d():
    rng = random.Random(7)
    worst = {}
    for trial in range(200):
        b = rng.randint(2, 7); MOD = 2**8
        coeffs = [rng.randrange(MOD) for _ in range(1 << b)]
        z = []
        for x in range(1 << b):
            z.append(sum(c * (-1)**bin(x & g).count('1') for g, c in enumerate(coeffs)) % MOD)
        for j in range(6):
            deg = anf_degree([(v >> j) & 1 for v in z], b)
            worst[j] = max(worst.get(j, -1), deg)
            assert deg <= 2**j, (j, deg)
    # existence: r tests precision M with b > r(2^M-1)
    fails = 0
    for trial in range(200):
        r = rng.randint(1, 2); M = rng.randint(1, 2); b = r*(2**M - 1) + 1
        if b > 10: continue
        MOD = 2**M
        zs = []
        for _ in range(r):
            coeffs = [rng.randrange(2**10) for _ in range(1 << b)]
            zs.append([sum(c*(-1)**bin(x & g).count('1') for g, c in enumerate(coeffs)) for x in range(1 << b)])
        ok = any(all((zz[x]-zz[0]) % MOD == 0 for zz in zs) for x in range(1, 1 << b))
        if not ok: fails += 1
    print(f"(d) cy:binary: max ANF degree of digit j over trials {worst} (claim <= 2^j); existence failures {fails}")

# (e) N(1+zeta+zeta^-1) = -1 and odd number of negative embeddings
def check_e():
    for n in range(1, 6):
        m = 2**(n+2)
        vals = [1 + 2*math.cos(2*math.pi*k/m) for k in range(1, m, 2) if k < m//2]
        prod = math.prod(vals); neg = sum(v < 0 for v in vals)
        print(f"(e) n={n}: degree {len(vals)}, norm {prod:.6f}, #negative embeddings {neg}")

# (f) rd:residual-euler symbolic mod 2 (Laurent polys in z with integer a)
def check_f():
    import itertools
    bad = 0
    for a in range(-6, 7):
        # represent Laurent poly as dict exponent->coef
        def mul(P, Q):
            R = {}
            for i, c in P.items():
                for j, d in Q.items(): R[i+j] = R.get(i+j, 0) + c*d
            return R
        for q in (3, 5, 7, 9, 11):
            # det(1 - q^{-1} z Fr) = 1 - (a/q) z + (1/q) z^2 ; multiply by q: q - a z + z^2 (q odd unit)
            L1 = {0: q, 1: -a, 2: 1}; L2 = {0: q, -1: -a, -2: 1}
            l = mul(L1, L2)
            e = {0: a, 1: -1, -1: -1}; e2 = mul(e, e)
            # l/q^2 vs e^2 : q^2 odd so compare l and e2 mod 2 (since q^2 = 1 mod 2... need q^2 l' ) check l - q^2 e2? mod 2 only parity
            keys = set(l) | set(e2)
            if any((l.get(k, 0) - e2.get(k, 0)) % 2 for k in keys): bad += 1
    print(f"(f) l_q = e_q^2 mod 2 identity: failures {bad}")

# (g) Tamagawa at twisted good prime: E^{(q)} with E: y^2=(x-e1)(x-e2)(x-e3), q odd good -> type I0*,
#     c_q = 1 + #{roots of cubic mod q} ; 2-adic valuation should equal t_q = dim E[2]^{Fr_q}
def check_g():
    bad = 0; tested = 0
    for (a, b) in [(1, -1), (3, 1), (2, 3), (-1, 5)]:   # one rational 2-torsion: x(x^2+ax+b)
        D = a*a - 4*b
        for q in PR[:100]:
            if (b*D) % q == 0: continue
            roots = 1 + (2 if leg(D, q) == 1 else 0)   # roots of x(x^2+ax+b) mod q
            c = 1 + roots
            t = 1 if roots == 1 else 2
            tested += 1
            if int(math.log2(c)) != t or c & (c-1): bad += 1
    print(f"(g) I0* Tamagawa c_q=1+#roots: v_2(c_q)=t_q: tested {tested}, failures {bad}")

# (h) rd:exact-frobenius and finite formula in Z[F]/(F^2 - aF + l)
def check_h():
    bad = 0
    for a in range(-20, 21):
        for l in range(3, 60, 2):
            # represent u + vF ; F^2 = aF - l
            def mul(x, y):
                u1, v1 = x; u2, v2 = y
                return (u1*u2 - l*v1*v2, u1*v2 + u2*v1 + a*v1*v2)
            lhs = mul((-a, 1), (-(l+1), a))
            if lhs != (a, -(l+1)): bad += 1
    print(f"(h) (F-a)(aF-(l+1)) = a-(l+1)F in Z[F]/(F^2-aF+l): failures {bad}")

if __name__ == "__main__":
    check_a(); check_b(); check_c(); check_d(); check_e(); check_f(); check_g(); check_h()
