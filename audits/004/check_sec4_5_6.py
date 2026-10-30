"""Independent checks for 004 (H10 over Q), Sections 4.1, 5, 6.1.
Run: .venv/bin/python -I check_sec4_5_6.py
"""
from sympy import sqrt, expand, simplify, factorint, isprime, primerange, Rational, nsimplify
from math import gcd

print("== Section 5: F=Q(sqrt2), d=75-53*sqrt2 ==")
s2 = sqrt(2)
d = 75 - 53*s2
ds = 75 + 53*s2
print("N(d) =", expand(d*ds))
print("d embeddings >0:", float(d) > 0, float(ds) > 0, float(d), float(ds))
w = -5 + 4*s2
print("w == d(1+sqrt2)^3 :", simplify(expand(d*(1+s2)**3) - w) == 0)
print("w^3 - w + 8d == 0 :", simplify(expand(w**3 - w + 8*d)) == 0)
Qx, Qy = -d*w, 2*s2*d**2
print("Q on E: y^2-(x^3-d^2 x) =", simplify(expand(Qy**2 - (Qx**3 - d**2*Qx))))
wsig = -5 - 4*s2
print("signs of -dw at two embeddings:", float(-d*w) > 0, float(-ds*wsig) > 0)
# -4 nonsquare mod 7
print("-4 mod 7 =", (-4) % 7, "squares mod 7:", sorted({x*x % 7 for x in range(1, 7)}))

print("\n== Section 6.1: reduction orders at split primes r=1 mod 8, r=3 mod 7 ==")
def count_points(r, c):
    # y^2 = x^3 - c^2 x over F_r, c in F_r
    sq = [0]*r
    for y in range(r):
        sq[y*y % r] += 1
    n = 1
    for x in range(r):
        n += sq[(x**3 - c*c*x) % r]
    return n

bad = 0; tested = 0
for r in primerange(3, 4000):
    if r % 8 != 1 or r % 7 != 3:
        continue
    # sqrt2 mod r
    b = next(x for x in range(r) if x*x % r == 2)
    for B in (b, r-b):
        dr = (75 - 53*B) % r
        dsr = (75 + 53*B) % r
        assert dr*dsr % r == 7
        d1, d2 = count_points(r, dr), count_points(r, dsr)
        tested += 1
        g = gcd(d1, d2)
        ok = (d1 + d2 == 2*(r+1)) and d1 % 4 == 0 and d2 % 4 == 0 and (g % 8 != 0) and (d1 % r) and (d2 % r)
        # odd common factor l must satisfy l|r+1 and l=1 mod 4
        oddg = g
        while oddg % 2 == 0:
            oddg //= 2
        for l in factorint(oddg):
            if (r+1) % l or l % 4 != 1:
                ok = False
        if not ok:
            bad += 1
            print("FAIL", r, B, d1, d2)
print(f"tested {tested} (prime,place) pairs; failures {bad}")
print("note: gcd 2-part is exactly 4 in all cases" if bad == 0 else "")

print("\n== Section 4.1: Atkin-Lehner fixed points on X^210 ==")
D = 210
Dp = [2, 3, 5, 7]
def kron(dd, p):
    # Kronecker symbol (dd/p) for prime p
    if p == 2:
        if dd % 2 == 0: return 0
        return 1 if dd % 8 in (1, 7) else -1
    r = dd % p
    if r == 0: return 0
    return 1 if pow(r, (p-1)//2, p) == 1 else -1

def classno(disc):
    # number of reduced primitive forms of discriminant disc<0
    h = 0
    a = 1
    while 3*a*a <= -disc:
        for b in range(-a+1, a+1):
            if (b*b - disc) % (4*a): continue
            c = (b*b - disc)//(4*a)
            if c < a: continue
            if c == a and b < 0: continue
            if gcd(gcd(a, abs(b)), c) != 1: continue
            h += 1
        a += 1
    return h

def fund_disc(m):  # disc of Q(sqrt(-m)), m squarefree
    return -m if (-m) % 4 == 1 else -4*m

divs = [m for m in range(2, D+1) if D % m == 0]
total = 0
for m in divs:
    fields = [fund_disc(m)] + ([-4] if m == 2 else [])
    nfix = 0
    for disc in fields:
        # field must be nonsplit at every p | D
        if any(kron(disc, p) == 1 for p in Dp):
            continue
        # orders: for fixed points of w_m: Z[sqrt(-m)] (disc -4m) and, if m=3 mod 4, maximal order disc -m;
        # for m=2 also Z[i] (disc -4).  Eichler: optimal embeddings into maximal order of disc D:
        # count h(O) * prod_{p|D} (1 - (O/p)) / ... ; for p|m (ramified) factor 1; for p|D/m factor 1-(K/p).
        if m == 2 and disc == -4:
            orders = [-4]
        else:
            orders = [-4*m] + ([-m] if m % 4 == 3 else [])
        for od in orders:
            # optimal embedding exists only if p|D does not divide conductor
            f2 = od // disc
            if any(f2 % (p*p) == 0 for p in Dp):
                continue
            fac = 1
            for p in Dp:
                if m % p == 0 and not (m == 2 and disc == -4):
                    fac *= 1  # ramified in K
                else:
                    fac *= (1 - kron(disc, p))
            nfix += classno(od) * fac
    total += nfix
    if nfix:
        print(f"w_{m}: fixed points = {nfix}")
print("sum of fixed points over nontrivial W:", total, "(Riemann-Hurwitz needs 2g-2 + 2*16 = 8 + 32 = 40)")
