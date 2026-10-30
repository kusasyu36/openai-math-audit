"""Numerical test of Lemma 3.2 (H10 over Q, Section 3.2) and Lemma 3.1.
For integer b, pick m, m_d satisfying (3.3),(3.4); sample rational t satisfying the local
preconditions used in the proof (H'>0, H' square at S(m), v_q(d t^2)>0 at Supp(m_d)\\S(m));
then check: B,J disjoint, avoid 2,3; theta square at q in B; beta square at q in J;
beta, theta squares in Q_2; theta square in Q_3; beta unit nonsquare at 3.
"""
import random
from fractions import Fraction
from math import isqrt

def factor(n):
    n = abs(n); f = {}; p = 2
    while p*p <= n:
        while n % p == 0: f[p] = f.get(p, 0)+1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f

def v(x, p):
    x = Fraction(x)
    if x == 0: return 10**9
    def vi(n):
        k = 0
        while n % p == 0: n //= p; k += 1
        return k
    return vi(abs(x.numerator)) - vi(x.denominator)

def unit(x, p, mod):
    x = Fraction(x); k = v(x, p)
    y = x / Fraction(p)**k
    return (y.numerator * pow(y.denominator, -1, mod)) % mod

def is_sq_Qp(x, p):
    x = Fraction(x)
    if x == 0: return False
    if v(x, p) % 2: return False
    if p == 2: return unit(x, 2, 8) == 1
    u = unit(x, p, p)
    return pow(u, (p-1)//2, p) == 1

def sqfree_class(x):
    x = Fraction(x)
    n = x.numerator * x.denominator  # same square class
    s = 1 if n > 0 else -1
    r = 1
    for p, e in factor(n).items():
        if e % 2: r *= p
    return s*r

def cond33(bp):
    return bp > 0 and is_sq_Qp(bp, 2) and v(bp, 3) % 2 == 0 and not is_sq_Qp(bp, 3)

rng = random.Random(2026)
M_list = [s*a for a in [1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,26,30,33,35,39,42,51,55,57,65,66,70,78] for s in (1,-1)]
tested = 0; fails = 0
for b in list(range(1, 400)) + list(range(-200, 0)):
    m = next((m for m in M_list for _ in [0] if cond33(m*b)), None)
    if m is None:
        # m with 3-adic and 2-adic adjustments
        cand = [s*2**i*3**j*u for s in (1,-1) for i in range(2) for j in range(2) for u in (1,5,7,11,13,17,19,23)]
        m = next((m for m in cand if cond33(m*b)), None)
    assert m is not None, b
    Sm = {2, 3} | set(factor(m))
    md = next((x for x in [s*2**i*3**j*u for u in (1,5,7,11,13,17,19,23,29,31,37,41,43) for s in (1,-1) for i in range(2) for j in range(2)]
               if x*b > 0 and all(is_sq_Qp(x*b, q) for q in Sm)), None)
    if md is None: continue
    d = md*b; bp = m*b
    extra = set(factor(md)) - Sm
    for attempt in range(4000):
        t = Fraction(rng.randint(-400, 400), rng.randint(1, 400))
        if 1 + d*t*t == 0: continue
        if any(v(d*t*t, q) <= 0 for q in extra): continue
        U = b*(1 - d*t*t)/(1 + d*t*t)
        Hp = bp*(bp - U*U)
        if Hp <= 0: continue
        if not all(is_sq_Qp(Hp, q) for q in Sm): continue
        beta = sqfree_class(bp); theta = sqfree_class(Hp)
        B = set(factor(beta)); J = set(factor(theta))
        ok = beta > 0 and theta > 0 and not (B & J) and not ({2,3} & (B|J))
        ok &= all(is_sq_Qp(theta, q) for q in B) and all(is_sq_Qp(beta, q) for q in J)
        ok &= is_sq_Qp(beta, 2) and is_sq_Qp(theta, 2) and is_sq_Qp(theta, 3)
        ok &= (v(beta, 3) == 0 and not is_sq_Qp(beta, 3)) and len(B) > 0
        tested += 1
        if not ok:
            fails += 1
            print("FAIL", b, m, md, t, beta, theta)
        break
print("Lemma 3.2 instances:", tested, "failures:", fails)

# Lemma 3.1 spot check: on Y^2 = X(X-L)(X+3L), L a q-unit (q>3), every nonzero X has even valuation.
bad = 0
for q in [5, 7, 11, 13]:
    for L in [1, 2, -1, 4, 6]:
        if L % q == 0: continue
        for _ in range(20000):
            X = Fraction(rng.randint(1, 3000), 1) * Fraction(q) ** rng.randint(-3, 3) * rng.choice([1, -1])
            f = X*(X-L)*(X+3*L)
            if f > 0 and is_sq_Qp(f, q) and v(X, q) % 2:
                bad += 1
print("Lemma 3.1 local parity violations found:", bad)
