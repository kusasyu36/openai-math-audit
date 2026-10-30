"""Independent full 2-descent for y^2=(x-e1)(x-e2)(x-e3), e_i integers.
Computes dim_F2 Sel_2(E/Q) via sampled local Kummer images at p in S and infinity.
Written for audit 004b (pointwise 2-converse).  No repository code used.
"""
import random, itertools, sys
from fractions import Fraction

def factor(n):
    n = abs(n); f = []; d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f

def vp(n, p):
    if n == 0: return 10**9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v

def loc_class(q, p):
    """class of nonzero rational q in Q_p^x/(Q_p^x)^2 as tuple of bits; p=0 means R."""
    q = Fraction(q)
    if p == 0:
        return (1 if q < 0 else 0,)
    a, b = q.numerator, q.denominator
    v = vp(a, p) - vp(b, p)
    u = Fraction(a, p ** vp(a, p)) / Fraction(b, p ** vp(b, p))
    if p == 2:
        um = (u.numerator * pow(u.denominator, -1, 8)) % 8
        return (v % 2, 1 if um in (3, 7) else 0, 1 if um in (5, 7) else 0)  # bits: val, -1, 5
    um = (u.numerator * pow(u.denominator, -1, p)) % p
    leg = pow(um, (p - 1) // 2, p)
    return (v % 2, 0 if leg == 1 else 1)

def is_square_loc(q, p):
    return all(b == 0 for b in loc_class(q, p))

def kummer(x, es, p):
    e1, e2, e3 = es
    x = Fraction(x)
    out = []
    for i in (0, 1):
        ei = es[i]
        if x == ei:
            others = [e for j, e in enumerate(es) if j != i]
            val = (ei - others[0]) * (ei - others[1])
        else:
            val = x - ei
        out.append(val)
    return loc_class(out[0], p) + loc_class(out[1], p)

def span(vecs):
    basis = []
    for v in vecs:
        v = list(v)
        for b in basis:
            piv = next(i for i, c in enumerate(b) if c)
            if v[piv]: v = [(a + c) % 2 for a, c in zip(v, b)]
        if any(v): basis.append(v)
    return basis

def in_span(basis, v):
    v = list(v)
    for b in basis:
        piv = next(i for i, c in enumerate(b) if c)
        if v[piv]: v = [(a + c) % 2 for a, c in zip(v, b)]
    return not any(v)

def local_image(es, p, rng):
    f = lambda x: (x - es[0]) * (x - es[1]) * (x - es[2])
    target = 1 if p == 0 else (3 if p == 2 else 2)
    vecs = [kummer(e, es, p) for e in es]  # 2-torsion
    basis = span(vecs)
    tries = 0
    while len(basis) < target and tries < 200000:
        tries += 1
        if p == 0:
            x = Fraction(rng.randint(-10**6, 10**6), rng.randint(1, 1000))
        else:
            k = rng.randint(-3, 6)
            x = Fraction(rng.randint(-p**8, p**8)) * Fraction(p) ** (-2 * k) if k > 0 else Fraction(rng.randint(-p**8, p**8))
            x += Fraction(rng.randint(0, p**3)) * rng.choice([1, p, p*p])
        fx = f(x)
        if fx == 0: continue
        if not is_square_loc(fx, p): continue
        basis = span(basis + [kummer(x, es, p)])
    assert len(basis) == target, (p, len(basis), target)
    return basis

def selmer_dim(es, seed=1):
    rng = random.Random(seed)
    e1, e2, e3 = es
    disc = 2 * (e1 - e2) * (e1 - e3) * (e2 - e3)
    S = sorted(set(factor(disc)) | {2})
    gens = [-1] + S
    places = [0] + S
    imgs = {p: local_image(es, p, rng) for p in places}
    count = 0
    elems = []
    for bits in itertools.product((0, 1), repeat=2 * len(gens)):
        d1 = 1; d2 = 1
        for g, b in zip(gens, bits[:len(gens)]):
            if b: d1 *= g
        for g, b in zip(gens, bits[len(gens):]):
            if b: d2 *= g
        ok = all(in_span(imgs[p], loc_class(d1, p) + loc_class(d2, p)) for p in places)
        if ok:
            count += 1; elems.append((d1, d2))
    import math
    return int(round(math.log2(count))), elems

if __name__ == "__main__":
    cases = []
    for n in (1, 2, 3, 5, 6, 7, 14, 34, 41):
        cases.append((f"congruent n={n}: y^2=x^3-{n}^2x", (0, n, -n)))
    for l in (1751, 5963, 1271, 23579, 323, 1763, 14663):
        cases.append((f"E_l l={l}: y^2=x(x-l)(x+3l)", (0, l, -3 * l)))
    for name, es in cases:
        d, el = selmer_dim(es)
        print(f"{name}: dim Sel2 = {d}  => s_2 <= {d-2}")
        sys.stdout.flush()
