"""Brute-force 2-Selmer check of Lemma 3.3 / Lemma 3.4 (Section 3) for E_l: y^2 = x(x-l)(x+3l).
Kummer coordinates (x+3l, x-l).  Pure python, run with -I.
"""
import itertools, random, sys
from math import gcd, isqrt
from fractions import Fraction

def factor(n):
    f = {}; p = 2
    while p*p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0)+1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f

def isprime(n):
    return n > 1 and factor(n) == {n: 1}

def vp(n, p):
    if n == 0: return 10**9
    k = 0
    while n % p == 0: n //= p; k += 1
    return k

def sqclass(x, p):
    """square class of nonzero Fraction x in Q_p^*/Q_p^*2 (p prime) or R (p=0)."""
    if p == 0:
        return 1 if x > 0 else -1
    a, b = x.numerator, x.denominator
    v = vp(abs(a), p) - vp(b, p)
    u_num = a // p**vp(abs(a), p); u_den = b // p**vp(b, p)
    if p == 2:
        u = (u_num * pow(u_den, -1, 8)) % 8
        return (v % 2, u)
    u = (u_num * pow(u_den, -1, p)) % p
    leg = pow(u, (p-1)//2, p)
    return (v % 2, 1 if leg == 1 else -1)

def is_local_square(x, p):
    c = sqclass(x, p)
    if p == 0: return c == 1
    if p == 2: return c == (0, 1)
    return c == (0, 1)

def local_image(l, p, target_size):
    """sample local points of E_l over Q_p (or R if p=0), return set of Kummer classes."""
    img = set()
    # two-torsion replacements
    e = [0, l, -3*l]  # roots
    # coordinates: first = x+3l (root -3l), second = x-l (root l)
    def cls(a1, a2):
        return (sqclass(Fraction(a1), p), sqclass(Fraction(a2), p))
    img.add(cls(1, 1))
    # (0,0): x=0 -> (3l, -l)
    img.add(cls(3*l, -l))
    # (l,0): first = 4l, second replaced by (l - 0)(l+3l) = 4l^2 ; class of second = prod_{j!=i}(e_i-e_j)
    img.add(cls(4*l, (l-0)*(l+3*l)))
    # (-3l,0): first replaced by (-3l-0)(-3l-l)=12 l^2 ; second = -4l
    img.add(cls((-3*l)*(-4*l), -4*l))
    pp = p if p else 2
    rng = random.Random(1234 + p)
    tries = 0
    while len(img) < target_size and tries < 200000:
        tries += 1
        j = rng.randint(-4, 4)
        n = rng.randint(-3000, 3000)
        if n == 0: continue
        x = Fraction(n) * (Fraction(pp)**j if p else Fraction(1, rng.randint(1, 50)))
        if x in (0, l, -3*l): continue
        fx = x*(x-l)*(x+3*l)
        if is_local_square(fx, p):
            img.add(cls(x+3*l, x-l))
    return img

def selmer(l):
    primes = sorted(set([2, 3] + list(factor(l).keys())))
    places = [0] + primes
    sizes = {0: 2, 2: 8}
    imgs = {}
    for p in places:
        size = sizes.get(p, 4)
        imgs[p] = local_image(l, p, size)
        if len(imgs[p]) != size:
            raise RuntimeError(f"local image at {p} incomplete: {len(imgs[p])}/{size}")
    gens = [-1] + primes
    sel = []
    for e1 in itertools.product([0, 1], repeat=len(gens)):
        a1 = 1
        for g, e in zip(gens, e1):
            if e: a1 *= g
        for e2 in itertools.product([0, 1], repeat=len(gens)):
            a2 = 1
            for g, e in zip(gens, e2):
                if e: a2 *= g
            ok = all((sqclass(Fraction(a1), p), sqclass(Fraction(a2), p)) in imgs[p] for p in places)
            if ok: sel.append((a1, a2))
    return sel

def bit_sq(a, q):  # nonsquare bit [a/q]
    a %= q
    return 0 if pow(a, (q-1)//2, q) == 1 else 1

def matrices(Q):
    n = len(Q)
    s = [bit_sq(-1, q) for q in Q]
    lam = [bit_sq(3, q) for q in Q]
    K = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                b = bit_sq(Q[j], Q[i])
                K[i][j] ^= b
                K[i][i] ^= b
    return K, s, lam

def matvec(M, x):
    return [sum(M[i][j]*x[j] for j in range(len(x))) % 2 for i in range(len(M))]

def kernel(M):
    n = len(M[0])
    return [x for x in itertools.product([0, 1], repeat=n) if not any(matvec(M, list(x)))]

def image(M):
    n = len(M[0])
    return {tuple(matvec(M, list(x))) for x in itertools.product([0, 1], repeat=n)}

def hyp_311(Q, B):
    K, s, lam = matrices(Q)
    ker = set(kernel(K))
    one = tuple([1]*len(Q)); oneB = tuple(1 if q in B else 0 for q in Q)
    span = {tuple([0]*len(Q)), one, oneB, tuple((a+b) % 2 for a, b in zip(one, oneB))}
    return ker == span and tuple(lam) not in image(K)

def find_point(l, beta, bound=400):
    for u in range(1, bound):
        for w in range(1, bound):
            if gcd(u, w) != 1: continue
            X = beta*u*u; W2 = w*w
            val = X*(X - l*W2)*(X + 3*l*W2)
            if val > 0 and isqrt(val)**2 == val:
                return Fraction(X, W2)
    return None

if __name__ == "__main__":
    random.seed(7)
    # candidate primes != 2,3
    P = [p for p in range(5, 120) if isprime(p)]
    checked = 0; agree = 0; pts = 0; nopts = []
    other = {}
    for trial in range(4000):
        k = random.choice([2, 3, 3, 4])
        Q = sorted(random.sample(P, k))
        l = 1
        for q in Q: l *= q
        if l % 4 != 3 or l % 3 != 2: continue
        for r in range(1, len(Q)):
            for B in itertools.combinations(Q, r):
                beta = 1
                for q in B: beta *= q
                if beta % 8 != 1 or beta % 3 != 2: continue
                if not hyp_311(Q, B): continue
                checked += 1
                sel = selmer(l)
                dim = len(sel).bit_length()-1
                ok = (len(sel) == 8) and ((beta, 1) in sel)
                agree += ok
                if not ok:
                    print("COUNTEREXAMPLE?", l, Q, B, len(sel))
                if checked <= 25:
                    pt = find_point(l, beta, 300)
                    if pt is not None: pts += 1
                    else: nopts.append((l, beta))
                    print(f"l={l} Q={Q} B={B} beta={beta} |Sel2|={len(sel)} (beta,1) in Sel: {(beta,1) in sel}  point x in beta*Q^2: {pt}")
        if checked >= 60: break
    print(f"\nLemma 3.3 instances checked: {checked}, consistent: {agree}")
    print(f"points found among first 25: {pts}; not found within search bound: {nopts}")
