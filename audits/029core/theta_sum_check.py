# -*- coding: utf-8 -*-
"""Sanity experiment over K = Q(omega) (the setting of result 003), self-written.

Purpose: test numerically the premise behind eq:reflection-angular-cancellation and
Lemma 6.1/6.3 of 029 in the one classical case where everything is explicit:
  completed cubic theta sums  T(X) = X^{-1/2} sum_{A = n b^3 primary} |b| gamma_2(n) * twist(A) * U(N(A)/X)
  (Patterson's coefficients, cube factors retained).
Prediction of the mechanism:
  - twist = 1 (no angular character): the constant term of theta produces a main term ~ c X^{1/3};
  - twist = conj(A/|A|) (angular type 1, as forced by Lambda in 029): no main term, T(X) -> 0 rapidly.
This does NOT test the quadratic-twist reflection formula itself nor any S-adic normalisation of 029.
"""
import cmath, math, sys
import numpy as np

W3 = cmath.exp(2j * cmath.pi / 3)
def norm(a, b): return a * a - a * b + b * b
def mul(x, y):
    a, b = x; c, d = y
    # (a+b w)(c+d w), w^2 = -1-w
    return (a * c - b * d, a * d + b * c - b * d)
def conj(x):
    a, b = x
    return (a - b, -b)
def tocomplex(x): return x[0] + x[1] * W3
def divides(x, y):
    """does x | y in Z[w]?  return quotient or None"""
    n = norm(*x)
    a, b = mul(y, conj(x))
    if a % n == 0 and b % n == 0:
        return (a // n, b // n)
    return None
def is_primary(x): return x[0] % 3 == 1 and x[1] % 3 == 0
UNITS = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]
def primary_assoc(x):
    for u in UNITS:
        y = mul(x, u)
        if is_primary(y):
            return y
    raise ValueError(x)

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

class Prime:
    """a primary prime pi of Z[w] prime to 3, with its cubic residue symbol and normalised Gauss sum"""
    def __init__(self, elt, p, inert):
        self.elt, self.p, self.inert = elt, p, inert
        self.N = p * p if inert else p
        if not inert:
            a, b = elt
            self.r = (-a * pow(b, -1, p)) % p            # omega = r mod pi
            assert (self.r * self.r + self.r + 1) % p == 0
            self.rp = [1, self.r, self.r * self.r % p]
    def symbol(self, x):
        """cubic residue symbol (x/pi)_3 as exponent k (value W3^k), or None if pi | x"""
        if self.inert:
            q = self.p
            a, b = x[0] % q, x[1] % q
            if a == 0 and b == 0: return None
            e = (q * q - 1) // 3
            res = (1, 0); base = (a, b)
            while e:
                if e & 1:
                    res = mul(res, base); res = (res[0] % q, res[1] % q)
                base = mul(base, base); base = (base[0] % q, base[1] % q)
                e >>= 1
            return {(1, 0): 0, (0, 1): 1, (q - 1, q - 1): 2}[res]
        p = self.p
        v = (x[0] + x[1] * self.r) % p
        if v == 0: return None
        t = pow(v, (p - 1) // 3, p)
        return self.rp.index(t)
    def gauss(self):
        """gamma_2(pi) = N^{-1/2} sum_{x mod pi} (x/pi)_3 e(x/pi),  e(z)=exp(4 pi i Im z / sqrt 3)"""
        if self.inert:
            q = self.p; tot = 0
            for i in range(q):
                for j in range(q):
                    k = self.symbol((i, j))
                    if k is None: continue
                    tot += W3 ** k * cmath.exp(2j * cmath.pi * j / q)
            return tot / q
        p = self.p; b = self.elt[1]
        e3 = (p - 1) // 3
        xs = np.arange(1, p, dtype=np.int64)
        res = np.ones(p - 1, dtype=np.int64); base = xs.copy(); e = e3
        while e:
            if e & 1: res = res * base % p
            base = base * base % p
            e >>= 1
        ks = np.where(res == 1, 0, np.where(res == self.r, 1, 2))
        assert np.all((res == 1) | (res == self.r) | (res == self.rp[2]))
        return np.sum(np.exp(2j * np.pi * ks / 3) * np.exp(-2j * np.pi * ((xs * b) % p) / p)) / math.sqrt(p)

def build_primes(Nmax):
    out = {}
    for p in primes_upto(Nmax):
        if p == 3: continue
        if p % 3 == 2:
            if p * p <= Nmax:
                out[("i", p)] = Prime(primary_assoc((p, 0)), p, True)
        else:
            # find a + b w of norm p
            found = None
            for b in range(1, int(math.sqrt(4 * p / 3)) + 2):
                # a^2 - a b + b^2 = p
                disc = b * b - 4 * (b * b - p)
                if disc < 0: continue
                s = math.isqrt(disc)
                if s * s == disc and (b + s) % 2 == 0:
                    found = ((b + s) // 2, b); break
            assert found and norm(*found) == p
            pi = primary_assoc(found)
            out[("s", p, 0)] = Prime(pi, p, False)
            out[("s", p, 1)] = Prime(primary_assoc(conj(pi)), p, False)
    return out

def brute_gauss(n):
    """direct normalised cubic Gauss sum mod a squarefree primary n (for testing multiplicativity)"""
    a, b = n; N = norm(a, b); g = math.gcd(a, b)
    zc = tocomplex(n)
    tot = 0
    for i in range(N // g):
        for j in range(g):
            x = (i, j)
            val = 1
            for P in FACT[n]:
                k = P.symbol(x)
                if k is None: val = 0; break
                val *= W3 ** k
            if val == 0: continue
            z = tocomplex(x) / zc
            tot += val * cmath.exp(4j * cmath.pi * z.imag / math.sqrt(3))
    return tot / math.sqrt(N)

def main():
    XS = [125, 250, 500, 1000, 2000, 4000, 8000, 16000, 32000]
    Nmax = 2 * XS[-1]
    PR = build_primes(Nmax)
    print("primary primes up to norm %d: %d" % (Nmax, len(PR)))
    gam = {key: P.gauss() for key, P in PR.items()}
    worst = max(abs(abs(g) - 1) for g in gam.values())
    print("max | |gamma_2(pi)| - 1 | over primes = %.2e" % worst)
    # classical check: g(pi)^3 = -N(pi) * pi  for primary pi  <=> gamma_2^3 = -pi/|pi|
    w2 = max(abs(gam[k] ** 3 + tocomplex(P.elt) / abs(tocomplex(P.elt))) for k, P in PR.items())
    w2c = max(abs(gam[k] ** 3 + (tocomplex(P.elt) / abs(tocomplex(P.elt))).conjugate()) for k, P in PR.items())
    print("max |gamma_2(pi)^3 + pi/|pi||       = %.2e   (orientation A)" % w2)
    print("max |gamma_2(pi)^3 + conj(pi)/|pi|| = %.2e   (orientation B)" % w2c)
    # enumerate primary elements and factor
    by_p = {}
    for key, P in PR.items():
        by_p.setdefault(P.p, []).append((key, P))
    R = int(math.sqrt(4 * Nmax / 3)) + 3
    terms = []          # (N(A), plain coefficient |b| gamma_2(n), alpha(A)=A/|A|, lam-symbol exponent)
    global FACT
    FACT = {}
    lam = (1, 2)
    checked_mult = 0
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            if a % 3 != 1 or b % 3 != 0: continue
            N = norm(a, b)
            if N > Nmax or N == 0: continue
            A = (a, b)
            # factor
            m = N; fac = []
            rest = A
            ok = True
            pp = 2
            mm = m
            plist = []
            d = 2
            while d * d <= mm:
                if mm % d == 0:
                    plist.append(d)
                    while mm % d == 0: mm //= d
                d += 1
            if mm > 1: plist.append(mm)
            for p in plist:
                for key, P in by_p[p]:
                    v = 0
                    while True:
                        qn = divides(P.elt, rest)
                        if qn is None: break
                        rest = qn; v += 1
                    if v: fac.append((key, P, v))
            assert rest == (1, 0), (A, rest)     # primary elements factor into primary primes with unit 1
            if any(v % 3 == 2 for _, _, v in fac):
                continue
            sq = [(key, P) for key, P, v in fac if v % 3 == 1]
            babs = 1.0
            for key, P, v in fac:
                babs *= P.N ** ((v // 3) / 2)
            g = 1.0 + 0j
            for key, P in sq: g *= gam[key]
            for i in range(len(sq)):
                for j in range(i + 1, len(sq)):
                    k1 = sq[i][1].symbol(sq[j][1].elt); k2 = sq[j][1].symbol(sq[i][1].elt)
                    g *= W3 ** (k1 + k2)
            # (lambda/n)_3 for the squarefree part
            kl = sum(P.symbol(lam) for key, P in sq) % 3
            n_elt = (1, 0)
            for key, P in sq: n_elt = mul(n_elt, P.elt)
            if len(sq) == 2 and checked_mult < 4 and norm(*n_elt) < 400:
                FACT[n_elt] = [P for _, P in sq]
                bg = brute_gauss(n_elt)
                print("multiplicativity test n=%s N=%d: brute %.6f%+.6fi  formula %.6f%+.6fi  diff=%.1e" %
                      (n_elt, norm(*n_elt), bg.real, bg.imag, g.real, g.imag, abs(bg - g)))
                checked_mult += 1
            z = tocomplex(A)
            terms.append((N, babs * g, z / abs(z), kl))
    terms.sort(key=lambda t: t[0])
    Ns = np.array([t[0] for t in terms], dtype=float)
    co = np.array([t[1] for t in terms])
    al = np.array([t[2] for t in terms])
    kl = np.array([t[3] for t in terms])
    lamf = np.exp(-2j * np.pi * kl / 3)         # conj((lambda/n)_3), Patterson's extra factor
    print("number of theta-supported primary A with N(A)<=%d: %d" % (Nmax, len(terms)))
    def U(x):
        t = (x - 1.25) / 0.75
        out = np.zeros_like(x)
        m = np.abs(t) < 1
        out[m] = np.exp(-1.0 / (1 - t[m] ** 2))
        return out
    print("\n  X      |T_plain|   T_plain/X^(1/3)       |T_ang|     |T_ang3|   |T_plain,lam|/X^(1/3)  |T_ang,lam|   |T_angbar|")
    for X in XS:
        w = U(Ns / X) / math.sqrt(X)
        Tp = np.sum(co * w)
        Ta = np.sum(co * np.conj(al) * w)
        Ta3 = np.sum(co * np.conj(al) ** 3 * w)
        Tpl = np.sum(co * lamf * w)
        Tal = np.sum(co * lamf * np.conj(al) * w)
        Tab = np.sum(co * al * w)
        r = Tp / X ** (1 / 3)
        print("%6d  %9.5f   %8.5f%+8.5fi   %10.3e  %10.3e   %10.5f          %10.3e   %10.3e   Tangbar=%+.4f%+.4fi  Tang=%+.4f%+.4fi" %
              (X, abs(Tp), r.real, r.imag, abs(Ta), abs(Ta3), abs(Tpl) / X ** (1 / 3), abs(Tal), abs(Tab), Tab.real, Tab.imag, Ta.real, Ta.imag))
    # Second weight: narrow log-Gaussian  U_w(t) = exp(-(log t)^2/(4w)); its Mellin transform is
    # sqrt(4 pi w) exp(w s^2), entire with Gaussian vertical decay, so an entire Dirichlet series
    # gives decay faster than any power of X (beyond the conductor).  Truncation at N<=Nmax: only
    # X with exp(-(log(Nmax/X))^2/(4w)) < 1e-9 are reported.
    for wdt in (0.05, 0.1):
        print("\n  log-Gaussian weight, w=%.2f" % wdt)
        print("  X      |T_plain|   T_plain/X^(1/3)    |T_ang|      |T_ang3|     |T_angbar|   |T_ang,lam|   truncation")
        for X in [125, 250, 500, 1000, 2000, 4000, 8000, 16000]:
            trunc = math.exp(-(math.log(Nmax / X)) ** 2 / (4 * wdt))
            if trunc > 1e-9: continue
            wv = np.exp(-(np.log(Ns / X)) ** 2 / (4 * wdt)) / math.sqrt(X)
            Tp = np.sum(co * wv); Ta = np.sum(co * np.conj(al) * wv); Ta3 = np.sum(co * np.conj(al) ** 3 * wv)
            Tab = np.sum(co * al * wv); Tal = np.sum(co * lamf * np.conj(al) * wv)
            print("%6d  %9.5f   %9.5f          %10.3e   %10.3e   %10.3e   %10.3e   %.1e" %
                  (X, abs(Tp), abs(Tp) / X ** (1 / 3), abs(Ta), abs(Ta3), abs(Tab), abs(Tal), trunc))
if __name__ == "__main__":
    main()
