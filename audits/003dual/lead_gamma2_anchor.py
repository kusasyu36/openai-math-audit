# -*- coding: utf-8 -*-
"""
Lead's anchor for the DIRECT side of layer E (written by the lead after the sub-auditor's report).

stage3_d_mellin.py takes  aW = conj(td.G) * W3**td.kl  from kubota.ThetaData as "gamma_2(n) |b|" for A = n b^3.
Nothing in the stage scripts ties this to the paper's definition.  Here gamma_2(n) is computed by brute force
from paper2.tex:310-313,
    gamma_2(n) = N(n)^-1/2 sum_{x mod n} chi_n(x)^2 e(x/n),  e(z) = exp(4 pi i Im z / sqrt 3),
with chi_n^2 = cubic residue symbol, computed HERE from x^((N(pi)-1)/3) mod pi (own code; kubota.py is used
only for ring arithmetic, for the factorisation of A, and to read the table).

CRITERION (fixed before running): for EVERY table entry A = n b^3 prime to 2 with N(A) <= 700,
    | aW(A) - gamma_2(n) |b| | <= 1e-10 ;  the set must contain n with 1, 2 and 3 prime factors and some b != 1.
Negative control: gamma_1 in place of gamma_2 (sextic symbol instead of cubic is not available for all n; use
    conj of the cubic symbol, i.e. gamma_4): must differ by > 1e-3 for at least one A.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy python -B -u audits/003dual/lead_gamma2_anchor.py > audits/003dual/lead_gamma2_anchor.out
"""
import sys, os, math, cmath
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "029reflect"))
import kubota as kb
W3 = kb.W3

def powmod(x, e, m):
    r = (1, 0); b = kb.divmod_e(x, m)[1]
    while e:
        if e & 1: r = kb.divmod_e(kb.mul(r, b), m)[1]
        b = kb.divmod_e(kb.mul(b, b), m)[1]; e >>= 1
    return r
def cub_own(x, pi):
    """exponent k of (x/pi)_3 = omega^k, or None if pi | x"""
    if kb.exact_div(x, pi) is not None: return None
    t = powmod(x, (kb.norm(pi) - 1) // 3, pi)
    for k, w in enumerate(((1, 0), (0, 1), (-1, -1))):
        if kb.exact_div(kb.sub(t, w), pi) is not None: return k
    raise AssertionError
def gamma_brute(n, primes, power):
    N = kb.norm(n); g = math.gcd(n[0], n[1]); nc = kb.conj(n); tot = 0
    for i in range(N // g):
        for j in range(g):
            x = (i, j); k = 0
            for pi in primes:
                c = cub_own(x, pi)
                if c is None: k = None; break
                k += c
            if k is None: continue
            tot += W3 ** ((power * k) % 3) * cmath.exp(2j * cmath.pi * kb.mul(x, nc)[1] / N)
    return tot / math.sqrt(N)

td = kb.ThetaData(200000, os.path.join(here, "cache"), verbose=False)
worst = 0.0; cnt = {}; nb = 0; ctl = 0.0; ntot = 0
for idx in range(len(td.A)):
    A = (int(td.A[idx, 0]), int(td.A[idx, 1])); NA = kb.norm(A)
    if NA > 700 or (A[0] % 2 == 0 and A[1] % 2 == 0): continue
    u, k, fac = kb.factor(A)
    primes = [pr.elt for pr, e in fac if e % 3 == 1]
    n = (1, 0)
    for p in primes: n = kb.mul(n, p)
    babs = math.sqrt(NA / kb.norm(n)) ** (1.0 / 3)
    aW = complex(td.G[idx]).conjugate() * W3 ** int(td.kl[idx])
    g2 = gamma_brute(n, primes, 1) if primes else 1.0       # chi^2 = cubic symbol
    g4 = gamma_brute(n, primes, 2) if primes else 1.0
    d = abs(aW - g2 * babs); worst = max(worst, d); ctl = max(ctl, abs(aW - g4 * babs))
    cnt[len(primes)] = cnt.get(len(primes), 0) + 1; nb += babs > 1.0001; ntot += 1
ok = worst <= 1e-10 and all(cnt.get(k, 0) > 0 for k in (1, 2, 3)) and nb > 0
print("entries A = n b^3 prime to 2 with N(A) <= 700: %d ; by number of prime factors of n: %s ; with b != 1: %d" % (ntot, dict(sorted(cnt.items())), nb))
print("G2   %s  max |aW - gamma_2(n)|b|| = %.2e (tolerance 1e-10)" % ("PASS" if ok else "FAIL", worst))
print("NG2  %s  control gamma_4: max difference %.2e (must be > 1e-3)" % ("PASS" if ctl > 1e-3 else "FAIL", ctl))
