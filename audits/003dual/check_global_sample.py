# -*- coding: utf-8 -*-
"""
003 (The Quasi-Riemann Hypothesis), paper2.tex -- audits/003dual, supplementary SAMPLE.

NOT part of the stage-1 gate (that is check_local.py).  check_local.py proves the
local identity eq:ray-local-transform for abstract sigma, eps in (O/p)^x.  This
script checks, on a finite sample of actual configurations, that the a, c,
delta', sigma_p, eps_p constructed in the proof (paper2.tex:3039-3230) are of
that form, i.e. the glue between the global construction and the local identity.

Construction followed (tex lines):
  z_h = lambda^2 (h_0/L + sum_p h_p/p) = a/c in lowest terms            :2970, :3046
  normalisation: a = 1 mod 3 if lambda | c, else c = 1 mod 3            :3047
  M = lambda^12 L^4, representatives h_p = 0 mod M^2                    :3039-3041
  delta' by CRT: a delta' = 1 mod q^{v_q(M c_0)} (q | c_0),
                 delta' = 0 mod q^{v_q(M)} (q | M, q not | c_0),
                 a delta' = 1 mod r                                     :3057-3059
  b_g = (a delta' - 1)/c, g = [[a, b_g], [c, delta']]                   :3063-3064
  H (three cases), g_1 = g H^-1 = I mod 3                               :3076-3089
  sigma_p = lambda^2 c/p, eps_p = -((lambda^3 c/p) sigma_p)^-1 in O/p   :3211-3212
  D_0 = lambda^3 c_0, psi(x) = e(-delta'_0 r^-1 x / D_0)                :3209, :3225

Sample: L = 18 (so that all three cases 3 | c, v_lambda(c) = 1, (c, lambda) = 1
occur), every h_0 in O/(18) (324 residues, (h_0, L) != 1 and h_0 = 0 included),
r in { p7, p13, -5, p7*p13 }, three tuples (h_p) per (h_0, r), and x = lambda^4 l
over the box [-3, 3]^2 of O.  Here c_0 := c / r (the tex says c = c_0 r "up to a
unit"; the unit is put into c_0).

CRITERIA (fixed before running).  For every configuration:
  S1  r | c in O, and c/r is an associate of the reduced denominator of lambda^2 h_0/L
  S2  a = sigma_p h_p (mod p) for every p | r
  S3  b_g in O (so det g = 1), and b_g = 0 mod q^{v_q(M)} (q | c_0),
      b_g = -c^-1 mod q^{v_q(M)} (q | M, q not | c_0)                    :3068-3075
  S4  eq:ray-additive-crt:  e~(-delta' l / c) = psi(lambda^4 l) prod_p e(eps_p h_p^-1 lambda^4 l / p)
      for every x = lambda^4 l in the box (exact, as rational phases mod 1)
  S5  for fixed (h_0, r): c and delta' mod D_0 do not change when the h_p vary
      (this is what lets the sums over h_p factor; "fixed by the choices above", :3228)
  S6  g_1 = g H^-1 = I (mod 3) with H as in the tex
  C   coverage: all three cases of H occur in the sample
  Negative controls (must FAIL in every configuration):
  NS2 sigma_p replaced by lambda^3 c/p in S2
  NS4 eps_p replaced by +((lambda^3 c/p) sigma_p)^-1 in S4
  Added AFTER the first run (the first run had S1-S6, C, NS2, NS4 and met all of
  them; nothing above was changed):
  NS5 control for S5: take the representatives h_p WITHOUT the factor M^2.  Then
      (c, delta' mod D_0) must change with h_p for at least one (h_0, r) for each
      of the four r.  (It cannot change when c_0 = 1 and (c, lambda) = 1 ... so
      "every (h_0, r)" is not the criterion.)  S6 has no negative control.
Assumed, not tested: eq:ray-multiplier, everything about theta.  S4 is close to a
tautology once a delta' = 1 (mod r); the discriminating checks are S2, S5, S6.

Run (standard library only, about a minute):
  python3 -u audits/003dual/check_global_sample.py > audits/003dual/check_global_sample.out
"""
import sys, time
from fractions import Fraction as Fr

UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
ZERO, ONE = (0, 0), (1, 0)
LAM, TWO = (1, 2), (2, 0)

def mul(x, y):
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def neg(x): return (-x[0], -x[1])
def conj(x): return (x[0] - x[1], -x[1])
def norm(x): return x[0] * x[0] - x[0] * x[1] + x[1] * x[1]
def powe(x, n):
    r = ONE
    for _ in range(n): r = mul(r, x)
    return r
def prod(xs):
    r = ONE
    for x in xs: r = mul(r, x)
    return r
def exact_div(y, x):
    n = norm(x); a, b = mul(y, conj(x))
    if a % n == 0 and b % n == 0: return (a // n, b // n)
    return None
def divides(x, y): return exact_div(y, x) is not None
def divmod_e(y, x):
    n = norm(x); a, b = mul(y, conj(x))
    q = ((2 * a + n) // (2 * n), (2 * b + n) // (2 * n))
    return q, sub(y, mul(q, x))                      # N(r) <= 3/4 N(x)
def xgcd(a, c):
    r0, r1 = a, c; s0, s1 = ONE, ZERO; t0, t1 = ZERO, ONE
    while r1 != ZERO:
        q, r = divmod_e(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, sub(s0, mul(q, s1))
        t0, t1 = t1, sub(t0, mul(q, t1))
    return r0, s0, t0
def gcd(a, c): return xgcd(a, c)[0]
def inv_mod(a, m):
    g, s, _ = xgcd(a, m)
    assert norm(g) == 1, "not invertible"
    return mul(s, conj(g))                           # g^-1 = conj(g) for a unit
def val(x, q):
    v = 0
    while True:
        y = exact_div(x, q)
        if y is None: return v
        x = y; v += 1
def is_primary(x): return x[0] % 3 == 1 and x[1] % 3 == 0
def mod3(x): return (x[0] % 3, x[1] % 3)
def T(num, den):
    """Tr(num / (den * lambda)) mod 1, so that e(num/den) = exp(2 pi i T)"""
    return Fr(mul(num, conj(den))[1], norm(den)) % 1
def crt(pairs):
    """x with x = t_i mod m_i (pairwise coprime moduli)"""
    x, Mod = pairs[0]
    for t, m in pairs[1:]:
        k = mul(sub(t, x), inv_mod(Mod, m))
        x = add(x, mul(Mod, k)); Mod = mul(Mod, m)
    return x


P7, P13, M5 = (1, 3), (4, 3), (-5, 0)
assert norm(P7) == 7 and norm(P13) == 13 and all(is_primary(p) for p in (P7, P13, M5))
L = (18, 0)
LAM2, LAM3 = powe(LAM, 2), powe(LAM, 3)
Mbig = mul(powe(LAM, 12), powe(L, 4)); M2 = mul(Mbig, Mbig)
VM = {LAM: val(Mbig, LAM), TWO: val(Mbig, TWO)}
BOX = [(a, b) for a in range(-3, 4) for b in range(-3, 4)]
KTUPLES = {1: [[(1, 0)], [(2, 0)], [(1, 1)]],
           2: [[(1, 0), (1, 0)], [(2, 0), (1, 1)], [(0, 1), (3, 0)]]}

NAMES = ["S1", "S2", "S3", "S4", "S5", "S6", "NS2", "NS4"]
bad = {n: 0 for n in NAMES}; tot = {n: 0 for n in NAMES}
cases = {"3|c": 0, "v=1": 0, "coprime": 0}


def build(h0, primes, ks, raw=False):
    r = prod(primes)
    hp = list(ks) if raw else [mul(M2, k) for k in ks]
    num = mul(LAM2, add(mul(h0, r), mul(L, _sum([mul(h, exact_div(r, p)) for h, p in zip(hp, primes)]))))
    den = mul(L, r)
    g = gcd(num, den)
    a, c = exact_div(num, g), exact_div(den, g)
    if divides(LAM, c):
        u = [u for u in UNITS if is_primary(mul(a, u))]
    else:
        u = [u for u in UNITS if is_primary(mul(c, u))]
    assert len(u) == 1
    a, c = mul(a, u[0]), mul(c, u[0])
    return r, hp, a, c

def _sum(xs):
    s = ZERO
    for x in xs: s = add(s, x)
    return s


def one_config(h0, primes, ks, raw=False):
    r, hp, a, c = build(h0, primes, ks, raw)
    out = {}
    # S1
    c0 = exact_div(c, r)
    n0 = mul(LAM2, h0)
    den0 = exact_div(L, gcd(n0, L)) if h0 != ZERO else ONE
    out["S1"] = c0 is not None and divides(c0, den0) and divides(den0, c0)
    if c0 is None: return out, None
    # delta'
    pairs = []
    for q in (LAM, TWO):
        if divides(q, c0):
            m = powe(q, VM[q] + val(c0, q)); pairs.append((inv_mod(a, m), m))
        else:
            pairs.append((ZERO, powe(q, VM[q])))
    pairs.append((inv_mod(a, r), r))
    dl = crt(pairs)
    # S3
    bg = exact_div(sub(mul(a, dl), ONE), c)
    ok3 = bg is not None
    if ok3:
        for q in (LAM, TWO):
            m = powe(q, VM[q])
            if divides(q, c0): ok3 &= divides(m, bg)
            else: ok3 &= divides(m, add(mul(bg, c), ONE))
    out["S3"] = ok3
    # S6
    ok6 = False
    if bg is not None:
        vl = val(c, LAM)
        if vl >= 2:
            g1 = (a, bg, c, dl); cases["3|c"] += 1
        elif vl == 1:
            u0 = [u for u in (LAM, neg(LAM)) if mod3(sub(u, c)) == (0, 0)]
            assert len(u0) == 1; u0 = u0[0]; cases["v=1"] += 1
            g1 = (sub(a, mul(bg, u0)), bg, sub(c, mul(dl, u0)), dl)
        else:
            u0 = mod3(a); cases["coprime"] += 1
            g1 = (neg(bg), add(a, mul(bg, u0)), neg(dl), add(c, mul(dl, u0)))
        ok6 = [mod3(e) for e in g1] == [(1, 0), (0, 0), (0, 0), (1, 0)]
    out["S6"] = ok6
    # S2, eps
    ok2 = True; nok2 = True; eps = []; hinv = []
    for p, h in zip(primes, hp):
        cp = exact_div(c, p)
        sig = mul(LAM2, cp)
        ok2 &= divides(p, sub(a, mul(sig, h)))
        nok2 &= divides(p, sub(a, mul(mul(LAM3, cp), h)))
        eps.append(neg(inv_mod(mul(mul(LAM3, cp), sig), p)))
        hinv.append(inv_mod(h, p))
    out["S2"] = ok2; out["NS2"] = nok2
    # S4
    D0 = mul(LAM3, c0); rinv = inv_mod(r, D0); lc = mul(LAM3, c)
    ok4 = True; nok4 = True
    for x in BOX:
        lhs = T(neg(mul(dl, x)), lc)
        psi = T(neg(mul(mul(dl, rinv), x)), D0)
        loc = [T(mul(mul(e, hi), x), p) for e, hi, p in zip(eps, hinv, primes)]
        ok4 &= (lhs - psi - sum(loc)) % 1 == 0
        nok4 &= (lhs - psi + sum(loc)) % 1 == 0
    out["S4"] = ok4; out["NS4"] = nok4
    return out, (c, dl, D0)


def main():
    t0 = time.time()
    print("003 audits/003dual: SAMPLE of the global construction behind eq:ray-local-transform")
    print("not part of the stage-1 gate; criteria in the header (fixed before the first run)\n")
    nconf = 0; ns5 = []
    for primes in ([P7], [P13], [M5], [P7, P13]):
        nchanged = 0
        for a0 in range(18):
            for b0 in range(18):
                h0 = (a0, b0); fixed = []
                for ks in KTUPLES[len(primes)]:
                    out, data = one_config(h0, primes, ks)
                    nconf += 1
                    for n, v in out.items():
                        tot[n] += 1
                        if n.startswith("N"): bad[n] += (1 if v else 0)     # control "passes" = bad
                        else: bad[n] += (0 if v else 1)
                    fixed.append(data)
                tot["S5"] += 1
                ok5 = all(d is not None for d in fixed)
                if ok5:
                    c1, d1, D1 = fixed[0]
                    ok5 = all(d[0] == c1 and divides(D1, sub(d[1], d1)) for d in fixed[1:])
                bad["S5"] += (0 if ok5 else 1)
                keep = dict(cases)                              # the raw runs must not count in C
                raw = [one_config(h0, primes, ks, raw=True)[1] for ks in KTUPLES[len(primes)]]
                cases.clear(); cases.update(keep)
                changed = any(d is None for d in raw) or any(
                    d[0] != raw[0][0] or not divides(raw[0][2], sub(d[1], raw[0][1])) for d in raw[1:])
                nchanged += 1 if changed else 0
        ns5.append((primes, nchanged))
        print("r = %s done, %d configurations so far [%.0fs]" % (primes, nconf, time.time() - t0)); sys.stdout.flush()
    print()
    allok = True
    for n in NAMES:
        if n.startswith("N"):
            ok = bad[n] == 0
            print("  %-4s %s  control held (i.e. did NOT fail) in %d / %d configurations (must be 0)" % (n, "PASS" if ok else "FAIL", bad[n], tot[n]))
        else:
            ok = bad[n] == 0
            print("  %-4s %s  %d failures / %d" % (n, "PASS" if ok else "FAIL", bad[n], tot[n]))
        allok &= ok
    ok5c = all(n > 0 for _, n in ns5); allok &= ok5c
    print("  %-4s %s  without M^2: (c, delta' mod D_0) changed in %s of 324 (h_0) per r (must be > 0 for each r)"
          % ("NS5", "PASS" if ok5c else "FAIL", [n for _, n in ns5]))
    okc = all(v > 0 for v in cases.values()); allok &= okc
    print("  %-4s %s  cases of H: %s" % ("C", "PASS" if okc else "FAIL", cases))
    print("\nconfigurations: %d ; x per configuration: %d" % (nconf, len(BOX)))
    print("SAMPLE (not a gate): %s" % ("all criteria met" if allok else "SOME CRITERIA NOT MET"))
    print("elapsed %.0fs" % (time.time() - t0))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
