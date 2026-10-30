"""Exact-arithmetic checks for 074 (Kakeya maximal, 3D), Section 10.1 closing cases,
Lemma 5.1 derivative bookkeeping, (5.5) outer-model exponent, and Prop 10.2 reduction exponents.
Only uses Fractions / sympy. Does not touch repository code."""
from fractions import Fraction as Fr
import random
import sympy as sp

random.seed(74)

def D(q, e, eta):
    # (1/eta) * int_{1/2-eta}^{1/2} min(q, r e) dr, exact
    lo, hi = Fr(1, 2) - eta, Fr(1, 2)
    if e == 0:
        return Fr(0)
    r0 = q / e
    if r0 >= hi:
        return e * (hi**2 - lo**2) / 2 / eta
    if r0 <= lo:
        return q
    return (e * (r0**2 - lo**2) / 2 + q * (hi - r0)) / eta

def Dq(q, e, eta):
    lo, hi = Fr(1, 2) - eta, Fr(1, 2)
    if e == 0:
        return Fr(0)
    r0 = q / e
    return max(Fr(0), min(Fr(1), (hi - max(r0, lo)) / eta))

def rnd(a, b, den=997):
    return a + (b - a) * Fr(random.randint(1, den - 1), den)

# ---------- symbolic identities ----------
h, q, tau, beta, p, Dv, w, eta = sp.symbols('h q tau beta p D w eta', real=True)
ell = 1 - tau
d = 2 - q
Bs = h + ell * (p * beta - Dv)
Ts = tau + Bs
bp = beta - q  # case q<beta
# case 1
U1 = tau / 2 * d + 2 * ell * bp
assert sp.simplify((Ts - U1) - (h + tau * q / 2 + ell * (p * beta - Dv - 2 * bp))) == 0
# case 2
U2 = tau * (1 - q) + (1 + ell) * bp
assert sp.simplify((Ts - U2) - (h + (ell * (p - 1) - 1) * beta - ell * Dv + 2 * q)) == 0
# case 3
U3 = tau / 2 * d + ell * beta
assert sp.simplify((Ts - U3) - (h + tau * q / 2 + ell * ((p - 1) * beta - Dv))) == 0
# case 4 (tau=ell=1/2)
sub = {tau: sp.Rational(1, 2)}
U4 = sp.Rational(1, 2) + sp.Rational(3, 2) * bp
gap4 = sp.expand((Ts - U4).subs(sub))
print('case4 gap (q<beta):', gap4)
# case 5 (tau=1/2, D=q)
U0 = sp.Rational(1, 2) + (beta - q) / 2
gap5 = sp.simplify((Ts - U0).subs(sub).subs(Dv, q))
print('case5 K0:', gap5)
assert sp.simplify(gap5 - (h + (p - 1) * beta / 2)) == 0
# (10.11) identity
a1, y1, u1, v1 = sp.symbols('a1 y1 u1 v1')
g = a1 + y1 + u1 + v1
M = tau * y1 + ell * u1 + v1
assert sp.simplify(M - tau * g / 2 - (-tau / 2 * (a1 - y1) + (ell - tau / 2) * u1 + (1 - tau / 2) * v1)) == 0
assert sp.simplify((ell - tau / 2) + (1 - tau / 2) - 2 * ell) == 0
# case5: d/2 - s0/2 == U0
assert sp.simplify(d / 2 - (1 - beta) / 2 - U0) == 0
print('symbolic identities OK')

# ---------- Lemma 5.1: (5.4) linearisation ----------
eps, pd, zd, hp, hz, Fk, IDq, a_, u_, Pi, tol = sp.symbols('eps pdot zdot hp hz Fk IDq a u Pi tol')
v = -hp
chi = a_ + u_
kappa = sp.symbols('kappa')
# log_N of (5.3) threshold <= log_N of (5.1) local bound ; b = kappa - u
hh, dd = sp.symbols('hh dd')
lhs = -tol - (Pi + eps * (pd * Fk + zd * IDq)) - (hh + eps * (hp * pd + hz * zd)) * kappa \
      + ((dd + eps * zd) - 1 - (hh + eps * (hp * pd + hz * zd))) * a_ + (1 + hh + eps * (hp * pd + hz * zd)) * (kappa - u_)
rhs = -Pi - hh * kappa + (dd - 1 - hh) * a_ + (1 + hh) * (kappa - u_)
expr = sp.expand(rhs - lhs)  # >= 0
claimed = tol + eps * (pd * (Fk - v * chi) + zd * (hz * chi + IDq - a_))
assert sp.simplify(expr - claimed) == 0
print('Lemma 5.1 linearisation (5.2) OK')
# (5.5): T*a/b exponent
T = -Pi - hh * kappa + (dd - 1 - hh) * a_ + (1 + hh) * (kappa - u_)
assert sp.simplify(T + a_ - (kappa - u_) - (-Pi - hh * chi + dd * a_)) == 0
print('(5.5) exponent OK')

# ---------- numeric sweep over admissible parameters (exact rationals) ----------
# For each case, maximise the stated upper bound of Q over rate vectors allowed by the
# constraints actually invoked, and compare to T*. Using LP vertices by brute grid.
worst = {}
N = 20000
for _ in range(N):
    pv = rnd(Fr(2), Fr(3))
    pmin = rnd(Fr(2), pv)
    etav = min(Fr(1, 32), (pmin - 2) / 16) * rnd(Fr(0), Fr(1))
    if etav == 0:
        continue
    qv = rnd(Fr(0), Fr(1))
    hv = rnd(Fr(0), 1 - qv)
    tv = rnd(Fr(0), Fr(1)) if random.random() < 0.8 else Fr(1, 2)
    lv = 1 - tv
    bv = rnd(Fr(0), Fr(1))
    Dval = D(qv, bv, etav)
    Dqv = Dq(qv, bv, etav)
    hzv = rnd(Fr(0), Fr(1))
    wv = hzv + lv * Dqv
    if not (0 < wv <= 1):
        continue
    dv = 2 - qv
    Bsv = hv + lv * (pv * bv - Dval)
    Tsv = tv + Bsv
    bplus = max(Fr(0), bv - qv)
    s0 = 1 - bv
    if tv != lv:
        if wv >= tv / 2:
            case = 1; U = tv / 2 * dv + 2 * lv * bplus
        elif tv < 2 * lv:
            if qv < (Fr(1, 2) - etav) * bv:
                raise SystemExit('case2 regime inconsistent: w<tau/2 but q small')
            case = 2; U = tv * (1 - qv) + (1 + lv) * bplus
        else:
            case = 3
            vstar = (2 * bv - qv) / 4 if 2 * bv > qv else Fr(0)  # from gbar<=d
            vstar = bv / 2  # paper uses v* <= beta/2
            U = tv / 2 * dv + 2 * lv * vstar
    else:
        if qv >= (Fr(1, 2) - etav) * bv:
            case = 4; U = Fr(1, 2) + Fr(3, 2) * bplus
        else:
            case = 5; U = Fr(1, 2) + (bv - qv) / 2
    gap = Tsv - U
    if gap <= 0:
        print('NEGATIVE GAP', case, dict(p=pv, eta=etav, q=qv, h=hv, tau=tv, beta=bv, w=wv), gap)
        raise SystemExit(1)
    rel = gap - hv  # gap beyond h
    worst.setdefault(case, rel)
    worst[case] = min(worst[case], rel)
print('numeric sweep: all gaps positive; min (gap - h) per case:', {k: float(v_) for k, v_ in sorted(worst.items())})

# ---------- case 4 coefficient: paper claims (p-2)/2 - 3eta/2 ; check bound with D<=q ----------
for _ in range(5000):
    pv = rnd(Fr(2), Fr(3)); etav = min(Fr(1, 32), (pv - 2) / 16) * rnd(Fr(0), Fr(1))
    bv = rnd(Fr(0), Fr(1)); qv = rnd((Fr(1, 2) - etav) * bv, bv)
    Dval = D(qv, bv, etav)
    gap_minus_h = (pv * bv - Dval - 3 * bv + 3 * qv) / 2
    assert gap_minus_h >= ((pv - 2) / 2 - 3 * etav / 2) * bv
    assert gap_minus_h >= ((pv - 2) / 2 - etav) * bv  # sharper
print('case 4 coefficient bound OK (paper constant is valid but not sharp)')

# ---------- Prop 10.2 exponent chain ----------
lam, Nn, xi, pp, nn, E = sp.symbols('lambda N xi p n E', positive=True)
# |A| >= N^-xi n N lam^(p+1) ; sigma <= C n N^-2 ; N^3|E| >= N^-xi n N lam^(p+1)
sigma_bound = sp.simplify((Nn**3 * E * Nn**xi / (Nn * lam**(pp + 1))) * Nn**-2)
print('sigma(S_lambda) bound:', sigma_bound)
assert sp.simplify(lam**3 * sigma_bound - Nn**xi * E * lam**(2 - pp)) == 0
print('lambda^3 sigma <= N^xi |E| lambda^(2-p) <= N^(p-2+xi)|E| for lambda>=1/N: OK')
