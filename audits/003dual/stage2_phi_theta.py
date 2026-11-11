# -*- coding: utf-8 -*-
"""
Stage 2: phi exactly as in 003 (paper2.tex:2909-2913, with the condition (n, S) = 1 at the prime 2), L (:2928-2933),
phi-hat (:2936), and the finite Fourier claim (:2964-2979)
    Theta_Psi(z, v) = sum_h c_F(h) conj theta(z + z_h, v),   z_h = lambda^2 (h_0/L + sum_p h_p/p),
    c_F(h) = phi-hat(h_0) prod_p C_{p,j_p}(h_p).

Choices (made before running).  S = {(2), (lambda)}.  L = 18.  Psi_0 in {A, B, C, D, E} (see stage2_lib.psi0_exp):
A trivial; B = conj (lambda/n)_3; C = B * omega^((n-1)/3 mod lambda); D = cubic character mod 2; E = C * D.
B and C are not in the paper; they are admissible only if they are ray class characters with conductor supported
on S -- P0 and P1 below check complete multiplicativity on primary elements and periodicity mod 18, nothing more.
They were chosen so that reduced denominators with (c, lambda) = 1 and v_lambda(c) = 1 occur with phi-hat != 0
(with A alone only 3 | c is expected); whether that happens is PRINTED in table P2 before any theta computation.

CRITERIA (fixed before the first run):
  P0  Psi_0(n1 n2) = Psi_0(n1) Psi_0(n2) for all primary n1, n2 prime to 2 in the box [-7, 7]^2.  0 violations.
  P1  phi(n) = phi(n mod 18) for every n in the box [-45, 45]^2 (phi = 0 off the primary elements prime to 2).
      0 violations.  (Informative only: whether smaller moduli 6, 9, 3*lambda*2 would do.)
  P2  sum_h phi-hat(h) = phi(0) = 0 (:2978) and Parseval, to 1e-12.  Table of the support of phi-hat with
      c_0, case of H, sigma, kappa_0 (active set empty).  [information for the later stages]
  S2a Theta_Psi (direct series, :2917-2921) against sum_h c_F(h) conj theta(z + z_h, v), where conj theta at
      the translated point is evaluated from its own Fourier series (exact rational phases for z_h).
      |lhs - rhs| <= 1e-9 * scale,  scale = sum_h |c_F(h)| |theta(z + z_h, v)|   (ABSOLUTE difference
      against the sum of absolute values of the translate terms; the difference relative to |lhs| is printed too).
      C_{p,j}(h) is taken from its definition (:2945), not from eq:ray-fourier.
  S2b the same right side with theta(z + z_h) evaluated by kubota.translate (cusp expansions and Kubota
      character of audits/029reflect, no formula of 003): only for the empty set of primes.  Same tolerance.
  Points: 2 random points per configuration, v in [0.05, 0.10], |Re z|, |Im z| <= 0.3 (seed fixed);
          for S2b v in [0.20, 0.35].
  Configurations: every Psi_0 with P = {} ; Psi_0 = A, B, C, E with P = {p7: j} for j = 0..5 ;
          A, C with {p13: 2}, {-5: 1}, {p7: 1, p13: 5}, {p7: 0, p13: 4}.
  Negative controls (each must give |lhs - rhs| > 1e-4 * scale, at every point where it is run):
  N2a the translate side built from phi WITHOUT the condition at 2 (the earlier phi_W), Psi_0 = A, B, C, P = {}
      and P = {p7: 1}
  N2b shifts -z_h instead of +z_h, same configurations
  The truncation is 4 pi |l| v <= 40 on both sides; a rigorous bound for the neglected tail is in stage4_bounds.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage2_phi_theta.py > audits/003dual/stage2_phi_theta.out 2>/dev/null
"""
import sys, time, random, itertools
from stage2_lib import *

t0 = time.time()
print(__doc__)
RES = []
def report(name, ok, detail):
    RES.append((name, ok)); print("  %-4s %s  %s" % (name, "PASS" if ok else "FAIL", detail)); sys.stdout.flush()
random.seed(20261112)
LINT = 18; setup = Setup((LINT, 0)); NAMES = ["A", "B", "C", "D", "E"]
P7, P13, M5 = (1, 3), (4, 3), (-5, 0)

# ---- P0
box = [(a, b) for a in range(-7, 8) for b in range(-7, 8) if is_primary((a, b)) and not (a % 2 == 0 and b % 2 == 0)]
for nm in NAMES:
    bad = sum(1 for n1 in box for n2 in box if psi0_exp(nm, mul(n1, n2)) != (psi0_exp(nm, n1) + psi0_exp(nm, n2)) % 3)
    report("P0", bad == 0, "Psi_0 = %s: %d violations / %d pairs ; values taken: %s" % (nm, bad, len(box) ** 2, sorted(set(psi0_exp(nm, n) for n in box))))

# ---- P1
TAB = {nm: phi_table(nm, LINT) for nm in NAMES}
TABno2 = {nm: phi_table(nm, LINT, with2=False) for nm in "ABC"}
for nm in NAMES:
    bad = tot = 0; smaller = {6: True, 9: True}
    seen = {6: {}, 9: {}}
    for a in range(-45, 46):
        for b in range(-45, 46):
            k = phi_exp(nm, (a, b)); v = 0 if k is None else W3 ** k
            tot += 1
            if abs(v - TAB[nm][a % LINT, b % LINT]) > 1e-12: bad += 1
            for q in (6, 9):
                key = (a % q, b % q)
                if seen[q].setdefault(key, k) != k: smaller[q] = False
    report("P1", bad == 0, "Psi_0 = %s: %d violations / %d ; nonzero classes mod 18: %d ; (info) periodic mod 6: %s, mod 9: %s"
           % (nm, bad, tot, int(np.sum(np.abs(TAB[nm]) > 0)), smaller[6], smaller[9]))

# ---- P2
PH = {nm: phi_hat(TAB[nm], LINT) for nm in NAMES}
PHno2 = {nm: phi_hat(TABno2[nm], LINT) for nm in "ABC"}
for nm in NAMES:
    s = sum(PH[nm].values()); pars = sum(abs(x) ** 2 for x in PH[nm].values()) - np.sum(np.abs(TAB[nm]) ** 2) / LINT ** 2
    report("P2", abs(s) < 1e-12 and abs(pars) < 1e-12, "Psi_0 = %s: |sum phi-hat| = %.1e, Parseval defect = %.1e" % (nm, abs(s), abs(pars)))
    rows = {}
    for h0, val_ in PH[nm].items():
        if abs(val_) < 1e-13: continue
        cfg = setup.build(h0, [], [])
        key = (norm(cfg["c0"]), cfg["case"], cfg["sigma"], kappa_formula(cfg)[1])
        d = rows.setdefault(key, [0, 0.0, set()]); d[0] += 1; d[1] += abs(val_) ** 2; d[2].add(round(abs(val_), 6))
    print("       support of phi-hat: %d of %d classes h_0.   N(c_0)  case     sigma  kappa_0-exp   #h_0   sum|phi-hat|^2   |phi-hat| values" % (sum(d[0] for d in rows.values()), LINT ** 2))
    for key in sorted(rows):
        print("                                              %5d   %-8s  %s      %d          %3d    %.6f      %s" % (key + (rows[key][0], rows[key][1], sorted(rows[key][2]))))
sys.stdout.flush()

# ---- theta
th = Theta(200000)
print("\n  theta tables: N(A) <= %d ; complete to |l| <= %.2f (theta), %.2f (cusp series)  [%.0fs]" % (th.Nmax, th.complete["0"], th.complete["+"], time.time() - t0))
d0 = th.d["0"]
LAMBAR = conj(LAM)
def n_of_m(m):
    """n = m / lambda (= lambda^3 l); returns (n array, mask of integrality)"""
    x0 = m[:, 0] * LAMBAR[0] - m[:, 1] * LAMBAR[1]; x1 = m[:, 0] * LAMBAR[1] + m[:, 1] * LAMBAR[0] - m[:, 1] * LAMBAR[1]
    ok = (x0 % 3 == 0) & (x1 % 3 == 0)
    return np.stack([x0 // 3, x1 // 3], axis=1), ok
NN, NOK = n_of_m(d0.m)
def twist_series(nm, jd, tab=None):
    tab = TAB[nm] if tab is None else tab
    f = np.where(NOK, tab[NN[:, 0] % LINT, NN[:, 1] % LINT], 0.0)
    for p, j in jd.items():
        lp = loc(p); f = f * lp.chitab(j)[lp.idx_arr(NN)]
    return d0.with_coef(f, const=0.0)

def translates(nm, jd, ph=None):
    """list of (c_F, A, Cden): z_h = A / Cden exactly"""
    ph = PH[nm] if ph is None else ph
    primes = list(jd); r = prod(primes); Cden = mul((LINT, 0), r); out = []
    locs = [loc(p) for p in primes]
    Cs = [lp.C(jd[p]) for lp, p in zip(locs, primes)]
    for h0, f0 in ph.items():
        if abs(f0) < 1e-13: continue
        for idxs in itertools.product(*[range(lp.N) for lp in locs]):
            cF = f0
            for Cc, i in zip(Cs, idxs): cF *= Cc[i]
            if abs(cF) < 1e-15: continue
            hs = [lp.L.reps[i] for lp, i in zip(locs, idxs)]
            A = mul(LAM2, add(mul(h0, r), mul((LINT, 0), esum([mul(h, exact_div(r, p)) for h, p in zip(hs, primes)]))))
            out.append((cF, A, Cden))
    return out

def rhs_series(tr, z, v, sign=1):
    tot = 0; sc = 0.0
    for cF, A, Cden in tr:
        val_ = th.thetabar_at(z, v, shift=(A if sign == 1 else neg(A), Cden))
        tot += cF * val_; sc += abs(cF) * abs(val_)
    return tot, sc

def rhs_kubota(tr, z, v):
    tot = 0
    for cF, A, Cden in tr:
        a, cc = kb.reduce_frac(A, Cden)
        tot += cF * np.conj(th.K.translate(a, cc, z, v))
    return tot

TOL, NEG = 1e-9, 1e-4
worst = {"S2a": 0.0, "S2b": 0.0}; nfail = {"S2a": 0, "S2b": 0}; cnt = {"S2a": 0, "S2b": 0}
negmin = {"N2a": 1e9, "N2b": 1e9}; negcnt = {"N2a": 0, "N2b": 0}
CONF = [(nm, {}) for nm in NAMES] + [(nm, {P7: j}) for nm in "ABCE" for j in range(6)] \
     + [(nm, jd) for nm in "AC" for jd in ({P13: 2}, {M5: 1}, {P7: 1, P13: 5}, {P7: 0, P13: 4})]
print("\n  Psi_0  twists                 #translates  point (z, v)                      |lhs|        |lhs-rhs|   scale      diff/scale  diff/|lhs|")
for nm, jd in CONF:
    S = twist_series(nm, jd); tr = translates(nm, jd)
    for it in range(2):
        v = random.uniform(0.05, 0.10); z = complex(random.uniform(-.3, .3), random.uniform(-.3, .3))
        lhs = S(z, v, complete_to=th.complete["0"]); rhs, sc = rhs_series(tr, z, v)
        e = abs(lhs - rhs) / sc; cnt["S2a"] += 1; worst["S2a"] = max(worst["S2a"], e)
        if e > TOL: nfail["S2a"] += 1
        print("  %s      %-22s %6d       (%+.3f%+.3fi, %.4f)   %.4e   %.2e   %.3e   %.1e     %.1e  %s"
              % (nm, {norm(p): j for p, j in jd.items()}, len(tr), z.real, z.imag, v, abs(lhs), abs(lhs - rhs), sc, e, abs(lhs - rhs) / abs(lhs), "ok" if e <= TOL else "FAIL"))
        if nm in "ABC" and (jd == {} or jd == {P7: 1}):
            tr2 = translates(nm, jd, PHno2[nm]); r2, _ = rhs_series(tr2, z, v)
            e2 = abs(lhs - r2) / sc; negmin["N2a"] = min(negmin["N2a"], e2); negcnt["N2a"] += 1
            r3, _ = rhs_series(tr, z, v, sign=-1)
            e3 = abs(lhs - r3) / sc; negmin["N2b"] = min(negmin["N2b"], e3); negcnt["N2b"] += 1
            print("         controls at this point: N2a (phi without the condition at 2, %d translates) diff/scale = %.2e ; N2b (shift -z_h) diff/scale = %.2e" % (len(tr2), e2, e3))
    if jd == {}:
        for it in range(2):
            v = random.uniform(0.20, 0.35); z = complex(random.uniform(-.3, .3), random.uniform(-.3, .3))
            lhs = S(z, v, complete_to=th.complete["0"]); _, sc = rhs_series(tr, z, v); rk = rhs_kubota(tr, z, v)
            e = abs(lhs - rk) / sc; cnt["S2b"] += 1; worst["S2b"] = max(worst["S2b"], e)
            if e > TOL: nfail["S2b"] += 1
            print("  %s      via kubota.translate   %6d       (%+.3f%+.3fi, %.4f)   %.4e   %.2e   %.3e   %.1e     %.1e  %s"
                  % (nm, len(tr), z.real, z.imag, v, abs(lhs), abs(lhs - rk), sc, e, abs(lhs - rk) / abs(lhs), "ok" if e <= TOL else "FAIL"))
    sys.stdout.flush()
print()
report("S2a", nfail["S2a"] == 0, "%d failures / %d points ; worst |lhs-rhs|/scale = %.2e (tolerance %.0e)" % (nfail["S2a"], cnt["S2a"], worst["S2a"], TOL))
report("S2b", nfail["S2b"] == 0, "%d failures / %d points ; worst |lhs-rhs|/scale = %.2e" % (nfail["S2b"], cnt["S2b"], worst["S2b"]))
report("N2a", negcnt["N2a"] > 0 and negmin["N2a"] > NEG, "control without the condition at 2: smallest diff/scale over %d points = %.2e (must be > %.0e)" % (negcnt["N2a"], negmin["N2a"], NEG))
report("N2b", negcnt["N2b"] > 0 and negmin["N2b"] > NEG, "control with shifts -z_h: smallest diff/scale over %d points = %.2e" % (negcnt["N2b"], negmin["N2b"]))
print("\nSTAGE 2: %s   [%.0fs]" % ("all criteria met" if all(ok for _, ok in RES) else "SOME CRITERIA NOT MET: " + ",".join(n for n, ok in RES if not ok), time.time() - t0))
