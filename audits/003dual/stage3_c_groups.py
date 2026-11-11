# -*- coding: utf-8 -*-
"""
Stage 3, layer D (floating point, one GROUP (h_0, active set) at a time): the closed dual form of
paper2.tex:3272-3282 with the scalar of :3420-3424 except its archimedean part.

For a fixed Psi_0, primes P with exponents j_p, a class h_0 with phi-hat(h_0) != 0 and an admissible active
set A ({p : j_p != 0} <= A <= P):

  LEFT   sum over h_p != 0 (p in A), h_p = 0 (p not in A) of  c_F(h) conj theta(z + z_h, v),
         c_F(h) = phi-hat(h_0) prod_p C_{p,j_p}(h_p)  with C_{p,j} from its definition (:2945),
         conj theta(z + z_h, v) from its own Fourier series (exact rational phases for z_h).
  RIGHT  C' * [ (3^{5/2}/2) 1_{sigma=0} v'^{2/3} prod_{p in A} B_{p,j_p}(0)
               + sum_l d_sigma(l) psi(lambda^4 l) prod_{p in A} B_{p,j_p}(lambda^4 l) v' K_{1/3}(4 pi |l| v') e~(l z'') ],
         C' = phi-hat(h_0) conj(kappa_0) prod_{p not in A}(1 - 1/N(p)) prod_{p in A} chi_p(sigma_p)^-2 omega_{p,j_p},
         z'' = -conj(z)/(c^2 (v^2 + |z|^2)),  v' = v/(N(c)(v^2 + |z|^2)),
         psi(x) = e(-delta'_0 r^-1 x / D_0), D_0 = lambda^3 c_0, c_0 = c/r,  kappa_0 from eq:ray-multiplier,
         sigma_p, eps_p, omega_{p,j} from :3209-3246, B_{p,j} from :1735-1742.
  So C of :3420 is  -(i/81) conj(alpha(c))^2 C'.  The factor -(i/81) conj(alpha(c))^2 and V^sharp are NOT
  tested here (stage3_d).  The constant mode is the l = 0 instance of the same local identity; the paper
  removes it by differentiating (:3363).  It matters only when every active prime has j = 4.
  The representatives h_p = M^2 k_p of :3039 are used to build a, c, delta', kappa_0 (k_p = 1); on the left
  the h_p run over small representatives of O/p (theta is 3O-periodic).
  NOT independent: chi_p, gamma_j on both sides come from check_local.Local; tau, d_sigma from kubota.py.

CRITERIA (fixed before the first run):
  D1  |LEFT - RIGHT| <= 1e-9 * scale,  scale = (sum_h |c_F(h)|) * (sum of absolute values of the terms of the
      series of conj theta at (z, v))   [ABSOLUTE difference against a sum of absolute values].
      Also printed: rel = |LEFT - RIGHT| / max(|LEFT|, |RIGHT|).  If rel > 1e-6 in any group although D1 holds,
      the group is listed (it would mean heavy cancellation, where D1 is weak).
      One random point per group: v = t/|c|, z = zeta/|c| with t in [0.9, 1.1], |Re zeta|, |Im zeta| <= 0.3 when
      the tables allow it (slack >= 1.27, slack = (complete range of the shorter table)/(T |c|/(4 pi))), otherwise
      t in [0.97, 1.03], |Re zeta|, |Im zeta| <= 0.1 (needs slack >= 1.05; else the script stops).
      [The first official run stopped at the first two-prime configuration with N(c) = 3276 because this
      threshold was 1.06 (miscomputed: the narrow rule needs 1.0494, the largest |c| = 57.24 has slack 1.0563).
      Its output, complete for the configurations with no prime or one prime, is kept as
      stage3_c_groups_run1.out.  Criteria unchanged.  This file is a full rerun with the same seed.]
      Truncation 4 pi |l| v <= T = 40 on both sides.  [CHANGED after the debug runs: T was 36 with the wide
      rule only; see the note in stage3_b_translate.py.  Tolerances unchanged.]
  COV printed table of what was exercised: (case of H, sigma, j, active/inactive).
  Negative controls; each is counted only on groups where it is live, and must give rel > 1e-3 on EVERY live group:
  ND1 kappa_0 in place of conj(kappa_0)                    live: kappa_0 != 1
  ND2 psi without r^-1                                     live: D_0 does not divide delta'(r^-1 - 1)
  ND3 chi_p^{+(j+2)} in B_{p,j}                            live: some active p with j in {2, 3, 5}
  ND4 chi_p(sigma_p)^-2 omitted                            live: prod chi_p(sigma_p)^2 != 1
  ND5 eps_p replaced by -eps_p in omega_{p,j}              live: prod over active p with j != 4 of chi_p(-1)^(j+2) = -1
  ND6 d_+ <-> d_-                                          live: sigma != 0
  ND7 c_0 replaced by c_0 * zeta6 (other unit) in psi and kappa_0, same c
      live: kappa_0 or the character psi actually changes.  [CHANGED after a debug run on the first 36
      configurations with the small table: the header said 'live: always'; in the 72 groups with c_0 = 1
      (psi trivial, (a/c_0)_3 = 1) the replacement changes nothing, so the control was dead there.  D1 and the
      other controls were not touched.]
  ND8 (1 - 1/N) for an inactive prime replaced by 1        live: some inactive prime
Configurations: see CONF below (written before running).
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage3_c_groups.py > audits/003dual/stage3_c_groups.out 2>/dev/null
"""
import sys, time, random, itertools
from scipy.special import kv
from stage2_lib import *

t0 = time.time()
print(__doc__)
RES = []
def report(name, ok, detail):
    RES.append((name, ok)); print("  %-4s %s  %s" % (name, "PASS" if ok else "FAIL", detail)); sys.stdout.flush()
random.seed(2026111204)
TOL, NEGT, T = 1e-9, 1e-3, 40.0
LINT = int(os.environ.get("STAGE_L", "18")); setup = Setup((LINT, 0))
P7, P13, M5, P19 = (1, 3), (4, 3), (-5, 0), cl.split_prime(19)
CONF = [(nm, {}) for nm in "ABCDEF"] \
     + [(nm, {P7: j}) for nm in "ABCEF" for j in range(6)] \
     + [(nm, jd) for nm in "ABF" for jd in ({P13: 1}, {P13: 4}, {P13: 0}, {M5: 1}, {M5: 2}, {P19: 3}, {P19: 5})] \
     + [(nm, jd) for nm in "ABC" for jd in ({P7: 1, P13: 5}, {P7: 0, P13: 4}, {P7: 2, P13: 0}, {P7: 4, P13: 4}, {P7: 1, P13: 1})]
if len(sys.argv) > 1: CONF = CONF[:int(sys.argv[1])]
if LINT != 18:
    # L-VARIATION RUN (added after the debug runs, before any official run of this layer): the same criteria with
    # L = 36 or 54 ("sufficiently high powers of the primes above 2 and 3", :2930), on a reduced list.
    CONF = [(nm, jd) for nm in "ABF" for jd in ({}, {P7: 1}, {P7: 0}, {P7: 4})]
print("  L = %d ; %d configurations (Psi_0, twists)" % (LINT, len(CONF)))
th = Theta(int(os.environ.get("STAGE_NMAX", "1000000")))
print("  theta tables: N(A) <= %d ; complete to |l| <= %.2f (theta), %.2f (cusp series)  [%.0fs]" % (th.Nmax, th.complete["0"], th.complete["+"], time.time() - t0))
PH = {nm: phi_hat(phi_table(nm, LINT), LINT) for nm in "ABCDEF"}
d0 = th.d["0"]

def left(nm, jd, h0, act, z, v):
    primes = list(jd); r = prod(primes); Cden = mul((LINT, 0), r)
    locs = [loc(p) for p in primes]
    n = d0.cut(v, T)
    assert T / (4 * math.pi * v) <= th.complete["0"] * (1 + 1e-12), "theta table too short on the left"
    base = d0.coef[:n] * v * kv(1.0 / 3, 4 * math.pi * d0.aell[:n] * v) * np.exp(4j * math.pi * (d0.ell[:n] * z).real)
    m = d0.m[:n]; F = np.zeros(n, dtype=complex); sF = 0; aF = 0.0; ntr = 0
    rng = [range(1, lp.N) if p in act else [0] for lp, p in zip(locs, primes)]
    Cs = [lp.C(jd[p]) for lp, p in zip(locs, primes)]
    for idxs in itertools.product(*rng):
        cF = PH[nm][h0]
        for Cc, i in zip(Cs, idxs): cF *= Cc[i]
        hs = [lp.L.reps[i] for lp, i in zip(locs, idxs)]
        A = mul(LAM2, add(mul(h0, r), mul((LINT, 0), esum([mul(h, exact_div(r, p)) for h, p in zip(hs, primes)]))))
        F += cF * AddChar(A, mul(LAM3, Cden))(m); sF += cF; aF += abs(cF); ntr += 1
    val_ = sF * SIGMA0 * v ** (2.0 / 3) + complex(np.sum(base * F))
    scale = aF * (SIGMA0 * v ** (2.0 / 3) + float(np.sum(np.abs(base))))
    return val_, scale, ntr

def rel(x, y): return abs(x - y) / max(abs(x), abs(y), 1e-300)
ctl = {k: [0, 0, 1e9] for k in ("ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "ND7", "ND8")}
def control(name, live, lhs, rhs):
    if not live: return
    e = rel(lhs, rhs); ctl[name][0] += 1; ctl[name][1] += e > NEGT; ctl[name][2] = min(ctl[name][2], e)
nnarrow = 0
cov = {}; nfail = 0; ngr = 0; worst = 0.0; worstrel = 0.0; flagged = []
print("\n  Psi_0  twists {N(p): j}        #groups  #translates/group   worst diff/scale   worst rel   N(c) values")
for nm, jd in CONF:
    primes = list(jd); must = [p for p in primes if jd[p] != 0]; opt = [p for p in primes if jd[p] == 0]
    w1 = w2 = 0.0; g1 = 0; ncs = set(); ntrs = set()
    for h0, f0 in PH[nm].items():
        if abs(f0) < 1e-13: continue
        for kopt in range(len(opt) + 1):
            for extra in itertools.combinations(opt, kopt):
                act = [p for p in primes if p in must or p in extra]; inact = [p for p in primes if p not in act]
                cfg = setup.build(h0, act, [ONE] * len(act))
                ca = abs(cxf(cfg["c"])); ncs.add(norm(cfg["c"]))
                slack = min(th.complete["0"], th.complete[cfg["sigma"]]) / (T / (4 * math.pi) * ca)
                if slack >= 1.27: tr_, zr_ = 0.1, 0.3
                elif slack >= 1.05: tr_, zr_ = 0.03, 0.1; nnarrow += 1
                else: raise RuntimeError("tables too short for |c| = %.2f" % ca)
                v = random.uniform(1 - tr_, 1 + tr_) / ca; z = complex(random.uniform(-zr_, zr_), random.uniform(-zr_, zr_)) / ca
                lhs, scale, ntr = left(nm, jd, h0, act, z, v); ntrs.add(ntr)
                z2, v2 = cusp_point(cfg["c"], z, v)
                Cp = Cprime(cfg, f0, jd, inact); dv = dual_value(th, cfg, jd, z2, v2, T=T)
                rhs = Cp * dv
                e = abs(lhs - rhs) / scale; er = rel(lhs, rhs); ngr += 1; g1 += 1
                w1 = max(w1, e); w2 = max(w2, er)
                if e > TOL: nfail += 1
                if er > 1e-6: flagged.append((nm, {norm(p): j for p, j in jd.items()}, h0, [norm(p) for p in act], e, er, abs(lhs)))
                for p in primes:
                    key = (cfg["case"], cfg["sigma"], jd[p], "active" if p in act else "inactive"); cov[key] = cov.get(key, 0) + 1
                if not primes: cov[(cfg["case"], cfg["sigma"], "-", "no prime")] = cov.get((cfg["case"], cfg["sigma"], "-", "no prime"), 0) + 1
                # ---- controls
                _, k0 = kappa_formula(cfg)
                control("ND1", k0 != 0, lhs, Cp * W3 ** ((2 * k0) % 3) * dv)
                live2 = not divides(cfg["D0"], mul(cfg["dl"], sub(cfg["rinv"], ONE)))
                if live2: control("ND2", True, lhs, Cp * dual_value(th, cfg, jd, z2, v2, T=T, variant="psi_no_rinv"))
                if any(jd[p] in (2, 3, 5) for p in act): control("ND3", True, lhs, Cp * dual_value(th, cfg, jd, z2, v2, T=T, variant="Bsign"))
                ks = sum(2 * loc(p).k(s) for p, s in zip(cfg["primes"], cfg["sig"])) % 6
                control("ND4", ks != 0, lhs, Cp * Z6 ** ks * dv)
                sgn = 1
                for p in act:
                    if jd[p] != 4 and loc(p).L.km1 == 3 and (jd[p] + 2) % 2 == 1: sgn = -sgn
                control("ND5", sgn == -1, lhs, -rhs)
                if cfg["sigma"] != "0": control("ND6", True, lhs, Cp * dual_value(th, cfg, jd, z2, v2, T=T, variant="swap"))
                cfg7 = setup.build(h0, act, [ONE] * len(act), c0_unit=UNITS[1])
                pa, pb = psi_char(cfg), psi_char(cfg7); bas = np.array([[1, 0], [0, 1]])
                live7 = kappa_formula(cfg7)[1] != k0 or float(np.max(np.abs(pa(bas) - pb(bas)))) > 1e-9
                control("ND7", live7, lhs, Cprime(cfg7, f0, jd, inact) * dual_value(th, cfg7, jd, z2, v2, T=T))
                if inact:
                    f = 1.0
                    for p in inact: f *= 1 - 1.0 / norm(p)
                    control("ND8", True, lhs, rhs / f)
    worst = max(worst, w1); worstrel = max(worstrel, w2)
    print("  %s      %-24s %5d    %-18s  %.2e           %.2e    %s   [%.0fs]" % (nm, {norm(p): j for p, j in jd.items()}, g1, sorted(ntrs), w1, w2, sorted(ncs), time.time() - t0)); sys.stdout.flush()

print("\n  coverage (case of H, sigma, j, role): number of (group, prime) incidences")
for key in sorted(cov, key=str): print("    %-8s sigma=%s  j=%s  %-9s : %d" % (key + (cov[key],)))
print()
report("D1", nfail == 0, "%d failures / %d groups ; worst |L-R|/scale = %.2e (tolerance %.0e) ; worst |L-R|/max(|L|,|R|) = %.2e ; narrow point rule in %d groups" % (nfail, ngr, worst, TOL, worstrel, nnarrow))
if flagged:
    print("       groups with rel > 1e-6 (heavy cancellation, D1 is weak there): %d" % len(flagged))
    for f in flagged[:20]: print("         ", f)
else:
    print("       no group has rel > 1e-6")
for k in sorted(ctl):
    report(k, ctl[k][0] > 0 and ctl[k][1] == ctl[k][0], "live in %d groups, rejected (rel > %.0e) in %d, smallest rel %.2e" % (ctl[k][0], NEGT, ctl[k][1], ctl[k][2]))
print("\nLAYER D: %s   [%.0fs]" % ("all criteria met" if all(ok for _, ok in RES) else "SOME CRITERIA NOT MET: " + ",".join(n for n, ok in RES if not ok), time.time() - t0))
