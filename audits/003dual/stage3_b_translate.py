# -*- coding: utf-8 -*-
"""
Stage 3, layer B (floating point, ONE translate at a time): eq:theta-cusp-automorphy (paper2.tex:3161-3170)
with kappa(g_1) taken from the right side of eq:ray-multiplier (:3175-3183) and the cusp expansion (:3134-3139).

  conj theta(z + a/c, v)  ==  conj kappa(g_1) * [ (3^{5/2}/2) 1_{sigma=0} v'^{2/3} + sum_l d_sigma(l) v' K_{1/3}(4 pi |l| v') e~(l z') ],
  (z', v') = g^-1(z + a/c, v) = (-delta'/c - conj(z)/(c^2 (v^2+|z|^2)),  v/(N(c)(v^2+|z|^2))).

Left: the Fourier series of conj theta (DR (5.7) as transcribed in kubota.py) at the translated point; the
rational shift a/c enters through exact phases.  Right: kappa from stage2_lib.kappa_formula (the paper's three
cases), sigma from H by the rule  [3 | c: 0;  v_lambda(c) = 1: u_0 = lambda -> '-', u_0 = -lambda -> '+';
(c, lambda) = 1: u_0 = x + y omega with y = 0, 1, 2 mod 3 -> 0, '-', '+'], d_+ = conj-series of
theta(gamma_10 w), d_- of theta(gamma_19 w) (kubota.py; gamma_+ = [[1,0],[omega,1]], :1773).  The phase
e~(-delta' l/c) is exact.  The rule H -> sigma is my reading of :3090-3101; it is tested separately in B0.

CRITERIA (fixed before the first run):
  B0  theta(H w) == theta(gamma_sigma w) for H = [[1,0],[u_0,1]], u_0 = +-lambda, and H = [[u_0,-1],[1,0]],
      u_0 in the 9 representatives {x + y omega : 0 <= x, y <= 2}; 2 random points each, relative difference
      <= 1e-9.  Control: with '+' and '-' exchanged the relative difference must exceed 1e-4 whenever sigma != 0.
  B1  [CRITERION CHANGED AFTER THE FIRST OFFICIAL RUN -- see "HISTORY" below.]
      v2 (in force):  |lhs - rhs| <= 1e-9 * S,  S = sum of the absolute values of the terms of the series of
      conj theta at the left point (constant mode included).  ABSOLUTE difference against a sum of absolute values.
      Reported besides: cond = |lhs| / S for every configuration; the number of configurations with cond < 1e-3
      ("ill-conditioned": theta is small there and the check says little about phases), the worst relative
      difference |lhs - rhs| / max(|lhs|, |rhs|) among the well-conditioned ones (cond >= 1e-3) and among all.
      v1 (as first registered, NOT met):  |lhs - rhs| / max(|lhs|, |rhs|) <= 1e-9.
      One random point per configuration: v = t/|c|, z = zeta/|c| with t in [0.9, 1.1], |Re zeta|, |Im zeta| <= 0.3
      when the tables allow it (slack >= 1.27, slack = (complete range of the shorter table)/(T |c|/(4 pi))),
      otherwise t in [0.97, 1.03], |Re zeta|, |Im zeta| <= 0.1 (needs slack >= 1.05; else the script stops).
      Configurations: L = 18, EVERY h_0 in O/(18) (324 classes, including (h_0, L) != 1 and phi-hat = 0),
      r in {1, p7, p13, -5, p7*p13}, one tuple (k_p) per (h_0, r) chosen at random from 4 (3) tuples.
      Truncation 4 pi |l| v <= T = 40 on both sides (bound for the tail: stage4_bounds).
  HISTORY (the outputs of the first official run are kept as stage3_b_translate_run1.out / _run1.err.txt):
      (a) Debug run on the small table (r in {1, p7}): T was 36; worst relative difference 4.3e-10.  I attributed
          it to the truncation and raised T to 40.  That diagnosis was WRONG (see (b)); T = 40 is kept.
      (b) First official run (criterion v1): worst relative difference 3.3e-10 (r = p7), 6.9e-10 (r = p13),
          7.0e-4 (r = -5): v1 NOT met at r = -5.  Cause, found by printing the worst configurations: at those
          points |theta| itself is tiny (e.g. 4.3e-8 with S about 30) and |lhs - rhs| is about 6e-15, i.e.
          rounding of a sum of terms of size S.  A criterion relative to the value was the defect (the same
          defect that the Codex review found in audits/029reflect).  v2 measures against S, as layers D, E0 and
          stage 2 did from the start.  This is a change of criterion AFTER seeing data.
      (c) The same run then stopped at r = p7*p13 (|c| = 57.24): my threshold 'slack >= 1.06' was miscomputed
          (the narrow rule needs slack >= 1.0494; the largest |c| has slack 1.0563).  Threshold now 1.05.  The
          series evaluator itself refuses to run beyond the complete range of a table, independently of this.
  COV every (case, sigma) that the construction can produce occurs, with kappa != 1 somewhere in each case.
  Negative controls (a control is counted only on configurations where it is live AND well-conditioned,
  cond >= 1e-3 [this restriction belongs to the change v1 -> v2]; the thresholds, relative difference > 1e-4,
  are unchanged):
  NB1 kappa instead of conj kappa: live where kappa != 1; relative difference must be > 1e-4 on all of them.
  NB2 d_+ <-> d_-: live where sigma != 0; must be > 1e-4 on all of them.
  NB3 delta' replaced by delta' + 1 (another solution of nothing): live everywhere; > 1e-4 on at least 95%
      (the phase e~(-l/c) could be nearly invisible at a few points; the fraction is reported).
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage3_b_translate.py > audits/003dual/stage3_b_translate.out 2>/dev/null
"""
import sys, time, random
from stage2_lib import *

t0 = time.time()
print(__doc__)
RES = []
def report(name, ok, detail):
    RES.append((name, ok)); print("  %-4s %s  %s" % (name, "PASS" if ok else "FAIL", detail)); sys.stdout.flush()
random.seed(2026111203)
TOL, NEG, T = 1e-9, 1e-4, 40.0
th = Theta(int(os.environ.get("STAGE_NMAX", "1000000")))
print("  theta tables: N(A) <= %d ; complete to |l| <= %.2f (theta), %.2f (cusp series)  [%.0fs]" % (th.Nmax, th.complete["0"], th.complete["+"], time.time() - t0))
def rel(x, y): return abs(x - y) / max(abs(x), abs(y), 1e-300)
def dser(sg, z, v, phase=None):
    return th.d[sg](z, v, T=T, complete_to=th.complete[sg], phase=phase)

# ---- B0
worst = 0.0; negmin = 1e9; nneg = 0; n0 = 0
Hs = [("v=1", u0, "-" if u0 == LAM else "+") for u0 in (LAM, neg(LAM))] + \
     [("coprime", (x, y), {0: "0", 1: "-", 2: "+"}[y]) for x in range(3) for y in range(3)]
for case, u0, sg in Hs:
    Hm = ((1, 0), (cxf(u0), 1)) if case == "v=1" else ((cxf(u0), -1), (1, 0))
    for it in range(2):
        z = complex(random.uniform(-.4, .4), random.uniform(-.4, .4)); v = random.uniform(0.7, 1.0)
        zz, vv = kb.mobius(Hm, z, v)
        lhs = np.conj(th.d["0"](zz, vv, T=T, complete_to=th.complete["0"]))        # theta(H w)
        rhs = np.conj(dser(sg, z, v))                                                # theta(gamma_sigma w)
        e = rel(lhs, rhs); worst = max(worst, e); n0 += 1
        if sg != "0":
            e2 = rel(lhs, np.conj(dser({"+": "-", "-": "+"}[sg], z, v))); negmin = min(negmin, e2); nneg += 1
    print("    H: %-8s u_0 = %-8s -> sigma = %s" % (case, u0, sg))
report("B0", worst <= TOL and negmin > NEG, "worst relative difference %.2e over %d points ; control (+ <-> -): smallest %.2e over %d points" % (worst, n0, negmin, nneg))

# ---- B1
P7, P13, M5 = (1, 3), (4, 3), (-5, 0)
RS = [[], [P7], [P13], [M5], [P7, P13]]
if os.environ.get('STAGE_DEBUG'): RS = RS[:2]
KT = {0: [[]], 1: [[(1, 0)], [(2, 0)], [(1, 1)], [(3, -2)]], 2: [[(1, 0), (1, 0)], [(2, 0), (1, 1)], [(0, 1), (3, 0)]]}
setup = Setup((18, 0))
nnarrow = 0; nill = 0; ill = []; worstrel_all = 0.0; worstrel_well = 0.0
stat = {}; nb = {"NB1": [0, 0, 1e9], "NB2": [0, 0, 1e9], "NB3": [0, 0, 1e9]}; nfail = 0; ntot = 0; worst1 = 0.0
def right(cfg, z, v, kexp, sg=None, dshift=ZERO):
    z2, v2 = cusp_point(cfg["c"], z, v)
    ph = AddChar(neg(add(cfg["dl"], dshift)), mul(LAM3, cfg["c"]))
    return W3 ** ((-kexp) % 3) * dser(cfg["sigma"] if sg is None else sg, z2, v2, phase=ph)
for primes in RS:
    for x0 in range(18):
        for x1 in range(18):
            h0 = (x0, x1); ks = random.choice(KT[len(primes)])
            cfg = setup.build(h0, primes, ks)
            kf, _ = kappa_formula(cfg); ca = abs(cxf(cfg["c"]))
            slack = min(th.complete["0"], th.complete[cfg["sigma"]]) / (T / (4 * math.pi) * ca)
            if slack >= 1.27: tr_, zr_ = 0.1, 0.3
            elif slack >= 1.05: tr_, zr_ = 0.03, 0.1; nnarrow += 1
            else: raise RuntimeError("tables too short for |c| = %.2f" % ca)
            v = random.uniform(1 - tr_, 1 + tr_) / ca; z = complex(random.uniform(-zr_, zr_), random.uniform(-zr_, zr_)) / ca
            lhs = th.thetabar_at(z, v, shift=(cfg["a"], cfg["c"]), T=T)
            rhs = right(cfg, z, v, kf)
            S = th.d["0"](z, v, T=T, complete_to=th.complete["0"], absum=True)
            ea = abs(lhs - rhs) / S; cond = abs(lhs) / S; e = rel(lhs, rhs); ntot += 1; worst1 = max(worst1, ea)
            worstrel_all = max(worstrel_all, e)
            if ea > TOL: nfail += 1
            key = (cfg["case"], cfg["sigma"], kf); d = stat.setdefault(key, [0, 0.0, 0.0, 0, 1e9]); d[0] += 1; d[1] = max(d[1], ea); d[4] = min(d[4], cond)
            if cond < 1e-3:
                nill += 1; d[3] += 1
                ill.append((cond, e, abs(lhs - rhs), str([norm(p) for p in primes]), str(h0), cfg["case"], cfg["sigma"], kf))
                continue
            worstrel_well = max(worstrel_well, e); d[2] = max(d[2], e)
            if kf != 0:
                e2 = rel(lhs, right(cfg, z, v, -kf)); nb["NB1"][0] += 1; nb["NB1"][1] += e2 > NEG; nb["NB1"][2] = min(nb["NB1"][2], e2)
            if cfg["sigma"] != "0":
                e2 = rel(lhs, right(cfg, z, v, kf, sg={"+": "-", "-": "+"}[cfg["sigma"]])); nb["NB2"][0] += 1; nb["NB2"][1] += e2 > NEG; nb["NB2"][2] = min(nb["NB2"][2], e2)
            e2 = rel(lhs, right(cfg, z, v, kf, dshift=ONE)); nb["NB3"][0] += 1; nb["NB3"][1] += e2 > NEG; nb["NB3"][2] = min(nb["NB3"][2], e2)
    print("    r = %s done: %d configurations so far, worst %.2e  [%.0fs]" % ([norm(p) for p in primes], ntot, worst1, time.time() - t0)); sys.stdout.flush()
print("\n    case      sigma  kappa-exp   #configurations   worst |lhs-rhs|/S   worst rel (well-conditioned)   #ill-conditioned   smallest cond")
for key in sorted(stat): print("    %-8s  %s      %d           %5d             %.2e            %.2e                       %4d               %.1e" % (key + tuple(stat[key])))
print("\n    ill-conditioned configurations (cond = |lhs|/S < 1e-3): %d of %d.  The 10 with the smallest cond:" % (nill, ntot))
print("      cond       rel        |lhs-rhs|   r         h_0        case     sigma kappa-exp")
for row in sorted(ill)[:10]: print("      %.2e   %.2e   %.2e    %-8s  %-9s  %-8s %s     %d" % row)
report("B1", nfail == 0, "[criterion v2] %d failures / %d configurations ; worst |lhs-rhs|/S = %.2e (tolerance %.0e) ; worst relative difference: %.2e among the %d well-conditioned, %.2e among all ; narrow point rule used in %d configurations"
       % (nfail, ntot, worst1, TOL, worstrel_well, ntot - nill, worstrel_all, nnarrow))
print("       criterion v1 (relative difference <= 1e-9 everywhere) would be: %s" % ("met" if worstrel_all <= TOL else "NOT met"))
cases = {}
for (case, sg, kf), d in stat.items(): cases.setdefault(case, [set(), set()]); cases[case][0].add(sg); cases[case][1].add(kf)
cov = cases.get("3|c", [set()])[0] == {"0"} and cases.get("v=1", [set()])[0] == {"+", "-"} and cases.get("coprime", [set()])[0] == {"0", "+", "-"} and all(len(c[1] - {0}) > 0 for c in cases.values())
report("COV", cov, "sigma and kappa exponents seen per case: %s" % {c: (sorted(v[0]), sorted(v[1])) for c, v in cases.items()})
report("NB1", nb["NB1"][0] > 0 and nb["NB1"][1] == nb["NB1"][0], "live in %d configurations, rejected in %d, smallest relative difference %.2e" % tuple(nb["NB1"]))
report("NB2", nb["NB2"][0] > 0 and nb["NB2"][1] == nb["NB2"][0], "live in %d configurations, rejected in %d, smallest relative difference %.2e" % tuple(nb["NB2"]))
report("NB3", nb["NB3"][0] > 0 and nb["NB3"][1] >= 0.95 * nb["NB3"][0], "live in %d configurations, rejected in %d, smallest relative difference %.2e" % tuple(nb["NB3"]))
print("\nLAYER B: %s   [%.0fs]" % ("all criteria met" if all(ok for _, ok in RES) else "SOME CRITERIA NOT MET: " + ",".join(n for n, ok in RES if not ok), time.time() - t0))
