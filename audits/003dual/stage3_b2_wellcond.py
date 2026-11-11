# -*- coding: utf-8 -*-
"""
Stage 3, layer B, SECOND PASS (written after the official run of stage3_b_translate; criteria below were fixed
before this script was run for the first time).  Reason: in layer B, 63 of 1620 configurations were evaluated at
points where theta itself is tiny (cond = |lhs|/S < 1e-3), so nothing was learned about the PHASE (kappa) there.
Also the 324 classes h_0 mod 18 are only 108 different translates (h_0 -> h_0 + 6 shifts z by -1, a period of
theta; stage2_h0_period.out), so this pass uses one representative per orbit.

Same identity as stage3_b_translate (eq:theta-cusp-automorphy :3161-3170 with kappa from eq:ray-multiplier).
Configurations: L = 18, h_0 = (x0, x1) with 0 <= x0 < 6, 0 <= x1 < 18 (108 orbit representatives),
  r in {1, p7, p13, -5, p7*p13}: 540 configurations, one tuple (k_p) each (random choice as in layer B).
POINT SELECTION (uses the LEFT side only, never the difference): at most 12 draws per configuration.
  Draws 1-6: v = t/|c|, z = zeta/|c|, t uniform in [0.6, 1.6], Re zeta, Im zeta uniform in [-0.8, 0.8]; a draw is
  skipped if the tables do not cover it (T/(4 pi v) > complete range of theta, or T/(4 pi v') > complete range of
  the table of d_sigma).  Draws 7-12: the narrow rule of layer B (t in [0.97, 1.03], |Re zeta|, |Im zeta| <= 0.1).
  The first covered draw with cond = |lhs|/S >= 1e-2 is SELECTED (S = sum of |terms| of the left series).
  A configuration with no selected point is reported as 'not resolved' and is not counted.
CRITERIA:
  B2  at the selected point, |lhs - rhs| / max(|lhs|, |rhs|) <= 1e-9  (RELATIVE to the value).  0 failures.
  NB1 (kappa instead of conj kappa; live where kappa != 1): relative difference > 1e-4 at every selected live point.
  NB2 (d_+ <-> d_-; live where sigma != 0): relative difference > 1e-4 at >= 99% of the selected live points
      (layer B showed that this control can be powerless at a point dominated by a nearby cusp).
  Reported: number of configurations not resolved, per case; distribution of the number of draws.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage3_b2_wellcond.py > audits/003dual/stage3_b2_wellcond.out 2>/dev/null
"""
import sys, time, random
from stage2_lib import *

t0 = time.time()
print(__doc__)
RES = []
def report(name, ok, detail):
    RES.append((name, ok)); print("  %-4s %s  %s" % (name, "PASS" if ok else "FAIL", detail)); sys.stdout.flush()
random.seed(2026111206)
TOL, NEG, T = 1e-9, 1e-4, 40.0
th = Theta(int(os.environ.get("STAGE_NMAX", "1000000")))
print("  theta tables: N(A) <= %d ; complete to |l| <= %.2f (theta), %.2f (cusp series)  [%.0fs]" % (th.Nmax, th.complete["0"], th.complete["+"], time.time() - t0))
def rel(x, y): return abs(x - y) / max(abs(x), abs(y), 1e-300)
P7, P13, M5 = (1, 3), (4, 3), (-5, 0)
RS = [[], [P7], [P13], [M5], [P7, P13]]
KT = {0: [[]], 1: [[(1, 0)], [(2, 0)], [(1, 1)], [(3, -2)]], 2: [[(1, 0), (1, 0)], [(2, 0), (1, 1)], [(0, 1), (3, 0)]]}
setup = Setup((18, 0))
def right(cfg, z, v, kexp, sg=None):
    z2, v2 = cusp_point(cfg["c"], z, v); sg = cfg["sigma"] if sg is None else sg
    ph = AddChar(neg(cfg["dl"]), mul(LAM3, cfg["c"]))
    return W3 ** ((-kexp) % 3) * th.d[sg](z2, v2, T=T, complete_to=th.complete[sg], phase=ph)
stat = {}; unresolved = {}; ndraw = {}; nfail = 0; nsel = 0; worst = 0.0; mincond = 1e9
nb = {"NB1": [0, 0, 1e9], "NB2": [0, 0, 1e9]}
for primes in RS:
    for x0 in range(6):
        for x1 in range(18):
            h0 = (x0, x1); ks = random.choice(KT[len(primes)])
            cfg = setup.build(h0, primes, ks); kf, _ = kappa_formula(cfg); ca = abs(cxf(cfg["c"])); sg = cfg["sigma"]
            sel = None
            for it in range(12):
                if it < 6:
                    t = random.uniform(0.6, 1.6); zeta = complex(random.uniform(-.8, .8), random.uniform(-.8, .8))
                else:
                    t = random.uniform(0.97, 1.03); zeta = complex(random.uniform(-.1, .1), random.uniform(-.1, .1))
                v = t / ca; z = zeta / ca; _, v2 = cusp_point(cfg["c"], z, v)
                if T / (4 * math.pi * v) > th.complete["0"] or T / (4 * math.pi * v2) > min(th.complete["0"], th.complete[sg]): continue
                lhs = th.thetabar_at(z, v, shift=(cfg["a"], cfg["c"]), T=T)
                S = th.d["0"](z, v, T=T, complete_to=th.complete["0"], absum=True)
                if abs(lhs) / S >= 1e-2: sel = (z, v, lhs, S, it + 1); break
            key = (cfg["case"], sg, kf)
            if sel is None:
                unresolved[key] = unresolved.get(key, 0) + 1; continue
            z, v, lhs, S, it = sel; ndraw[it] = ndraw.get(it, 0) + 1; nsel += 1; mincond = min(mincond, abs(lhs) / S)
            rhs = right(cfg, z, v, kf); e = rel(lhs, rhs); worst = max(worst, e); nfail += e > TOL
            d = stat.setdefault(key, [0, 0.0]); d[0] += 1; d[1] = max(d[1], e)
            if kf != 0:
                e2 = rel(lhs, right(cfg, z, v, -kf)); nb["NB1"][0] += 1; nb["NB1"][1] += e2 > NEG; nb["NB1"][2] = min(nb["NB1"][2], e2)
            if sg != "0":
                e2 = rel(lhs, right(cfg, z, v, kf, sg={"+": "-", "-": "+"}[sg])); nb["NB2"][0] += 1; nb["NB2"][1] += e2 > NEG; nb["NB2"][2] = min(nb["NB2"][2], e2)
    print("    r = %s done: %d selected so far, worst relative difference %.2e  [%.0fs]" % ([norm(p) for p in primes], nsel, worst, time.time() - t0)); sys.stdout.flush()
print("\n    case      sigma  kappa-exp   #selected   worst relative difference   #not resolved")
for key in sorted(set(stat) | set(unresolved)):
    d = stat.get(key, [0, 0.0]); print("    %-8s  %s      %d           %4d        %.2e                    %d" % (key + (d[0], d[1], unresolved.get(key, 0))))
print("    number of draws needed: %s ; smallest cond among selected points: %.2e" % (sorted(ndraw.items()), mincond))
nun = sum(unresolved.values())
report("B2", nfail == 0 and nsel > 0, "%d failures / %d selected configurations (of 540 ; %d not resolved) ; worst RELATIVE difference %.2e (tolerance %.0e)" % (nfail, nsel, nun, worst, TOL))
report("NB1", nb["NB1"][0] > 0 and nb["NB1"][1] == nb["NB1"][0], "live at %d selected points, rejected at %d, smallest relative difference %.2e" % tuple(nb["NB1"]))
report("NB2", nb["NB2"][0] > 0 and nb["NB2"][1] >= 0.99 * nb["NB2"][0], "live at %d selected points, rejected at %d, smallest relative difference %.2e" % tuple(nb["NB2"]))
print("\nLAYER B, SECOND PASS: %s   [%.0fs]" % ("all criteria met" if all(ok for _, ok in RES) else "SOME CRITERIA NOT MET: " + ",".join(n for n, ok in RES if not ok), time.time() - t0))
