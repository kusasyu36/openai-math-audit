# -*- coding: utf-8 -*-
"""
Stage 5 (b): the translates for which layer B / B2 found no point where theta is not tiny (so the PHASE kappa of
eq:theta-cusp-automorphy, paper2.tex:3161-3183, was not tested by a relative comparison).

Same identity as stage3_b_translate:  conj theta(z + a/c, v) = conj kappa(g_1) [d_sigma-expansion at g^-1(z + a/c, v)].
Configurations: the 540 (orbit representative h_0 = (x0, x1), 0 <= x0 < 6, 0 <= x1 < 18; r in {1, p7, p13, -5, p7*p13}),
k_p = 1 for every active prime (deterministic; B2 drew k_p at random, so the set of hard translates need not be the same 50).

PROCEDURE (fixed before the first run).  Points: v = t/|c|, z = zeta/|c|, (t, zeta) from the fixed list
  t in {1, 0.6, 1.6, 2.5},  zeta in {0} + {rho e^{i(0.2 + k pi/2)} : rho in {0.5, 1.2, 2.5}, k = 0..3}   (52 points, this order: zeta outer).
Truncation levels T in (40, 34, 28, 22), tried in this order.  At level T a point is ADMISSIBLE if both series are inside
the complete part of the tables (T/(4 pi v) and T/(4 pi v')).  For an admissible point:
  cond = |lhs|/S (S = sum of |terms| of the left series, truncated at T),
  tol  = 10 * (tb(v, T) + tb(v', T)) / |lhs|,   tb = rigorous tail bound (T1) of stage4_bounds (absolute, per series).
The first point (in list order, at the largest T) with cond >= 1e-2 and tol <= 1e-2 is SELECTED.
CRITERIA at the selected point:
  G1  |lhs - rhs| / |lhs| <= tol + 1e-9        (RELATIVE; tol is an upper bound for the truncation, 1e-9 covers rounding)
  NG1 (kappa in place of conj kappa; live where kappa != 1): |lhs - rhs_ctl| / |lhs| > 0.5 at every selected live point
      (the true deviation is |1 - omega^(+-2)| = 1.73).
A configuration with no selected point at any T is UNRESOLVED; for it the script prints N(c), the best cond per T and
the reason in numbers: (i) the reach limits, i.e. the largest value of (v^2 + |z|^2)/v * N(c) = 1/v' and of 1/v that the
tables allow at each T; (ii) the nearest cusps p/q with N(q) <= 13 and the height v_q = v/(N(q)(v^2 + |a/c - p/q|^2)) of
the base point (zeta = 0, t = 1) in their horoballs.  (ii) is a diagnostic; the interpretation is an inference.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage5_b_unresolved.py > audits/003dual/stage5_b_unresolved.out 2>/dev/null
"""
import sys, time, cmath
from scipy.special import kv
from stage2_lib import *
import stage4_bounds as sb

t0 = time.time()
print(__doc__)
th = Theta(1000000)
P7, P13, M5 = (1, 3), (4, 3), (-5, 0)
RS = [[], [P7], [P13], [M5], [P7, P13]]
setup = Setup((18, 0))
TS = (40.0, 34.0, 28.0, 22.0)
ZETAS = [0j] + [rho * cmath.exp(1j * (0.2 + k * math.pi / 2)) for rho in (0.5, 1.2, 2.5) for k in range(4)]
PTS = [(t, ze) for ze in ZETAS for t in (1.0, 0.6, 1.6, 2.5)]
_tb = {}
def tb(v, T):
    key = (round(math.log(v), 2), T)          # bound evaluated at a slightly smaller v (monotone: smaller v, larger bound)
    if key not in _tb: _tb[key] = sb.theta_tail_bound(math.exp(key[0] - 0.005) , T)
    return _tb[key]
def left(cfg, z, v, T):
    s = th.d["0"]; n = s.cut(v, T)
    terms = s.coef[:n] * v * kv(1.0 / 3, 4 * math.pi * s.aell[:n] * v) * np.exp(4j * math.pi * (s.ell[:n] * z).real) * AddChar(cfg["a"], mul(LAM3, cfg["c"]))(s.m[:n])
    c0 = SIGMA0 * v ** (2.0 / 3)
    return c0 + complex(np.sum(terms)), c0 + float(np.sum(np.abs(terms)))
def right(cfg, z, v, kexp, T):
    z2, v2 = cusp_point(cfg["c"], z, v); sg = cfg["sigma"]
    ph = AddChar(neg(cfg["dl"]), mul(LAM3, cfg["c"]))
    return W3 ** ((-kexp) % 3) * th.d[sg](z2, v2, T=T, complete_to=th.complete[sg], phase=ph)
def horoballs(cfg):
    """the three largest v_q over q with N(q) <= 13 (q up to units), base point z = 0, v = 1/|c|"""
    cc = cxf(cfg["c"]); Z = cxf(emod(cfg["a"], mul((3, 0), cfg["c"]))) / cc; v = 1 / abs(cc); out = []
    for q0 in range(0, 5):
        for q1 in range(-4, 5):
            q = (q0, q1)
            if q == (0, 0) or norm(q) > 13 or (q0 == 0 and q1 < 0): continue
            qc = cxf(q); w = Z * qc
            # nearest lattice point p to w = Z q
            y = w.imag / (math.sqrt(3) / 2); x = w.real + y / 2
            best = None
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    p = (round(x) + dx, round(y) + dy)
                    if norm(gcd(p, q)) != 1: continue
                    d = abs(Z - cxf(p) / qc)
                    if best is None or d < best[0]: best = (d, p)
            if best: out.append((v / (norm(q) * (v * v + best[0] ** 2)), norm(q), q, best[1], best[0] * abs(cc)))
    out.sort(reverse=True)
    seen = set(); res = []
    for row in out:
        key = round(row[0], 6)
        if key in seen: continue
        seen.add(key); res.append(row)
    return res[:3]

sel_stat = {T: 0 for T in TS}; nfail = 0; ng = [0, 0]; worst = {T: 0.0 for T in TS}; worst_tol = {T: 0.0 for T in TS}; unres = []
for primes in RS:
    for x0 in range(6):
        for x1 in range(18):
            h0 = (x0, x1); cfg = setup.build(h0, primes, [ONE] * len(primes))
            kf, _ = kappa_formula(cfg); ca = abs(cxf(cfg["c"])); sg = cfg["sigma"]; done = False; bestc = {}
            for T in TS:
                lim0 = 4 * math.pi * th.complete["0"] / T; lims = 4 * math.pi * min(th.complete["0"], th.complete[sg]) / T
                bc = 0.0
                for t, ze in PTS:
                    v = t / ca; z = ze / ca; _, v2 = cusp_point(cfg["c"], z, v)
                    if 1 / v > lim0 or 1 / v2 > lims: continue
                    lhs, S = left(cfg, z, v, T); cond = abs(lhs) / S; bc = max(bc, cond)
                    if cond < 1e-2: continue
                    tol = 10 * (tb(v, T) + tb(v2, T)) / abs(lhs)
                    if tol > 1e-2: continue
                    rhs = right(cfg, z, v, kf, T); e = abs(lhs - rhs) / abs(lhs)
                    sel_stat[T] += 1; worst[T] = max(worst[T], e); worst_tol[T] = max(worst_tol[T], tol)
                    if e > tol + 1e-9:
                        nfail += 1; print("    FAIL: r=%s h0=%s T=%g rel=%.2e tol=%.2e" % ([norm(p) for p in primes], h0, T, e, tol))
                    if kf != 0:
                        e2 = abs(lhs - right(cfg, z, v, -kf, T)) / abs(lhs); ng[0] += 1; ng[1] += e2 > 0.5
                    done = True; break
                bestc[T] = bc
                if done: break
            if not done: unres.append((cfg, h0, primes, bestc, kf))
    print("    r = %s done: selected per T so far %s, unresolved %d  [%.0fs]" % ([norm(p) for p in primes], {int(T): n for T, n in sel_stat.items()}, len(unres), time.time() - t0)); sys.stdout.flush()

print("\n  selected at T = 40/34/28/22: %s ; worst relative difference per T: %s ; largest tol per T: %s"
      % ([sel_stat[T] for T in TS], ["%.1e" % worst[T] for T in TS], ["%.1e" % worst_tol[T] for T in TS]))
print("  G1 %s: %d failures / %d selected configurations (of 540) ; %d unresolved" % ("PASS" if nfail == 0 else "FAIL", nfail, sum(sel_stat.values()), len(unres)))
print("  NG1 %s: live at %d selected points, rejected (> 0.5) at %d" % ("PASS" if ng[0] > 0 and ng[0] == ng[1] else "FAIL", ng[0], ng[1]))
print("\n  reach limits: at level T the tables allow 1/v <= %s (theta), 1/v' <= the same for sigma = 0 and <= %s for sigma = +-"
      % (["%.1f" % (4 * math.pi * th.complete["0"] / T) for T in TS], ["%.1f" % (4 * math.pi * th.complete["+"] / T) for T in TS]))
print("  since (1/v)(1/v') = N(c)(1 + |z|^2/v^2) >= N(c), both points fit only if N(c)(1 + |z/v|^2) <= (product of the two limits).")
print("\n  UNRESOLVED configurations (best cond per T; three deepest horoballs at the base point: v_q, N(q), |a/c - p/q| |c|):")
print("    r        h_0       case     sigma kappa  N(c)   best cond at T=40/34/28/22            horoballs")
for cfg, h0, primes, bestc, kf in unres:
    hb = horoballs(cfg)
    print("    %-8s %-9s %-8s %s     %d      %5d  %s   %s" % ([norm(p) for p in primes], h0, cfg["case"], cfg["sigma"], kf, norm(cfg["c"]),
          " ".join("%.1e" % bestc.get(T, float("nan")) for T in TS), " ; ".join("v_q=%.1f N(q)=%d d|c|=%.2f" % (r_[0], r_[1], r_[4]) for r_ in hb)))
print("\nSTAGE 5 (b): %s   [%.0fs]" % ("criteria met on the selected configurations" if nfail == 0 and ng[0] == ng[1] else "SOME CRITERIA NOT MET", time.time() - t0))
