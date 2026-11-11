# -*- coding: utf-8 -*-
"""
Stage 3, layer E (end to end, with the weight): T(X; Psi) of eq:T (paper2.tex:1287-1293) against the sum over the
groups (h_0, active set) of the terms eq:reflection (:1806-1812)

    C * sum_{l in lambda^-4 O, l != 0}  d(l) alpha(l) N(l)^-1/2  psi(lambda^4 l)  prod_{p active} B_{p,j_p}(lambda^4 l)
          * V_*^sharp( N(l) X / N(c)^2 ),
    C = -(i/81) conj(alpha(c))^2 phi-hat(h_0) conj(kappa_0) prod_{inactive}(1 - 1/N(p)) prod_{active} chi_p(sigma_p)^-2 omega_{p,j_p}
                                                                                            (:3420-3424, as printed)
    V_*^sharp from eq:theta-weight (:1752-1757, as printed), alpha(x) = x/|x| (:323).

LEFT (direct):  X^-1/2 sum over primary A = n b^3 prime to 2 (n squarefree) of
                conj(alpha(A)) gamma_2(n) |b| Psi(A) W(N(A)/X),   Psi = Psi_0 prod_p chi_p^{j_p}   [= eq:T with V_*(y) = sqrt(y) W(y)]
                gamma_2(n) = (lambda/n)_3 g~(n) from kubota.ThetaData (brute-force cubic Gauss sums), all A with N(A) <= Nmax.
RIGHT (dual):   the display above, all table entries with N(l) <= R_g (R_g printed), d_sigma from kubota.py.
WEIGHT: W(y) = exp(-(log y)^2/(4w)).  NOT compactly supported: a deviation from the hypothesis of the proposition
        (see stage4_bounds.py); what is tested is the identity for this W.

E0 (intermediate gate, pointwise, per group): d/dzbar at z = 0 of the left side of layer D equals
        -(c v)^-2 * C' * d/dz' [dual series] at (z', v') = (0-part, 1/(N(c) v))   (:3347-3358),
        i.e. the factor -1/c^2 = -conj(alpha(c))^2/N(c) before any Mellin transform.

CRITERIA (fixed before the first run; the bounds are those of stage4_bounds.py and do not use the data):
  E0  |lhs - rhs| <= 1e-9 * scale, scale = (sum_h |c_F|) * (sum of |terms| of the differentiated series of conj theta);
      Psi_0 in {A, B, F}, P = {p7: 1}, every group, v = 1/|c| (T = 40).
  E1  for each run (Psi_0, P, j, w, X):  |T_direct - T_dual| <= Btot, where
        Btot = (T2) direct tail + sum over groups of (T3) dual tails + (Q) quadrature bound * sum|coefficients|
               + rounding allowance (ESTIMATE, not a bound): 1e-13 * (sum |direct terms| + sum |dual terms|) + 2e-15 * sum |dual coefficients|.
      The run counts only if it is INFORMATIVE:  Btot <= 1e-7 * |T_direct|  (otherwise it is listed as 'not informative',
      which is a failure of my choice of (X, w), not a statement about the paper).
      Reported for every run: |diff| (ABSOLUTE), |diff|/|T_direct| (RELATIVE), each part of Btot.
  Negative controls, run on every informative run where live; each must give |T_direct - T_ctl| > max(100 Btot, 1e-3 |T_direct|):
  NE1 +i/81 in place of -i/81                                   live: always
  NE2 alpha(c)^2 in place of conj(alpha(c))^2                   live: alpha(c)^4 != 1 in some group
  NE3 conj(alpha(l)) in place of alpha(l)                       live: always
  NE4 V^sharp(N(l) X / N(c)) (N(c) not squared)                 live: N(c) != 1 in some group
  NE5 phi-hat from phi WITHOUT the condition at 2               live: always (Psi_0 in A, B, C, F)
  NE6 kappa_0 in place of conj(kappa_0)                         live: kappa_0 != 1 in some group
  NE7 the gamma quotient of eq:theta-weight inverted            live: always
Runs: RUNS below, chosen from the tables of stage4_bounds.out (bounds only) before any comparison.
Two primes are NOT run here: with tables up to N(A) <= 10^6 the bounds cannot be met for N(c) >= 364 (stage4_bounds.out).
[Added later: with STAGE_EXTRA=1 the script runs six two-prime configurations {p7, conj p7} with Psi_0 = B instead
 (output stage3_d_mellin_extra.out); see the comment at RUNS.]
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage3_d_mellin.py > audits/003dual/stage3_d_mellin.out 2>/dev/null
"""
import sys, time, itertools
from scipy.special import kv
from stage2_lib import *
import stage4_bounds as sb

t0 = time.time()
print(__doc__)
RES = []
def report(name, ok, detail):
    RES.append((name, ok)); print("  %-4s %s  %s" % (name, "PASS" if ok else "FAIL", detail)); sys.stdout.flush()
NMAX = int(os.environ.get("STAGE_NMAX", "1000000"))
LINT = 18; setup = Setup((LINT, 0))
P7, P13, M5 = (1, 3), (4, 3), (-5, 0)
RUNS = [(nm, {}, w, X) for nm in "ABCDEF" for (w, X) in ((0.1, 500.0), (0.1, 2000.0))] \
     + [("A", {P7: j}, w, X) for j in (1, 0, 4) for (w, X) in ((0.05, 64000.0), (0.1, 22000.0))] \
     + [("B", {P7: j}, w, X) for j in range(6) for (w, X) in ((0.1, 2000.0), (0.1, 8000.0))] \
     + [("F", {P7: j}, w, X) for j in (1, 0, 4) for (w, X) in ((0.1, 10000.0), (0.1, 16000.0))] \
     + [("C", {P7: 1}, 0.1, 16000.0), ("E", {P7: 1}, 0.1, 16000.0), ("E", {P7: 2}, 0.1, 16000.0)] \
     + [("B", {P13: 1}, 0.1, 8000.0), ("B", {M5: 1}, 0.1, 16000.0), ("F", {P13: 5}, 0.05, 64000.0)]
if os.environ.get("STAGE_EXTRA"):
    # EXTRA RUNS, added after the first official run of this layer (which had no two-prime run).  Criteria unchanged.
    # Two primes of norm 7 (p7 and its conjugate): N(c) in {49, 196} for Psi_0 = B, the only two-prime family for
    # which the bounds of stage4_bounds can be met with the tables N(A) <= 10^6.  Chosen from the bounds alone.
    P7b = (-2, -3)
    assert norm(P7b) == 7 and is_primary(P7b)
    RUNS = [("B", jd, 0.05, 64000.0) for jd in ({P7: 1, P7b: 1}, {P7: 1, P7b: 5}, {P7: 0, P7b: 4}, {P7: 2, P7b: 0}, {P7: 0, P7b: 0}, {P7: 3, P7b: 4})]
if os.environ.get("STAGE_DEBUG"):
    RUNS = [("B", {}, 0.1, 500.0), ("A", {}, 0.1, 500.0), ("B", {P7: 1}, 0.1, 2000.0), ("F", {P7: 0}, 0.1, 3000.0), ("B", {P7: 4}, 0.1, 2000.0)]
XMAX = {0.05: 3.0e4, 0.1: 1.5e4, 0.2: 1.0e4}       # dual terms are kept up to x = N(l) X/N(c)^2 <= XMAX (and within the tables)

th = Theta(NMAX)
print("  theta tables: N(A) <= %d ; complete to |l| <= %.2f (theta), %.2f (cusp series)  [%.0fs]" % (th.Nmax, th.complete["0"], th.complete["+"], time.time() - t0))
PH = {nm: phi_hat(phi_table(nm, LINT), LINT) for nm in "ABCDEF"}
PHno2 = {nm: phi_hat(phi_table(nm, LINT, with2=False), LINT) for nm in "ABCF"}
d0 = th.d["0"]
# (D) check of the coefficient bound on the tables
for sg in ("0", "+", "-"):
    s = th.d[sg]; nn = s.cut(1.0, 4 * math.pi * th.complete[sg])
    print("  bound (D): max |d_%s(l)| / N(l)^(1/6) over the complete part of the table = %.4f  (bound 27)" % (sg, float(np.max(np.abs(s.coef[:nn]) / s.aell[:nn] ** (1.0 / 3)))))

# ---------------------------------------------------------------- E0
import random
random.seed(2026111205)
T = 40.0
def left_dzbar(nm, jd, h0, act, v):
    primes = list(jd); r = prod(primes); Cden = mul((LINT, 0), r); locs = [loc(p) for p in primes]
    n = d0.cut(v, T)
    assert T / (4 * math.pi * v) <= th.complete["0"] * (1 + 1e-12)
    base = d0.coef[:n] * v * kv(1.0 / 3, 4 * math.pi * d0.aell[:n] * v) * 2j * math.pi * np.conj(d0.ell[:n])
    m = d0.m[:n]; F = np.zeros(n, dtype=complex); aF = 0.0
    rng = [range(1, lp.N) if p in act else [0] for lp, p in zip(locs, primes)]
    Cs = [lp.C(jd[p]) for lp, p in zip(locs, primes)]
    for idxs in itertools.product(*rng):
        cF = PH[nm][h0]
        for Cc, i in zip(Cs, idxs): cF *= Cc[i]
        hs = [lp.L.reps[i] for lp, i in zip(locs, idxs)]
        A = mul(LAM2, add(mul(h0, r), mul((LINT, 0), esum([mul(h, exact_div(r, p)) for h, p in zip(hs, primes)]))))
        F += cF * AddChar(A, mul(LAM3, Cden))(m); aF += abs(cF)
    return complex(np.sum(base * F)), aF * float(np.sum(np.abs(base)))
w0 = 0.0; n0 = 0; f0n = 0; negmin = 1e9
for nm in ("A", "B", "F"):
    jd = {P7: 1}
    for h0, f0 in PH[nm].items():
        if abs(f0) < 1e-13: continue
        cfg = setup.build(h0, [P7], [ONE]); cc = cxf(cfg["c"]); ca = abs(cc)
        if min(th.complete["0"], th.complete[cfg["sigma"]]) < T / (4 * math.pi) * ca: continue
        v = 1.0 / ca
        lhs, sc = left_dzbar(nm, jd, h0, [P7], v)
        rhs = -1.0 / (cc * cc * v * v) * Cprime(cfg, f0, jd) * dual_value(th, cfg, jd, 0.0, 1.0 / (ca * ca * v), T=T, deriv="z")
        e = abs(lhs - rhs) / sc; w0 = max(w0, e); n0 += 1; f0n += e > 1e-9
        negmin = min(negmin, abs(lhs + rhs) / sc)
report("E0", f0n == 0 and n0 > 0, "%d failures / %d groups ; worst |lhs-rhs|/scale = %.2e ; (control: sign of -(cv)^-2 flipped, smallest diff/scale %.2e)" % (f0n, n0, w0, negmin))

# ---------------------------------------------------------------- E1
td = th.K.td
A_ = td.A; NA = (A_[:, 0] ** 2 - A_[:, 0] * A_[:, 1] + A_[:, 1] ** 2).astype(float)
Ac = A_[:, 0] + A_[:, 1] * W3; albar = np.conj(Ac / np.abs(Ac))
aW = np.conj(td.G) * W3 ** td.kl                      # gamma_2(n)|b|
def psi0_tab(nm):
    tab = np.zeros((LINT, LINT), dtype=complex)
    for x0 in range(LINT):
        for x1 in range(LINT):
            n = (x0, x1)
            if is_primary(n) and not (x0 % 2 == 0 and x1 % 2 == 0): tab[x0, x1] = W3 ** psi0_exp(nm, n)
    return tab
P0T = {nm: psi0_tab(nm) for nm in "ABCDEF"}

def direct(nm, jd, w, X):
    f = aW * albar * P0T[nm][A_[:, 0] % LINT, A_[:, 1] % LINT]
    for p, j in jd.items():
        lp = loc(p); f = f * lp.chitab(j)[lp.idx_arr(A_)]
    t = f * np.exp(-(np.log(NA / X)) ** 2 / (4 * w)) / math.sqrt(X)
    return complex(np.sum(t)), float(np.sum(np.abs(t)))

def groups(nm, jd, w, X, ph):
    """per group: dict(cfg, Cp (C'), inact, key, arrays).  arr = d(l) psi prod B / conj(l)  [alpha(l)/sqrt N(l) = 1/conj(l)]"""
    primes = list(jd); must = [p for p in primes if jd[p] != 0]; opt = [p for p in primes if jd[p] == 0]
    out = []
    for h0, f0 in ph.items():
        if abs(f0) < 1e-13: continue
        for kopt in range(len(opt) + 1):
            for extra in itertools.combinations(opt, kopt):
                act = [p for p in primes if p in must or p in extra]; inact = [p for p in primes if p not in act]
                cfg = setup.build(h0, act, [ONE] * len(act)); Nc = norm(cfg["c"]); sg = cfg["sigma"]; s = th.d[sg]
                R = min(th.complete[sg] ** 2, XMAX[w] * Nc * Nc / X)
                n = int(np.searchsorted(s.aell, math.sqrt(R), side="right")); m = s.m[:n]
                arr = s.coef[:n] * psi_char(cfg)(m) / np.conj(s.ell[:n]); Bmax = 1.0
                for p in act:
                    lp = loc(p); j = jd[p]; arr = arr * lp.Btab(j)[lp.idx_arr(m)]
                    Bmax *= {4: (lp.N - 1) / math.sqrt(lp.N), 0: 1 / math.sqrt(lp.N)}.get(j, 1.0)
                out.append(dict(cfg=cfg, Cp=Cprime(cfg, f0, jd, inact), key=(sg, Nc, n), arr=arr, R=R, Bmax=Bmax,
                                k0=kappa_formula(cfg)[1], ac=alpha(cfg["c"])))
    return out

_VS = {}
def vs_cached(Nm, Nc, X, w, variant):
    """V^sharp at x = (Nm/81) X / Nc^2 for an integer array Nm (unique), cached"""
    key = (Nc, X, w, variant, len(Nm), int(Nm[-1]) if len(Nm) else 0)
    if key not in _VS:
        x = Nm / 81.0 * X / (Nc * Nc if variant != "NE4" else Nc)
        if variant == "NE7":
            saved = sb.logQ; sb.logQ = lambda t: -saved(t)
            try: _VS[key] = sb.vsharp(x, w)
            finally: sb.logQ = saved
        else:
            _VS[key] = sb.vsharp(x, w)
    return _VS[key]

def dual(G, w, X, variant=None):
    agg = {}
    for g in G:
        C = -1j / 81 * np.conj(g["ac"]) ** 2 * g["Cp"]
        if variant == "NE1": C = -C
        if variant == "NE2": C = -1j / 81 * g["ac"] ** 2 * g["Cp"]
        if variant == "NE6": C = C * W3 ** ((2 * g["k0"]) % 3)
        arr = g["arr"]
        if variant == "NE3": arr = arr * (np.conj(th.d[g["key"][0]].ell[:g["key"][2]]) / th.d[g["key"][0]].ell[:g["key"][2]])
        if g["key"] in agg: agg[g["key"]][0] += C * arr; agg[g["key"]][1] += abs(C) * np.abs(arr)
        else: agg[g["key"]] = [C * arr, abs(C) * np.abs(arr)]
    tot = 0; sabs = 0.0; cabs = 0.0; qb = 0.0
    for (sg, Nc, n), (E, Ea) in agg.items():
        m = th.d[sg].m[:n]; Nm = m[:, 0] ** 2 - m[:, 0] * m[:, 1] + m[:, 1] ** 2
        un, inv = np.unique(Nm, return_inverse=True)
        V = vs_cached(un, Nc, X, w, variant if variant in ("NE4", "NE7") else None)[inv]
        tot += complex(np.sum(E * V)); sabs += float(np.sum(Ea * np.abs(V))); cabs += float(np.sum(Ea))
        if len(un):
            xs = un[[0, -1]] / 81.0 * X / (Nc * Nc)
            qb += float(np.sum(Ea)) * max(sum(sb.quad_error_bound(float(x), w)) for x in xs)
    return tot, sabs, cabs, qb

print("\n  columns: |T_dir|, |T_dir - T_dual| (ABSOLUTE), the same divided by |T_dir| (RELATIVE), Btot and its parts")
ninf = 0; nfail = 0; nrun = 0; worst_abs = 0.0; worst_rel = 0.0; worst_ratio = 0.0
ctl = {k: [0, 0, 1e9] for k in ("NE1", "NE2", "NE3", "NE4", "NE5", "NE6", "NE7")}
covE = {}
for nm, jd, w, X in RUNS:
    t1 = time.time()
    Td, Sd = direct(nm, jd, w, X)
    G = groups(nm, jd, w, X, PH[nm])
    Tu, Su, Cu, qb = dual(G, w, X)
    b_dir = sb.direct_tail_bound(X, w, NMAX)
    b_dual = sum(sb.dual_tail_bound(X, w, g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"])[0] for g in G)
    b_round = 1e-13 * (Sd + Su) + 2e-15 * Cu
    Btot = b_dir + b_dual + qb + b_round
    diff = abs(Td - Tu); nrun += 1
    inf_ = Btot <= 1e-7 * abs(Td); ok = diff <= Btot
    print("  Psi_0=%s  P=%-10s w=%.2f X=%-7g #groups=%3d N(c) in %s cases %s" % (nm, {norm(p): j for p, j in jd.items()}, w, X, len(G),
          sorted(set(g["key"][1] for g in G)), sorted(set((g["cfg"]["case"], g["cfg"]["sigma"]) for g in G))))
    print("      T_dir = % .12e %+.12ei   T_dual = % .12e %+.12ei" % (Td.real, Td.imag, Tu.real, Tu.imag))
    print("      |T_dir| = %.4e  |diff| = %.2e (abs)  %.2e (rel)   Btot = %.2e = direct tail %.1e + dual tail %.1e + quadrature %.1e + rounding allowance (estimate) %.1e ; sum|direct terms| = %.2e, sum|dual terms| = %.2e   %s%s  [%.0fs]"
          % (abs(Td), diff, diff / abs(Td), Btot, b_dir, b_dual, qb, b_round, Sd, Su, "PASS" if ok else "FAIL", "" if inf_ else "  (NOT INFORMATIVE: Btot > 1e-7 |T|)", time.time() - t1))
    if not inf_: ninf += 1; continue
    if not ok: nfail += 1
    worst_abs = max(worst_abs, diff); worst_rel = max(worst_rel, diff / abs(Td)); worst_ratio = max(worst_ratio, diff / Btot)
    for g in G:
        for p in jd:
            key = (g["cfg"]["case"], g["cfg"]["sigma"], jd[p], "active" if p in g["cfg"]["primes"] else "inactive"); covE[key] = covE.get(key, 0) + 1
        if not jd: covE[(g["cfg"]["case"], g["cfg"]["sigma"], "-", "no prime")] = covE.get((g["cfg"]["case"], g["cfg"]["sigma"], "-", "no prime"), 0) + 1
    thr = max(100 * Btot, 1e-3 * abs(Td)); line = []
    lives = {"NE1": True, "NE2": any(abs(g["ac"] ** 4 - 1) > 1e-9 for g in G), "NE3": True, "NE4": any(g["key"][1] != 1 for g in G),
             "NE5": nm in PHno2, "NE6": any(g["k0"] != 0 for g in G), "NE7": True}
    for k in sorted(ctl):
        if not lives[k]: line.append("%s dead" % k); continue
        if k == "NE5": Tc = dual(groups(nm, jd, w, X, PHno2[nm]), w, X)[0]
        else: Tc = dual(G, w, X, variant=k)[0]
        d = abs(Td - Tc); ctl[k][0] += 1; ctl[k][1] += d > thr; ctl[k][2] = min(ctl[k][2], d / abs(Td))
        line.append("%s %.1e" % (k, d / abs(Td)))
    print("      controls, |T_dir - T_ctl|/|T_dir|: " + "  ".join(line)); sys.stdout.flush()

print("\n  coverage of the informative runs (case of H, sigma, j, role): number of (group, prime) incidences")
for key in sorted(covE, key=str): print("    %-8s sigma=%s  j=%s  %-9s : %d" % (key + (covE[key],)))
print()
report("E1", nfail == 0 and nrun - ninf > 0, "%d failures / %d informative runs (%d runs not informative, %d runs in all) ; worst |diff| = %.2e (ABSOLUTE), worst |diff|/|T| = %.2e (RELATIVE), worst |diff|/Btot = %.2e"
       % (nfail, nrun - ninf, ninf, nrun, worst_abs, worst_rel, worst_ratio))
for k in sorted(ctl):
    report(k, ctl[k][0] > 0 and ctl[k][1] == ctl[k][0], "live in %d runs, rejected in %d, smallest |T_dir - T_ctl|/|T_dir| = %.2e" % tuple(ctl[k]))
print("\nLAYER E: %s   [%.0fs]" % ("all criteria met" if all(ok for _, ok in RES) else "SOME CRITERIA NOT MET: " + ",".join(n for n, ok in RES if not ok), time.time() - t0))
