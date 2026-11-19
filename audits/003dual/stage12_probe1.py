# -*- coding: utf-8 -*-
"""
Stage 12, probe 1 (PLANNING ONLY: bounds, the direct sum |T_dir| as the scale, and checks of the derivation; the dual sum
is NOT computed and no comparison is made).  Written 2026-11-19 after stage 5 had been seen.

 (1) conversion check of stage12_lib (EQUALITY expected): for each sigma and windows R1 < N(l) <= R2 inside the table,
     sum |d_sigma(l)| N(l)^-s computed from the d-tables  ==  the same sum computed through the family/A decomposition of
     stage12_lib with the table values |G(A)|.  Then the INEQUALITY P-bound(s, M) >= sum_{M < N(A) <= Nmax} |G(A)| N(A)^-s.
     These test the derivation; they are not the basis of the bound.
 (2) |Q(A + i tau)| monotone non-decreasing on the tau grid (used by Mt_upper).
 (3) Mt_upper (one-sided Riemann) vs Mt (log-trapezoid of stage5_bump) and vs Mt with 1600 Gauss-Legendre nodes.
 (4) planning table: for candidate weights (p, a), X = floor(10^6 e^-a), and candidate runs, the registered stage-5 tail
     bound (T3'), the structured tail bound of stage12_lib, and |T_dir|.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage12_probe1.py > audits/003dual/stage12_probe1.out
"""
from stage5_lib import *
import stage5_bump as s5
import stage12_lib as L
t0 = time.time()
print(__doc__)
E = EndToEnd(); th = E.th; td = th.K.td
NA_all = E.NA; absG = np.abs(td.G)

# ---------------- (1) conversion check
def table_sum(sg, s, R1, R2):
    S = th.d[sg]; N = S.aell ** 2
    m = (N > R1) & (N <= R2)
    return float(np.sum(np.abs(S.coef[m]) * N[m] ** (-s)))
def formula_sum(sg, s, R1, R2):
    """same sum via families, using table |G(A)|"""
    tot = 0.0
    def win(scale, mult, units):
        Nl = NA_all * scale; m = (Nl > R1) & (Nl <= R2)
        return units * mult * float(np.sum(absG[m] * Nl[m] ** (-s)))
    if sg == "0":
        for k in range(0, 12): tot += win(3.0 ** (3 * k - 3), 3 ** (k / 2 + 2.5), 2)
        for k in range(1, 12): tot += win(3.0 ** (3 * k - 4), 3 ** (k / 2 + 2), 6)
    else:
        tot = win(1 / 81.0, 9.0, 3)
    return tot
print("(1) conversion check, sum |d| N^-s over windows (table vs family decomposition):")
worst = 0.0
for sg in ("0", "+", "-"):
    lim = th.complete[sg] ** 2
    for (R1, R2) in ((0.0, lim / 100), (lim / 100, lim / 4), (lim / 4, lim)):
        for s in (0.0, 1.5, 3.0):
            a_ = table_sum(sg, s, R1, R2); b_ = formula_sum(sg, s, R1, R2); rel = abs(a_ - b_) / max(a_, 1e-300)
            worst = max(worst, rel)
            print("   sigma %s  (%.4g, %.4g]  s=%.1f : table %.10e  formula %.10e  rel diff %.1e" % (sg, R1, R2, s, a_, b_, rel))
print("   worst relative difference: %.1e  (%s)" % (worst, "OK" if worst < 1e-10 else "MISMATCH"))
print("   P-bound vs table partial sums (bound must be >= the partial sum):")
for s in (1.5, 2.5, 4.0):
    for M in (0.5, 10.0, 1e3, 1e5):
        msk = NA_all > M; part = float(np.sum(absG[msk] * NA_all[msk] ** (-s))); bd = L.P(s, M)
        print("     s=%.1f M=%g : partial sum (N(A) <= 10^6) %.6e   bound %.6e   ratio bound/partial %.3f %s" % (s, M, part, bd, bd / part, "" if bd >= part else "VIOLATION"))
print("   max |G(A)|/N(A)^(1/6) on table = %.6f (must be <= 1)" % float(np.max(absG / NA_all ** (1 / 6))))
for p in ((1, 3), (-2, -3), (4, 3)):
    lp = loc(p); print("   max |B_{p,1}| for p = %s (N = %d): %.6f" % (p, lp.N, float(np.max(np.abs(lp.Btab(1))))))
sys.stdout.flush()

# ---------------- (2), (3)
CANDW = [(4, 3.0), (4, 2.5), (3, 2.5), (5, 2.5), (4, 2.0), (6, 2.5)]
bumps = {}
for p, a in CANDW:
    t1 = time.time(); bp = s5.Bump(a, p); bp.setup_quadrature(); bumps[(p, a)] = (bp, L.MtUp(bp))
    print("weight p=%d a=%g: built [%.0fs]" % (p, a, time.time() - t1))
print("(2) |Q(A+i tau)| monotone on the tau grid:")
bp0 = bumps[(4, 3.0)][0]
for A in (1.0, 2.0, 4.0, 8.0, 12.0):
    lq = np.array([sb.logQ(complex(A, x)).real for x in bp0.taus])
    print("   A=%g: min diff of log|Q| along the grid = %.3e (>= 0 means non-decreasing)" % (A, float(np.min(np.diff(lq)))))
print("(3) Mt variants for p=4, a=3:")
s5._GL = np.polynomial.legendre.leggauss(1600)
bp16 = s5.Bump(3.0, 4)
s5._GL = np.polynomial.legendre.leggauss(800)
mu = bumps[(4, 3.0)][1]
for A in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0):
    print("   A=%g: Mt(trapezoid) %.6e  Mt_upper %.6e  Mt(GL 1600) %.6e  Mt_upper(GL 1600) %.6e" % (A, bp0.Mt(A), mu(A), bp16.Mt(A), L.Mt_upper(bp16, A)))
sys.stdout.flush()

# ---------------- (4) planning table
P7, P7b, P13 = (1, 3), (-2, -3), (4, 3)
RUNS = [("B", {P7: 1}), ("B", {P7b: 1}), ("B", {P13: 1}), ("F", {P7: 1}), ("A", {}), ("A", {P7: 1}), ("B", {P7: 2})]
print("(4) planning table: tail bounds (sum over groups) and |T_dir|; ratio = struct tail / |T_dir|")
for (p, a) in CANDW:
    bp, mu = bumps[(p, a)]
    X = float(math.floor(1.0e6 * math.exp(-a)))
    assert X * math.exp(a) <= 1.0e6
    wt = WBump(bp, 3.0e7)
    for nm, jd in RUNS:
        t1 = time.time()
        Td, Sd = E.direct(nm, jd, wt, X)
        G = E.groups(nm, jd, wt, X)
        old = sum(wt.dual_tail(X, g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"]) for g in G)
        new = 0.0; As = set(); per = {}
        for g in G:
            b, Ab = L.struct_tail(mu, X, g["cfg"]["sigma"], g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"])
            new += b; As.add(Ab); k = (g["cfg"]["sigma"], g["key"][1]); per[k] = per.get(k, 0) + b
        top = sorted(per.items(), key=lambda kv: -kv[1])[:3]
        print("   p=%d a=%g X=%d  %s %-10s |T_dir|=%.3e  old tail %.2e (%.1f digits)  struct tail %.2e (%.1f digits)  A* %s  largest (sigma,N(c)): %s  [%.0fs]"
              % (p, a, X, nm, {norm(q): j for q, j in jd.items()}, abs(Td), old, -math.log10(old / abs(Td)), new, -math.log10(new / abs(Td)),
                 sorted(As), " ".join("%s,%d:%.1e" % (k[0], k[1], v) for k, v in top), time.time() - t1))
        del G
        sys.stdout.flush()
print("\nprobe 1 done [%.0fs]" % (time.time() - t0))
