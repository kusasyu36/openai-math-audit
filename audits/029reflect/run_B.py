# -*- coding: utf-8 -*-
"""Task B: why the "paper's combination" showed no visible rapid decay in audits/029core/theta_sum_check.py.

Objects (K = Q(omega), A = c d^3 primary, c squarefree):
  Patterson's own :  a_P(A) = g~(c)|d|                      (= 3^{-5/2} conj tau(lam^-3 A),  DR (5.7))
  paper's (003/029 normalisation of the Gauss sum):  a_W(A) = gamma_2(c)|d| = (lam/c)_3 g~(c)|d|
  T_W(X) = X^{-1/2} sum_A a(A) conj(alpha(A)) W(N(A)/X),   alpha(A) = A/|A|.
Both are sums of coefficients of  d/dzbar Theta_phi(z,v) at z = 0,  Theta_phi = sum conj(tau(mu)) phi(lam^3 mu) v K e^(mu z),
  phi_P(n) = 1[n = 1 mod 3],    phi_W(n) = (lam/n)_3 1[n = 1 mod 3].
Steps: B0 period of (lam/n)_3;  B1 finite Fourier expansion, reduced denominators c_0, coset types;
       B2 pointwise identity direct = dual (function and d/dzbar);  B3 decay with the Bessel kernel;
       B4 log-Gaussian weights: dual formula reproduces the direct sums (and the numbers printed by the old script);
       B5 phi_W combined with a moving sextic twist chi_p (brute-force layer).
Tolerance for identities: relative error <= 1e-8.
"""
import os, math, cmath, random
import numpy as np
from scipy.special import loggamma
from kubota import *

TOL = 1e-8
random.seed(2026111102)
CACHE = os.environ.get("TMPDIR", "/tmp")
NMAX = 200000
K = Kubota(NMAX, math.sqrt(NMAX / 27.0), CACHE)
thb = K.theta.conj()
nfail = 0; worst = 0.0
def report(name, lhs, rhs, expect_pass=True, scale=None):
    global nfail, worst
    err = abs(lhs - rhs) / (scale if scale else max(abs(lhs), abs(rhs), 1e-300))
    ok = (err <= TOL) == expect_pass
    if not ok: nfail += 1
    if expect_pass: worst = max(worst, err)
    print("  %-50s lhs=% .10e%+.10ei  rhs=% .10e%+.10ei  err=%.1e  %s" % (name, lhs.real, lhs.imag, rhs.real, rhs.imag, err,
          ("OK" if ok else "FAIL") + ("" if expect_pass else " (negative control)")))

# ---------------------------------------------------------------- B0
print("=== B0: period of (lam/n)_3 on primary n (brute force) ===")
vals = {}; bad9 = 0; cnt = 0
R = 80
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        n = (a, b)
        if not is_primary(n) or norm(n) > 5000: continue
        k = cubic_symbol(LAM, n); cnt += 1
        key = (a % 9, b % 9)
        if key in vals and vals[key] != k: bad9 += 1
        vals.setdefault(key, k)
print("  primary n with N(n) <= 5000: %d ; residue classes mod 9 seen: %d (expected 9) ; inconsistencies mod 9: %d" % (cnt, len(vals), bad9))
# is a smaller modulus enough?  test modulus 3*lam: n ~ n' iff (n - n')/(3 lam) integral
smaller = True
keys = list(vals)
for k1 in keys:
    for k2 in keys:
        if exact_div(sub(k1, k2), mul((3, 0), LAM)) is not None and vals[k1] != vals[k2]: smaller = False
print("  depends only on n mod 3*lam ? %s   => exact period ideal is (9) [conductor of the twist]" % smaller)
# DR (1.5): d = 1 + a2 lam^2 + a3 lam^3 mod 9  =>  (lam/d)_3 = omega^{-a3}
ok15 = True
lam2, lam3 = mul(LAM, LAM), powe(LAM, 3)
for a2 in (-1, 0, 1):
    for a3 in (-1, 0, 1):
        d = add(add(ONE, mul((a2, 0), lam2)), mul((a3, 0), lam3)); key = (d[0] % 9, d[1] % 9)
        if vals[key] != (-a3) % 3: ok15 = False
print("  agrees with the supplementary law DR (1.5): %s" % ok15)
if bad9 or smaller or not ok15: nfail += 1

def phiW(n):
    key = (n[0] % 9, n[1] % 9)
    return W3 ** vals[key] if key in vals else 0.0
def phiP(n):
    return 1.0 if (n[0] % 3 == 1 and n[1] % 3 == 0) else 0.0

# gamma_2 (003 normalisation) = (lam/pi)_3 g~(pi)
for elt in [(-2, -3), (4, 3), (-5, 0), (7, 3)]:
    pr = get_prime(elt)
    d = abs(gauss_gamma(pr, 2) - W3 ** pr.cub(LAM) * gauss_tilde(pr))
    print("  gamma_2^{003}(%s) - (lam/pi)_3 g~(pi) = %.1e" % (elt, d))
    if d > 1e-10: nfail += 1

# ---------------------------------------------------------------- B1
def fourier(phi, q):
    """phi-hat(h) = N(q)^-1 sum_{n mod q} phi(n) e^(-n h/(lam q)),  q rational integer"""
    out = {}
    for h0 in range(q):
        for h1 in range(q):
            h = (h0, h1); s = 0
            for n0 in range(q):
                for n1 in range(q):
                    f = phi((n0, n1))
                    if f != 0: s += f * cmath.exp(-2j * cmath.pi * omcoef(mul((n0, n1), h)) / q)
            out[h] = s / q ** 2
    return out

lam2 = mul(LAM, LAM)
def expansion(phi, q, label):
    fh = fourier(phi, q); rows = []
    for h, c in fh.items():
        if abs(c) < 1e-13: continue
        a, cc = reduce_frac(mul(lam2, h), (q, 0))
        g, kexp, rep, t = K.translate_parts(a, cc)
        rows.append((h, c, a, cc, rep, kexp, t))
    print("  %s: period %d, nonzero phi-hat(h): %d of %d" % (label, q, len(rows), q * q))
    summ = {}
    for h, c, a, cc, rep, kexp, t in rows:
        key = (norm(cc), rep); summ.setdefault(key, [0, 0.0]); summ[key][0] += 1; summ[key][1] += abs(c) ** 2
    for key in sorted(summ):
        print("     N(c_0) = %2d  cusp expansion %-4s : %2d translates, sum |phi-hat|^2 = %.5f" % (key[0], key[1], summ[key][0], summ[key][1]))
    return rows

print("\n=== B1: finite Fourier expansion at 3; translates z + lam^2 h/q, reduced denominators and cusp types ===")
rowsP = expansion(phiP, 3, "phi_P (Patterson's own, primary indices only)")
rowsW = expansion(phiW, 9, "phi_W (paper's combination: extra twist (lam/n)_3)")

def dual_value(rows, z, v):
    return sum(c * np.conj(K.translate(a, cc, z, v)) for h, c, a, cc, rep, kexp, t in rows)

def dzbar_translate(a, cc, v):
    """d/dzbar [ conj theta(z + a/c, v) ] at z = 0, through the cusp expansion"""
    g, kexp, rep, t = K.translate_parts(a, cc)
    c = cx(cc); d = cx(g[1][1])
    zz = -d / c + cx(t); vv = 1.0 / (abs(c) ** 2 * v)
    if rep == "I": F = K.theta(zz, vv, deriv="zbar")
    elif rep == "P": F = W3 * K.theta(W3 ** 2 * zz, vv, deriv="zbar")
    elif rep == "P2": F = W3 ** 2 * K.theta(W3 * zz, vv, deriv="zbar")
    elif rep == "g10": F = K.F10(zz, vv, deriv="zbar")
    else: F = K.F19(zz, vv, deriv="zbar")
    return -1.0 / (c * c * v * v) * np.conj(W3 ** kexp * F)

def J_dual(rows, v): return sum(c * dzbar_translate(a, cc, v) for h, c, a, cc, rep, kexp, t in rows)
def series_phi(phi):
    f = np.array([phi((int(n[0]), int(n[1]))) for n in thb.num], dtype=complex)
    return thb.twisted(f, const=0.0)
SP, SW = series_phi(phiP), series_phi(phiW)

# ---------------------------------------------------------------- B2
TD = 44.0     # truncation for derivative sums
def dzbar_T(a, cc, v):
    g, kexp, rep, t = K.translate_parts(a, cc)
    c = cx(cc); d = cx(g[1][1])
    zz = -d / c + cx(t); vv = 1.0 / (abs(c) ** 2 * v)
    if rep == "I": F = K.theta(zz, vv, T=TD, deriv="zbar")
    elif rep == "P": F = W3 * K.theta(W3 ** 2 * zz, vv, T=TD, deriv="zbar")
    elif rep == "P2": F = W3 ** 2 * K.theta(W3 * zz, vv, T=TD, deriv="zbar")
    elif rep == "g10": F = K.F10(zz, vv, T=TD, deriv="zbar")
    else: F = K.F19(zz, vv, T=TD, deriv="zbar")
    return -1.0 / (c * c * v * v) * np.conj(W3 ** kexp * F)
def J_dual(rows, v, conj_phihat=False):
    return sum((np.conj(c) if conj_phihat else c) * dzbar_T(a, cc, v) for h, c, a, cc, rep, kexp, t in rows)

print("\n=== B2: pointwise identity, direct twisted series = sum of translates through cusp expansions ===")
print("  Kubota exponents of the translates: phi_P %s ; phi_W %s" % (sorted(set(r[5] for r in rowsP)), sorted(set(r[5] for r in rowsW))))
print("  function values: error relative to the value.  d/dzbar at z=0: the value is small by cancellation, so the error is")
print("  measured against the sum of the absolute values of the direct terms (the first run used the value itself as the")
print("  denominator and flagged |diff| = 2.5e-13 on a value 1.0e-8 as a failure; that was a defect of my criterion).")
for label, S, rows in (("phi_P", SP, rowsP), ("phi_W", SW, rowsW)):
    for i in range(3):
        v = random.uniform(0.05, 0.12); z = complex(random.uniform(-.05, .05), random.uniform(-.05, .05))
        report("%s function at w=(%.3f%+.3fi, %.3f)" % (label, z.real, z.imag, v), S(z, v), dual_value(rows, z, v))
        report("%s d/dzbar at z=0, v=%.3f" % (label, v), S(0.0, v, T=TD, deriv="zbar"), J_dual(rows, v), scale=S.abs_sum(v, T=TD, deriv="zbar"))
v = 0.08
report("phi_W d/dzbar with conj(phi-hat)  [= twist by conj (lam/n)]", SW(0.0, v, T=TD, deriv="zbar"), J_dual(rowsW, v, conj_phihat=True), expect_pass=False, scale=SW.abs_sum(v, T=TD, deriv="zbar"))
report("phi_W series against phi_P translates", SW(0.0, v, T=TD, deriv="zbar"), J_dual(rowsP, v), expect_pass=False, scale=SW.abs_sum(v, T=TD, deriv="zbar"))
z = 0.02 + 0.03j
K.conj_kappa = True
report("phi_P function with conj(chi_Kubota)", SP(z, v), dual_value(rowsP, z, v), expect_pass=(set(r[5] for r in rowsP) == {0}))
K.conj_kappa = False

# ---------------------------------------------------------------- B3
print("\n=== B3: decay with the Bessel kernel.  J(v) = d/dzbar Theta_phi(0,v) = 3^{5/2} 2 pi i sum_A a(A) conj(mu_A) v K_{1/3}(4 pi |mu_A| v) ===")
print("  N_eff = 27/(4 pi v)^2 is the norm N(A) at which 4 pi |mu_A| v = 1.  'abs' = sum of absolute values of the direct terms.")
print("  The untwisted column is Theta_phi(0,v) (no angular factor): it grows, the angular sums collapse.")
print("        v     N_eff     abs        |J_P| direct   |J_P| dual     |J_W| direct   |J_W| dual     err_P/abs  err_W/abs   |Theta_P(0,v)|  |Theta_W(0,v)|")
tabJ = []
for v in (0.30, 0.20, 0.14, 0.10, 0.075, 0.060, 0.050, 0.0425):
    jp, jpd, jw, jwd = SP(0.0, v, T=TD, deriv="zbar"), J_dual(rowsP, v), SW(0.0, v, T=TD, deriv="zbar"), J_dual(rowsW, v)
    sc = SP.abs_sum(v, T=TD, deriv="zbar")
    ep, ew = abs(jp - jpd) / sc, abs(jw - jwd) / sc
    if max(ep, ew) > TOL: nfail += 1
    worst = max(worst, ep, ew); tabJ.append((v, abs(jpd), abs(jwd)))
    print("  %.4f  %7.1f  %.3e   %.6e   %.6e   %.6e   %.6e   %.1e    %.1e     %.4e      %.4e" % (v, 27 / (4 * math.pi * v) ** 2, sc, abs(jp), abs(jpd), abs(jw), abs(jwd), ep, ew, abs(SP(0.0, v)), abs(SW(0.0, v))))
print("  empirical exponent  -d log|J| / d(1/v)  between successive rows (dual values), compared with 4 pi |ell_min| / N(c_0):")
for (v1, p1, w1), (v2, p2, w2) in zip(tabJ[:-1], tabJ[1:]):
    print("     v %.4f -> %.4f :  P %.3f   W %.3f" % (v1, v2, -math.log(p2 / p1) / (1 / v2 - 1 / v1), -math.log(w2 / w1) / (1 / v2 - 1 / v1)))

# ---------------------------------------------------------------- B4
print("\n=== B4: log-Gaussian weights W(t) = exp(-(log t)^2/(4w)) : dual formula against direct sums ===")
A = K.td.A; NA = (A[:, 0] ** 2 - A[:, 0] * A[:, 1] + A[:, 1] ** 2).astype(float)
Ac = A[:, 0] + A[:, 1] * W3; alb = np.conj(Ac / np.abs(Ac))
aP = np.conj(K.td.G)                               # g~(c)|d|
aW = aP * W3 ** K.td.kl                            # (lam/c)_3 g~(c)|d| = gamma_2(c)|d|
def T_direct(a, X, w): return np.sum(a * alb * np.exp(-(np.log(NA / X)) ** 2 / (4 * w))) / math.sqrt(X)

def dual_terms(rows):
    """J(v) = sum_k e_k v^-3 K_{1/3}(beta_k / v).  Terms are aggregated over the translates by the frequency ell.
       returns (E_k = e_k beta_k^-3 summed over translates, beta_k^2, N(ell_k), N(c))"""
    agg = {}
    for h, c, a, cc, rep, kexp, t in rows:
        g, _, _, _ = K.translate_parts(a, cc)
        cc_c = cx(cc); d = cx(g[1][1]); Nc = norm(cc); shift = -d / cc_c + cx(t)
        if rep in ("I", "P", "P2"):
            S = K.theta; rot = {"I": ONE, "P": OMEGA2, "P2": OMEGA}[rep]
        else:
            S = K.F10 if rep == "g10" else K.F19; rot = ONE
        num = np.stack([S.num[:, 0] * rot[0] - S.num[:, 1] * rot[1], S.num[:, 0] * rot[1] + S.num[:, 1] * rot[0] - S.num[:, 1] * rot[1]], axis=1)
        ell = S.mu * cx(rot)
        Fz = S.coef * 2j * math.pi * np.conj(ell) * np.exp(4j * math.pi * (ell * shift).real)
        ek = c * (-1.0 / (cc_c ** 2)) * np.conj(W3 ** kexp * Fz) / Nc
        beta = 4 * math.pi * np.abs(ell) / Nc
        key = (Nc, S.den)
        d_ = agg.setdefault(key, {})
        for n0, n1, val in zip(num[:, 0].tolist(), num[:, 1].tolist(), (ek / beta ** 3).tolist()):
            d_[(n0, n1)] = d_.get((n0, n1), 0.0) + val
    E, B2, NL, NC = [], [], [], []
    for (Nc, den), d_ in agg.items():
        for (n0, n1), val in d_.items():
            if abs(val) < 1e-13: continue
            nl = (n0 * n0 - n0 * n1 + n1 * n1) / 3.0 ** den
            E.append(val); NL.append(nl); NC.append(Nc); B2.append((4 * math.pi) ** 2 * nl / Nc ** 2)
    return np.array(E), np.array(B2), np.array(NL), np.array(NC)

def H_kernel(u, w, sigma=-0.5):
    """H(u) = (1/2 pi i) int_(sigma) W^(s) u^s R(s) ds,
       R(s) = 2^{1-2s} Gamma(5/3-s) Gamma(4/3-s) * 4 * 3^{-5/2} ((2 pi)^2/27)^s / (Gamma(s+1/3) Gamma(s+2/3))"""
    Tm = math.sqrt(75.0 / w); n = int(2 * Tm / 0.004) + 1
    tau = np.linspace(-Tm, Tm, n); s = sigma + 1j * tau; dt = tau[1] - tau[0]
    logR = (1 - 2 * s) * math.log(2) + loggamma(5.0 / 3 - s) + loggamma(4.0 / 3 - s) + math.log(4 * 3 ** -2.5) \
           + s * math.log((2 * math.pi) ** 2 / 27) - loggamma(s + 1.0 / 3) - loggamma(s + 2.0 / 3)
    base = np.exp(0.5 * math.log(4 * math.pi * w) + w * s * s + logR)
    out = np.empty(len(u), dtype=complex); lu = np.log(u)
    for i0 in range(0, len(u), 2000):
        blk = lu[i0:i0 + 2000]
        out[i0:i0 + 2000] = (np.exp(np.outer(blk, s)) * base[None, :]).sum(axis=1) * dt / (2 * math.pi)
    return out

DT = {"P": dual_terms(rowsP), "W": dual_terms(rowsW)}
TABLE_N = {"P": min(K.theta.amu[-1], K.F10.amu[-1]) ** 2, "W": K.theta.amu[-1] ** 2}     # the dual tables are complete up to this N(ell)
for lab in ("P", "W"):
    E, B2, NL, NC = DT[lab]
    i = np.argmin(NL)
    print("  dual support of %s after aggregation over translates: %d frequencies, N(c) in %s, smallest N(ell) = %.6f (= 1/%.0f),"
          % (lab, len(E), sorted(set(NC.tolist())), NL[i], 1 / NL[i]))
    print("      so x = (2 pi)^4 N(ell) X / (27 N(c)^2) >= X / %.3f" % (27 * NC[i] ** 2 / ((2 * math.pi) ** 4 * NL[i])))
def T_dual(lab, X, w, xcut):
    E, B2, NL, NC = DT[lab]
    u = X * B2; x = u * math.pi ** 2 / 27
    sel = x <= xcut
    need = xcut * 27 * NC.max() ** 2 / ((2 * math.pi) ** 4 * X)
    if need > TABLE_N[lab]: return None, 0, float(x.min())
    uu, inv = np.unique(np.round(u[sel], 9), return_inverse=True)
    Hh = H_kernel(uu, w)
    return np.sum(E[sel] * Hh[inv]) / math.sqrt(X), int(sel.sum()), float(x.min())

# envelope of the dual kernel
print("  size of the dual kernel |H|(x), x = (2 pi)^4 N(ell) X / (27 N(c)^2):")
xs = np.array([1, 10, 100, 1e3, 1e4, 1e5, 1e6])
for w in (0.05, 0.1, 1.0):
    Hh = np.abs(H_kernel(xs * 27 / math.pi ** 2, w))
    print("    w=%.2f : " % w + "  ".join("x=%g: %.1e" % (x, h) for x, h in zip(xs, Hh)))

old = {0.1: {"W": {125: 3.415e-01, 250: 3.445e-02, 500: 1.234e-01, 1000: 8.578e-02, 2000: 2.891e-02},
             "P": {125: 8.162e-02, 250: 1.658e-02, 500: 2.310e-02, 1000: 6.368e-03, 2000: 2.441e-03}},
       0.05: {"W": {125: 3.562e-01, 250: 5.446e-03, 500: 2.260e-01, 1000: 1.726e-01, 2000: 3.350e-03, 4000: 3.342e-02, 8000: 3.968e-02},
              "P": {125: 2.796e-02, 250: 6.957e-02, 500: 2.181e-03, 1000: 2.863e-02, 2000: 1.153e-02, 4000: 3.297e-02, 8000: 1.006e-02}}}
for w, xcut in ((0.1, 3.0e5), (0.05, 1.2e6)):
    print("  w = %.2f   (dual terms kept: x <= %.1e)" % (w, xcut))
    print("       X    combo   T direct (re, im)               T dual (re, im)                 |diff|    |T|        old script |T|   #dual freq   x_min")
    for X in (125, 250, 500, 1000, 2000, 4000, 8000, 16000):
        trunc = math.exp(-(math.log(NMAX / X)) ** 2 / (4 * w))
        if trunc > 1e-11: continue
        for lab, a in (("P", aP), ("W", aW)):
            td_ = T_direct(a, X, w)
            tdu, nterms, xmin = T_dual(lab, X, w, xcut)
            o = old[w][lab].get(X)
            if tdu is None:
                print("  %6d    %s     % .8e % .8e   (dual table too short for this X; skipped)       %.4e   %s" % (X, lab, td_.real, td_.imag, abs(td_), ("%.3e" % o) if o else "")); continue
            diff = abs(td_ - tdu)
            if diff > 1e-8: nfail += 1
            worst = max(worst, diff)
            print("  %6d    %s     % .8e % .8e   % .8e % .8e   %.1e   %.4e   %s   %7d   %.3g" % (X, lab, td_.real, td_.imag, tdu.real, tdu.imag, diff, abs(td_),
                  ("%.3e" % o) if o else "    -    ", nterms, xmin))
    print("  beyond the range of the direct sums (N(A) <= %d): values from the dual formula only (a prediction, not a direct check)" % NMAX)
    print("       X        |T_P| dual     |T_W| dual     x_min(P)   x_min(W)")
    for X in (3.2e4, 1e5, 3.2e5, 1e6, 3.2e6, 1e7):
        tp, _, xp = T_dual("P", X, w, xcut); tw_, _, xw = T_dual("W", X, w, xcut)
        print("  %9.0f   %.3e      %.3e      %.3g    %.3g" % (X, abs(tp), abs(tw_), xp, xw))

# ---------------------------------------------------------------- B5
print("\n=== B5 (rung P3): phi_W combined with a moving sextic twist chi_p^j, brute-force layer ===")
def sext_tab(pr, j): return sextic_table(pr, j)
for elt, j in (((-2, -3), 1), ((4, 3), 1), ((-2, -3), 5), ((-5, 0), 1)):
    pr = get_prime(elt); P = pr.rat; tab = sextic_table(pr, j)
    S = SW.twisted(tab[thb.num[:, 0] % P, thb.num[:, 1] % P])
    res = [(a, b) for a in range(P) for b in range(P)] if pr.inert else [(a, 0) for a in range(P)]
    def chi(y):
        k = pr.sext(y); return 0 if k is None else Z6 ** ((k * j) % 6)
    Cj = {h: sum(chi(y) * cmath.exp(-2j * cmath.pi * omcoef(mul(mul(h, y), conj(elt))) / pr.N) for y in res) / pr.N for h in res}
    a_abs = 3 * abs(cx(elt))
    z = complex(random.uniform(-.3, .3), random.uniform(-.3, .3)) / a_abs; v = random.uniform(0.7, 0.9) / a_abs
    tot = 0; ntr = 0
    for h0, c0, _, _, _, _, _ in rowsW:
        for hp, cp in Cj.items():
            if abs(cp) < 1e-14: continue
            # shift lam^2 (h0/9 + hp/p) = lam^2 (h0 p + 9 hp) / (9 p)
            a, cc = reduce_frac(mul(lam2, add(mul(h0, elt), mul((9, 0), hp))), mul((9, 0), elt))
            tot += c0 * cp * np.conj(K.translate(a, cc, z, v)); ntr += 1
    report("p=%s j=%d, %d translates, v=%.4f" % (elt, j, ntr, v), S(z, v), tot)

print("\nworst error among identities expected to hold: %.2e" % worst)
print("TASK B checks: %s  (failures: %d)" % ("PASS" if nfail == 0 else "FAIL", nfail))
