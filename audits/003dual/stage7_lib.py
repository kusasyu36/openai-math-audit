# -*- coding: utf-8 -*-
"""
003dual stage 7 library (lead, Claude Opus 5.5, 2026-11-12): solve the Fourier coefficients of the cubic theta
function from automorphy alone (Hejhal-type linear system), WITHOUT using the coefficient formula DR (5.7)-(5.8).

Inputs used (and only these):
  * the constant term sigma = 3^{5/2}/2           (paper2.tex:1650-1656 = DR (5.6)/(sigma))
  * the action of SL2(C) on H^3                    (DR (5.1), kubota.mobius)
  * the multiplier: theta(g w) = chi(g1) theta(rep(w + t)) for g = g1 sigma rep T_t   (DR (5.4), kubota.decompose)
  * for rep in {I, P, P2}: F_I = theta, F_P(w) = theta(omega^2 z, v), F_P2(w) = theta(omega z, v)
  * the period lattice L = Z + 3 omega Z (T_1 in SL2(Z), T_{3 omega} = I mod 3)
NOT used: kubota.ThetaData / build_series (the DR coefficient table) -- except in the comparison step, which is
done in a separate script after the solve parameters were frozen.

Unknowns: tau(mu) for every mu in the dual lattice L* (inside lambda^{-3} O), 0 < |mu| <= M.  No support pattern,
no parity, no rotation symmetry is imposed.
Equation for a sample point w = (z, v0) and h in SL2(O) with w'' = h w, decompose(h^{-1}) = (k, rep, t):
   sum_mu tau(mu) [phi_mu(z, v0) - omega^k phi_mu(R_rep(z'' + t), v'')] = sigma (omega^k v''^{2/3} - v0^{2/3}),
   phi_mu(z, v) = v K_{1/3}(4 pi |mu| v) exp(2 pi i (mu z + conj(mu z))).
"""
import sys, os, math, cmath
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "029reflect"))
import numpy as np
from scipy.special import kv
import kubota as kb

W3 = kb.W3
SIGMA = 3 ** 2.5 / 2
LAMC3 = kb.LAMC ** 3
ROT = {"I": 1.0, "P": W3 ** 2, "P2": W3}


def dual_indices(M):
    """numerators m = lambda^3 mu (pairs) with mu in L*, 0 < |mu| <= M, sorted by |mu|"""
    R = int(M * 3 ** 1.5 * 2 / math.sqrt(3)) + 3
    out = []
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            if a == 0 and b == 0: continue
            mu = (a + b * W3) / LAMC3
            if abs(mu) > M * (1 + 1e-12): continue
            t1 = 2 * mu.real; t2 = 2 * (3 * W3 * mu).real
            if abs(t1 - round(t1)) < 1e-9 and abs(t2 - round(t2)) < 1e-9:
                out.append((a, b))
    out.sort(key=lambda m: (abs((m[0] + m[1] * W3) / LAMC3), m))
    return out


def mu_of(nums):
    nums = np.asarray(nums)
    return (nums[:, 0] + nums[:, 1] * W3) / LAMC3


def phi_row(mu, z, v):
    return v * kv(1.0 / 3, 4 * math.pi * np.abs(mu) * v) * np.exp(4j * math.pi * (mu * z).real)


def complete(c, d):
    """h = ((a, b), (c, d)) in SL2(O), or None if (c, d) is not coprime"""
    g0, s, t = kb.xgcd(d, c)                 # s d + t c = g0
    u = kb.exact_div(kb.ONE, g0)
    if u is None: return None
    a, b = kb.mul(s, u), kb.neg(kb.mul(t, u))  # a d - b c = 1
    h = ((a, b), (c, d))
    assert kb.mat_det(h) == kb.ONE
    return h


_DEC = {}
def dec_inv(h):
    if h not in _DEC: _DEC[h] = kb.decompose(kb.mat_inv(h))
    return _DEC[h]


def pullbacks(z, v, Cmax, reps=("I", "P", "P2"), gain=1.0 + 1e-6):
    """all h (c != 0, N(c) <= Cmax) raising the height above gain*v; returns list of
       (height'', kexp, rep, z_eval, h) with z_eval = R_rep(z'' + t)"""
    res = []; best_any = (v, "none")
    R = int(math.sqrt(4 * Cmax / 3)) + 2
    for c0 in range(-R, R + 1):
        for c1 in range(-R, R + 1):
            c = (c0, c1); Nc = kb.norm(c)
            if Nc == 0 or Nc > Cmax: continue
            cc = kb.cx(c)
            if Nc * v * v >= 1.0 / gain: continue   # fixed 11-12 (Codex w10 (9)): was ">= 1", wrong only for gain < 1 (stage 8)
            # d near -c z: need |c z + d|^2 < 1/gain*... - |c|^2 v^2  (height condition)
            rad2 = 1.0 / gain - Nc * v * v
            if rad2 <= 0: continue
            ctr = -cc * z
            rad = math.sqrt(rad2)
            # enumerate d = d0 + d1 omega with |d - ctr| <= rad
            y1 = ctr.imag / (math.sqrt(3) / 2)
            for d1 in range(math.floor(y1 - rad / 0.866) - 1, math.ceil(y1 + rad / 0.866) + 2):
                x0 = ctr.real + d1 / 2.0
                for d0 in range(math.floor(x0 - rad) - 1, math.ceil(x0 + rad) + 2):
                    d = (d0, d1); dc = kb.cx(d)
                    den = abs(cc * z + dc) ** 2 + Nc * v * v
                    if den >= 1.0 / gain: continue
                    h = complete(c, d)
                    if h is None: continue
                    k, rep, t = dec_inv(h)
                    z2, v2 = kb.mobius(kb.mat_cx(h), z, v)
                    if v2 > best_any[0]: best_any = (v2, rep)
                    if rep in reps:
                        res.append((v2, k % 3, rep, ROT[rep] * (z2 + kb.cx(t)), h))
    res.sort(key=lambda r: -r[0])
    return res, best_any


def sample_points(n_side, v0, shift=(0.0, 0.0)):
    """grid on the fundamental parallelogram of L = Z + 3 omega Z"""
    pts = []
    for i in range(n_side):
        for j in range(3 * n_side):
            x = (i + 0.5 + shift[0]) / n_side; y = (j + 0.5 + shift[1]) / (3 * n_side)
            pts.append(x + y * 3 * W3)
    return pts


def build_system(nums, v0, pts, Cmax, per_point=2, vmin_gain=1.05):
    """rows of the linear system; returns A (scaled columns), b, scale, info"""
    mu = mu_of(nums)
    scale = v0 * kv(1.0 / 3, 4 * math.pi * np.abs(mu) * v0)
    rows, rhs, info = [], [], {"chi": [0, 0, 0], "rep": {}, "skipped": 0, "skipped_rep": {}, "heights": []}
    for z in pts:
        cand, best_any = pullbacks(z, v0, Cmax)
        cand = [c for c in cand if c[0] >= vmin_gain * v0]
        if not cand:
            info["skipped"] += 1
            info["skipped_rep"][best_any[1]] = info["skipped_rep"].get(best_any[1], 0) + 1
            continue
        # take the highest, then the highest with a different chi exponent if any (to cover chi != 1)
        chosen = [cand[0]]
        for c in cand[1:]:
            if len(chosen) >= per_point: break
            if all(c[1] != x[1] or c[2] != x[2] for x in chosen): chosen.append(c)
        for v2, k, rep, z2, h in chosen:
            f = W3 ** k
            row = phi_row(mu, z, v0) - f * phi_row(mu, z2, v2)
            rows.append(row / scale); rhs.append(SIGMA * (f * v2 ** (2.0 / 3) - v0 ** (2.0 / 3)))
            info["chi"][k] += 1; info["rep"][rep] = info["rep"].get(rep, 0) + 1; info["heights"].append(v2)
    return np.array(rows), np.array(rhs), scale, info


def solve(A, b, scale, rcond=1e-14):
    x, res, rank, sv = np.linalg.lstsq(A, b, rcond=rcond)
    tau = x / scale
    resid = np.linalg.norm(A @ x - b) / max(np.linalg.norm(b), 1e-300)
    return tau, sv, rank, resid
