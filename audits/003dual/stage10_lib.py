# -*- coding: utf-8 -*-
"""
003dual stage 10 library (lead, Claude Opus 5.5, 2026-11-13): extend the stage 7 solve (stage7_lib.py) to the
cusp functions  F_10(w) = theta(G10 w),  F_19(w) = theta(G19 w)  (G10 = [[1,0],[omega,1]], G19 = [[1,0],[omega^2,1]],
kubota.COSET_REPS).  Their coefficients d_10, d_19 become unknowns too; no DR (5.13)/(5.14)/App. A input.

Index sets are derived from the group alone (kubota.decompose, i.e. DR (5.4)):
  decompose(R T_t) = (k(t), R, 0) for t in 3O  =>  F_R(w + t) = omega^{k(t)} F_R(w),
  so d_R(mu) can be nonzero only where 2 Re(mu t) = k(t)/3 mod 1 for t = 3 and t = 3 omega, and F_R has no
  constant term (nontrivial character).  We enumerate mu in lambda^{-6} O (a superset) and keep that coset.
Equations: for X in {I, g10, g19} and a sample point w = (z, v0):
  F_X(w) = theta(R_X w) = theta(g^{-1} w''),  w'' = g R_X w,  decompose(g^{-1}) = (k, rep, t)
         = omega^k F_rep(w'' + t)   with F_P(w) = theta(omega^2 z, v), F_P2 = theta(omega z, v)  (as stage 7).
  Rows use only g with height(w'') >= gain * v0; the identity row (w'' = w) is excluded by gain > 1.
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage7_lib as S7
kb = S7.kb
W3 = kb.W3; SIGMA = S7.SIGMA
LAM2 = kb.mul(kb.LAM, kb.LAM)
FAM = {"I": "th", "P": "th", "P2": "th", "g10": "g10", "g19": "g19"}


def char_of(rep):
    """k(3), k(3 omega) for F_rep (group only)"""
    out = []
    for t in ((3, 0), (0, 3)):
        k, r, tp = kb.decompose(kb.mat_mul(kb.COSET_REPS[rep], kb.Tmat(t)))
        assert r == rep and tp == (0, 0)
        out.append(k % 3)
    return tuple(out)


def cusp_indices(rep, M):
    """numerators m = lambda^6 mu (pairs), mu in the character coset, |mu| <= M, sorted by |mu|"""
    k3, k3w = char_of(rep)
    L6 = kb.LAMC ** 6
    R = int(M * abs(L6) * 2 / math.sqrt(3)) + 3
    out = []
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            mu = (a + b * W3) / L6
            if abs(mu) > M * (1 + 1e-12) or (a == 0 and b == 0): continue
            s1 = 2 * (3 * mu).real - k3 / 3; s2 = 2 * (3 * W3 * mu).real - k3w / 3
            if abs(s1 - round(s1)) < 1e-9 and abs(s2 - round(s2)) < 1e-9:
                out.append((a, b))
    out.sort(key=lambda m: (abs((m[0] + m[1] * W3) / L6), m))
    return out


def mu6(nums):
    nums = np.asarray(nums); return (nums[:, 0] + nums[:, 1] * W3) / kb.LAMC ** 6


def reduce_pt(z, v):
    """exact h0 in SL2(O) moving (z, v) high (translate to the O-parallelogram, invert while |w| < 1)"""
    h = kb.I2
    for _ in range(200):
        # translate: z -> z - n, n nearest in O
        y1 = round(z.imag / (math.sqrt(3) / 2)); x0 = round(z.real + y1 / 2.0)
        n = (x0, y1); z = z - kb.cx(n); h = kb.mat_mul(kb.Tmat(kb.neg(n)), h)
        if abs(z) ** 2 + v * v >= 1 - 1e-12: return h, z, v
        h = kb.mat_mul(kb.EMAT, h)
        z, v = kb.mobius(kb.mat_cx(kb.EMAT), z, v)
    raise RuntimeError("reduction did not stop")


def candidates(z, v, vmin, Cmax=60):
    """list of g in SL2(O) with height(g (z, v)) >= vmin (via reduction, then raise/lower around the top point)"""
    h0, z0, v0 = reduce_pt(z, v)
    out = []
    if v0 >= vmin: out.append(h0)
    # other top points: h1 (z0, v0) with height >= vmin, c != 0
    R = int(math.sqrt(4 * Cmax / 3)) + 2
    for c0 in range(-R, R + 1):
        for c1 in range(-R, R + 1):
            c = (c0, c1); Nc = kb.norm(c)
            if Nc == 0 or Nc > Cmax: continue
            rad2 = v0 / vmin - Nc * v0 * v0
            if rad2 <= 0: continue
            cc = kb.cx(c); ctr = -cc * z0; rad = math.sqrt(rad2)
            y1 = ctr.imag / (math.sqrt(3) / 2)
            for d1 in range(math.floor(y1 - rad / 0.866) - 1, math.ceil(y1 + rad / 0.866) + 2):
                x0 = ctr.real + d1 / 2.0
                for d0 in range(math.floor(x0 - rad) - 1, math.ceil(x0 + rad) + 2):
                    d = (d0, d1)
                    if abs(cc * z0 + kb.cx(d)) ** 2 + Nc * v0 * v0 > v0 / vmin: continue
                    h1 = S7.complete(c, d)
                    if h1 is None: continue
                    out.append(kb.mat_mul(h1, h0))
    return out


_DEC = {}
def dec(g):
    if g not in _DEC: _DEC[g] = kb.decompose(g)
    return _DEC[g]


def sample_points(X, n_side, shift=(0.0, 0.0)):
    if X == "I":
        return S7.sample_points(n_side, None, shift)
    # period lattice 3O: parallelogram 3 x 3 omega, centred near the cusp -1/omega (g10) or -1/omega^2 (g19)
    ctr = -1 / W3 if X == "g10" else -1 / W3 ** 2
    pts = []
    for i in range(n_side):
        for j in range(n_side):
            x = (i + 0.5 + shift[0]) / n_side - 0.5; y = (j + 0.5 + shift[1]) / n_side - 0.5
            pts.append(ctr + 3 * x + 3 * W3 * y)
    return pts


class System:
    def __init__(self, M_th, M_cu, v0):
        self.v0 = v0
        self.nums = {"th": S7.dual_indices(M_th), "g10": cusp_indices("g10", M_cu), "g19": cusp_indices("g19", M_cu)}
        self.mu = {"th": S7.mu_of(self.nums["th"]), "g10": mu6(self.nums["g10"]), "g19": mu6(self.nums["g19"])}
        self.off = {}; o = 0
        for f in ("th", "g10", "g19"): self.off[f] = o; o += len(self.nums[f])
        self.n = o
        self.scale = np.concatenate([v0 * S7.kv(1.0 / 3, 4 * math.pi * np.abs(self.mu[f]) * v0) for f in ("th", "g10", "g19")])

    def block(self, fam, z, v):
        row = np.zeros(self.n, dtype=complex)
        o = self.off[fam]; row[o:o + len(self.nums[fam])] = S7.phi_row(self.mu[fam], z, v)
        return row

    def rows(self, X, pts, per_point=4, gain=1.05, info=None):
        A, b = [], []
        RX = kb.COSET_REPS[X]
        for z in pts:
            zp, vp = kb.mobius(kb.mat_cx(RX), z, self.v0)
            cands = []
            for g in candidates(zp, vp, gain * self.v0):
                z2, v2 = kb.mobius(kb.mat_cx(g), zp, vp)
                k, rep, t = dec(kb.mat_inv(g))
                cands.append((v2, k % 3, rep, S7.ROT.get(rep, 1.0) * (z2 + kb.cx(t))))
            if not cands:
                if info is not None: info["skipped"] += 1
                continue
            cands.sort(key=lambda c: -c[0])
            chosen = [cands[0]]
            for c in cands[1:]:
                if len(chosen) >= per_point: break
                if all(c[1] != x[1] or c[2] != x[2] for x in chosen): chosen.append(c)
            for v2, k, rep, z2 in chosen:
                f = W3 ** k
                row = self.block(FAM[X], z, self.v0) - f * self.block(FAM[rep], z2, v2)
                rhs = (f * SIGMA * v2 ** (2.0 / 3) if FAM[rep] == "th" else 0) - (SIGMA * self.v0 ** (2.0 / 3) if X == "I" else 0)
                A.append(row / self.scale); b.append(rhs)
                if info is not None:
                    info["chi"][k] += 1; info["rep"][(X, rep)] = info["rep"].get((X, rep), 0) + 1
        return A, b

    def build(self, n_th, n_cu, shift=(0.0, 0.0), per_point=4, pp_cu=None):
        info = {"chi": [0, 0, 0], "rep": {}, "skipped": 0}
        A, b = [], []
        for X, n, pp in (("I", n_th, per_point), ("g10", n_cu, pp_cu or per_point), ("g19", n_cu, pp_cu or per_point)):
            a_, b_ = self.rows(X, sample_points(X, n, shift), pp, info=info); A += a_; b += b_
        return np.array(A), np.array(b), info

    def split(self, x):
        """x in scaled variables -> dict fam -> {num: coefficient}"""
        c = x / self.scale
        return {f: dict(zip(self.nums[f], c[self.off[f]:self.off[f] + len(self.nums[f])])) for f in ("th", "g10", "g19")}


def table_vector(sysm, Mmax=9.0):
    """DR table (kubota.build_series) on the unknowns of sysm -- for code checks and the final comparison only"""
    td = kb.ThetaData(200000, os.path.join(here, "cache"), verbose=False)
    th, F10, F19 = kb.build_series(td, Mmax)
    out = np.zeros(sysm.n, dtype=complex); extra = {}
    for fam, ser, conv in (("th", th, lambda m: m), ("g10", F10, lambda m: kb.mul(m, LAM2)), ("g19", F19, lambda m: kb.mul(m, LAM2))):
        d = {}
        for m, cf in zip(map(tuple, ser.num), ser.coef):
            mm = conv(m); d[mm] = d.get(mm, 0) + cf
        idx = {m: i for i, m in enumerate(sysm.nums[fam])}
        Mf = np.abs(sysm.mu[fam]).max()
        miss = 0
        for mm, cf in d.items():
            if mm in idx: out[sysm.off[fam] + idx[mm]] = cf
            else:
                amu = abs((mm[0] + mm[1] * W3) / kb.LAMC ** (3 if fam == "th" else 6))
                if amu <= Mf * (1 - 1e-9) and abs(cf) > 0: miss += 1
        extra[fam] = miss          # table entries inside the range but outside the group-derived index set
    return out, extra
