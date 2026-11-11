# -*- coding: utf-8 -*-
"""
Stage 5, common code for the end-to-end comparison (layer E of stage 3) with an arbitrary weight object.
The direct sum, the groups, the constant C and the dual sum are the SAME formulas as in stage3_d_mellin.py (copied, since
that script cannot be imported without running it); only the weight is a parameter.  Importing runs nothing.
  direct: X^-1/2 sum_A conj(alpha(A)) gamma_2(n)|b| Psi(A) W(N(A)/X)                                   (eq:T, :1287)
  dual:   sum_groups C sum_l d(l) alpha(l) N(l)^-1/2 psi(lambda^4 l) prod B_{p,j}(lambda^4 l) V^sharp(N(l) X/N(c)^2)   (:1806, :3420, :1752)
"""
import sys, os, math, itertools, time
import numpy as np
from stage2_lib import *
import stage4_bounds as sb

LINT = 18

class WLogGauss:
    """W(y) = exp(-(log y)^2/(4w)); everything from stage4_bounds"""
    def __init__(self, w, xmax):
        self.w = w; self.xmax = xmax; self.name = "logGauss w=%g" % w
    def W(self, y): return np.exp(-(np.log(y)) ** 2 / (4 * self.w))
    def vsharp(self, x, invert=False):
        if not invert: return sb.vsharp(x, self.w)
        saved = sb.logQ; sb.logQ = lambda t: -saved(t)
        try: return sb.vsharp(x, self.w)
        finally: sb.logQ = saved
    def quad_bound(self, x): return sum(sb.quad_error_bound(x, self.w))
    def dual_tail(self, X, Nc, R, Cabs, Bmax): return sb.dual_tail_bound(X, self.w, Nc, R, Cabs, Bmax)[0]
    def direct_tail(self, X, Nmax): return sb.direct_tail_bound(X, self.w, Nmax)

class WBump:
    """compactly supported weight of stage5_bump.Bump"""
    def __init__(self, bump, xmax):
        self.b = bump; self.xmax = xmax; self.name = "bump p=%d a=%g" % (bump.p, bump.a)
    def W(self, y): return self.b.W(y)
    def vsharp(self, x, invert=False): return self.b.vsharp(x, invert=invert)
    def quad_bound(self, x): return sum(self.b.quad_error_bound(x))
    def dual_tail(self, X, Nc, R, Cabs, Bmax): return self.b.dual_tail_bound(X, Nc, R, Cabs, Bmax)
    def direct_tail(self, X, Nmax):
        assert X * math.exp(self.b.a) <= Nmax, "support of W(N/X) exceeds the table"
        return 0.0


class EndToEnd:
    def __init__(self, Nmax=1000000):
        self.th = Theta(Nmax); self.Nmax = Nmax
        self.setup = Setup((LINT, 0))
        self.PH = {nm: phi_hat(phi_table(nm, LINT), LINT) for nm in "ABCDEF"}
        self.PHno2 = {nm: phi_hat(phi_table(nm, LINT, with2=False), LINT) for nm in "ABCF"}
        td = self.th.K.td
        self.A = td.A; A_ = td.A
        self.NA = (A_[:, 0] ** 2 - A_[:, 0] * A_[:, 1] + A_[:, 1] ** 2).astype(float)
        Ac = A_[:, 0] + A_[:, 1] * W3; self.albar = np.conj(Ac / np.abs(Ac))
        self.aW = np.conj(td.G) * W3 ** td.kl
        self.P0T = {}
        for nm in "ABCDEF":
            tab = np.zeros((LINT, LINT), dtype=complex)
            for x0 in range(LINT):
                for x1 in range(LINT):
                    n = (x0, x1)
                    if is_primary(n) and not (x0 % 2 == 0 and x1 % 2 == 0): tab[x0, x1] = W3 ** psi0_exp(nm, n)
            self.P0T[nm] = tab
        self._VS = {}

    def direct(self, nm, jd, wt, X):
        A_ = self.A
        f = self.aW * self.albar * self.P0T[nm][A_[:, 0] % LINT, A_[:, 1] % LINT]
        for p, j in jd.items():
            lp = loc(p); f = f * lp.chitab(j)[lp.idx_arr(A_)]
        t = f * wt.W(self.NA / X) / math.sqrt(X)
        return complex(np.sum(t)), float(np.sum(np.abs(t)))

    def groups(self, nm, jd, wt, X, ph=None):
        th = self.th; ph = self.PH[nm] if ph is None else ph
        primes = list(jd); must = [p for p in primes if jd[p] != 0]; opt = [p for p in primes if jd[p] == 0]
        out = []
        for h0, f0 in ph.items():
            if abs(f0) < 1e-13: continue
            for kopt in range(len(opt) + 1):
                for extra in itertools.combinations(opt, kopt):
                    act = [p for p in primes if p in must or p in extra]; inact = [p for p in primes if p not in act]
                    cfg = self.setup.build(h0, act, [ONE] * len(act)); Nc = norm(cfg["c"]); sg = cfg["sigma"]; s = th.d[sg]
                    R = min(th.complete[sg] ** 2, wt.xmax * Nc * Nc / X)
                    n = int(np.searchsorted(s.aell, math.sqrt(R), side="right")); m = s.m[:n]
                    arr = s.coef[:n] * psi_char(cfg)(m) / np.conj(s.ell[:n]); Bmax = 1.0
                    for p in act:
                        lp = loc(p); j = jd[p]; arr = arr * lp.Btab(j)[lp.idx_arr(m)]
                        Bmax *= {4: (lp.N - 1) / math.sqrt(lp.N), 0: 1 / math.sqrt(lp.N)}.get(j, 1.0)
                    out.append(dict(cfg=cfg, Cp=Cprime(cfg, f0, jd, inact), key=(sg, Nc, n), arr=arr, R=R, Bmax=Bmax,
                                    k0=kappa_formula(cfg)[1], ac=alpha(cfg["c"])))
        return out

    def _vs(self, wt, Nm, Nc, X, variant):
        key = (wt.name, Nc, X, variant, len(Nm), int(Nm[-1]) if len(Nm) else 0)
        if key not in self._VS:
            x = Nm / 81.0 * X / (Nc * Nc if variant != "NE4" else Nc)
            self._VS[key] = wt.vsharp(x, invert=(variant == "NE7"))
        return self._VS[key]

    def dual(self, G, wt, X, variant=None):
        th = self.th; agg = {}
        for g in G:
            C = -1j / 81 * np.conj(g["ac"]) ** 2 * g["Cp"]
            if variant == "NE1": C = -C
            if variant == "NE2": C = -1j / 81 * g["ac"] ** 2 * g["Cp"]
            if variant == "NE6": C = C * W3 ** ((2 * g["k0"]) % 3)
            arr = g["arr"]
            if variant == "NE3":
                el = th.d[g["key"][0]].ell[:g["key"][2]]; arr = arr * (np.conj(el) / el)
            if g["key"] in agg: agg[g["key"]][0] += C * arr; agg[g["key"]][1] += abs(C) * np.abs(arr)
            else: agg[g["key"]] = [C * arr, abs(C) * np.abs(arr)]
        tot = 0; sabs = 0.0; cabs = 0.0; qb = 0.0
        for (sg, Nc, n), (E, Ea) in agg.items():
            m = th.d[sg].m[:n]; Nm = m[:, 0] ** 2 - m[:, 0] * m[:, 1] + m[:, 1] ** 2
            un, inv = np.unique(Nm, return_inverse=True)
            V = self._vs(wt, un, Nc, X, variant if variant in ("NE4", "NE7") else None)[inv]
            tot += complex(np.sum(E * V)); sabs += float(np.sum(Ea * np.abs(V))); cabs += float(np.sum(Ea))
            if len(un):
                xs = un[[0, -1]] / 81.0 * X / (Nc * Nc)
                qb += float(np.sum(Ea)) * max(wt.quad_bound(float(x)) for x in xs)
        return tot, sabs, cabs, qb

    def run(self, nm, jd, wt, X):
        """returns dict with T_dir, T_dual, bounds, groups"""
        Td, Sd = self.direct(nm, jd, wt, X)
        G = self.groups(nm, jd, wt, X)
        Tu, Su, Cu, qb = self.dual(G, wt, X)
        b_dir = wt.direct_tail(X, self.Nmax)
        b_dual = sum(wt.dual_tail(X, g["key"][1], g["R"], abs(g["Cp"]) / 81, g["Bmax"]) for g in G)
        b_round = 1e-13 * (Sd + Su) + 2e-15 * Cu
        return dict(Td=Td, Tu=Tu, Sd=Sd, Su=Su, G=G, b_dir=b_dir, b_dual=b_dual, qb=qb, b_round=b_round,
                    Btot=b_dir + b_dual + qb + b_round)

    def control(self, k, nm, jd, wt, X, G):
        """T_ctl for control k, or None if the control is dead in this run"""
        lives = {"NE1": True, "NE2": any(abs(g["ac"] ** 4 - 1) > 1e-9 for g in G), "NE3": True, "NE4": any(g["key"][1] != 1 for g in G),
                 "NE5": nm in self.PHno2, "NE6": any(g["k0"] != 0 for g in G), "NE7": True}
        if not lives[k]: return None
        if k == "NE5": return self.dual(self.groups(nm, jd, wt, X, self.PHno2[nm]), wt, X)[0]
        return self.dual(G, wt, X, variant=k)[0]
