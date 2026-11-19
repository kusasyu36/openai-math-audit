# -*- coding: utf-8 -*-
"""Stage 12: control helpers shared by stage12_probe2 and stage12_main2 (copied from stage12_main.py, unchanged in
substance; stage12_main.py keeps its own copy so that its output stays reproducible).  Importing runs nothing."""
import numpy as np
from stage5_lib import *

def groups_jdual(E, nm, jd, jd_dual, wt, X):
    """E.groups with B_{p, jd_dual[p]} in the dual coefficients, everything else (C included) from jd"""
    th = E.th
    G = E.groups(nm, jd, wt, X)
    for g in G:
        cfg = g["cfg"]; sg, Nc, n = g["key"]; m = th.d[sg].m[:n]
        arr = g["arr"].copy()
        for p in jd:
            if p not in list(cfg["primes"]): continue
            lp = loc(p); i = lp.idx_arr(m)
            old = lp.Btab(jd[p])[i]; new = lp.Btab(jd_dual[p])[i]
            nz = old != 0
            arr[nz] = arr[nz] / old[nz] * new[nz]
            assert np.all(new[~nz] == 0)
        g["arr_orig"] = g["arr"]; g["arr"] = arr
    return G

def liveness(E, k, nm, jd, wt, X, G):
    """size of the change made by the switch of control k to its own ingredient (0 means dead)"""
    th = E.th
    if k == "NE1": return max(float(np.max(np.abs(2 * (-1j / 81) * np.conj(g["ac"]) ** 2 * g["Cp"]))) for g in G)
    if k == "NE2": return max(abs(-1j / 81 * (g["ac"] ** 2 - np.conj(g["ac"]) ** 2) * g["Cp"]) for g in G)
    if k == "NE6": return max(abs(-1j / 81 * np.conj(g["ac"]) ** 2 * g["Cp"] * (W3 ** ((2 * g["k0"]) % 3) - 1)) for g in G)
    if k == "NE3":
        return max(float(np.max(np.abs(g["arr"] * (np.conj(th.d[g["key"][0]].ell[:g["key"][2]]) / th.d[g["key"][0]].ell[:g["key"][2]] - 1)))) if g["key"][2] else 0.0 for g in G)
    if k == "NE5":
        if nm not in E.PHno2: return 0.0
        a_, b_ = E.PH[nm], E.PHno2[nm]
        return max(abs(a_.get(h, 0) - b_.get(h, 0)) for h in set(a_) | set(b_))
    if k in ("NE4", "NE7"):
        best = 0.0
        for g in G:
            sg, Nc, n = g["key"]
            if n == 0: continue
            mm = th.d[sg].m[:min(n, 50)]; Nm = (mm[:, 0] ** 2 - mm[:, 0] * mm[:, 1] + mm[:, 1] ** 2).astype(float)
            x0 = Nm / 81.0 * X / (Nc * Nc)
            v0 = wt.vsharp(x0)
            v1 = wt.vsharp(Nm / 81.0 * X / Nc) if k == "NE4" else wt.vsharp(x0, invert=True)
            best = max(best, float(np.max(np.abs(v1 - v0))))
        return best
    if k == "J":
        Gj = groups_jdual(E, nm, jd, {p: 5 for p in jd}, wt, X)
        val = max(float(np.max(np.abs(g["arr"] - g["arr_orig"]))) if g["key"][2] else 0.0 for g in Gj)
        return val, Gj
    raise KeyError(k)
