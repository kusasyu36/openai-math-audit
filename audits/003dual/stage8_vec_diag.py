# -*- coding: utf-8 -*-
"""
003dual stage 8, TABLE-BLIND diagnostic (lead, Claude Opus 5.5, 2026-11-12): shape of the right singular vector of
the isolated smallest singular value (system of stage8_tune with g_low = 0.85, v0 = 0.30, M = 9).
Reports, in unscaled coordinates delta tau(mu) (normalised to max 1 on |mu| <= 3):
  parity: || delta(mu) + delta(-mu) || vs || delta(mu) - delta(-mu) ||;
  rotation by omega: delta(omega mu) vs delta(mu);  reflection: delta(conj mu)...;
  weight per lambda-valuation class;  and how big the function it represents is at a high point (v = 1, 2).
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage7_lib as L
from stage8_tune import build
kb = L.kb
print(__doc__)
nums, mu, A, b, sc, chi, sk, low = build(0.30, 9.0, 30, (0, 0), 0.85)
U, s, Vh = np.linalg.svd(A, full_matrices=False)
print("smallest 4 sv", s[-4:])
for idx in (-1, -2):
    v = Vh[idx].conj() / sc
    a = np.abs(mu); sel = a <= 3 * (1 + 1e-12)
    v = v / np.abs(v[sel]).max()
    D = dict(zip(nums, v))
    def pair_norm(f):
        p = [(D[m], D.get(f(m))) for m in nums if abs(L.mu_of([m])[0]) <= 3 and f(m) in D]
        x = np.array([q[0] for q in p]); y = np.array([q[1] for q in p])
        return np.linalg.norm(x - y), np.linalg.norm(x + y)
    print(f"\nsingular vector {idx} (sv {s[idx]:.3e}):")
    print("  parity  |d(mu)-d(-mu)|, |d(mu)+d(-mu)| =", ["%.3e" % t for t in pair_norm(kb.neg)])
    print("  omega-rotation |d(mu)-d(omega mu)|, |+| =", ["%.3e" % t for t in pair_norm(lambda m: kb.mul(m, kb.OMEGA))])
    by = {}
    for m, x in zip(nums, v):
        if abs(L.mu_of([m])[0]) > 3: continue
        val = 0; y = m
        while kb.exact_div(y, kb.LAM) is not None: y = kb.exact_div(y, kb.LAM); val += 1
        by[val] = by.get(val, 0) + abs(x) ** 2
    print("  weight by lambda-valuation of lambda^3 mu (|mu|<=3):", {k: "%.3f" % w for k, w in sorted(by.items())})
    top = np.argsort(-np.abs(v) * sel)[:6]
    print("  largest entries:", ", ".join(f"m={nums[i]} |mu|={a[i]:.3f} |d|={abs(v[i]):.3f}" for i in top))
