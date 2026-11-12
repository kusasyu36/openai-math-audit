# -*- coding: utf-8 -*-
"""
003dual stage 8, TABLE-BLIND (lead, Claude Opus 5.5, 2026-11-12): split the stage-8 system by parity.
iota: (z, v) -> (-z, v) is conjugation by diag(i, -i); it maps [[a,b],[c,d]] to [[a,-b],[-c,d]], preserves Gamma(3),
SL2(Z) and chi ((-c/a)_3 = (c/a)_3), so the solution space splits into even (delta(-mu) = delta(mu)) and odd parts.
An odd automorphic form has zero constant term and cannot change the even (= table-type) coefficients.
(i) even subsystem: columns phi_mu + phi_{-mu}; (ii) odd subsystem: phi_mu - phi_{-mu}, homogeneous.
Smallest singular values vs truncation M at v0 = 0.30 (g_low = 0.85).  If the odd minimum decreases like the
truncation error as M grows, it points to a genuine odd form; if it stays put, the odd direction is weak but not null.
"""
import sys, os, math, time
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import numpy as np
import stage7_lib as L
from stage8_tune import build
kb = L.kb
print(__doc__)
for M in (7.0, 8.0, 9.0, 10.0, 11.0):
    t = time.time()
    nums, mu, A, b, sc, chi, sk, low = build(0.30, M, 30, (0, 0), 0.85)
    idx = {m: i for i, m in enumerate(nums)}
    pairs = [(i, idx[kb.neg(m)]) for i, m in enumerate(nums) if i < idx[kb.neg(m)]]
    P = np.array(pairs)
    Ae = A[:, P[:, 0]] + A[:, P[:, 1]]; Ao = A[:, P[:, 0]] - A[:, P[:, 1]]
    se = np.linalg.svd(Ae, compute_uv=False); so = np.linalg.svd(Ao, compute_uv=False)
    xe, *_ = np.linalg.lstsq(Ae, b, rcond=1e-14); re = np.linalg.norm(Ae @ xe - b) / np.linalg.norm(b)
    print(f"M={M}: unknown pairs {len(pairs)}; EVEN sv max {se[0]:.3g} min3 {np.array2string(se[-3:], precision=3)} resid {re:.1e}; "
          f"ODD sv max {so[0]:.3g} min3 {np.array2string(so[-3:], precision=3)}  ({time.time()-t:.0f}s)", flush=True)
