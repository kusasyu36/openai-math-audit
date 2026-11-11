# -*- coding: utf-8 -*-
"""Control G4b: the X^{1/3} constant 0.075 that the old script found for the UNTWISTED paper's combination
   X^{-1/2} sum_A gamma_2(c)|d| U(N(A)/X)   (no angular factor; gamma_2 in the 003/029 normalisation).
   Prediction from the cusp structure: the 9 translates z + a/3 are all equivalent to infinity, so
   Theta_W(0,v) ~ sigma * (9 v)^{-2/3} * sum_h phi-hat(h) conj(chi_h), against sigma/3 * v^{-2/3} for Patterson's own.
   Hence constant_W = constant_P * 3 * 9^{-2/3} * sum_h phi-hat(h) conj(chi_h)."""
import os, math, cmath
import numpy as np
from scipy.special import gamma as Gam
from scipy.integrate import quad
from kubota import *
K = Kubota(200000, math.sqrt(200000 / 27.0), os.environ.get("TMPDIR", "/tmp"))
lam2 = mul(LAM, LAM); vals = {}
for a in range(-40, 41):
    for b in range(-40, 41):
        if is_primary((a, b)) and norm((a, b)) <= 1500: vals.setdefault((a % 9, b % 9), cubic_symbol(LAM, (a, b)))
tot = 0
for h0 in range(9):
    for h1 in range(9):
        s = sum(W3 ** vals[(n0, n1)] * cmath.exp(-2j * cmath.pi * omcoef(mul((n0, n1), (h0, h1))) / 9) for (n0, n1) in vals) / 81
        if abs(s) > 1e-13:
            a, cc = reduce_frac(mul(lam2, (h0, h1)), (9, 0)); g, kexp, rep, t = K.translate_parts(a, cc)
            assert rep == "I" and norm(cc) == 9
            tot += s * W3 ** (-kexp)
print("sum_h phi-hat(h) conj(chi_h) = %.12f%+.12fi" % (tot.real, tot.imag))
def U(x):
    t = (np.asarray(x, dtype=float) - 1.25) / 0.75; out = np.zeros_like(t); m = np.abs(t) < 1; out[m] = np.exp(-1.0 / (1 - t[m] ** 2)); return out
Ut = quad(lambda x: float(U(np.array([x]))[0]) * x ** (-1.0 / 6), 0.5, 2.0, epsabs=1e-13)[0]
cP = (2 * math.pi) ** (5.0 / 3) / (3 * Gam(2.0 / 3) * 27 ** (5.0 / 6)) * Ut
cW = cP * 3 * 9 ** (-2.0 / 3) * tot
print("predicted constants: Patterson's own %.6f ; paper's combination %.6f%+.6fi" % (cP, cW.real, cW.imag))
A = K.td.A; NA = (A[:, 0] ** 2 - A[:, 0] * A[:, 1] + A[:, 1] ** 2).astype(float)
aW = np.conj(K.td.G) * W3 ** K.td.kl
for X in (1000, 4000, 16000, 32000, 64000, 100000):
    S = np.sum(aW * U(NA / X)) / math.sqrt(X) / X ** (1 / 3)
    print("  X=%6d  measured %.6f%+.6fi   ratio to prediction %.4f     (old script: 0.0758 at 1000, 0.0751 at 4000, 0.0743 at 16000, 0.0753 at 32000)" % (X, S.real, S.imag, S.real / cW.real))
