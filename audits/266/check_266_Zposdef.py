"""Companion paper (Exact Fourier certificates), Section 05 / Appendix A:
re-check positive definiteness of the six integer matrices Z_lambda with exact
rational Schur complements (all pivots > 40000), using the literal `data` list
read from verify.py with ast.literal_eval (NO code from verify.py is executed)."""
import ast, re, sys
from fractions import Fraction as Fr
src = open("data/openai_math/preprints/Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026/verification/verify.py").read()
i = src.index("\ndata=[")
j = src.index("\nEpair=pairtable()")
data = ast.literal_eval(src[i+len("\ndata="):j].strip())
print("shapes (lambda, m, e_lambda):", [(lam, len(ws), e) for lam, ws, e, a in data])
expected = [((1,1,1,2),5,14),((1,1,3),21,13),((1,2,2),11,13),((1,4),27,14),((2,3),22,14),((5,),12,15)]
print("matches table in 04-mixed.tex:234-247:", [(lam,len(ws),e) for lam,ws,e,a in data]==expected)
total_pivots = 0; minpiv = None; allok = True
for lam, ws, e, a in data:
    m = len(ws)
    assert len(a) == m*(m+1)//2, "lower triangle length"
    assert sum(lam) == 5 and all(t in (0,1) and len(v) == 6 and sum(v) == 0 for t, v in ws)
    Z = [[Fr(0)]*m for _ in range(m)]
    it = iter(a)
    for r in range(m):
        for c in range(r+1):
            Z[r][c] = Z[c][r] = Fr(next(it))
    # Schur complement elimination from the last index (as described in 05-verification.tex:326-344)
    B = [row[:] for row in Z]
    piv = []
    while B:
        q = B.pop(); p = q.pop()
        piv.append(p)
        B = [[x - y*aa/p for x, aa in zip(b, q)] for b, y in zip(B, q)]
    ok = all(p > 40000 for p in piv)
    allok &= ok
    total_pivots += len(piv); minpiv = min(piv) if minpiv is None else min(minpiv, min(piv))
    print(f"lambda={lam}: m={m}, pivots>40000: {ok}, min pivot={float(min(piv)):.1f}")
print("total pivots:", total_pivots, "(paper: 98); min pivot:", float(minpiv), "; all positive-definite:", allok)
import numpy as np
mins = []
for lam, ws, e, a in data:
    m = len(ws); Z = np.zeros((m, m)); it = iter(a)
    for r in range(m):
        for c in range(r+1):
            Z[r, c] = Z[c, r] = next(it)
    mins.append(float(np.linalg.eigvalsh(Z).min()))
print("numpy smallest eigenvalues:", [f"{x:.1f}" for x in mins])
