# Lem 3.1 (completion.tex:55-75) の数値確認（主担当、2026-11-11）。
# M=(1 p), p は絶対値 1 の成分 4 個。M^*M=[[4,S],[conj S,4]]、特異値^2 = 4±|S|、
# M(M^*M)^{-1} = (4*1 - p conj(S), 4p - 1*S)/(16-|S|^2) を乱数 2000 例で確認。
import numpy as np
rng = np.random.default_rng(266)
worst = [0.0, 0.0, 0.0]
for _ in range(2000):
    p = np.exp(2j*np.pi*rng.random(4)); one = np.ones(4, dtype=complex)
    S = p.sum(); M = np.column_stack([one, p]); G = M.conj().T @ M
    worst[0] = max(worst[0], np.abs(G - np.array([[4, S], [np.conj(S), 4]])).max())
    sv2 = np.sort(np.linalg.svd(M, compute_uv=False)**2)
    worst[1] = max(worst[1], np.abs(sv2 - np.array([4-abs(S), 4+abs(S)])).max())
    if abs(S) < 3.99:
        lhs = M @ np.linalg.inv(G)
        rhs = np.column_stack([4*one - p*np.conj(S), 4*p - one*S])/(16-abs(S)**2)
        worst[2] = max(worst[2], np.abs(lhs-rhs).max()*(16-abs(S)**2))
names = ["M^*M = [[4,S],[conj S,4]]", "squared singular values = 4 +- |S|", "inverse-frame formula (scaled residual)"]
for n, w in zip(names, worst):
    print(("OK  " if w < 1e-9 else "FAIL"), n, "max residual %.2e" % w)
