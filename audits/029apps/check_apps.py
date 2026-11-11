# -*- coding: utf-8 -*-
"""029 §7 (07-splitting.tex) と 142 §11 (10-analytic.tex) の、体が動く応用の節の算術を確かめる自作検算。
標準ライブラリのみ。確かめるのは不等式の帳尻だけで、零点なし帯そのものや明示公式の解析は確かめない。"""
from fractions import Fraction as Fr
from math import log, gcd
ok = 0; ng = 0
def chk(name, cond):
    global ok, ng
    print(("OK  " if cond else "NG  ") + name)
    if cond: ok += 1
    else: ng += 1
def primes(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p ** .5) + 1))]
def phi(m):
    r = m; x = m; p = 2
    while p * p <= x:
        if x % p == 0:
            while x % p == 0: x //= p
            r -= r // p
        p += 1
    if x > 1: r -= r // x
    return r
def log_cyc_disc(m):
    """log |disc Q(mu_m)| = phi(m) ( log m - sum_{p|m} log p/(p-1) )"""
    s = 0.0; x = m; p = 2; ps = []
    while p * p <= x:
        if x % p == 0:
            ps.append(p)
            while x % p == 0: x //= p
        p += 1
    if x > 1: ps.append(x)
    return phi(m) * (log(m) - sum(log(p) / (p - 1) for p in ps))

# ---- 1. 指数の帳尻（有理数で厳密）
delta = Fr(1, 10**6)
chk("142: x=T*B^(2e7) で x^(-delta) の B の指数が -20", -delta * 2 * 10**7 == -20)
chk("142: 誤差/主項 の B の指数が -20+3 = -17 (10-analytic.tex:130)", -delta * 2 * 10**7 + 3 == -17)
# ---- 2. 142: 次数と判別式 (10-analytic.tex:113-122)
allok = True; worst = 0.0
for q in primes(400):
    dK = phi(12 * q)
    if not dK <= 12 * q: allok = False
    lDK = log_cyc_disc(12 * q)
    if not lDK <= dK * log(12 * q) + 1e-9: allok = False
chk("142: [K:Q]=phi(12q)<=12q, log D_K<=[K:Q]log(12q)  (q<=400 の全素数。q=2,3 を含む)", allok)
# log D_M + [M:Q] <= C*B^3 が絶対定数 C で成り立つこと。B=20+(n+1)(L+1), q<=n, L=log p
allok = True; worst = 0.0
for q in primes(400):
    for n in (q, q + 1, 2 * q, 10 * q):
        for L in (1.0, 5.0, 50.0, 1e3, 1e6):
            B = 20 + (n + 1) * (L + 1)
            dK = phi(12 * q)
            lDM = q * dK * log(12 * q) + dK * (q * log(q) + (q - 1) * L)   # tex の上界そのもの
            val = (lDM + 12 * q * q) / B ** 3
            worst = max(worst, val)
chk("142: (log D_M + [M:Q])/B^3 が有界（走査範囲の最大 %.4f）" % worst, worst < 5)
# ---- 3. 029: Lemma 7.1 の判別式の上界 (07-splitting.tex:22-25)
allok = True; worst = 0.0
for a in (2, 3, 5, 10, -7, 1000):
    for q in primes(2000):
        if q < 5: continue
        n_q = q * (q - 1)
        # log of  q^{q(q-2)} (q^q |a|^{q-1})^{q-1}
        lD = q * (q - 2) * log(q) + (q - 1) * (q * log(q) + (q - 1) * log(abs(a)))
        # |d_{E_q}| = q^{q-2} と一般公式の照合
        if abs(log_cyc_disc(q) - (q - 2) * log(q)) > 1e-6: allok = False
        worst = max(worst, (lD + n_q) / (n_q * log(2 * q) * (1 + log(abs(a)))))
chk("029: |d_{Q(mu_q)}| = q^(q-2) が円分体の一般公式と一致 (5<=q<=2000)", allok)
chk("029: (log D_q 上界 + n_q)/(n_q log 2q (1+log|a|)) が有界（最大 %.4f）" % worst, worst < 5)
# ---- 4. 029: 最終の和 (01a-reduction.tex:185-192)。 x^(1-d0) exp(L^0.3) = x exp(-d0 L + L^0.3) = o(x/L^2)
# -d0 L + L^0.3 + 2 log L -> -inf。符号が負に転じる L を数値で出す（閾値の大きさの記録用）
d0 = 1e-6
def g(L): return -d0 * L + L ** 0.3 + 2 * log(L)
lo, hi = 1.0, 1e30
for _ in range(300):
    mid = (lo * hi) ** .5
    if g(mid) > 0: lo = mid
    else: hi = mid
print("     029: -d0*L + L^0.3 + 2 log L が負になるのは L = log x > %.3e から（x > exp(%.2e)）" % (hi, hi))
chk("029: L=1e9, 1e12 で負", g(1e9) < 0 and g(1e12) < 0)
# 分解評価の第 2 項: x^(1-d0) log(2q)/L <= x^(1-d0) は q<=exp(L^0.3) で log(2q)/L<=1 から
chk("029: q<=exp(L^0.3) なら log(2q)/L <= 1 (L>=3)", all((log(2) + L ** 0.3) / L <= 1 for L in (3, 10, 1e3, 1e9)))
# ---- 5. 142: 相殺の組合せ。K の 1 次素イデアル（ノルム l）が M で完全分解すれば q 個、q^{-1} 倍で相殺
#        小さい例で「p が l を法として q 乗剰余 <=> X^q - p が F_l で q 個の相異なる根」を総当たり
allok = True; cnt = 0
for q in (2, 3, 5, 7):
    for l in primes(700):
        if l % (12 * q) != 1: continue
        for p in primes(60):
            if p in (2, 3, q, l): continue
            roots = sum(1 for t in range(l) if (pow(t, q, l) - p) % l == 0)
            isres = pow(p, (l - 1) // q, l) == 1
            if (roots == q) != isres or roots not in (0, q): allok = False
            cnt += 1
chk("142: l≡1 (12q) で X^q-p の F_l での根の数は 0 か q、q 個 <=> p^((l-1)/q)=1  (%d 組)" % cnt, allok)
# q=2: 補助素数 l≡1 (mod 24), p は l の非剰余 => 相互法則で l は p の非剰余
allok = True; cnt = 0
for l in primes(3000):
    if l % 24 != 1: continue
    for p in primes(300):
        if p in (2, 3, l): continue
        a = pow(p, (l - 1) // 2, l) == 1; b = pow(l % p, (p - 1) // 2, p) == 1
        if a != b: allok = False
        cnt += 1
chk("142: q=2 で (p/l)=(l/p)（l≡1 mod 4）。補助素数は p の平方非剰余 (%d 組)" % cnt, allok)
print("\n合計: OK %d, NG %d" % (ok, ng))
