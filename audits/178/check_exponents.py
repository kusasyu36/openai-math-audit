"""
成果178 監査: setup.tex / repair.tex / cleanup.tex / foundations.tex に現れる n の冪の算術を a=1e-4 で検算。
各項目: (出典行, 論文の主張する冪, 定義から計算した冪)
"""
import sys
from fractions import Fraction as Fr
out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)
a = Fr(1, 10000)
def E(x):  # 冪を a の倍数と有理数で表示
    return f"{float(x):+.6f}"
# 基本スケール (setup.tex:55-73, repair.tex:157-160, 287-289)
e   = Fr(-2,3) + a
l0  = Fr(2,3) + Fr(15,100)*a
b   = Fr(1,3) + Fr(345,1000)*a
qL2 = l0/2 - Fr(29,100)*a           # q^{L/2}
qR0 = Fr(4,100)*a
alpha = Fr(-1,3) + a/2
Wb  = Fr(-2,3) + 100*a
M   = Fr(1,3) + 10*a
JT  = Fr(3,2)*a + a/10
vstar = Fr(125,1000)*a
checks = []
def chk(src, claim, val, label):
    ok = (claim == val)
    P(f"{'OK ' if ok else 'NG '} {src:28s} {label:48s} claim={E(claim)} calc={E(val)}")
chk("foundations.tex:364", Fr(1,3)+Fr(34,100)*a, qR0 + 2*l0 - 1, "q^R0 l1^2/n = n^{1/3+.34a}")
chk("foundations.tex:365", b, b, "b = n^{1/3+.345a} (margin over above = .005a)")
chk("foundations.tex:592", -Fr(58,100)*a, 2*qL2 - l0, "q^L/l0 = n^{-.58a}")
chk("foundations.tex:611", Fr(29,100)*a - l0/2, -qL2, "q^{-L/2} = n^{.29a}/sqrt(l0)")
chk("foundations.tex:942", -Fr(45,1000)*a, Fr(53,100)*a - alpha - l0/2, "n^{.53a}/(α* sqrt l0) = n^{-.045a}")
chk("foundations.tex:934", (Fr(125,1000)-Fr(16,100))*a, -Fr(16,100)*a + vstar, "b'||S_h|| = n^{(.125-.16)a}")
fb = b - l0
chk("cleanup.tex:29", Fr(27,100)*a, fb + l0/2, "f_b sqrt(l0) = n^{.27a}")
chk("cleanup.tex:32", Fr(575,1000)*a, alpha + l0/2, "α* sqrt(l0) = n^{.575a}")
P(f"    cleanup.tex:30  f_b q^{{L/2}} = n^{E(fb+qL2)} (needs η < {E(-(fb+qL2))})")
P(f"    cleanup.tex:31  f_b n^v*/α* = n^{E(fb+vstar-alpha)}")
P(f"    cleanup.tex:28  b^{{3/2}}/l0 = n^{E(Fr(3,2)*b - l0)}")
chk("repair.tex:153", Fr(-1,3)+Fr(325,1000)*a, Fr(40,100)*a - l0/2, "B_B = n^{.40a}/sqrt(l2) with g+r=0 -> n^{-1/3+.325a}")
P(f"    repair.tex:45,153  B_B = o(sqrt e) requires (g+r)/2 < {E(e/2 - (Fr(40,100)*a - l0/2))}  i.e. g+r < .35a")
chk("repair.tex:507", Fr(-2,3)+Fr(84,10)*a, M - JT - 1, "M/(J_T n) = n^{-2/3+8.4a}")
chk("repair.tex:573", Fr(-5,6)+155*a, M/2 + Fr(3,2)*Wb, "sqrt(M) W_b^{3/2} = n^{-5/6+155a}")
chk("repair.tex:599", -15*a, M - Fr(1,2) - Wb/4, "M/(sqrt n W_b^{1/4}) = n^{-15a}")
chk("repair.tex:600", -9*a, e + 1 - M, "en/M = n^{-9a}")
chk("repair.tex:609", Fr(-1,3)+20*a, 2*M - 1, "M^2/n = n^{-1/3+20a}")
P("--- 表 repair.tex:626-631 (比の冪) ---")
chk("repair.tex:626", Fr(74,10)*a, (M-JT-1) - e, "mean bias / e")
chk("repair.tex:627", Fr(-1,10)+159*a, M - Fr(1,10) + Fr(3,2)*Wb - e, "conditioning error / e")
chk("repair.tex:628", Fr(-1,6)+154*a, M/2 + Fr(3,2)*Wb - e, "near fluctuation / e")
chk("repair.tex:629", -6*a, (M - Fr(1,2) - Wb/4) - (e + 1 - M), "bulk fluctuation / (en/M)")
chk("repair.tex:630", Fr(-2,100)+9*a, Fr(-2,100) - (e+1-M), "typical-block error / (en/M)")
chk("repair.tex:631", Fr(-1,3)+20*a, 2*M-1, "collision prob")
P("--- Lemma 8.2 の例外割合 (repair.tex:263-277) ---")
P(f"    l0 (1+L)^6 / n = n^{E(l0-1)} polylog; fraction with row^2 > n^-.03: n^{E(l0-1+Fr(3,100))} polylog  (claimed <= n^-.2)")
P(f"    r_n = .08 log_q n: q^(-r_n/2) = n^-.04 ; n^.08 return-bad vertices / n = n^-.92")
P("--- n_0 の規模感 ---")
import math
for target in [2, 10]:
    for exp_ in [Fr(74,10)*a, 6*a, Fr(5,1000)*a]:
        lg = math.log10(target)/float(exp_)
        P(f"    n^{{{float(exp_):.5f}}} >= {target}  needs  n >= 10^{lg:.0f}")
P("    (比較: 可観測宇宙の原子数 ~ 10^80)")
open(sys.argv[0].replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
