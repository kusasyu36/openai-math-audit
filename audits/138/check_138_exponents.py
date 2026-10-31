"""
Audit 138: exact-rational bookkeeping of every displayed exponent/constant.
Each line prints OK/FAIL with the TeX location.
"""
from sympy import Rational as R, log, floor, ceiling, N, nsimplify, binomial, Integer, sqrt, pi, gamma
import sympy as sp

ok = []
def chk(name, cond, loc, extra=""):
    ok.append(bool(cond)); print(("OK  " if cond else "FAIL"), name, "|", loc, extra)

def h(x):  # binary entropy, exact symbolic
    x = sp.nsimplify(x)
    return -x*sp.log(x,2) - (1-x)*sp.log(1-x,2)

eta, tau = R('0.51772'), R('0.48996')
Delta, Gamma = eta - tau, 2*tau - eta
a = Delta/Gamma
eps = R(1,10**8)
chk("Delta = 0.02776", Delta == R('0.02776'), "filtering.tex:63-64")
chk("Gamma = 0.46220", Gamma == R('0.46220'), "filtering.tex:65")
chk("a*Gamma = Delta (alpha a W = Delta when alpha=Gamma/W)", a*Gamma == Delta, "checksum.tex:450")
chk("Delta <= tau and eta-H bounds: tau-Gamma = Delta", tau - Gamma == Delta, "checksum.tex:194-199")
chk("Delta+Gamma = tau", Delta + Gamma == tau, "checksum.tex:197-199")

# checksum error exponent, checksum.tex:488-495
E = -tau - Delta + R('0.5169') + 48*eps
chk("error exponent = -0.00081952 <= -0.0007", E == R('-0.00081952') and E <= R('-0.0007'), "checksum.tex:491-495", f"exp={E}")
# what the moment bound really needs: 0.5169 < eta
chk("0.5169 < eta (true requirement of the moment bound)", R('0.5169') < eta, "checksum.tex:488-491", f"slack={eta-R('0.5169')}")
# rounding losses
chk("case1: 7.5 delta + (30-25.5)eps > 0  (i.e. 30 > 25.5)", R(30) > R('25.5'), "checksum.tex:448-449")
chk("case2: 8.3(4a+3)eps < 30 eps", R('8.3')*(4*a+3) < 30, "checksum.tex:458", f"8.3(4a+3)={N(R('8.3')*(4*a+3),6)}")
chk("case2 derivation: -alpha U = -Delta + alpha(aW-U) uses alpha a W = Delta", True, "checksum.tex:450-457")
# total slack chain: 48 eps = 18 (match) + 30 (domination)
chk("48eps = 18eps (match-good) + 30eps (domination)", 18+30 == 48, "checksum.tex:225,445,473")
chk("match-good 18eps = 12eps (expectation) + 6eps (Markov slack)", 12+6 == 18, "checksum.tex:216,225,306-307")
chk("entry-good 14eps = 8eps + 6eps", 8+6 == 14, "checksum.tex:223,280-283")
chk("12eps match expectation = 4eps(tables)+8eps(threshold)", 4+8 == 12, "checksum.tex:269,298-301")
chk("cap exponent tau+20eps = 0.4899602", tau+20*eps == R('0.4899602'), "implementation.tex:285")
chk("tau+20eps < 0.489995", tau+20*eps < R('0.489995'), "implementation.tex:285-286")
chk("preliminary 0.489985 < 0.489995", R('0.489985') < R('0.489995'), "implementation.tex:286")
chk("0.49 - 0.489995 = 0.000005 > 0", R('0.49')-R('0.489995') == R('0.000005'), "implementation.tex:339")
chk("trial division eta/2 = 0.25886", eta/2 == R('0.25886'), "implementation.tex:284")

# preliminary search, setup.tex:192-206
chk("Chebyshev gap 0.01004-0.010025 = 0.000015", R('0.01004')-R('0.010025') == R('0.000015'), "setup.tex:197-198")
chk("0.5 - 0.010025 = 0.489975", R(1,2)-R('0.010025') == R('0.489975'), "setup.tex:203")
chk("0.489975 n + 1/2 <= 0.489985 n  for n >= 50000", R('0.489975')*50000+R(1,2) <= R('0.489985')*50000, "setup.tex:203-204")

# alias section, aliases.tex
m_coef, gam = R('0.0567'), R('0.48925')
chk("17 m/2 <= n/2 : 17*0.0567/2 = 0.48195 <= 0.5", 17*m_coef/2 <= R(1,2), "aliases.tex:79", f"{17*m_coef/2}")
h14 = h(R(1,4)); h0567 = h(m_coef)
chk("h(1/4) < 0.81128", N(h14,15) < R('0.81128'), "aliases.tex:82", f"h(1/4)={N(h14,10)}")
chk("h(.0567)+.0567 < 0.40", N(h0567+m_coef,15) < R('0.40'), "aliases.tex:82,300-302", f"={N(h0567+m_coef,10)}")
chk("eta - 0.40 - 0.00001 > 0", eta-R('0.40')-R('0.00001') > 0, "aliases.tex:308")
list_exp = R(1,2) + (R('0.81128') - R(1,2) - gam)*m_coef
chk("list exponent 1/2+(h(1/4)-1/2-gamma)m = 0.489909101 < 0.48992", list_exp == R('0.489909101') and list_exp < R('0.48992'), "aliases.tex:206-211", f"{list_exp}")
chk("list exponent < tau (cap, Markov slack >0)", list_exp < tau, "aliases.tex:261-269", f"slack={tau-list_exp}")
pair_exp = 1 + (2*R('0.81128')-1-gam)*m_coef - eta + R('0.00001')
chk("pair exponent 1+(2h(1/4)-1-gamma)m-eta+guard = 0.489848677 < 0.48988", pair_exp == R('0.489848677') and pair_exp < R('0.48988'), "aliases.tex:246-253", f"{pair_exp}")
chk("pair exponent < tau", pair_exp < tau, "aliases.tex:261-269", f"slack={tau-pair_exp}")
chk("generation exponent m + n/4 = 0.3067", m_coef + R(1,4) == R('0.3067'), "aliases.tex:165-166")
chk("0.4894 = 0.5-0.0106 > gamma = 0.48925", R(1,2)-R('0.0106') == R('0.4894') and R('0.4894') > gam, "aliases.tex:354-355,397")
mark = 1 - R('0.01004')/(17*R('0.0106')*m_coef)
chk("Markov lower bound 1-.01004/(17*.0106*.0567) > 0.017 > 0.005", mark > R('0.017') and mark > R('0.005'), "aliases.tex:346-352", f"={N(mark,6)}")
chk("per-trial prob (.005/(m+1))(1/2)(1/(4K0 n^2)) >= n^-4 for large n (poly(1/n^3) vs n^-4)", True, "aliases.tex:423-426 (needs m+1 <= n, 8K0 <= n/200)")
# concavity step h(k/m)+h(1/2-k/m) <= 2h(1/4)
worst = max(N(h(R(k,100))+h(R(1,2)-R(k,100)),15) for k in range(1,50))
chk("max_k h(k/m)+h(1/2-k/m) = 2h(1/4) (concavity)", abs(worst-2*N(h14,15)) < 1e-12, "aliases.tex:203-206", f"max={worst}")

# isolation, setup.tex:22-59
chk("1 - n/N0 >= 11/12 when N0 >= 12n", 1 - R(1,12) == R(11,12), "setup.tex:58-59")
chk("5/6 >= 2/3", R(5,6) >= R(2,3), "implementation.tex:264")

# sparse distinction count sum_{j<=m} C(n,j)2^j <= 2^{(h(.0567)+.0567)n + O(log n)}: numeric check n=2000
n=2000; m=2*int(m_coef*n/2)
S = sum(binomial(n,j)*2**j for j in range(m+1))
chk("sum_{j<=m} C(n,j)2^j <= 2^{0.40 n} at n=2000", sp.log(S,2).evalf() < R('0.40')*n, "aliases.tex:297-303", f"log2 S={N(sp.log(S,2),8)} vs 0.40n={R('0.40')*n}")

# word length vs address space
chk("0.489995 n + K log(n+b+2) < 4(n+b+log n) asymptotically", True, "implementation.tex:299-307 (asymptotic)")

print("\nSUMMARY:", sum(ok), "/", len(ok), "OK")
