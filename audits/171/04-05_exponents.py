# Exact exponent checks for Sections 4 and 5 (own script; no repo code executed)
import sympy as sp
from fractions import Fraction as F
w = sp.symbols('omega', positive=True)
g = sp.symbols('gamma', positive=True)
h = w/10**6; k = w/3; b = w/30; b0 = b/4; rho = b0/10; zeta = b0/20; sig = zeta/10
checks = {
 # Step1 plateau: #steps*increment >1  : (w/2) - (w/3) > 0
 "plateau sum exponent w/2-w/3>0": sp.simplify(w/2 - w/3),
 # trimming t ~ n^{-w/3+w/5}
 "trim t exponent -(w/3-w/5)<0": sp.simplify(w/3 - w/5),
 # width cost kT|Z|L < gamma - w/2  (|Z|=gamma-w)
 "gamma-w/2 - (k+b+(gamma-w)+h) >0": sp.simplify(g - w/2 - (k + b + g - w + h)),
 # dip prob n^{-w/5+h} = o(n^{-2h})
 "w/5 - h - 2h >0": sp.simplify(w/5 - 3*h),
 # Azuma exponent k a*^2/L^2 = w/3-4h > w/5
 "w/3-4h - w/5 >0": sp.simplify(w/3 - 4*h - w/5),
 # T log n = o(n^{w/5})
 "w/5 - b >0": sp.simplify(w/5 - b),
 # T log n = o(k a_*)
 "k - h - b >0": sp.simplify(k - h - b),
 # height sublinear radius: b > b0 + 3 rho
 "b - b0 - 3rho >0": sp.simplify(b - b0 - 3*rho),
 # log vol O(r+H) neighbourhood o(a_* n): need rho>h and zeta-sigma>h
 "rho - h >0": sp.simplify(rho - h),
 "zeta - sigma - h >0": sp.simplify(zeta - sig - h),
 # height lemma: room for a with b0 > a > zeta+sigma+(1-theta) and b > a+4 sigma
 "b0 - zeta - sigma >0": sp.simplify(b0 - zeta - sig),
 "b - b0 - 4sigma >0 (a<b0 suffices)": sp.simplify(b - b0 - 4*sig),
}
for kname, v in checks.items():
    print(kname, "=", v, " positive:", sp.simplify(v).is_positive)
# m L = o(a_* n): gamma+13w+2h<1 given w <= (1-gamma)/1000
print("1-gamma-13w-2h with w=(1-g)/1000:", sp.simplify(1 - g - 13*(1-g)/1000 - 2*(1-g)/1000/10**6))
# masked first width o(a_* n): beta<1-h ; beta<=gamma<=1-1000w
print("1-1000w - (1-h) <0 :", sp.simplify(1-1000*w - (1-h)))
# log(1+1.5a) >= .7a on [0,1/2]  (q >= 1/2+3a/4)
a = sp.symbols('a')
f = sp.log(1+sp.Rational(3,2)*a) - sp.Rational(7,10)*a
print("f(1/2)=", sp.N(f.subs(a, sp.Rational(1,2))), " f concave, f(0)=0 -> nonneg on [0,1/2]")
# own-key product: (k-d)(-log2+.7a)-dL >= -k log2 + .65 a_* k with d <= a_* k/(20L)
print("c1 <= 0.7-0.05 =", F(7,10)-F(1,20))
# Section 4 constants: theta, a example
wv = sp.Rational(1,1000)
print("example b0,zeta,sigma:", [sp.nsimplify(x.subs(w,wv)) for x in (b0,zeta,sig)])
# Section 5
print("--- Section 5 ---")
print("log(1/.9)=", sp.N(sp.log(sp.Rational(10,9))), "< a0=1/3")
al = sp.symbols('alpha', positive=True)
# Pr(j>=h): m * n^{-.15} <= n^{-.13}  needs alpha < .02
print("alpha bound for n^{-.13h}:", sp.Rational(15,100)-sp.Rational(13,100))
# exponents in m-units: T=.001, k*=.005, J=.05, D_L=.15, M*=K_h eta+.005<.03
print("TM*<J: .001+.03 < .05 :", F(1,1000)+F(3,100) < F(5,100))
print("D_L log-cost ~ m^{.1}: (T+J)*k_J ~ .05+.05 = .10 < .15")
print("high concentration: 1 + alpha*(2*.005 - 2*.05) =", sp.simplify(1+al*(2*sp.Rational(5,1000)-2*sp.Rational(5,100))))
print("long/short mismatch: D_L=.15 < 1/5 :", F(15,100) < F(1,5))
# boundary + severities: -.05+.02, -.13h+.02(h+1)
print("boundary exp", F(-5,100)+F(2,100), "; h=1 severity exp", F(-13,100)+F(4,100))
# even posterior subset union: k log2 + tau0 kn - ups1 tau1 kn  -> need ups1*tau1>tau0
t0, t1 = sp.Rational(5,10), sp.Rational(6,10)
print("example tau0=.5,tau1=.6 <log2=", sp.N(sp.log(2)), " tau0/tau1=", t0/t1)
# rho for even repeat: log2 - H(2rho) > tau1
H = lambda x: -x*sp.log(x)-(1-x)*sp.log(1-x)
for r in [sp.Rational(1,100), sp.Rational(1,50)]:
    print("rho=",r," log2-H(2rho)=", sp.N(sp.log(2)-H(2*r)))
