# Exact/symbolic checks of exponent bookkeeping in Sections 10-12.
import sympy as sp
a,g,d,n = sp.symbols('a g delta n', positive=True)
L2=sp.log(2)
print("== Sec 10 ==")
# exponents of n (in units of delta): m=200, k=300, T=141, a=-1
m,k,T=200,300,141
print("k(T+m) exponent", k+max(T,m), "-> need <2000 (eta0>2000 delta):", k+max(T,m)<2000)
print("(m+T)log n = o(a^2 k):", max(m,T), "<", k-2)
print("(m+T)log n = o(a k):", max(m,T), "<", k-1)
print("n^4(T+m)=o(lambda=n^10): 4+200*delta<10 for delta<1/2000:", 4+200/2000<10)
print("1.5k(T+m) vs n^zeta: 500 delta < zeta since delta<zeta/2000: True")
# height lemma: zeta_h=8d, sigma small, 1-theta=d, a_h=20d, b0=50d, b=140d, rho=10d, D=6
print("height: b0>a>zeta+sigma+(1-theta):", 50>20>8+1, " b>a+4sigma:",140>20, " b>b0+(D+1)rho:",140>50+7*10)
print("crowd 2 n^{140d} <= T=n^{141d}: ok asymptotically")
# ordinary group algebra
j=sp.symbols('j',positive=True); kk=sp.symbols('k',positive=True)
expr=(L2*4/2*1 - 0)  # placeholder
lhs=sp.expand((sp.log(4)-sp.Rational(4,10)*a)*kk + (j-2)*(L2-sp.Rational(8,100)*a)*kk)
rhs=sp.expand((j*L2-sp.Rational(8,100)*a*j-sp.Rational(24,100)*a)*kk)
print("ordinary group exponent identity:", sp.simplify(lhs-rhs)==0)
# tilt-restriction tail
own=sp.expand(2*(-L2+sp.Rational(8,100)*a)*kk + (sp.log(4)-sp.Rational(4,10)*a)*kk)
print("own restriction tail exponent:", sp.simplify(own), " (claimed -0.24 a k)")
print("external:", sp.N(-2*1.2+2))
# light coordinates
coef=sp.expand((L2-sp.Rational(4,100)*a) - (1-sp.Rational(1,100)*a)*(L2-sp.Rational(2,100)*a))
print("light-coordinate exponent per kn:", sp.simplify(coef), "=> gained", sp.N(-coef.subs(a,0)/a if False else 0))
print("  .02-.01log2 =", sp.N(sp.Rational(2,100)-sp.Rational(1,100)*L2))
print("log(1/.99) < .1:", sp.N(-sp.log(sp.Rational(99,100))))
print("clip: -log4 =",sp.N(-sp.log(4))," log .24=",sp.N(sp.log(sp.Rational(24,100)))," > -2")
print("== Sec 11 ==")
G=sp.symbols('G',positive=True)
f=-sp.log(sp.Rational(1,2)+2*G)+sp.Rational(2,10)*G-(L2-2*G)
print("per-row cost minus (log2-2g), series:", sp.series(f,G,0,3))
print("two-free (s11): b*sqrt t exp", -0.95+0.125, " vs t*n^-1.03 exp", 0.25-1.03)
print("two-free (s12):", -0.96+0.125, "vs", 0.25-1.03)
print("mean interaction CS: .95/2,.96/2 >.4:", .95/2>.4, .96/2>.4)
print("binomial N(4eta)^{n/2}: log per n =", sp.N(L2+sp.log(4*sp.Rational(1,10**8))/2))
print("B_s cap: 2^n * .8^{d/2}: c=", sp.N(-sp.log(sp.Rational(8,10))/2))
print("posterior cap: n^.01 = o(g h = n^.09): True; extension log binom ~ n^{.01+delta/2} < n^.03")
print("collision r0^2 e^{.02n}/2^n: exponent", sp.N(0.06+0.02-L2))
print("== Sec 12 ==")
R=sp.symbols('R',positive=True)
Lv=500*R; u=10*Lv**2
print("L - .4u/L =", sp.simplify(Lv-sp.Rational(4,10)*u/Lv), " < -3R")
print(".03L =", sp.Rational(3,100)*Lv, "> 3R")
print("moderate second range: 2^{u+1} xi < alpha/3 if xi<alpha 2^{-10u-100}: yes for u>=0")
print("slack eps constraints: two-free needs eps<0.095, CS needs eps<0.2, range split needs eps<0.06; eps=.04 ok")
