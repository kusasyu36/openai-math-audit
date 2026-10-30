import sympy as sp, mpmath as mp
x=sp.symbols('x',positive=True)
I=sp.integrate(sp.sin(x)**4/x**2,(x,0,sp.oo))
print("int_0^inf sin^4 x/x^2 =",I)
val=64*sp.pi*2*sp.pi*I
print("||E1||_4^4 =",sp.simplify(val)," vs 32 pi^3:",sp.simplify(val-32*sp.pi**3)==0)
# numeric cross-check with large cutoff + tail mean (sin^4 avg 3/8)
mp.mp.dps=30
X=200
head=mp.quad(lambda r: mp.sin(2*mp.pi*r)**4/r**2, mp.linspace(0,X,4*X+1))
tail=mp.mpf(3)/8/X
print("numeric:",64*mp.pi*(head+tail), float(32*mp.pi**3))
