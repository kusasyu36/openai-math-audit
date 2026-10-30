"""Numerical test of Theorem 1.2 (m=1,2,3): gap = h - m - 1/2 log det Cov >= 0, own code."""
import mpmath as mp
mp.mp.dps=20
def gap2(g, logg, xr, ylo, yhi):
    def I(fn):
        return mp.quad(lambda x: mp.quad(lambda y: fn(x,y), [ylo(x), yhi(x)]), xr)
    Z=I(g); m1=I(lambda x,y:x*g(x,y))/Z; m2=I(lambda x,y:y*g(x,y))/Z
    cxx=I(lambda x,y:x*x*g(x,y))/Z-m1**2; cyy=I(lambda x,y:y*y*g(x,y))/Z-m2**2; cxy=I(lambda x,y:x*y*g(x,y))/Z-m1*m2
    ent=-I(lambda x,y: g(x,y)*logg(x,y))/Z+mp.log(Z)
    return ent-2-mp.log(cxx*cyy-cxy**2)/2
inf=mp.inf
tests=[]
tests.append(('exp product (equality)', lambda x,y: mp.e**(-x-y), lambda x,y:-x-y, [0,inf], lambda x:0, lambda x:inf))
tests.append(('e^{-x} on cone |y|<=x (cone over interval = simplex; equality)', lambda x,y: mp.e**(-x), lambda x,y:-x, [0,inf], lambda x:-x, lambda x:x))
tests.append(('e^{-x-y} on cone 0<=y<=2x (equality expected)', lambda x,y: mp.e**(-x-y), lambda x,y:-x-y, [0,inf], lambda x:0, lambda x:2*x))
for Lc in (0.5,1,2,4,8):
    tests.append((f'e^(-x-y) on x,y>=0, x+y<={Lc}', lambda x,y: mp.e**(-x-y), lambda x,y:-x-y, [0,Lc], lambda x:0, (lambda Lc: lambda x: Lc-x)(Lc)))
tests.append(('uniform triangle', lambda x,y: mp.mpf(1), lambda x,y:mp.mpf(0), [0,1], lambda x:0, lambda x:1-x))
tests.append(('uniform square', lambda x,y: mp.mpf(1), lambda x,y:mp.mpf(0), [0,1], lambda x:0, lambda x:1))
tests.append(('uniform disk', lambda x,y: mp.mpf(1), lambda x,y:mp.mpf(0), [-1,1], lambda x:-mp.sqrt(1-x*x), lambda x:mp.sqrt(1-x*x)))
tests.append(('e^{-max(x,y)} on orthant', lambda x,y: mp.e**(-max(x,y)), lambda x,y:-max(x,y), [0,inf], lambda x:0, lambda x:x))  # half by symmetry -> use region y<=x
tests.append(('e^{-x} 1{|y|<=sqrt x}', lambda x,y: mp.e**(-x), lambda x,y:-x, [0,inf], lambda x:-mp.sqrt(x), lambda x:mp.sqrt(x)))
tests.append(('x e^{-x} on |y|<=x', lambda x,y: x*mp.e**(-x), lambda x,y:mp.log(x)-x if x>0 else mp.mpf(0), [0,inf], lambda x:-x, lambda x:x))
tests.append(('e^{-x-y^2} on x>=0', lambda x,y: mp.e**(-x-y*y), lambda x,y:-x-y*y, [0,inf], lambda x:-inf, lambda x:inf))
tests.append(('e^{-x-xy?}: e^{-(x+|y|)} on x>=0,|y|<=x', lambda x,y: mp.e**(-x-abs(y)), lambda x,y:-x-abs(y), [0,inf], lambda x:-x, lambda x:x))
tests.append(('e^{-sqrt(x^2+y^2)} on quarter plane 0<=y<=x', lambda x,y: mp.e**(-mp.sqrt(x*x+y*y)), lambda x,y:-mp.sqrt(x*x+y*y), [0,inf], lambda x:0, lambda x:x))
tests.append(('e^{-x} on cone 0<=y<=x minus x>=1 cap? (e^{-x} on {0<=y<=x<=3})', lambda x,y: mp.e**(-x), lambda x,y:-x, [0,3], lambda x:0, lambda x:x))
for name,g,lg,xr,lo,hi in tests:
    try:
        print(f'{name:60s} gap = {mp.nstr(gap2(g,lg,xr,lo,hi),10)}')
    except Exception as e:
        print(name,'ERR',e)
# 1D sanity: Gamma(k) gaps  (MNR theorem)
for k in (1,1.5,2,5):
    k=mp.mpf(k); h=k+mp.loggamma(k)+(1-k)*mp.digamma(k)
    print(f'1D Gamma({k}) gap = {mp.nstr(h-1-mp.log(k)/2,10)}')
# simplex in R^n via uniform: gap = -n(1+log L)
for n in range(1,9):
    L=mp.factorial(n)**(mp.mpf(1)/n)/((n+1)**((n+1)/(2*mp.mpf(n)))*mp.sqrt(n+2))
    print(f'uniform simplex n={n}: gap = {mp.nstr(-n*(1+mp.log(L)),8)}')
