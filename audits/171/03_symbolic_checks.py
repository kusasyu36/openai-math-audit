# Section 3 audit: exact/symbolic checks of identities and exponent bookkeeping
import sympy as sp
from fractions import Fraction as F
import itertools, random

print("== Clock Step 4: walk-sum identity")
th,T=sp.symbols('theta T',positive=True)
j=sp.symbols('j',integer=True,nonnegative=True)
k=sp.symbols('k',integer=True,nonnegative=True)
even=sp.summation(th**k*T**(2*k)/sp.factorial(2*k),(k,0,sp.oo))
odd=sp.summation(th**k*T**(2*k+1)/sp.factorial(2*k+1),(k,0,sp.oo))
rhs=sp.cosh(sp.sqrt(th)*T)+sp.sinh(sp.sqrt(th)*T)/sp.sqrt(th)
print(sp.simplify((even+odd-rhs).rewrite(sp.exp)))

print("== Clock Step 6: ODE z'=-qz, q'=-theta z q; q=1-theta+theta z; (log zq)'=-(q+theta z); int zq = 1-z")
t=sp.symbols('t',positive=True)
z=sp.Function('z')(t)
q=1-th+th*z
qprime=sp.diff(q,t).subs(sp.Derivative(z,t),-q*z)
print("q' + theta z q =",sp.simplify(qprime+th*z*q))
dlog=sp.diff(sp.log(z*q),t).subs(sp.Derivative(z,t),-q*z)
print("(log zq)' + (q+theta z) =",sp.simplify(dlog+(q+th*z)))
print("zq = -z' so int_0^T zq = z(0)-z(T) = 1-z(T): trivially since z'=-qz")

print("== Clock Step 5: branching exponential-moment closure (exact integration)")
# Suppose Hl-1<=4 d E(s), Hr-1<=4 d E(s)/sqrt(th), E(s)=exp(2 sqrt(th) s)
d_=sp.symbols('delta',positive=True)
s=sp.symbols('s',positive=True)
E=sp.exp(2*sp.sqrt(th)*s)
int_l=sp.integrate(4*d_*E,(s,0,t))
int_r=sp.integrate(4*d_*E/sp.sqrt(th),(s,0,t))
Et=sp.exp(2*sp.sqrt(th)*t)
xr=d_+int_l          # exponent for H_row
xl=d_+th*int_r       # exponent for H_label
print("x_row <= 2 d E/sqrt(th)?  slack =",sp.simplify(2*d_*Et/sp.sqrt(th)-xr))
print("x_lab <= 2 d E?           slack =",sp.simplify(2*d_*Et-xl))
# numeric with theta=1e-6, K0 s.t. small exponents
for K0 in [100,1000]:
  for n in [10**3,10**6]:
    thv=1e-6; T0=K0*sp.log(n); dv=n**(-0.01*K0)
    xv=float(2*dv*sp.exp(2*sp.sqrt(thv)*T0)/sp.sqrt(thv))
    print(f"K0={K0} n={n}: max exponent x_row={xv:.3e} (need small for e^x-1<=2x, i.e. x<=1.256)")

print("== Height lemma exponent bookkeeping (linear constraints)")
def feasible(z,sg,thv,a,b0,b):
    c=[sg<z, z<1, 0<thv<1, b0>a, a>z+sg+(1-thv), b>a+4*sg, b0<b<1]
    return all(c)
ex=dict(z=F(2,100),sg=F(1,100),thv=F(95,100),a=F(1,10),b0=F(12,100),b=F(2,10))
print("example feasible:",feasible(**ex),ex)
# random search for counterexamples to derived claims
random.seed(1)
bad=0;cnt=0
for _ in range(200000):
    zv,sg,thv,a,b0,b=[F(random.randint(1,999),1000) for _ in range(6)]
    if not feasible(zv,sg,thv,a,b0,b): continue
    cnt+=1
    claims={
     'top: a+(1-zeta)theta>1': a+(1-zv)*thv>1,
     'step4: a-sigma-(1-theta)>0': a-sg-(1-thv)>0,
     'step5(ii): a-sigma-(1-theta)(1-zeta)>0': a-sg-(1-thv)*(1-zv)>0,
     'overlap: b-a-2sigma>0': b-a-2*sg>0,
     'zeta+sigma<theta': zv+sg<thv,
    }
    if thv>F(9,10):
        for al in [F(1,100),F(1,2),F(1),2*(1-zv+sg)]:
            claims[f'(iii) alpha={al}']= a-sg*thv+al*(thv/2-F(1,5))>0
    for kk,v in claims.items():
        if not v: bad+=1; print("FAIL",kk,zv,sg,thv,a,b0,b); break
print("feasible samples",cnt,"failures",bad)

print("== near-product Step 4 gain identity (exact rationals)")
random.seed(2)
for trial in range(2000):
    m=random.randint(2,7)
    p=[F(random.randint(1,20)) for _ in range(m)]; S=sum(p); p=[x/S for x in p]
    c=[F(random.randint(-10,10)) for _ in range(m)]
    cb=sum(pi*ci for pi,ci in zip(p,c))
    dl=F(1,random.randint(2,50))
    sg=[(1 if ci>cb else -1 if ci<cb else 0) for ci in c]
    w=[pi*(1+dl*s_) for pi,s_ in zip(p,sg)]; Z=sum(w); pp=[x/Z for x in w]
    gain=sum(x*(ci-cb) for x,ci in zip(pp,c))
    A=sum(pi*abs(ci-cb) for pi,ci in zip(p,c))
    assert gain==dl*A/(1+dl*sum(pi*s_ for pi,s_ in zip(p,sg)))
    assert gain>=dl*A/2
print("gain identity and >= delta/2 * mean abs deviation: OK on 2000 random cases")

print("== near-product: harmonic sum  sum_{j<i} (1/d)/(1-s_j) vs -log(1-s_i)")
for d in [100,1000,10000]:
    i=(2*d)//3
    lhs=sum(F(1,d)/(1-F(jj-1,d)) for jj in range(1,i))
    print(d, float(lhs), float(-sp.log(1-sp.Rational(i-1,d))), "diff",float(lhs)-float(-sp.log(1-sp.Rational(i-1,d))))
