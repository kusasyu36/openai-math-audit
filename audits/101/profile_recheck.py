"""Recheck profile bodies with exact piecewise-polynomial integration (no quadrature)."""
import random, itertools
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps=30
exec(open('bodies.py').read().split("print('n, L_simplex")[0])  # helper defs only
def pl_exact(n, sl, ints, K0):
    # rho(t)=max(0,min_i(ints[i]+sl[i] t)), piecewise linear; integrate exactly via breakpoints
    sl=[Fr(x) for x in sl]; ints=[Fr(x) for x in ints]
    bps={Fr(0),Fr(1)}
    for i,j in itertools.combinations(range(3),2):
        if sl[i]!=sl[j]:
            t=(ints[j]-ints[i])/(sl[i]-sl[j])
            if 0<t<1: bps.add(t)
    for i in range(3):
        if sl[i]!=0:
            t=-ints[i]/sl[i]
            if 0<t<1: bps.add(t)
    bps=sorted(bps)
    def mom(pw, tk):
        tot=Fr(0)
        for a,b in zip(bps,bps[1:]):
            mid=(a+b)/2
            i=min(range(3),key=lambda i: ints[i]+sl[i]*mid)
            if ints[i]+sl[i]*mid<=0: continue
            # integrate t^tk (ints+sl t)^pw over [a,b] exactly by expanding
            c0,c1=ints[i],sl[i]
            from math import comb
            for m in range(pw+1):
                coef=comb(pw,m)*c0**(pw-m)*c1**m
                e=m+tk
                tot+=coef*(b**(e+1)-a**(e+1))/(e+1)
        return tot
    In1=mom(n-1,0); In1p=mom(n+1,0)
    m1=mom(n-1,1)/In1; m2=mom(n-1,2)/In1
    v0,dc0=K0
    vol=v0*mp.mpf(In1.numerator)/In1.denominator
    det=mp.mpf((m2-m1**2).numerator)/(m2-m1**2).denominator*dc0*(mp.mpf(In1p.numerator)/In1p.denominator/(mp.mpf(In1.numerator)/In1.denominator))**(n-1)
    return (det/vol**2)**(1/(2*mp.mpf(n)))
def Lcross(n): n = M(n); return mp.factorial(n)**(1/n)/mp.sqrt(2*(n+1)*(n+2))
_src=open('bodies.py').read()
exec(_src[_src.index('def base('):_src.index('random.seed(1)')])
random.seed(1)
worst=(0,None)
for n in range(2,8):
    for kind in ('simplex','cube','ball','cross'):
        K0=base(kind,n-1)
        for _ in range(400):
            sl=sorted([random.uniform(-3,3) for _ in range(3)],reverse=True)
            ints=[random.uniform(0.2,1.5) for _ in range(3)]
            r=pl_exact(n,sl,ints,K0)/Ls[n]
            if r>worst[0]: worst=(r,(n,kind,sl,ints))
print('exact piecewise-linear profiles: max L_K/L_simplex =',mp.nstr(worst[0],12),' (ratio-1 =',mp.nstr(worst[0]-1,5),')',worst[1][:2])
