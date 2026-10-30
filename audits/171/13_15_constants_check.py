# Exact/symbolic checks of exponent arithmetic in Secs 13-15
import sympy as sp
from sympy import Rational as Q, log, exp, symbols, simplify, nsimplify
a=sp.symbols('a',positive=True)
L2=sp.log(2)
out=[]
# --- Sec 13 direct extraction 1/8,7/8 arithmetic
g,n=sp.symbols('g n',positive=True)
expr=Q(1,8)*(Q(1,2)+2*g/n)+Q(7,8)*(Q(1,2)+g/(2*n))
out.append(("direct extraction average = 1/2+11g/(16n)", sp.simplify(expr-(Q(1,2)+Q(11,16)*g/n))==0))
# --- Sec 13 cleaning codegree (1/4+3a-r)/(1-r) >= 1/4+a for 0<=r<a<1/4
r=sp.symbols('r',nonnegative=True)
diff=sp.together((Q(1,4)+3*a-r)/(1-r)-(Q(1,4)+a))
num=sp.factor(sp.numer(diff))
# numerator = 2a - 3r/4 + a r ; at r=a: 2a-3a/4+a^2>0 and linear in r
numexp=sp.expand(sp.numer(diff)*(1 if sp.denom(diff).subs({a:Q(1,100),r:Q(1,200)})>0 else -1))
out.append(("codegree numerator (sign-normalised)", numexp))
out.append(("  numerator at r=a positive", sp.simplify(numexp.subs(r,a))))
# --- Sec 14 deletion-test retained prob: rho^2 * A(F-c)/A(F) exponent
e1=2*(-L2+Q(8,100)*a)+(sp.log(4)-Q(4,10)*a)
out.append(("bin restriction failure exponent per c (should be -.24a)", sp.simplify(e1)))
# likelihood ratio exponent  (log4-.4a) + (j-2)(log2-.08a)  vs  j(log2-.08a)
j=sp.symbols('j',integer=True,positive=True)
lr=(sp.log(4)-Q(4,10)*a)+(j-2)*(L2-Q(8,100)*a)-j*(L2-Q(8,100)*a)
out.append(("LR exponent minus claimed (should be <=0)", sp.simplify(lr)))
# posterior B-set exponent per kh: (log2-.04a) + (1-.01a)(-log2+.02a)
pb=sp.expand((L2-Q(4,100)*a)+(1-Q(1,100)*a)*(-L2+Q(2,100)*a))
out.append(("posterior B-set exponent per kh", pb, "numeric at a=1e-3:", sp.N(pb.subs(a,Q(1,1000)))))
# clipped log-ratio mean via chord: L=1/4-eps, U=1, E R>=1/4+a => E log R >= log L + (1/4+a-L)*(-log L)/(1-L)
eps=sp.symbols('eps',nonnegative=True)
Lw=Q(1,4)-eps
chord=sp.log(Lw)+(Q(1,4)+a-Lw)*(-sp.log(Lw))/(1-Lw)
c0=sp.simplify(chord.subs(eps,0)+sp.log(4))
out.append(("chord bound E log R + log4 at eps=0 (should be >= .7a)", sp.simplify(c0), "coef", sp.N(c0/a)))
# with crude clip L=.24 only (no high-prob lower bound): fails for small a
Lc=Q(24,100)
chord2=sp.log(Lc)+(Q(1,4)+a-Lc)*(-sp.log(Lc))/(1-Lc)+sp.log(4)
out.append(("chord with L=.24 only, at a=1e-3 (negative => needs the 1/4-O(n^-eta0) w.h.p. bound)", sp.N(chord2.subs(a,Q(1,1000)))))
# --- Sec 15 high direct: D^{-L} with D>=1/2+g/(4n): log(2^L D^{-L}) = -L log(1+g/(2n)) <= -.4 g when L>=.9n, g<=n^{.005}
for nn in [10**4,10**6,10**9]:
    gg=sp.Float(nn)**sp.Float(0.005); LL=nn-int(nn**0.01)
    val=-LL*sp.log(1+gg/(2*nn))+Q(4,10)*gg
    out.append((f"high direct D^-L check n={nn} (should be <0)", sp.N(val)))
# capacity certificate: 2^l z^{l-1}/(l-1)! with l-1>=20z
import math
worst=0
for l1 in range(1,400):
    for zfrac in [1,0.5,0.1,0.01]:
        z=l1/20*zfrac
        v=math.log(2)*(l1+1)+(l1)*math.log(max(z,1e-300))-math.lgamma(l1+1)
        worst=max(worst, v/ l1)
out.append(("capacity: max over l of log(2^l z^{l-1}/(l-1)!)/(l-1) (neg => exp decay in l)", worst))
for t in out: print(*t)
