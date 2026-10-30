import sympy as sp
from math import comb, factorial
from fractions import Fraction as Fr
# (5.16) simplex integral of |u|_1 over {u>=0, sum<=b} in dim d equals d b^{d+1}/((d+1) d!)
b=sp.symbols('b',positive=True)
for d in range(1,6):
    # integral of sum u_i = d * integral of u_1 ; integral of u_1 over simplex = b^{d+1}/(d+1)!
    val = d*b**(d+1)/sp.factorial(d+1)
    assert sp.simplify(val - d*b**(d+1)/((d+1)*sp.factorial(d)))==0
    # brute-force check via iterated integration for d<=3
    if d<=3:
        us=sp.symbols('u1:%d'%(d+1),nonnegative=True)
        expr=sum(us); lim=b
        integ=expr
        for i in reversed(range(d)):
            upper=b-sum(us[:i])
            integ=sp.integrate(integ,(us[i],0,upper))
        assert sp.simplify(integ-val)==0
print("5.16 simplex integral OK")
# (3.3): mean total degree over lattice points of degree<=r in Z^n_{>=0} -> n r/(n+1)
for n in range(1,6):
    r=400
    num=sum(k*comb(k+n-1,n-1) for k in range(r+1)); den=comb(r+n,n)
    print(n, float(Fr(num,den)/r), n/(n+1))
# Lemma 5.2 / (5.3): d! * D0 * C(q-D0+d,d) >= (k(q-D0))^d with D0=k^d, q>=D0
bad=0
for d in range(1,5):
  for k in range(1,5):
    D0=k**d
    for q in range(D0,D0+60):
        lhs=factorial(d)*D0*comb(q-D0+d,d); rhs=(k*(q-D0))**d
        if lhs<rhs: bad+=1
print("5.3 violations:",bad)
# (3.7): p*P-(K_B+E)-H/h with K_B=p*K+(n-1)E, H=p*(h t' L)-n h E
n,tp,h=sp.symbols('n tp h',positive=True)
pK,pL,E=sp.symbols('pK pL E')
P=pK+(n+1)*pL; KB=pK+(n-1)*E; H=h*tp*pL-n*h*E
print("3.7 residual:",sp.simplify(P-(KB+E)-H/h-(n+1-tp)*pL))
# (6.3): scale derivative => b.q=1, m sum b = t
# (5.17)/(I_*): I_* > 1/(n+1) for a_*>0 and Riemann sum limit
s,a=sp.symbols('s a',positive=True)
for nn in range(1,6):
    I=sp.integrate(s**nn*sp.exp((1-s)*a/(nn+1)),(s,0,1))
    ser=sp.series(I,a,0,2).removeO()
    print(nn, sp.simplify(ser))
