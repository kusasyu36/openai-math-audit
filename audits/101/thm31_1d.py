"""End-to-end 1D check of Theorem 3.1 remainder (3.2): 2h-2-log Var >= E Y^2/500 + omega(B)/1000,
with A = q^{-1}(P J), P chosen so that E A = 0 (Prop 2.3), X=E[xA], Y=A-xX, B=X^2."""
import mpmath as mp
mp.mp.dps=30
SQ2=mp.sqrt(2)
def q(a): return mp.npdf(a)/(mp.erfc(a/SQ2)/2)
def qinv(lam):
    lam=mp.mpf(lam)
    a0 = lam if lam>2 else (-mp.sqrt(2*abs(mp.log(lam))) if lam<0.05 else mp.mpf(0))
    return mp.findroot(lambda a: mp.log(q(a))-mp.log(lam), a0)
def omega(l):
    return 3*(l-1)**2/(37+3*l)*(mp.mpf(3)/20+mp.mpf(49)/100*37/(37+3*l))*min(mp.mpf(13)/100, (mp.mpf(32)/25)/mp.sqrt(l) if l>0 else mp.inf)
# Gauss-Hermite (probabilists') nodes
N=40
xs,ws=[],[]
nodes=mp.polyroots(mp.taylor(lambda t: mp.hermite(N,t),0,N)[::-1], maxsteps=200, extraprec=400)
for r in nodes:
    r=mp.re(r); x=r*SQ2
    w=2**(N-1)*mp.factorial(N)*mp.sqrt(mp.pi)/(N**2*mp.hermite(N-1,r)**2)/mp.sqrt(mp.pi)
    xs.append(x); ws.append(w)
E=lambda vals: mp.fsum(w*v for w,v in zip(ws,vals))
print('GH check E x^2 =',mp.nstr(E([x*x for x in xs]),12),' E x^4=',mp.nstr(E([x**4 for x in xs]),12))
def run(V, lo, hi, label):
    Z=mp.quad(lambda y: mp.e**(-V(y)), [lo,0,hi])
    f=lambda y: mp.e**(-V(y))/Z
    mean=mp.quad(lambda y: y*f(y), [lo,0,hi]); var=mp.quad(lambda y: (y-mean)**2*f(y), [lo,0,hi])
    h=mp.quad(lambda y: f(y)*(V(y)+mp.log(Z)), [lo,0,hi])
    gap2=2*h-2-mp.log(var)
    # transport: F(x)=Q(Phi(x)); J = phi(x)/f(F(x))
    Js=[]
    for x in xs:
        if x<0:
            fn=lambda y: mp.log(mp.quad(f,[lo,y]))-mp.log(mp.ncdf(x))
        else:
            fn=lambda y: mp.log(mp.quad(f,[y,hi]))-mp.log(mp.ncdf(-x))
        # bracket by scanning
        a,b=mp.mpf(lo)+mp.mpf(10)**-6,mp.mpf(hi)-mp.mpf(10)**-6
        y=mp.findroot(fn,(a,b),solver='anderson') if False else None
        # bisection
        fa=fn(a)
        for _ in range(110):
            mid=(a+b)/2; fm=fn(mid)
            if (fm>0)==(fa>0): a,fa=mid,fm
            else: b=mid
        y=(a+b)/2
        Js.append(mp.npdf(x)/f(y))
    # choose P with E q^{-1}(P J) = 0
    g=lambda logP: E([qinv(mp.e**logP*J) for J in Js])
    logP=mp.findroot(g, -mp.log(E(Js))+mp.log(q(0)))
    A=[qinv(mp.e**logP*J) for J in Js]
    X=E([x*a for x,a in zip(xs,A)]); Y=[a-x*X for x,a in zip(xs,A)]
    B=X**2; EY2=E([y*y for y in Y])
    rem=EY2/500+omega(B)/1000
    print(f'{label}: 2*gap={mp.nstr(gap2,8)}  EA={mp.nstr(E(A),3)}  B={mp.nstr(B,6)}  EY^2={mp.nstr(EY2,6)}  remainder={mp.nstr(rem,6)}  OK={gap2>=rem}')
for tau in (1, mp.mpf(1)/2, mp.mpf(1)/4):
    eps=mp.mpf(1)/50
    run(lambda y, tau=tau, eps=eps: y+mp.e**(-y/tau)+eps*y*y/2, -6, 60, f'V=y+exp(-y/{float(tau)})+.01y^2')
run(lambda y: y*y/2+mp.log(1+mp.e**(3*y)), -12, 12, 'V=y^2/2+softplus(3y)')
run(lambda y: y**4/4+y*y/2, -6, 6, 'V=y^4/4+y^2/2 (not (2.1) at inf, smooth)')
