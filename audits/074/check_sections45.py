"""Exact checks of exponent bookkeeping in 074 (3D) Sections 4-5: (4.9), Lemma 4.3 fraction,
Lemma 4.5 (4.18), Lemma 4.6 cases, (4.31), (4.35), Lemma 5.3 inner equality, Prop 5.4 passage."""
import sympy as sp
H,h,h0,p,q,al,e,eps,Pi1,A,B,lam,kap,D,tau,v,Dp,s,mu,F,taus = sp.symbols(
 'H h h0 p q alpha e eps Pi1 A B lambda kappa Dstar tau v Dp s mu F taustar', real=True)
d=2-q; ps=p-sp.Rational(1,2)
# (4.9): log_N(k^2 m/n) >= 2(1-al) + H + Pi1 - d  with Pi1 >= p al - q
expr=2*(1-al)+H+(p*al-q)-d
assert sp.simplify(expr-(H+(p-2)*al))==0
# Lemma 4.3 bad-ball incidence fraction: K^C (n/m)^{3/2} * m / (n k) == K^C (n/(k^2 m))^{1/2}
n,m,k=sp.symbols('n m k',positive=True)
assert sp.simplify((n/m)**sp.Rational(3,2)*m/(n*k)-(n/(k**2*m))**sp.Rational(1,2))==0
# Lemma 4.5: density, Nikodym, multiplicity -> (4.18)
Fk=sp.Symbol('Fk'); 
rho_exp=-Fk+kap*(1+mu)-2*kap       # log_N rho
Eplus=-2*rho_exp                    # |E+|/|E| exponent
assert sp.simplify(Eplus-(2*Fk+2*kap*(1-mu)))==0
mplus=Fk-Eplus                      # m+/m exponent
assert sp.simplify(mplus-(-Fk-2*kap*(1-mu)))==0
# Lemma 4.6, e>0 at eps=0
ee=sp.Symbol('e',positive=True); Hh=sp.Symbol('H',positive=True)
s_=1-ee; t_=1+Hh+(ps-1)*ee
assert sp.simplify(t_-(s_+Hh+ps*ee))==0
half=(s_+t_)/2
assert sp.simplify(2*(1-half)-((2-ps)*ee-Hh))==0
gap=(ps-1)*ee-((2-ps)*ee-Hh)
assert sp.simplify(gap-((2*ps-3)*ee+Hh))==0
print('Lemma 4.6: gap =',sp.factor(sp.expand(gap)),'(>0 since p*>3/2)')
# (4.31)
lhs=(H+Pi1-sp.Symbol('e1')+A+B)-((d-1-h0)*A+(1+h0)*B)   # <=0
claimed=(q+h0)*A+h0*(1-B)+(H-h0)+Pi1-sp.Symbol('e1')
assert sp.simplify(lhs-claimed)==0
# (4.35): Delta theta^{-(1+h0)} (N theta)^H N^Pi  vs m = Delta N^{H+Pi}; theta=N^{-lam}
up=-(-(1+h0)*lam)+0; up=(1+h0)*lam - H*lam   # exponent of N beyond m
low=lam
print('(4.35) upper-lower exponent =',sp.simplify(up-low),' = -(H-h0)*lambda <0')
# Lemma 5.3 inner equality: parameter m N^{-h kap - Pi(kap) + y D}, times N^{h(kap-D)+Pi(kap)} -> m N^{(y-h)D}
y=1+h
assert sp.simplify((-h*kap+y*D)+(h*(kap-D))-D)==0
# Prop 5.4 passage
lo=1+Dp-Dp/taus; hi=Dp+(taus-Dp)/taus
assert sp.simplify(lo-hi)==0
sv=1-(1-taus)/taus*Dp
Fs=v*(sv-Dp)
assert sp.simplify(Fs-v/(1-taus)*(sv-taus))==0
# tau*>0 argument: s=1-h/4, gamma=0 => t=1+3h/4, mu=1
hp_=sp.Symbol('h',positive=True)
s2=1-hp_/4; t2=s2+hp_
assert sp.simplify(t2-(1+3*hp_/4))==0 and sp.simplify((s2+t2)/2-(1+hp_/4))==0
print('Sections 4-5 bookkeeping identities OK')
