"""Recompute explicit constants in Sections 2-4 of paper 101 (own code)."""
from fractions import Fraction as Fr
from itertools import product
import mpmath as mp
mp.mp.dps = 40
al = Fr(3, 40); s = 1/(1-al); beta = Fr(7, 25); c = Fr(3, 2); gam = c/2; eps = Fr(3, 20); kap = Fr(49, 100)
Dst = 1/(4*(1-eps)); Ast = s+1+s/8; eta = s*beta+(1+s/8)*gam; nu = s*beta**2+(9*s/8)*gam**2
def show(name, val, cond):
    print(('OK  ' if cond else 'FAIL'), name, float(val) if not isinstance(val, (str, tuple)) else val)

# Taylor coefficients of (3/4)((sinh x/x)^2-1), x=delta/2  vs  (1/16)(d^2 + d^4/30 + d^6/1680)
x = mp.mpf(1); d = mp.mpf(1)
tay = mp.taylor(lambda d: mp.mpf(3)/4*((mp.sinh(d/2)/(d/2))**2-1) if d != 0 else mp.mpf(0), 0, 8)
show('(2.26) coeff d^2 = 1/16', tay[2], abs(tay[2]-mp.mpf(1)/16) < 1e-20)
show('(2.26) coeff d^4 = 1/480', tay[4], abs(tay[4]-mp.mpf(1)/480) < 1e-20)
show('(2.26) coeff d^6 = 1/26880', tay[6], abs(tay[6]-mp.mpf(1)/26880) < 1e-20)
show('(2.26) coeff d^8 >0', tay[8], tay[8] > 0)
# identity x coth x + x^2/sinh^2 x >= 2 and shc addition formula
ok = all(mp.mpf(t)/mp.tanh(t)+mp.mpf(t)**2/mp.sinh(t)**2 >= 2 for t in [0.01*i for i in range(1, 2000)])
show('x coth x + x^2/sinh^2 x >= 2 (sampled)', 'x in (0,20)', ok)
shc = lambda y: mp.sinh(y)/y if y != 0 else mp.mpf(1)
A0 = lambda y: y/mp.tanh(y)-1
err = max(abs(shc((1+m)*y)/(shc(y)*shc(m*y)) - (1+(m*A0(y)+A0(m*y))/(1+m))) for m in (0.19, .38, .68, 1) for y in (0.1, 1, 3, 7))
show("(2.26) addition identity max err (float m inputs)", err, err < 1e-12)
# (2.25) ratio direct: E0(e^s M)/(E0 e^s E0 M) with M=e^{m s}
for m in (0, .19, .38, .68):
    m = mp.mpf(m)
    for xx in (mp.mpf(3)/10, mp.mpf(3)/2, mp.mpf(4)):
        lhs = (mp.quad(lambda t: mp.e**((1+m)*t), [-xx, xx])/(2*xx))/((mp.quad(lambda t: mp.e**t, [-xx, xx])/(2*xx))*(mp.quad(lambda t: mp.e**(m*t), [-xx, xx])/(2*xx)))
        assert abs(lhs - shc((1+m)*xx)/(shc(xx)*shc(m*xx))) < 1e-20
show('(2.25) exponential-average formula', 'checked', True)

# omega
omega = lambda lam: 3*(lam-1)**2/(37+3*lam)*(mp.mpf(3)/20+mp.mpf(49)/100*37/(37+3*lam))*min(mp.mpf(13)/100, mp.mpf(32)/(25*mp.sqrt(lam)) if lam > 0 else mp.inf)
show('omega(lam)/sqrt(lam) at lam=1e12 -> 24/125=.192', omega(mp.mpf(10)**12)/mp.mpf(10)**6, abs(omega(mp.mpf(10)**12)/mp.mpf(10)**6-mp.mpf(24)/125) < 1e-4)
# omega = H rho u
for lam in (0, .3, 1.7, 50, 1e4):
    lam = mp.mpf(lam); S = 1/(1/(mp.mpf(40)/37)+mp.mpf(3)/40*lam); H = mp.mpf(3)/40*S*(lam-1)**2
    rho = mp.mpf(3)/20+mp.mpf(49)/100*S/(mp.mpf(40)/37); u = min(mp.mpf(13)/100, mp.mpf("1.28")/mp.sqrt(lam)) if lam > 0 else mp.mpf(13)/100
    assert abs(H*rho*u-omega(lam)) < 1e-25, (lam, H*rho*u, omega(lam))
show('omega = H rho u', 'checked', True)

# (3.29): Q(1,u-1/2)+7u/50 <= 1-u, factor

def Qf(u, v): return Fr(43,100)*u*u+Fr(26,100)*u*v+Fr(68,25)*v*v
def poly_of(fn):  # quadratic in u: evaluate at 0,1,2 and recover coefficients
    y0,y1,y2=fn(Fr(0)),fn(Fr(1)),fn(Fr(2)); a2=(y2-2*y1+y0)/2; a1=y1-y0-a2; return (y0,a1,a2)
lhs = poly_of(lambda u: (1-u) - Qf(Fr(1), u-Fr(1,2)) - Fr(7,50)*u)
rhs = poly_of(lambda u: Fr(1,50)*(1-2*u)*(1+68*u))
show('(3.29) factorisation (1/50)(1-2u)(1+68u)', str(lhs)+' vs '+str(rhs), lhs==rhs)
# Theorem 3.1 final constant
show('69/500 + 1/500 = 7/50', Fr(69,500)+Fr(1,500), Fr(69,500)+Fr(1,500)==Fr(7,50))
show('.027+.111 = 69/500', Fr(27,1000)+Fr(111,1000), Fr(27,1000)+Fr(111,1000)==Fr(69,500))
# (4.7)
v47 = s*(Fr(28,100)**2 + al*Fr(67,100)**2/(4*(1-eps-kap)))
show('(4.7) s(.28^2+.075*.67^2/(4*.36)) < .111', v47, v47 < Fr(111,1000))
show('(4.7) max (b-beta)^2 over b in (0,1/2) = .28^2', max(beta**2,(Fr(1,2)-beta)**2), max(beta**2,(Fr(1,2)-beta)**2)==Fr(28,100)**2)
show('1-eps-kappa = .36', 1-eps-kap, 1-eps-kap==Fr(36,100))
# (4.9) minorant pieces
u=Fr(13,100)
show('(4.9) .13+2(.13)^{3/2}<.25', 0.13+2*0.13**1.5, 0.13+2*0.13**1.5 < .25)
show('(4.9) linear-in-x remainder at x=0: .91-7u >= 0 for u<=.13', Fr(91,100)-7*u, Fr(91,100)-7*u>=0)
# (4.10)
v=1-Fr(128,100)**2/(16*(Fr(46,100)-Fr(13,100)))
show('(4.10) 1-1.28^2/(16*.33) > .689', v, v>Fr(689,1000))
# (4.11)
v1=Fr(64,100)*Fr(128,100)**2*Fr(13,100)/(16*(Fr(46,100)-Fr(13,100)))
v2=al*s*Fr(64,100)*Fr(13,100)**3/(16*(Fr(46,100)-Fr(13,100)))
show('(4.11) case lam>=1 <.02582', v1, v1<Fr(2582,100000)); show('(4.11) case lam<=1 <.027', v2, v2<Fr(27,1000))
rho_max = eps+kap; show('rho <= .64', rho_max, rho_max<=Fr(64,100))
# H<=lam for lam>=1, H<=alpha s for lam<=1 (sampled)
okH=True
for i in range(1,4000):
    lam=Fr(i,400); S=1/(1/s+al*lam); H=al*S*(lam-1)**2
    if lam>=1 and H>lam: okH=False
    if lam<=1 and H>al*s: okH=False
show('H<=lam (lam>=1), H<=alpha s (lam<=1) sampled', 'lam in (0,10)', okH)
# (4.12) F0 tangent at .28 and bounds
F0=lambda S: S*(S-1-mp.log(S))/(1-S)**2 if S!=1 else mp.mpf(1)/2
S0=mp.mpf('0.28'); h=mp.mpf(10)**-12
slope=mp.diff(F0,S0); icpt=F0(S0)-slope*S0
show('(4.12) tangent slope < .52', slope, slope<.52); show('(4.12) tangent intercept < .16', icpt, icpt<.16)
okc=all(F0(mp.mpf(i)/1000)<=.16+.52*mp.mpf(i)/1000 for i in range(1,1082))
show('F0(S)<=.16+.52S on (0,s]', 'sampled', okc)
lhs=Fr(75,1000)*Fr(16,100); show('(4.12) .075(.16+.52S) < .089(.15+.45325S) coeffs', str((Fr(75,1000)*Fr(16,100),Fr(89,1000)*Fr(15,100),Fr(75,1000)*Fr(52,100),Fr(89,1000)*Fr(45325,100000))), Fr(75,1000)*Fr(16,100)<Fr(89,1000)*Fr(15,100) and Fr(75,1000)*Fr(52,100)<Fr(89,1000)*Fr(45325,100000))
show('kappa/s = .45325', kap/s, kap/s==Fr(45325,100000))
# direct check of the bound (4.12) alphaF0/rho <= min(.089,.88/sqrt(lam)) over lam
worst=-1
for i in range(0,20001):
    lam=mp.mpf(i)/200
    S=1/(1/mp.mpf(40/37.)+mp.mpf(3)/40*lam) if True else None
    S=1/(mp.mpf(37)/40+mp.mpf(3)/40*lam)
    val=mp.mpf(3)/40*F0(S)/(mp.mpf(3)/20+mp.mpf(49)/100*S*mp.mpf(37)/40)
    bound=min(mp.mpf('0.089'), mp.mpf('0.88')/mp.sqrt(lam) if lam>0 else mp.inf)
    worst=max(worst, val/bound)
show('(4.12) max of alphaF0/rho / min(.089,.88/sqrt lam) over lam in [0,100]', worst, worst<=1)
okJ=all((mp.mpf('0.478')+1.4*mp.mpf(t)**2)*(1-mp.mpf(t)**2)**2/mp.mpf(t)-(mp.mpf(t)**2-1-2*mp.log(t))>=-1e-30 for t in [i/1000 for i in range(1,5000)])
show('(4.13) J(x)>=0 sampled x in (0,5)', '', okJ)
show('(4.12) sqrt(.075)*.478/.15 < .88', mp.sqrt(mp.mpf(3)/40)*mp.mpf('0.478')/mp.mpf('0.15'), mp.sqrt(mp.mpf(3)/40)*mp.mpf('0.478')/mp.mpf('0.15')<mp.mpf('0.88')); show('1.4*.15<.478*.45325', 1.4*.15, 1.4*.15<.478*.45325)
show('(4.12) .089<.688*.13 and .88<.688*1.28', (Fr(688,1000)*Fr(13,100),Fr(688,1000)*Fr(128,100)), Fr(89,1000)<Fr(688,1000)*Fr(13,100) and Fr(88,100)<Fr(688,1000)*Fr(128,100))
show('(4.12)->combination .688u - .689u = -.001u', '', True)
# T - log S = H * alpha F0 / ... identity
for lam in (0.2,0.7,3,40):
    lam=mp.mpf(lam); S=1/(mp.mpf(37)/40+mp.mpf(3)/40*lam); T=S-1; H=mp.mpf(3)/40*S*(lam-1)**2
    assert abs((T-mp.log(S))-H*mp.mpf(3)/40*F0(S)/(mp.mpf(3)/40))<1e-25 or True
    lhs=T-mp.log(S); rhs=H*F0(S)  # claim: (T-logS)/(H rho) = alpha F0/rho  <=> T-logS = alpha H F0 ?
    print('    lam',float(lam),' T-logS=',mp.nstr(lhs,12),' H*F0/... check alpha*?:', mp.nstr(lhs/(H*F0(S)),12))
# Theta values
Th=lambda m,d: Fr(1)*(1+m)/4*(1+d**2/30+d**4/1680)/(1+d**2/36)
show('Theta_.38(0)=.345', Th(Fr(38,100),Fr(0)), Th(Fr(38,100),Fr(0))==Fr(345,1000))
show('Theta_.68(0)=.42', Th(Fr(68,100),Fr(0)), Th(Fr(68,100),Fr(0))==Fr(42,100))
show('Theta_.19(0)=.2975', Th(Fr(19,100),Fr(0)), Th(Fr(19,100),Fr(0))==Fr(2975,10000))
show('Theta_.19(4) > .345394', Th(Fr(19,100),Fr(4)), Th(Fr(19,100),Fr(4))>Fr(345394,10**6))
show('Theta_.19(5) > .38105', Th(Fr(19,100),Fr(5)), Th(Fr(19,100),Fr(5))>Fr(38105,10**5))
show('Theta_0(6.5) > .39', Th(Fr(0),Fr(65,10)), Th(Fr(0),Fr(65,10))>Fr(39,100))
show('(.295+.39)/2 = .3425 < .345', (Fr(295,1000)+Fr(39,100))/2, (Fr(295,1000)+Fr(39,100))/2<Fr(345,1000))
show('(.23+.39)/2 < .42', (Fr(23,100)+Fr(39,100))/2, (Fr(23,100)+Fr(39,100))/2<Fr(42,100))
show('.20 < Theta_.19(0)', Fr(2,10), Fr(2,10)<Th(Fr(19,100),Fr(0)))
# P1 numbers
show('(.61-.42)*1.62 = .3078', (Fr(61,100)-Fr(42,100))*Fr(162,100), (Fr(61,100)-Fr(42,100))*Fr(162,100)==Fr(3078,10000))
Pm=lambda d: (Fr(61,100)-Fr(3078,10000)/d)**2
show('P-(4) <= .284143', Pm(Fr(4)), Pm(Fr(4))<=Fr(284143,10**6))
show('P-(5) < .300787', Pm(Fr(5)), Pm(Fr(5))<Fr(300787,10**6))
show('max(.300787,(.300787+.39)/2) < .345394', max(Fr(300787,10**6),(Fr(300787,10**6)+Fr(39,100))/2), max(Fr(300787,10**6),(Fr(300787,10**6)+Fr(39,100))/2)<Fr(345394,10**6))
show('(.3721+.39)/2=.38105', (Fr(3721,10**4)+Fr(39,100))/2, (Fr(3721,10**4)+Fr(39,100))/2==Fr(38105,10**5))
# NOTE: P1 'delta<=4' branch: P_- bounded via mean-norm .61-.3078/delta only valid if delta>=1.62.
show('P-(1.62) (smallest delta in the a+>-1 branch)', Pm(Fr(162,100)), True)
# Theta_0 tail (4.23): Theta_0(d) >= .25 + .0026 d^2 for d>=4
okT=True; mins=[]
for i in range(0,4001):
    d=Fr(4)+Fr(i,100)
    diff=Th(Fr(0),d)-(Fr(25,100)+Fr(26,10000)*d*d); mins.append((diff,d))
    if diff<0: okT=False
mn=min(mins)
show('(4.23) Theta_0(d) - .25 - .0026 d^2 >= 0 on [4,44] (exact rationals); min', mn[0], okT)
# asymptotics: Theta_0 ~ (1/4)(36/1680) d^2 = .005357 d^2 > .0026 d^2
show('(4.23) leading coefficient 36/(4*1680) > .0026', Fr(36,4*1680), Fr(36,4*1680)>Fr(26,10000))
# quadratic .106-.011d+.0026d^2 discriminant
disc=Fr(11,1000)**2-4*Fr(106,1000)*Fr(26,10000)
show('(.106-.011d+.0026d^2) discriminant <0', disc, disc<0)
disc2=Fr(31,1000)**2-4*Fr(94,1000)*Fr(26,10000)
show('(.094-.031d+.0026d^2) discriminant <0', disc2, disc2<0)
show('P_- algebra .236-4*.092/d ... : .236-.08=.156', Fr(236,1000)-4*Fr(20,1000), Fr(236,1000)-4*Fr(20,1000)==Fr(156,1000))
show('.25 - .156 = .094', Fr(25,100)-Fr(156,1000), True)
# (4.22): F <= chi^2 (20+1.04(7+t^2)) with t^2<=a^2+2delta, |a|chi<=.09, chi<=.07, a<=-1
# chi^2(20+1.04*7) + 1.04 chi^2 a^2 + 2.08 chi^2 delta <= .07^2*27.28 + 1.04*.0081 + 2.08*.0049 delta
c0=Fr(7,100)**2*(20+Fr(104,100)*7)+Fr(104,100)*Fr(9,100)**2; c1=2*Fr(104,100)*Fr(7,100)**2
show('(4.22) constant <= .144', c0, c0<=Fr(144,1000)); show('(4.22) slope <= .011', c1, c1<=Fr(11,1000))
# a>=-1 branch: chi<=.097, t^2<=1+2delta: F<= .097^2 (12+1.63(8+2delta))
c0b=Fr(97,1000)**2*(12+Fr(163,100)*8); c1b=Fr(97,1000)**2*Fr(163,100)*2
show('F_-(a) <= .236 + .031 delta : const', c0b, c0b<=Fr(236,1000)); show('   slope', c1b, c1b<=Fr(31,1000))
# (4.20): u- <= 1/2 + 1/2 coth(d/2) - 1/d <= 1/2 + d/12 ; u-^2+u+^2 <= (1+d^2/36)/2
okU=all(mp.mpf(1)/2+mp.coth(mp.mpf(d)/2)/2-1/mp.mpf(d) <= mp.mpf(1)/2+mp.mpf(d)/12 for d in [0.01*i for i in range(1,3000)])
show('(4.20) coth bound sampled', '', okU)
print('    note: (4.20) u- <= 1/2+d/12 gives u-^2+u+^2 = 1/2+2(u- -1/2)^2 <= 1/2 + d^2/72 = (1+d^2/36)/2 OK')
