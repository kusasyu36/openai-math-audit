"""Independent high-precision sampling check of Lemma 2.1, 2.2, 4.1 (paper 101).
Own code; does not import repository code.  Run: .venv/bin/python -I scalar_check.py
"""
import mpmath as mp
mp.mp.dps = 50
F = mp.mpf
al = F(3)/40; s = 1/(1-al); beta = F(7)/25; c = F(3)/2; gam = c/2
eps = F(3)/20; kap = F(49)/100
Dst = 1/(4*(1-eps)); Ast = s+1+s/8; eta = s*beta+(1+s/8)*gam; nu = s*beta**2+(9*s/8)*gam**2
SQ2 = mp.sqrt(2)

def funcs(a):
    a = F(a)
    phi = mp.npdf(a)
    Pm = mp.erfc(a/SQ2)/2           # Phi(-a)
    q = phi/Pm
    if a < 0:
        t = -mp.log1p(-mp.erfc(-a/SQ2)/2)
    else:
        t = -mp.log(Pm)
    d = q - a
    k = t/q
    b = 1 - k*d
    l = k - q + a*b
    lp = 2*b - q*d + a*l
    p = b + k*l - l**2/al
    E = p + b*(F(3)/4 - b) - l**2/al
    h0 = k*(b - F(1)/2) - 2*(b - F(7)/25)*l/al
    slope = 2*(1-q*d)/d**2 - 2*b*l/(d*(1-b**2))
    return dict(a=a,q=q,d=d,k=k,b=b,l=l,lp=lp,p=p,E=E,h0=h0,slope=slope,t=t)

BINS = [  # name, lo, hi, bmin,bmax,Lmax,Dmax,Kmin,Kmax,Jmax,ell1,m0
 ('O', -15, -3, 0,.09,.60,None,.37,.56,.72,.13,0),
 ('P1', -3, -2, .083,.15,1.01,4.3,.29,.46,.72,.13,.19),
 ('R1', -2, -1, .14,.23,1.19,3.00,.29,.371,.68,.13,.38),
 ('U1', -1, .5, .22,.342,1.19,1.72,.29,.398,.61,.123,.68),
 ('V1', .5, 4, .337,.45,.83,.87,.31,.4,.50,.060,.68),
 ('T1', 4, 30, .445,.5,.20,.640,.32,.4,.48,.022,.68)]

def grid(lo, hi, h):
    n = int(round((hi-lo)/h))
    return [F(lo) + (F(hi)-F(lo))*i/n for i in range(n+1)]

results = {}
def upd(key, val):
    if key not in results or val < results[key][0]:
        results[key] = (val, cur)

gl = []
for lo, hi, h in [(-15,-6,.01),(-6,-3.5,.002),(-3.5,-0.5,.001),(-0.5,1,.001),(1,5,.002),(5,30,.01)]:
    gl += grid(lo, hi, h)
# Lemma 2.1/2.2 + 4.1
cur = None
for a in gl:
    f = funcs(a); cur = float(a)
    q,d,k,b,l,lp,p,E,h0 = (f[x] for x in 'q d k b l lp p E h0'.split())
    upd('L2.1 d>0', d); upd('L2.1 1-qd>0 (q\'<1)', 1-q*d); upd('L2.1 d^2-(1-qd)>=0', d**2-(1-q*d))
    upd('L2.1 qd>0', q*d)
    upd('L2.2 b>0', b); upd('L2.2 1/2-b>0', F(1)/2-b); upd('L2.2 l>0', l)
    upd('4.2 p>0', p)
    upd('4.2 .67sqrt(p)-|h0|', F(67)/100*mp.sqrt(p)-abs(h0))
    if a <= 0: upd('4.2 p-1/(7+a^2) (a<=0)', p-1/(7+a**2))
    if a >= 0: upd('4.2 p-.25 (a>=0)', p-F(1)/4)
    if a >= 4: upd('4.2 p-.46 (a>=4)', p-F(46)/100)
    z = lp/l; hst = 1+3*b+k*z; chi = l/d
    if a <= -1:
        upd('4.6 z+.08 (a<=-1)', z+F(8)/100); upd('4.6 .70-z (a<=-1)', F(70)/100-z)
        upd('4.6 h*-.8 (a<=-1)', hst-F(8)/10); upd('4.6 1.77-h* (a<=-1)', F(177)/100-hst)
        upd('4.6 .070-chi (a<=-1)', F(7)/100-chi)
    if a >= -1:
        upd('4.6 z+.475 (a>=-1)', z+F(475)/1000); upd('4.6 -z (a>=-1)', -z)
        upd('4.6 h*-.8 (a>=-1)', hst-F(8)/10); upd('4.6 1.77-h* (a>=-1)', F(177)/100-hst)
        upd('4.6 .097-chi (a>=-1)', F(97)/1000-chi)
    for (nm, lo, hi, bmin,bmax,Lmax,Dmax,Kmin,Kmax,Jmax,ell1,m0) in BINS:
        if not (lo <= a <= hi): continue
        D = Dst/p; L = l/al; K0 = D*E; J0 = D*c*b
        upd(f'{nm} b-bmin', b-F(bmin)); upd(f'{nm} bmax-b', F(bmax)-b)
        upd(f'{nm} Lmax-L', F(Lmax)-L)
        if Dmax is None: upd(f'{nm} .295(7+a^2)-D', F(295)/1000*(7+a**2)-D)
        else: upd(f'{nm} Dmax-D', F(Dmax)-D)
        upd(f'{nm} K0-Kmin', K0-F(Kmin)); upd(f'{nm} Kmax-K0', F(Kmax)-K0)
        upd(f'{nm} Jmax-J0', F(Jmax)-J0)
        upd(f'{nm} slope-m0', f['slope']-F(m0))
        # matrix (4.4)
        N11 = Ast*b**2-eta*b+nu/4+l**2/al + D*E*E
        N12 = nu/2-eta*b + D*E*c*b
        N22 = nu + D*(c*b)**2
        R11 = F(43)/100 - N11 - F(ell1); R12 = F(13)/100 - N12; R22 = F(272)/100 - N22 - F(3)/4
        upd(f'{nm} (4.4) R11', R11); upd(f'{nm} (4.4) det', R11*R22-R12**2)

for kkey in sorted(results):
    v, at = results[kkey]
    flag = 'OK ' if v >= 0 else 'VIOLATION'
    print(f'{flag} {kkey:32s} min={mp.nstr(v,6):>12s} at a={at}')

# endpoint facts
qf = lambda a: funcs(a)['q']
print('q(-1)=', mp.nstr(qf(-1),12), ' >.287'); print('q(-2)=', mp.nstr(qf(-2),12), ' <.056')
print('q(-3)=', mp.nstr(qf(-3),12), ' <.005'); print('q(4)=', mp.nstr(qf(4),12), ' >4')
for (x,y,thr) in [(-1,-2,1.62),(-1,-3,4),(4,-2,4),(4,-3,6.5)]:
    print(f'log q({x})/q({y}) = {mp.nstr(mp.log(qf(x)/qf(y)),10)} > {thr}')
print('q(0)=',qf(0),' sqrt(2/pi)=',mp.sqrt(2/mp.pi))
print('q(10)=',mp.nstr(qf(10),15),' paper 10.098093234')
print('k(-6)=',mp.nstr(funcs(-6)['k'],15),' paper .16237766')
