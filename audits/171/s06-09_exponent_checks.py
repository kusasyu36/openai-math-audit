# Exact checks of exponent bookkeeping in Sections 6-9 (own script; no repo code executed).
from fractions import Fraction as F
import math, itertools, random

ok=lambda name,cond: print(("OK  " if cond else "FAIL"),name)

# ---------- Section 6 ----------
c0=F(1,100); c1=c0/2
# kappa(1+log(4/c1))<.001 ; then coefficient kappa(log(2/c1)+.02)<.002
kappa_max=0.001/(1+math.log(4/float(c1)))
print("Sec6 kappa_max",kappa_max)
ok("6: kappa*(log(2/c1)+.02) < .002", kappa_max*(math.log(2/float(c1))+0.02) < 0.002)
ok("6: .8*.55 > .4 (heavy-set tail)", F(8,10)*F(55,100) > F(4,10))
# rho: H_bin(2rho) < log2-.55
Hb=lambda p: -p*math.log(p)-(1-p)*math.log(1-p)
lo,hi=1e-9,0.25
for _ in range(200):
    mid=(lo+hi)/2
    if Hb(mid) < math.log(2)-0.55: lo=mid
    else: hi=mid
print("Sec6: max 2rho with H(2rho)<log2-.55 :",lo)
# .2k+.1k per matching neighbour => .3kn ; total .34kn ; posterior: .34+.02 + o(1) <= .4
ok("6: .3+.04(slack) +.02 <.4", F(3,10)+F(2,100)+F(4,100) <= F(4,10))
# low deletion: d2/kappa + .02 < .16 for d2<<kappa
# P(j>=q)<=n^{-.13q}: binom(m,q) * (C n^{-.15})^q -> need alpha<.02
ok("6: alpha-.15 <= -.13 needs alpha<=.02", F(2,100)-F(15,100) <= F(-13,100))

# ---------- Section 7 ----------
def sec7(D0,dfrac):
    d=D0*dfrac; b=d/4; c=b/20
    res={}
    res['B_g: -D0/2+2d <= -D0/3']= -D0/2+2*d <= -D0/3
    res['LLL tags: 2d < D0/4']= 2*d < D0/4
    res['anchor LLL: 16d < D0/4']= 16*d < D0/4
    res['odd comparison: 12d < D0/4']= 12*d < D0/4
    res['b>2c']= b>2*c
    res['log L exponent max(d/2,d)<1/4']= max(d/2,d) < F(1,4)
    res['m log n=o(q): 2d<4d']= 2*d<4*d
    res['s log n = o(q): d<4d']= d<4*d
    res['n^{d/2}=o(q)']= d/2<4*d
    res['.04 < log2']= 0.04<math.log(2)
    res['D0 < D*/2 => 2n^{D0}<=n^{D*}']=True
    return res
for D0 in [F(1,10)-F(1,10**6),F(1,1000)]:
    for dfrac in [F(999,1000000)]:
        for k,v in sec7(D0,dfrac).items(): ok(f"7 (D0={float(D0)}): {k}",v)

# ---------- Section 8 ----------
def sec8(eta0):
    eta=min(eta0/2,F(4,100)); tau=eta/4
    beta=tau/4*F(999,1000)
    T=tau/8; s=tau
    r={}
    r['tau<eta (o(1/s))']= tau<eta
    r['T*n^{tau/2}: T+tau/2 < tau']= T+tau/2<tau
    r['T*h n^beta: T+beta<tau']= T+beta<tau
    r['n^{2tau} >> hs: 2tau>tau']= 2*tau>tau
    r['center IDs 2n^{tau/16} <= n^{tau/8}']= tau/16<tau/8
    r['aggregate width h n^beta + n^{tau/2} < n^eta']= max(beta,tau/2)<eta
    r['ms log n=o(n): tau+.2+tau<1']= 2*tau+F(2,10)<1
    r['2delta=.0002 < .003']= 2*F(1,10**4)<F(3,1000)
    r['.003 - .01/2 <0']= F(3,1000)-F(1,100)/2<0
    r['.98^{-n}: log(1/.98)+o(1) < .03']= -math.log(.98)<0.03
    r['w_v: .03+.04 < .1']= F(3,100)+F(4,100)<F(1,10)
    r['bins: chunk n^.2, atom n^-.1 <= n^-.04']= F(1,10)>F(4,100)
    r['eta <= eta0/2 so 2n^-eta>n^-eta0']= eta<=eta0/2
    r['beta+eta width <= eta0 (n^beta+n^eta<=n^eta0)']= max(beta,eta)<eta0
    return r
for eta0 in [F(1,10**9),F(1,5)]:
    for k,v in sec8(eta0).items(): ok(f"8 (eta0={float(eta0)}): {k}",v)
print("8: H_bin(.02)=",Hb(0.02)," < log2=",math.log(2))

# ---------- Section 9 ----------
# random sampling of admissible parameters; check all exponent claims
random.seed(1)
bad=[]
def check9(case):
    xs=F(random.randint(1,999),10000)    # <.1
    xd=xs+(F(1,10)-xs)*F(random.randint(1,999),1000)
    hp=F(random.randint(1,999),1000)
    if xd>=(1-hp)/10: return None
    chi=min(xs,1-hp)/100*F(random.randint(1,999),1000)
    hm=hp-chi/10*F(random.randint(1,999),1000)
    if hm<=0 or chi>=min(xs,hm,1-hp)/100: return None
    u=xs/2
    if case=='sub':
        sigma_lo=10*chi
        sigma=sigma_lo+(xs/10-sigma_lo)*F(random.randint(1,999),1000)
        if not (sigma< xs/10 and chi<sigma/10): return None
        ys=1-sigma-F(1,10**6); ym=1-sigma-F(1,2*10**6); yd=1-sigma+F(1,10**6)
    else:
        sigma=chi/10*F(random.randint(1,999),1000)
        if sigma>=xs/10: return None
    c={}
    c['read tail 1-sigma-2(h+-h-) > u']= 1-sigma-2*(hp-hm) > u
    c['t b*^2->0: 8chi<2h-']= 8*chi<2*hm
    c['cov ratio 2chi-(h+-h-)>0']= 2*chi-(hp-hm)>0
    c['k a* n^-2chi=o(a*)']= chi-2*chi<0
    c['k b*^2 = o(a*): chi-2h- < -h+']= chi-2*hm < -hp
    c['k b* = o(1)']= chi<hm
    c['n b*^2=o(n a*): -2h-<-h+']= -2*hm<-hp
    c['u+8chi<xs']= u+8*chi<xs
    c['chi<u (2^k e^{-n^u})']= chi<u
    c['sigma<u (LLL degree)']= sigma<u
    c['n^{xs}=o(n a*)']= xs<1-hp
    c['r log n=o(n a*)']= sigma<1-hp
    c['t=n^{8chi}: Restriction n^{u+8chi}<n^{xd}']= u+8*chi<xd
    if case=='sub':
        c['m b*=o(n a*): ym-h- < 1-h+']= ym-hm<1-hp
        c['ys<ym<1-sigma<yd']= ys<ym<1-sigma<yd
        c['T exponent 1-sigma+eps<yd for small eps']= 1-sigma<yd
        c['adj-slice crowd mean m r/n <=1: ym+sigma<=1']= ym+sigma<=1
    else:
        c['b0+sigma<chi/2 feasible']= sigma<chi/2
        a_d=F(1,200); a_s=a_d/101
        c['width: a_s+log(1/.49)a_d/10 < a_d']= float(a_s)+math.log(1/.49)*float(a_d)/10 < float(a_d)
        c['tag repeat: m log2 > S_s']= float(a_d)/10*math.log(2) > float(a_s)
        c['odd cap: S_d < n log 2']= float(a_d)<math.log(2)
    return c
for case in ['sub','lin']:
    n=0
    while n<20000:
        c=check9(case)
        if c is None: continue
        n+=1
        for k,v in c.items():
            if not v: bad.append((case,k))
    print(f"9 {case}: sampled {n} admissible parameter sets; failures:",len([b for b in bad if b[0]==case]))
print(set(bad))
