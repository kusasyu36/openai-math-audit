"""Corner enumeration for (4.19), and the consequences (4.17)-(4.18). Exact rationals."""
from fractions import Fraction as Fr
from itertools import product
al=Fr(3,40); s=1/(1-al); beta=Fr(7,25); c=Fr(3,2); gam=c/2; eps=Fr(3,20)
Ast=s+1+s/8; eta=s*beta+(1+s/8)*gam
D=lambda x: Fr(x)
BINS={'O':('0','.09','.60',None,'.37','.56','.72','.13'),'P1':('.083','.15','1.01','4.3','.29','.46','.72','.13'),
'R1':('.14','.23','1.19','3','.29','.371','.68','.13'),'U1':('.22','.342','1.19','1.72','.29','.398','.61','.123'),
'V1':('.337','.45','.83','.87','.31','.4','.50','.060'),'T1':('.445','.5','.20','.640','.32','.4','.48','.022')}
SIDES={'a<=-1':(Fr('-.08'),Fr('.70'),Fr('.8'),Fr('1.77'),Fr('.070')),'a>=-1':(Fr('-.475'),Fr(0),Fr('.8'),Fr('1.77'),Fr('.097'))}
PAPER={('O','a<=-1'):('1.867','.52','1.314'),('O','a>=-1'):('2.34','.55','1.155'),('P1','a<=-1'):('1.766','.49','1.813'),('P1','a>=-1'):('2.564','.433','1.155'),
('R1','a<=-1'):('1.681','.664','2.057'),('R1','a>=-1'):('2.621','.420','1.155'),('U1','a<=-1'):('1.55','.847','2.101'),('U1','a>=-1'):('2.461','.615','1.155'),
('V1','a<=-1'):('1.263','.912','1.787'),('V1','a>=-1'):('1.885','.770','1.205'),('T1','a<=-1'):('.913','.869','1.386'),('T1','a>=-1'):('1.071','.840','1.251')}
F_CLAIM={'R1':(Fr('.20'),Fr('.295')),'U1':(Fr('.23'),Fr('.23')),'V1':(Fr('.23'),Fr('.23')),'T1':(Fr('.39'),Fr('.39')),'P1':(Fr('.42')**2,Fr('.61')**2)}
print('Ast',float(Ast),'eta',float(eta))
for nm,(bmin,bmax,Lmax,Dmax,Kmin,Kmax,Jmax,ell1) in BINS.items():
    for side,(zmin,zmax,hmin,hmax,chimax) in SIDES.items():
        mW=mX=mY=0
        # J0 >= 0 (b>=0); L>=0
        for b,L,K,J,z,h in product((Fr(bmin),Fr(bmax)),(Fr(0),Fr(Lmax)),(Fr(Kmin),Fr(Kmax)),(Fr(0),Fr(Jmax)),(zmin,zmax),(hmin,hmax)):
            W=h-2*b-2*L*z
            X=Ast*b-eta/2+L*z+K*W
            Y=-eta+J*W
            mW=max(mW,abs(W)); mX=max(mX,abs(X)); mY=max(mY,abs(Y))
        pw,px,py=(Fr(v) for v in PAPER[(nm,side)])
        ok = mW<=pw and mX<=px and mY<=py
        zz=max(zmin**2,zmax**2)
        base=Ast+zz/al+mX**2/Fr(ell1)+mY**2/Fr('.75')
        if Dmax is None:
            msg=f'F/chi^2 <= {float(base):.4f} + {float(Fr(".295")*mW**2):.4f}(7+t^2)'
        else:
            Fb=chimax**2*(base+Fr(Dmax)*mW**2)
            idx=0 if side=='a<=-1' else 1
            cl=F_CLAIM[nm][idx]
            msg=f'F<= {float(Fb):.4f} (claim {float(cl):.4f}) {"OK" if Fb<=cl else "EXCEEDS"}'
        print(f'{nm:3s} {side}: max|W|={float(mW):.4f} max|X|={float(mX):.4f} max|Y|={float(mY):.4f}  paper=({pw},{px},{py}) {"OK" if ok else "MISMATCH"} | {msg}')
