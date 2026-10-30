# Joint satisfiability of Sec 13 step-1 constants and regime-table inequalities, with sample admissible values.
# Works in log-space (exact rationals for exponents; mpmath 60 digits for logs of thresholds).
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps=60
eta0=F(1,100); u=1000  # u: some fixed even moment order (actual value >10L^2 with L=500R; larger u only raises Q0)
theta=F(1,10**6); a=theta/100
aC=F(1,10**8)                        # < min(eta0,1)/1e6 = 1e-8 ... take just below
aC=F(9,10**9)
aB=aC/100                            # "much smaller"
Cb=100*aC/aB+101                     # > 100 aC/aB + 100
Mlo=Cb+101                           # > Cb+100
cq=F(1,20)/Mlo*F(99,100)             # < 1/(20 Mlo)
Mhi=(Mlo+10/cq)*10                   # >> Mlo + 10/cq
omega=aC/(500*Mhi)*F(99,100)         # 5 omega Mhi < aC/100
print("Cb",Cb,"Mlo",Mlo,"cq",float(cq),"Mhi",float(Mhi),"omega",float(omega))
chk=[]
chk.append(("cq*Mlo < 1/20 (=> 2q^Mlo <=(log n)^.1 eventually)", cq*Mlo<F(1,20)))
chk.append(("cq*Mhi > 5 (=> h>(log n)^5)", cq*Mhi>5, float(cq*Mhi)))
chk.append(("4 omega Mhi < aC/2 (kT << q^{aC/2})", 4*omega*Mhi<aC/2))
chk.append(("8 omega Mhi < 1/2 (kT << sqrt(log n) at small-bin cutoff q<=(log n)^2)", 8*omega*Mhi<F(1,2)))
chk.append(("Cb*aB > aC (bias-test width beats input width)", Cb*aB>aC, float(Cb*aB/aC)))
chk.append(("Mlo - Cb > 0 (10q^Cb < a h/(1e8 u))", Mlo-Cb>0))
for t in chk: print(*t)
# thresholds (log of the minimal scale) implied:
# (i) (2^aB-1) g^aB > log 8  => log g > (1/aB) log( log8/(2^aB-1) )
aBm=mp.mpf(aB.numerator)/aB.denominator
logg1=(1/aBm)*mp.log(mp.log(8)/(mp.power(2,aBm)-1))
# (ii) low direct budget: g^aB + log 8 <= g/(1e6 u) * (ln 2 for log2 in ell) ; solve roughly
f=lambda L: mp.e**(aBm*L)+mp.log(8)+1 - mp.e**L/(mp.mpf(10)**6*u)*mp.log(2)
lo,hi=mp.mpf(1),mp.mpf(10)**3
while f(hi)>0: hi*=2
for _ in range(200):
    mid=(lo+hi)/2
    if f(mid)>0: lo=mid
    else: hi=mid
print("log g needed for (2g)^aB-g^aB>log8:", mp.nstr(logg1,6))
print("log g needed for low-direct allocation budget:", mp.nstr(hi,6))
# (iii) low cluster: q^{Mlo-Cb} > 1e9 u / a
logq=mp.log(mp.mpf(10)**9*u*mp.mpf(a.denominator)/a.numerator)/(Mlo-Cb)
print("log q needed for 10q^Cb < a h/(1e8 u):", mp.nstr(logq,6))
# (iv) internal solver needs T log h = o(a^2 k): a^2 h^{2 omega} >> log h ; log h >= Mlo log q
om=mp.mpf(omega.numerator)/omega.denominator
# find log h with a^2 * exp(2 om logh) >= 100*logh
am=mp.mpf(a.numerator)/a.denominator
g2=lambda Lh: am**2*mp.e**(2*om*Lh)-100*Lh
lo,hi=mp.mpf(1),mp.mpf(10)
while g2(hi)<0: hi*=2
for _ in range(300):
    mid=(lo+hi)/2
    if g2(mid)<0: lo=mid
    else: hi=mid
print("log h needed for a^2 k >= 100 log h:", mp.nstr(hi,6), "-> log log n needed in low mode (h<=(log n)^.1):", mp.nstr(mp.log(hi*10),6))
