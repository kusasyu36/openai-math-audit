# 266 付録 A の誤差予算を、自分で組んだ伝播式で再評価する(厳密部分)。
# mpmath.iv は実数関数の外向き区間包含であり、binary64 の丸めそのものを模擬するものではない。
# ここで区間に入れるのは「自分の導出した閉じた上界式」と「スケールの主張」で、
# 上界式の各項は REPORT.md 2 節の導出に対応する。結果は、上端が論文の値を下回るかで判定する。
from mpmath import iv, mpf
iv.dps = 40

U = iv.mpf('2e-16')          # 論文の u(rounding.tex:61)。真の 2^-53 より大きい(保守的)
U53 = iv.mpf(2)**-53
def gam(m, u=U):             # 自分の gamma: (1+u)^(m+1)-1 の区間上界
    return (1+u)**(m+1)-1
def gpaper(m):
    return iv.mpf('2.1e-16')*(m+1)

rows = []
def chk(name, mine, paper, note=''):
    hi = mine.b if hasattr(mine, 'b') else mine
    ok = hi < mpf(paper)
    rows.append((name, hi, paper, ok, note))
    print(f"{'OK ' if ok else 'NG '} {name:<46} mine<{float(hi):.3e}  paper {paper}  {note}")
    return mine

print("== 0. 前提 ==")
for m in (2, 9, 10, 12, 16, 20, 25, 36, 41, 2000):
    chk(f"(1+u)^(m+1)-1 < gamma_m (m={m})", gam(m), str(float(gpaper(m).b)))
chk("2^-53 <= u", U53, '2.0000001e-16')

print("== 1. expair (Lemma A.2) ==")
pihat = iv.mpf(884279719003555)/iv.mpf(281474976710656)
pi = iv.pi
chk("|pihat - pi|", abs(pihat-pi), '1.225418e-16')
ymax = 2*pihat*iv.mpf('0.125')*(1+U)
chk("|y| bound", ymax, '0.786')
yerr = 2*abs(pihat-pi)*iv.mpf('0.125') + U*iv.mpf('0.786')
chk("angle error |y-2pi(x-k/4)|", yerr, '3e-16')
y = iv.mpf('0.786')
chk("cos tail y^20/20!", y**20/iv.factorial(20), '3.330e-21')
chk("sin tail y^21/21!", y**21/iv.factorial(21), '1.247e-22')
yv = iv.mpf('0.786')
chk("even nonconstant abs sum cosh(.786)-1", ((iv.exp(yv)+iv.exp(-yv))/2)-1, '0.4')
chk("even abs sum cosh(.786)", ((iv.exp(yv)+iv.exp(-yv))/2), '1.4')
chk("odd after-y abs sum sinh(.786)-.786", ((iv.exp(yv)-iv.exp(-yv))/2)-yv, '0.2')
chk("odd abs sum sinh(.786)", ((iv.exp(yv)-iv.exp(-yv))/2), '1')
rho = (1+U)**27-1; alpha = (1+U)**9-1
chk("rho=(1+u)^27-1", rho, '5.5e-15'); chk("alpha=(1+u)^9-1", alpha, '1.81e-15')
ec = iv.mpf('0.4')*rho + iv.mpf('1.4')*(1+rho)*alpha + iv.mpf('4e-20')
es = iv.mpf('0.2')*rho + (1+rho)*alpha + iv.mpf('4e-20')
chk("cos coordinate error", ec, '4.74e-15'); chk("sin coordinate error", es, '2.92e-15')
# 角度誤差は |e'|<=1 で効く。単位円上の Lipschitz 定数 1(角度は radian)
eexp = iv.sqrt(ec**2+es**2) + yerr + iv.mpf('1e-295')
chk("expair complex error (|x|<=16, exact arg)", eexp, '6e-15')
# 引数形成: (b+.5)/N 等で 5e-16 cycles、位相差の減算で計 2e-15 cycles
eexp_arg = eexp + 2*pi*iv.mpf('2e-15')
chk("expair + argument formation", eexp_arg, '6e-14', 'epsilon_exp')
EPS = iv.mpf('6e-14')

print("== 2. frame -> W -> K, D (rounding.tex:370-452) ==")
# stored m=sqrt(4-|S|) > 1.25  =>  stored |S| < 4-1.5625 = 2.4375; ideal |S| <= that + 2.6e-13
S_err = 4*EPS + 3*U*iv.mpf(4)          # 4 個の指数誤差 + 3 回の加算(部分和 <=4)
chk("S error", S_err, '2.5e-13')
Smax = iv.mpf('2.4375') + iv.mpf('2.6e-13')
S2_err = 2*Smax*S_err + S_err**2 + 3*U*Smax**2
chk("|S|^2 error", S2_err, '2.2e-12')
den = 16 - Smax**2
chk("16-|S|^2 lower (paper uses 6.1)", -den, '-6.1')
inv_err = (S2_err + U*16)/(den*(den - S2_err)) + U/den
chk("reciprocal 1/(16-|S|^2) error", inv_err, '7e-14')
num_err = EPS*Smax + S_err + iv.mpf('1e-15')*Smax + U*8
chk("pseudoinverse numerator error", num_err, '5.1e-13')
Gerr = num_err/den + 8*inv_err + U*iv.mpf('0.81')*2
chk("frame G entry error", Gerr, '1e-12')
Gmax = 1/(4-Smax)
chk("|G ideal| (max (4+s)/(16-s^2)=1/(4-s))", Gmax, '0.81')
# basis[1]
b1num = EPS + iv.mpf('0.25')*S_err + U*iv.mpf(2)
chk("basis[1] numerator error", b1num, '1.3e-13')
b1den = iv.sqrt(4 - Smax**2/4)
chk("basis[1] denominator lower (>1.5)", -b1den, '-1.5')
b1den_err = (S2_err/4 + 2*U*4)/(2*(b1den - iv.mpf('1e-12'))) + U*2
chk("basis[1] denominator error", b1den_err, '2e-13')
b1err = (b1num + 1*b1den_err)/(b1den - b1den_err) + U
chk("basis[1] entry error", b1err, '4e-13')
# 第 1 の補空間(k=2): q_i = e_j - sum_{l<2} b_l[i] conj(b_l[j])
cmul = iv.mpf('1e-15')
q2 = 2*b1err + cmul + 2*U*2          # 積の摂動 + 複素乗算丸め + 2 回の減算
chk("k=2 projection entry error", q2, '9e-13')
n2 = 2*q2 + 2.2e-16*10*2              # 4 成分の摂動ノルム + ノルム丸め(||q||<=2)
chk("k=2 norm error", n2, '2e-12')
low = iv.mpf('0.45')*(1-U) - n2      # stored sqrt(r) > .45(r>.45^2 から)
chk("k=2 ideal norm lower (>.449999)", -low, '-0.449999')
c2 = (q2 + n2)/low + 2*U
chk("k=2 normalized entry error", c2, '7e-12')
q3 = q2 + 2*c2 + cmul + U*2
chk("k=3 projection entry error", q3, '1.6e-11')
n3 = 2*q3 + 2.2e-16*10*2
chk("k=3 norm error", n3, '3.5e-11')
c3 = (q3 + n3)/(iv.mpf('0.45')*(1-U) - n3) + 2*U
chk("k=3 normalized entry error", c3, '1.2e-10')
Perr = 2*c2 + 2*c3 + 2*cmul + U
chk("kernel projector entry error", Perr, '3e-10')
s6 = iv.sqrt(6)
Werr = (c2 + c3)*s6*(1+U)**3 + cmul*s6 + iv.mpf('2e-15')
chk("W entry error (worst: c3 on both factors)", (2*c3)*s6 + cmul*s6 + iv.mpf('2e-15'), '8e-10', 'paper 8e-10 covers c3+c3')
Werr = (2*c3)*s6 + cmul*s6 + iv.mpf('2e-15')
basef = 2*EPS + cmul + U*2           # 1 + C conj(a)
chk("base factor error", basef, '1.3e-13')
baseerr = 2*(2*Gerr + basef*Gmax + cmul*2*Gmax) + U*4
chk("base entry error", baseerr, '1e-11')
Uerr = 2*EPS + U
chk("uc.U entry error", Uerr, '2e-13')
Kerr = 4*(Werr*1 + s6*Uerr + cmul*s6) + 3*U*iv.mpf('4.6')*2
chk("K entry error", Kerr, '4e-9')
Derr = Kerr + baseerr + U*iv.mpf(5)
chk("D entry error", Derr, '4e-9')
DERR = Derr
magerr = 2*iv.mpf('4.6')*Derr + Derr**2 + 2*U*iv.mpf('4.6')**2
chk("mag=|D|^2 error", magerr, '4e-8')
sqrterr = Derr + 2.2e-16*4*iv.mpf('4.6')
chk("sqrt(mag) vs |D| error", sqrterr, '5e-9')

Dop = iv.sqrt(4+Smax)*2/(iv.mpf('1.25')-iv.mpf('1e-12'))
chk("||D|| <= max(sqrt6, top*A0/mB) (completion.tex:388, imported)", Dop, '4.6')
chk("seed derivative modulus 1.62+.81*4.6", iv.mpf('1.62')+iv.mpf('0.81')*iv.mpf('4.6'), '6')
print("== 3. A flat norm (ub/patchx:422-424) ==")
dev = iv.mpf('3.6')
x_err = sqrterr + U*dev
flat = 16*(2*dev*x_err + x_err**2 + U*dev**2) + gam(16)*16*dev**2
chk("flat sum error (paper scale |dev|<=3.6)", flat, '1e-6', 'margin 1e-5')
flat_tight = 2*6*x_err + 16*x_err**2 + gam(16)*iv.mpf('2.26')
chk("flat sum error if ideal sum<=z^2<2.26", flat_tight, '1e-6', 'exclusion direction only')
FLAT = flat

print("== 4. A directional (ub.cpp:236-250) ==")
# jac: fixed[1+k] は projector 誤差 x 1.62; jac[1+k] は + D*bm; jac[5+k] は + K*cm; jac[9+k] は sum W dU
fixed_err = Perr*iv.mpf('1.62') + iv.mpf('3e-12')
chk("fixed entry error", fixed_err, '2e-9')
jac_err = fixed_err + Derr*iv.mpf('0.81') + iv.mpf('4.6')*iv.mpf('2e-12') + cmul*5
jac_err = max(jac_err.b, (4*(Werr + s6*iv.mpf('2e-13'))).b)
jac_err = iv.mpf(jac_err)
chk("jac entry error", jac_err, '1e-8')
Jerr = Derr*6 + iv.mpf('4.6')*jac_err + 2*U*28
chk("flat real J error", Jerr, '1e-7')
cflat = magerr/2 + U*11
chk("flat c error", cflat, '5e-8')
const_pair = 16*cflat + 4*iv.mpf('4e-13') + gam(20)*(16*iv.mpf('10.1')+4)
chk("constant pairing error", const_pair, '9e-7')
dpair = 16*Jerr + 4*EPS + gam(20)*(16*28+4)
chk("derivative pairing error", dpair, '1.7e-6')
lin = (iv.mpf('0.203')*3 + iv.mpf('0.197')*4)*dpair
chk("linear support undererror", lin, '2.38e-6')
sd = 4*Derr + 2.2e-16*18*iv.mpf('4.6')*4
chk("sqrt(sd) error", sd, '2e-8')
Adir = const_pair + lin + sd*1 + iv.mpf('1e-9')
chk("A directional total", Adir, '4e-6', 'margin 1e-4')
# 論文自身の小項目を足した値
chk("A directional total (paper sub-bounds)", iv.mpf('9e-7')+iv.mpf('2.38e-6')+iv.mpf('2e-8')+iv.mpf('1e-9'), '4e-6')

print("== 5. Candidate (rounding.tex:711-821) ==")
mlow = iv.mpf('0.15')*(1-U) - sqrterr
chk("ideal m_j lower (>.14999998)", -mlow, '-0.14999998')
uderr = 2*Derr/mlow + 4*U
chk("unit-direction error", uderr, '5.5e-8')
chord = uderr + iv.mpf('1e-7') + EPS + 4*U
chk("true chord error", chord, '1.56e-7')
Xerr = uderr*6 + jac_err + 4*U*6
chk("X entry error", Xerr, '1e-6')
chk("mi lower (>.148)", -(mlow - Xerr - iv.mpf('1e-6')*(1+U)), '-0.148')
chk("invmag upper (<6.8)", 1/(iv.mpf('0.15')-iv.mpf('1e-6')-iv.mpf('1e-15'))+iv.mpf('1e-9'), '6.8')
# スケール: q<3.8 を通った枝
q = iv.mpf('3.8')
ratio = 1/iv.sqrt(1-q/4) + iv.mpf('1e-7')
arc = iv.sqrt(q)*ratio + iv.mpf('1e-8')
chk("arc (q<3.8)", arc, '8.8')
chk("slack (<4.1) from z<1.5", iv.mpf('1.5')**2 + iv.mpf('1e-5') + iv.mpf('1e-10'), '4.1')
weighted = iv.sqrt(iv.mpf('1.5')**2 + iv.mpf('1e-5') + iv.mpf('2e-10'))*ratio + iv.mpf('1e-8')
chk("weighted (slack<z^2+1e-5, ratio at q<3.8)", weighted, '10')
h = pi/64 + iv.mpf('1e-10')
pnorm = iv.sqrt(iv.mpf('0.203')**2 + 4*h**2)
chk("parameter vector norm (N>=64)", pnorm, '0.226')
ABerr = iv.sqrt(208)*Xerr*iv.mpf('1.00001') + iv.sqrt(130)*gam(36)*iv.mpf('25.0003')
chk("||AB - T*X|| (Frobenius)", ABerr, '1.5e-5')
ABpaper = iv.mpf('1.5e-5')*iv.mpf('0.226')
chk("AB substitution (paper 1.5e-5 x .226)", ABpaper, '3.4e-6')
hor = iv.mpf('1e-5')*iv.sqrt(iv.mpf('0.203')**2 + iv.mpf('8.8')**2)
chk("hor*sqrt(Rs^2+arc^2)", hor, '8.81e-5', 'mathematical allowance')
cc = 2*pi*iv.sqrt(1+iv.mpf('1e-5'))*iv.mpf('2e-6')
chk("center correction 2pi sqrt(1+ort) 2e-6", cc, '1.27e-5', 'mathematical allowance')
offerr = 2*pi*iv.mpf('1.00001')*iv.mpf('1e-11') + gam(36)*iv.mpf('0.15')*7
chk("off error", offerr, '1e-9')
ext = hor + cc + ABpaper + iv.mpf('1e-9') + iv.mpf('1e-11') + iv.mpf('1e-10')
chk("external tangent total (paper sub-bounds)", ext, '0.00012', 'margin .00015')
ext_mine = hor + cc + ABerr*pnorm + offerr + iv.mpf('1e-11') + iv.mpf('1e-10')
chk("external tangent total (mine)", ext_mine, '0.00012')
extra = 2*2*pi*iv.mpf('2e-5')
print("   margin/paperX =", float(mpf('0.00015')/mpf('0.00012')),
      " margin/ext(paper sub) =", float((iv.mpf('0.00015')/ext).a),
      " margin/ext(mine) =", float((iv.mpf('0.00015')/ext_mine).a),
      " (margin+2x.00002cyc)/ext =", float(((iv.mpf('0.00015')+extra)/ext).a))
print("   rounding-only part of ext (mine):", float((ABerr*pnorm + offerr + iv.mpf('1.1e-10')).b))

print("== 6. covariance G entries (rounding.tex:823-846) ==")
z1 = gam(9)*9*iv.mpf('25.0003')**2
chk("first G: 9-term AB dot", z1, '1.3e-11')
Rs2 = iv.mpf('1.45')*iv.mpf('0.203')**2
chk("1.45 Rs^2", Rs2, '0.06')
G1 = Rs2*z1 + 4*U*Rs2*9*iv.mpf('25.0003')**2 + (1+1/iv.mpf('0.45'))*(gam(16)*64 + 4*U*64) + 2*U*10
chk("first G entry error", G1, '1e-11')
dw = gam(17)*110                      # 2 積/項 -> 1 余分の因子
chk("Dweighted dot error (abs sum 110)", dw, '5e-13')
coef2 = (1+1/iv.mpf('0.8'))*weighted**2
chk("second multiplier 2.25*weighted^2", coef2, '225')
G2 = coef2*dw + 6*U*coef2*110 + iv.mpf('1.8')*iv.mpf('0.203')**2*(gam(9)*iv.mpf('1.00002')+4*U) + 2*U*26000
chk("second G entry error", G2, '1e-9')
chk("k*second G (k<=10)", 10*G2, '1e-7', 'diag margin 1e-7')
chk("k*first G (k<=10)", 10*G1, '1e-7')
G2t = iv.mpf('2.25')*100*(gam(17)*iv.mpf('6.81')) + 6*U*iv.mpf('2.25')*100*iv.mpf('6.81') + 2*U*10
print("   (tighter, Cauchy-Schwarz sum|T_ih T_jh| <= 1.00001: second G entry <", float(G2t.b), ")")

print("== 7. spectralUpper / ell (rounding.tex:851-969) ==")
chk("row sum error", gam(10)*100, '3e-13', 'margin 1e-10')
chk("G^2 dot rel. error /m^2", gam(10)*10, '3e-14')
fro = 10*iv.mpf('3e-14') + 2.2e-16*102
chk("fl(sqrt fro4) rel error vs A", fro, '4e-13')
chk("fourth-root abs error at sqrt(A)<=100", fro*100 + U*100, '1e-10', 'margin 2e-7')
gamma_max = iv.mpf('100.1') + iv.sqrt(iv.mpf(1000)*iv.mpf('100.1'))
chk("gamma upper", gamma_max, '418')
res_c = gam(12)*528
chk("residual component /s_z", res_c, '2e-12')
res_n = iv.sqrt(10)*iv.mpf('2e-12') + 2.2e-16*12*iv.sqrt(10)*528
chk("residual norm /s_z", res_n, '2e-11', 'code adds 1e-10*s_z')
chk("U.z dot /s_z", gam(10)*100, '3e-13', 'code adds 1e-10*s_z')
chk("U norm", 2.2e-16*12*iv.sqrt(1000), '1e-13')
dsur = iv.mpf('1e-10') - gam(10)*100 - U*101 - U*iv.mpf('1e-10')
chk("dot surplus >= 9e-11 s_z", -dsur, '-9e-11')
f1 = (iv.mpf('1e-10') - iv.mpf('1e-13') - U*40)/iv.sqrt(1000)
chk("first factor rel surplus >= 2e-12", -f1, '-2e-12')
Bmax = iv.sqrt(1000)*(iv.sqrt(10)*528*iv.mpf('1e8')*1.001)/iv.mpf('0.999e-6')
chk("quotient B scale", Bmax, '6e18')
chk("final add err vs surplus (s_z=1,B=0 worst)", U*101 + iv.mpf('1e-300'), '9e-11')

print("== 8. Root identity (patchx.cpp:167-182) ==")
qerr = 2*EPS + cmul
chk("row-edge product error", qerr, '1.3e-13')
LPerr = 12*(qerr*iv.mpf(25)/36 + 1*iv.mpf('1e-15')) + gam(36)*iv.mpf('8.34')
chk("LP entry error", LPerr, '5e-12')
ent = 30*2*LPerr + iv.mpf('1e-15') + gam(41)*560
chk("residual entry error", ent, '1e-9')
chk("absolute sum (1+10*4+30*2*8.4)", iv.mpf(1)+40+30*2*iv.mpf('8.4'), '560')
fr = 36*ent
chk("Frobenius of entry errors", fr, '3.6e-8')
acc = 2.2e-16*1298*iv.mpf(20161)
chk("norm accumulation at norm 20161", acc, '1e-8')
chk("Root identity total", fr + acc, '5e-8', 'margin 2e-6')
chk("Root identity total (paper entry 1e-9)", 36*iv.mpf('1e-9') + acc, '5e-8')

print("== 8b. Root stored bounds used by rows 1 and 6 (rounding.tex:550-569, 651-658) ==")
chk("T*T-I dot (len 36, |T|<=2, initial addend)", gam(37)*145, '1.2e-12')
ortall = 10*iv.mpf('1.2e-12') + iv.mpf('2.24e-14')*iv.mpf('1e-5')
chk("ort aggregate+norm error (accepted path ort<1e-5)", ortall, '1e-10', 'margin 1e-8')
pc = gam(36)*72/36 + (gam(6)*12*2 + U*24)/6 + 3*U*4
chk("proj component error", pc, '1e-13')
chk("hor Frobenius error (360 comps)", iv.sqrt(360)*iv.mpf('1e-13') + iv.mpf('2.24e-14')*iv.mpf('1e-5'), '1e-11', 'margin 1e-8')
chk("HV/LV Gram dot (len<=36, |.|<=2)", gam(36)*144, '1.2e-12')
chk("Gram row sum (25 terms, scale 3600)", 25*iv.mpf('1.2e-12') + gam(25)*3600, '1e-10', 'inside +1e-8')
chk("M - HV LV^* entry (25-term, abs sum 100)", gam(25)*100, '6e-13')
chk("M - HV LV^* Frobenius (36x30)", iv.sqrt(1080)*iv.mpf('6e-13'), '2.1e-11')
chk("mnorm product/root error at scale 3601", 4*U*3601 + U*3601, '5e-12', 'final +1e-8 covers this + 2.1e-11')
print("== 9. C arrays and joint test (rounding.tex:1226-1262, 1385-1455) ==")
inv6 = 1/s6
qC = (iv.mpf('6e-15') + 2*pi*U*6)*inv6*(1+U) + U*inv6*3
chk("q error (expair internal 6e-15 + arg 6 cycles)", qC, '2.5e-14')
chk("q error (with lemma's 2e-14)", (iv.mpf('2e-14') + 2*pi*U*6)*inv6 + U*3, '1e-13')
SC = 6*iv.mpf('1e-13') + 5*U*s6
chk("S error", SC, '7e-13')
S2C = 2*s6*SC + SC**2 + 2*U*6
chk("|S|^2 error", S2C, '4e-12')
chk("c error", S2C/(4*pi) + 3*U, '1e-12')
JC = SC*inv6 + s6*iv.mpf('1e-13') + 2*U*1
chk("J, Dre error", JC, '1e-12')
BC = 6*iv.mpf('1e-12')*iv.mpf('1.00001') + gam(6)*6*iv.mpf('1.00001')
chk("B error", BC, '7e-12')
Qe = 2*EPS + cmul
chk("pair product Q error", Qe, '1.3e-13')
Se = 6*Qe + 5*U*6
chk("pair S error", Se, '8e-13')
chk("pair |S|^2 error", 2*6*Se + Se**2 + 2*U*36, '1.1e-11')
chk("U derivative error", (Se + 6*Qe)/6 + 4*U, '1e-12')
cpair = 14*iv.mpf('1e-12') + gam(14)*iv.mpf('6.7')
chk("joint constant pairing", cpair, '1.5e-11')
zB = 12*iv.mpf('7e-12') + gam(12)*12*iv.mpf('2.5')
chk("shared-Hadamard coefficient", zB, '9e-11')
zx = 8*iv.mpf('1e-12') + gam(8)*8
chk("vector-phase coefficient", zx, '9e-12')
jt = cpair + iv.mpf('0.030001')*iv.sqrt(10)*zB + 2*(iv.mpf('0.1')*s6*zx) \
   + 2*(pi*iv.mpf('0.01')*iv.mpf('1.1e-11') + 2*pi*iv.mpf('0.1')*iv.mpf('0.03')*iv.mpf('1e-11')/s6
        + pi*iv.mpf('0.0009')*iv.mpf('1e-11')/s6 + iv.mpf('0.019')*iv.mpf('1e-11')) \
   + 40*U*10
chk("joint total (mine)", jt, '1e-8', 'margin 3e-7')

print()
bad = [r for r in rows if not r[3]]
print("checks:", len(rows), " failed:", len(bad))
for r in bad: print("  FAILED:", r[0], float(r[1]), r[2])
