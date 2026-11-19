# 補助(観測であって上界ではない): ub.cpp の expair / frame::init / block::setup の一部 / makeD /
# testA の平坦和を Python の binary64(+,-,*,/,sqrt は IEEE の正しい丸め)で写し、
# mpmath 50 桁の「理想値」(同じ添字選択、厳密な単位指数)と比べる。
# 写しは ub.cpp:31-62, 119-131, 155-176, 422-423、seed.cpp:21-31 に対応する。
import math, random
from mpmath import mp, mpf, mpc, cos as mcos, sin as msin, sqrt as msqrt, pi as mpi
mp.dps = 50
PI = float.fromhex('0x1.921fb54442d18p+1')
assert PI == math.pi

def expair(x):
    assert abs(x) <= 16
    k = int(_nearbyint(4*x))
    y = (x - .25*k)*(2*PI); yy = -y*y
    a = 1.; b = y; ta = 1.; tb = y
    for j in range(1, 10):
        ta = (ta*yy)/((2*j-1.)*(2*j)); tb = (tb*yy)/((2*j)*(2*j+1.))
        a += ta; b += tb
    r = k % 4
    return [(a, b), (-b, a), (-a, -b), (b, -a)][r]
def _nearbyint(v):
    f = math.floor(v); d = v - f
    if d > .5: return f+1
    if d < .5: return f
    return f if f % 2 == 0 else f+1

class FZ:  # float complex with the same operation order as struct Z
    __slots__ = ('r', 'i')
    def __init__(s, r, i): s.r = r; s.i = i
    def __add__(s, x): return FZ(s.r+x.r, s.i+x.i)
    def __sub__(s, x): return FZ(s.r-x.r, s.i-x.i)
    def __mul__(s, x):
        if isinstance(x, FZ): return FZ(s.r*x.r - s.i*x.i, s.r*x.i + s.i*x.r)
        return FZ(s.r*x, s.i*x)
    def __truediv__(s, y): return FZ(s.r/y, s.i/y)
def conjf(x): return FZ(x.r, -x.i)
def nrf(x): return x.r*x.r + x.i*x.i
def abf(x): return math.sqrt(nrf(x))
def Ifn(x): return FZ(-x.i, x.r)
def exf(x): a, b = expair(x); return FZ(a, b)

# ideal side in mpmath complex
def exm(x): t = 2*mpi*mpf(x); return mpc(mcos(t), msin(t))

def frame_f(raw):
    S = FZ(0., 0.)
    row = list(raw)
    for j in range(4): S = S + row[j]
    m = math.sqrt(4 - abf(S))
    if not (m > 1.25): return None
    basis = [[None]*4 for _ in range(4)]
    one = FZ(1., 0.)
    for j in range(4):
        basis[0][j] = one*.5
        basis[1][j] = (row[j] - S*.25)/math.sqrt(4. - nrf(S)/4)
    Q = [[None, None] for _ in range(4)]; sel = []
    for k in (2, 3):
        best = -1; q = None; bj = None
        for j in range(4):
            s = 1.
            for l in range(k): s -= nrf(basis[l][j])
            if s < best: continue
            best = s; bj = j
            q = []
            for i in range(4):
                v = one if i == j else FZ(0., 0.)
                for l in range(k): v = v - basis[l][i]*conjf(basis[l][j])
                q.append(v)
        r = 0.
        for i in range(4): r += nrf(q[i])
        if not (r > .45*.45): return None
        sel.append(bj)
        for i in range(4):
            basis[k][i] = q[i]/math.sqrt(r); Q[i][k-2] = basis[k][i]
    return dict(S=S, m=m, Q=Q, sel=sel, r_last=r, row=row)

def frame_m(raw, sel):
    S = sum(raw, mpc(0))
    row = raw
    basis = [[None]*4 for _ in range(4)]
    for j in range(4):
        basis[0][j] = mpc(.5)
        basis[1][j] = (row[j] - S/4)/msqrt(4 - abs(S)**2/4)
    Q = [[None, None] for _ in range(4)]
    for idx, k in enumerate((2, 3)):
        j = sel[idx]
        q = []
        for i in range(4):
            v = mpc(1) if i == j else mpc(0)
            for l in range(k): v -= basis[l][i]*basis[l][j].conjugate()
            q.append(v)
        nrm = msqrt(sum(abs(x)**2 for x in q))
        for i in range(4):
            basis[k][i] = q[i]/nrm; Q[i][k-2] = basis[k][i]
    return dict(S=S, Q=Q)

def Gf(fr):  # FB.G / FC.G
    S = fr['S']; row = fr['row']; one = FZ(1., 0.)
    inv = 1./(16 - nrf(S))
    return [((one*4 - row[j]*conjf(S))*inv, (row[j]*4 - S)*inv) for j in range(4)]
def Gm(frm, raw):
    S = frm['S']; inv = 1/(16 - abs(S)**2)
    return [((4 - raw[j]*S.conjugate())*inv, (raw[j]*4 - S)*inv) for j in range(4)]

SQ6 = math.sqrt(6.)
def trial(phases, bins, N):
    # block::setup (ub.cpp:119-131) の a, B, C, frames, base, W
    a = exf(phases[0]); am = exm(phases[0])
    B = [exf(phases[1+i]) for i in range(4)]; C = [exf(phases[5+i]) for i in range(4)]
    Bm = [exm(phases[1+i]) for i in range(4)]; Cm = [exm(phases[5+i]) for i in range(4)]
    FB = frame_f([conjf(b) for b in B]); FC = frame_f(C)
    if FB is None or FC is None: return None
    if not (min(FB['m'], FC['m']) >= 1.): return None
    FBm = frame_m([b.conjugate() for b in Bm], FB['sel']); FCm = frame_m(Cm, FC['sel'])
    GB = Gf(FB); GBm = Gm(FBm, [b.conjugate() for b in Bm])
    one = FZ(1., 0.)
    # uc::init (ub.cpp:155-163) with exact table entries
    b0, b1, b2, b3 = bins
    def whole(i): return exf(i/N)
    def half(i): return exf((i+.5)/N)
    t = half(b1); dP = whole((b0+b2+1) % N); dM = whole((b0-b2+N) % N)
    oP = whole((b0+b3+1) % N); oM = whole((b0-b3+N) % N); oM = FZ(-oM.r, -oM.i)
    Uf = [dP*t.r, oP*t.i, oM*t.i, dM*t.r]
    tm = exm((b1+.5)/N)
    Um = [exm(mpf(b0+b2+1)/N)*tm.real, exm(mpf(b0+b3+1)/N)*tm.imag,
          -exm(mpf(b0-b3)/N)*tm.imag, exm(mpf(b0-b2)/N)*tm.real]
    maxD = 0.; maxW = 0.; flat_f = 0.; flat_m = mpf(0)
    for i in range(4):
        for j in range(4):
            base = (one + C[i])*conjf(GB[j][0])*(-1.) - (one + C[i]*conjf(a))*conjf(GB[j][1])
            # 注: C++ は -(one+C[i])*conj(...) で、単項マイナスは厳密なので同値
            basem = -(1 + Cm[i])*GBm[j][0].conjugate() - (1 + Cm[i]*am.conjugate())*GBm[j][1].conjugate()
            z = FZ(0., 0.); zm = mpc(0)
            for l in range(4):
                W = FC['Q'][i][l//2]*conjf(FB['Q'][j][l % 2])*SQ6
                Wm = FCm['Q'][i][l//2]*FBm['Q'][j][l % 2].conjugate()*msqrt(6)
                maxW = max(maxW, float(abs(mpc(W.r, W.i) - Wm)))
                z = z + W*Uf[l]; zm += Wm*Um[l]
            D = z + base; Dm = zm + basem
            maxD = max(maxD, float(abs(mpc(D.r, D.i) - Dm)))
            mag = nrf(D); x = math.sqrt(mag) - 1; flat_f += x*x
            flat_m += (abs(Dm) - 1)**2
    return maxW, maxD, abs(flat_f - float(flat_m)), float(flat_m), FB['m'], FC['m'], FB['r_last'], FC['r_last']

random.seed(266)
res = []
tries = 0
while len(res) < 1500 and tries < 200000:
    tries += 1
    spread = random.choice([.05, .1, .15, .2, .3, .5])
    ph = [random.random()] + [random.gauss(0, spread) % 1.0 for _ in range(8)]
    N = random.choice([64, 128, 256])
    bins = [random.randrange(N//2), random.randrange(N//4), random.randrange(N), random.randrange(N)]
    out = trial(ph, bins, N)
    if out is None: continue
    res.append(out)
print("samples", len(res), "tries", tries)
print("max observed W error   %.3e  (my bound 1.29e-10, paper 8e-10)" % max(r[0] for r in res))
print("max observed D error   %.3e  (my bound 5.18e-10, paper 4e-9)" % max(r[1] for r in res))
print("max observed flat-sum error %.3e (my bound 5.95e-8, paper 1e-6)" % max(r[2] for r in res))
print("min frame m observed   %.5f  (threshold 1.25)" % min(min(r[4], r[5]) for r in res))
print("min last complement r  %.5f  (threshold .2025)" % min(min(r[6], r[7]) for r in res))
near = [r for r in res if min(r[4], r[5]) < 1.3]
if near: print("near m<1.3: n=%d max D err %.3e" % (len(near), max(r[1] for r in near)))
# expair 単体
worst = 0.
for _ in range(20000):
    x = random.uniform(-16, 16)
    a, b = expair(x); e = exm(x)
    worst = max(worst, float(abs(mpc(a, b) - e)))
print("expair max observed error %.3e (lemma 2e-14 incl. none of arg formation; my internal bound 5.7e-15)" % worst)
print("NOTE: all numbers above are observations, not upper bounds.")
