"""
109 監査: 07-resampling.tex と 06-transforms.tex の変換恒等式を数値（mpmath 高精度）で確かめる。
 (a) Gaussian 恒等式 T P_s F_s = P_t F_t S（HvdH Thm 4.2 の本論文版の符号・置換）
 (b) N = C T D, ||N - I|| < 2.01 e^{-pi alpha^2 theta/2}（HvdH Lemma 4.6）
 (c) eq:no-sort-one-dimensional: P_s F_s = 2^{2 alpha^2} B_0 P_t F_t A
 (d) eq:permuted-chirp: Q F_t u = conj(a) . (a * (conj(a) . u)) / T
 (e) eq:source-permutation-cancel: (R F_s^+)[(R F_s u).(R F_s v)] = u*v/S^2
 (f) 06: R_r = C[y]/(y^r+1) 上の合成変換 F_t^+ F_t^- = I/t、DIF 式 (eq:synthetic-dif)、ビット反転順、正規化畳み込み
 (g) 08: CRT 写像 Φ が畳み込みを保つこと、(eq:triangular-crt) の座標
実行: uv run --with sympy --with mpmath --with numpy python check_resampling.py
"""
import mpmath as mp
import itertools, random
from math import gcd
mp.mp.dps = 40

ok_all = True
def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK' if cond else 'NG'}] {name} {detail}")

def F(q, u, sign=-1):
    return [sum(u[j] * mp.expj(sign * 2 * mp.pi * j * k / q) for j in range(q)) / q for k in range(q)]
def P_s(s, t, u):   # (P_s u)_j = u_{tj}
    return [u[(t * j) % s] for j in range(s)]
def P_t(s, t, v):   # (P_t v)_k = v_{-sk}
    return [v[(-s * k) % t] for k in range(t)]
def Smap(s, t, alpha, u, J=60):
    return [sum(mp.exp(-mp.pi * (j - mp.mpf(s) * k / t)**2 / alpha**2) * u[j % s] for j in range(-J, J + s)) / alpha for k in range(t)]
def Tmap(s, t, alpha, u, J=60):
    return [sum(mp.exp(-mp.pi * alpha**2 * (mp.mpf(t) * j / s - k)**2) * u[j % s] for j in range(-J, J + s)) for k in range(t)]
def norm(v): return max(abs(x) for x in v)
def sub(u, v): return [a - b for a, b in zip(u, v)]

rng = random.Random(3)
for (s, t, alpha) in ((7, 11, 2), (13, 16, 2), (11, 16, 3)):
    assert gcd(s, t) == 1
    u = [mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(s)]
    lhs = Tmap(s, t, alpha, P_s(s, t, F(s, u)))
    rhs = P_t(s, t, F(t, Smap(s, t, alpha, u)))
    err = norm(sub(lhs, rhs))
    check(f"(a) T P_s F_s = P_t F_t S  (s,t,alpha)=({s},{t},{alpha})", err < mp.mpf(10)**-25, mp.nstr(err, 5))
    # (b) N = C T D
    theta = mp.mpf(t) / s - 1
    q = [int(mp.floor(mp.mpf(t) * j / s + mp.mpf(1) / 2)) for j in range(s)]
    beta = [mp.mpf(t) * j / s - q[j] for j in range(s)]
    def Nmap(u):
        Du = [mp.exp(mp.pi * alpha**2 * beta[j]**2) * u[j] for j in range(s)]
        Tu = Tmap(s, t, alpha, Du)
        return [Tu[q[j]] for j in range(s)]
    # operator norm (max row sum) of N - I: compute columns
    cols = []
    for j in range(s):
        e = [mp.mpf(0)] * s; e[j] = mp.mpf(1)
        cols.append(Nmap(e))
    NmI = max(sum(abs(cols[j][i] - (1 if i == j else 0)) for j in range(s)) for i in range(s))
    bound = mp.mpf('2.01') * mp.exp(-mp.pi * alpha**2 * theta / 2)
    check(f"(b) ||N-I|| < 2.01 e^(-pi alpha^2 theta/2)", NmI < bound, f"{mp.nstr(NmI,5)} < {mp.nstr(bound,5)}")
    check(f"(b) q_0=0, q_{{s-1}}<t, q increasing", q[0] == 0 and q[-1] < t and all(q[j + 1] > q[j] for j in range(s - 1)))
    # (c) with A = S/2, B_0 = D' J' C, D' = 2^{-(2alpha^2-2)} D, J' = N^{-1}/2
    if NmI < 0.5:
        Nmat = mp.matrix(s, s)
        for j in range(s):
            for i in range(s):
                Nmat[i, j] = cols[j][i]
        Ninv = mp.inverse(Nmat)
        def B0(v):  # v length t
            Cv = [v[q[j]] for j in range(s)]
            Jv = [sum(Ninv[i, j] * Cv[j] for j in range(s)) / 2 for i in range(s)]
            return [mp.mpf(2)**(-(2 * alpha**2 - 2)) * mp.exp(mp.pi * alpha**2 * beta[i]**2) * Jv[i] for i in range(s)]
        A = [x / 2 for x in Smap(s, t, alpha, u)]
        rhs2 = [mp.mpf(2)**(2 * alpha**2) * x for x in B0(P_t(s, t, F(t, A)))]
        lhs2 = P_s(s, t, F(s, u))
        err2 = norm(sub(lhs2, rhs2))
        check(f"(c) P_s F_s = 2^(2alpha^2) B_0 P_t F_t A", err2 < mp.mpf(10)**-20, mp.nstr(err2, 5))
        # contraction norms (||S'||<3/4, ||J'||<7/8, ||D'||<1) as claimed (checked numerically here)
        Scols = [Smap(s, t, alpha, [mp.mpf(1) if i == j else mp.mpf(0) for i in range(s)]) for j in range(s)]
        Snorm = max(sum(abs(Scols[j][k]) for j in range(s)) for k in range(t)) / 2
        Jnorm = max(sum(abs(Ninv[i, j]) for j in range(s)) for i in range(s)) / 2
        Dnorm = max(mp.mpf(2)**(-(2 * alpha**2 - 2)) * mp.exp(mp.pi * alpha**2 * beta[i]**2) for i in range(s))
        check(f"(c) ||S'||<3/4, ||J'||<7/8, ||D'||<1", Snorm < 0.75 and Jnorm < 0.875 and Dnorm < 1, f"{mp.nstr(Snorm,4)},{mp.nstr(Jnorm,4)},{mp.nstr(Dnorm,4)}")

# (d) permuted chirp, d=2 axes
def tensorF(ts, u, sign=-1):
    # u dict over product index
    d = len(ts)
    out = {}
    for k in itertools.product(*[range(t) for t in ts]):
        acc = mp.mpc(0)
        for j in itertools.product(*[range(t) for t in ts]):
            acc += u[j] * mp.expj(sign * 2 * mp.pi * sum(mp.mpf(ji * ki) / ti for ji, ki, ti in zip(j, k, ts)))
        out[k] = acc / mp.fprod(ts)
    return out
ss, ts = (3, 5), (4, 8)
idxs = list(itertools.product(*[range(t) for t in ts]))
u = {j: mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for j in idxs}
Ftu = tensorF(ts, u)
Q = {k: Ftu[tuple((-si * ki) % ti for si, ki, ti in zip(ss, k, ts))] for k in idxs}
a = {j: mp.expj(-mp.pi * sum(mp.mpf(si * ji * ji) / ti for si, ji, ti in zip(ss, j, ts))) for j in idxs}
T = ts[0] * ts[1]
abu = {j: mp.conj(a[j]) * u[j] for j in idxs}
conv = {j: sum(a[tuple((ji - ki) % ti for ji, ki, ti in zip(j, k, ts))] * abu[k] for k in idxs) for j in idxs}
rhs = {j: mp.conj(a[j]) * conv[j] / T for j in idxs}
err = max(abs(Q[j] - rhs[j]) for j in idxs)
check("(d) Q F_t u = conj(a).(a*(conj(a).u))/T  (d=2, s=(3,5), t=(4,8))", err < mp.mpf(10)**-25, mp.nstr(err, 5))
check("(d) a has period t_i", all(abs(a[j] - mp.expj(-mp.pi * sum(mp.mpf(si * (ji + ti) * (ji + ti)) / ti if i == 0 else mp.mpf(si * ji * ji) / ti for i, (si, ji, ti) in enumerate(zip(ss, j, ts))))) < mp.mpf(10)**-30 for j in idxs))

# (e) source permutation cancellation, d=2, s=(3,5), t=(4,8) -> (P_s u)_j = u_{tj}
sidx = list(itertools.product(*[range(s) for s in ss]))
u = {j: mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for j in sidx}
v = {j: mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for j in sidx}
def R(z): return {j: z[tuple((ti * ji) % si for ti, ji, si in zip(ts, j, ss))] for j in sidx}
RFu = R(tensorF(ss, u)); RFv = R(tensorF(ss, v))
prod = {j: RFu[j] * RFv[j] for j in sidx}
lhs = R(tensorF(ss, prod, sign=+1))
S = ss[0] * ss[1]
conv = {j: sum(u[k] * v[tuple((ji - ki) % si for ji, ki, si in zip(j, k, ss))] for k in sidx) for j in sidx}
err = max(abs(lhs[j] - conv[j] / S**2) for j in sidx)
check("(e) (R F_s^+)[(R F_s u).(R F_s v)] = u*v/S^2", err < mp.mpf(10)**-25, mp.nstr(err, 5))
# F^+ R = R^-1 F^+ and F_s(u*v) = S (F_s u).(F_s v), F^+ F^- = I/S
Fconv = tensorF(ss, conv); Fu = tensorF(ss, u); Fv = tensorF(ss, v)
check("(e) F_s(u*v) = S (F_s u).(F_s v)", max(abs(Fconv[j] - S * Fu[j] * Fv[j]) for j in sidx) < mp.mpf(10)**-25)
FpFu = tensorF(ss, Fu, sign=+1)
check("(e) F_s^+ F_s = I/S", max(abs(FpFu[j] - u[j] / S) for j in sidx) < mp.mpf(10)**-25)

# (f) synthetic transforms over R_r = C[y]/(y^r+1): omega_t = y^{2r/t}
def mulmono(z, e, r):   # z * y^e mod y^r+1 ; z list length r
    e %= 2 * r
    out = [mp.mpc(0)] * r
    for k in range(r):
        idx = k + e
        sign = 1
        while idx >= r:
            idx -= r; sign = -sign
        out[idx] += sign * z[k]
    return out
def Ft(t, r, z, sign=-1):  # z: list of t ring elements (each list length r)
    out = []
    for j in range(t):
        acc = [mp.mpc(0)] * r
        for k in range(t):
            term = mulmono(z[k], sign * (2 * r // t) * j * k, r)
            acc = [x + y for x, y in zip(acc, term)]
        out.append([x / t for x in acc])
    return out
r, t = 8, 8
z = [[mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(r)] for _ in range(t)]
back = Ft(t, r, Ft(t, r, z, -1), +1)
err = max(abs(back[j][k] - z[j][k] / t) for j in range(t) for k in range(r))
check("(f) F_t^+ F_t^- = I/t in C[y]/(y^r+1), r=t=8", err < mp.mpf(10)**-25, mp.nstr(err, 5))
r, t = 8, 4
z = [[mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(r)] for _ in range(t)]
back = Ft(t, r, Ft(t, r, z, -1), +1)
check("(f) F_t^+ F_t^- = I/t, r=8, t=4", max(abs(back[j][k] - z[j][k] / t) for j in range(t) for k in range(r)) < mp.mpf(10)**-25)
# character sum: sum_k omega_t^{jk} = t if t|j else 0
r, t = 8, 8
for j in range(0, 2 * t):
    acc = [mp.mpc(0)] * r
    for k in range(t):
        one = [mp.mpc(0)] * r; one[0] = mp.mpc(1)
        acc = [x + y for x, y in zip(acc, mulmono(one, (2 * r // t) * j * k, r))]
    expected = t if j % t == 0 else 0
    ok = abs(acc[0] - expected) < 1e-25 and all(abs(x) < 1e-25 for x in acc[1:])
    if not ok: check(f"(f) character sum j={j}", ok)
check("(f) character cancellation sum_k omega^{jk}", True)
# DIF: (F_t z)_{2h+b} = (1/(t/2)) sum_{k<t/2} (omega_t^2)^{-hk} omega_t^{-bk} (z_k + (-1)^b z_{k+t/2})/2
r, t = 8, 8
z = [[mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(r)] for _ in range(t)]
full = Ft(t, r, z)
okdif = True
for h_ in range(t // 2):
    for bb in (0, 1):
        acc = [mp.mpc(0)] * r
        for k in range(t // 2):
            comb_ = [(x + (-1)**bb * y) / 2 for x, y in zip(z[k], z[k + t // 2])]
            term = mulmono(comb_, -(2 * r // t) * (2 * h_ * k + bb * k), r)
            acc = [x + y for x, y in zip(acc, term)]
        acc = [x / (t // 2) for x in acc]
        okdif &= max(abs(acc[i] - full[2 * h_ + bb][i]) for i in range(r)) < 1e-25
check("(f) decimation-in-frequency formula (eq:synthetic-dif)", okdif)
# recursive DIF procedure writes output j at slot rev(j): implement and compare
def dif_rec(z, t, r, a_):
    # returns list indexed by slot; after full recursion slot = rev(j)
    if t == 1: return [z[0]]
    half = t // 2
    top = [[(x + y) / 2 for x, y in zip(z[k], z[k + half])] for k in range(half)]
    bot = [mulmono([(x - y) / 2 for x, y in zip(z[k], z[k + half])], -(2 * r // t) * k, r) for k in range(half)]
    return dif_rec(top, half, r, a_ - 1) + dif_rec(bot, half, r, a_ - 1)
slots = dif_rec(z, t, r, 3)
def rev(j, a_): return int(format(j, f'0{a_}b')[::-1], 2)
okrev = all(max(abs(slots[rev(j, 3)][i] - full[j][i]) for i in range(r)) < 1e-25 for j in range(t))
check("(f) DIF recursion stores (F_t z)_j at slot rev(j)", okrev)
# normalized convolution: M F^+((F^- f)(F^- g)/r) = f*g/(rM) for one axis M=t, negacyclic product in R_r
def ringmul(x, y, r):
    out = [mp.mpc(0)] * r
    for i in range(r):
        for j in range(r):
            idx = i + j; sign = 1
            if idx >= r: idx -= r; sign = -1
            out[idx] += sign * x[i] * y[j]
    return out
t, r = 4, 8
f_ = [[mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(r)] for _ in range(t)]
g_ = [[mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(r)] for _ in range(t)]
Ff, Fg = Ft(t, r, f_), Ft(t, r, g_)
pw = [[x / r for x in ringmul(Ff[j], Fg[j], r)] for j in range(t)]
res = Ft(t, r, pw, +1)
res = [[t * x for x in row] for row in res]
conv = []
for j in range(t):
    acc = [mp.mpc(0)] * r
    for k in range(t):
        acc = [x + y for x, y in zip(acc, ringmul(f_[k], g_[(j - k) % t], r))]
    conv.append([x / (r * t) for x in acc])
check("(f) M F^+((F^- f)(F^- g)/r) = f*g/(rM)", max(abs(res[j][i] - conv[j][i]) for j in range(t) for i in range(r)) < 1e-25)
# coefficient twist: x_d -> zeta y turns cyclic convolution mod x^r-1 into negacyclic
zeta = mp.expj(mp.pi / r)
fa = [mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(r)]
ga = [mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(r)]
cyc = [sum(fa[i] * ga[(j - i) % r] for i in range(r)) for j in range(r)]
tw = ringmul([zeta**k * fa[k] for k in range(r)], [zeta**k * ga[k] for k in range(r)], r)
untw = [zeta**(-k) * tw[k] for k in range(r)]
check("(g) twist x->zeta y: cyclic conv = untwist(negacyclic conv of twisted)", max(abs(untw[k] - cyc[k]) for k in range(r)) < 1e-25)

# (g) CRT map Φ: x -> prod x_i^{mu_i}; coordinates b_i = a_i + mu_i sum_{j<i} a_j P_j mod s_i
sl = (3, 5, 7)
Sprod = 3 * 5 * 7
P = [1, 3, 15]
mu = [pow(P[i], -1, sl[i]) for i in range(3)]
def phi_coord(k):
    # mixed radix digits a_j
    a_ = []; kk = k
    for i in range(3):
        a_.append(kk % sl[i]); kk //= sl[i]
    b_ = [(a_[i] + mu[i] * sum(a_[j] * P[j] for j in range(i))) % sl[i] for i in range(3)]
    direct = [(k * mu[i]) % sl[i] for i in range(3)]   # exponent of x_i in x^k = prod x_i^{k mu_i}
    return tuple(b_), tuple(direct)
check("(g) (eq:triangular-crt) coordinates = k*mu_i mod s_i", all(phi_coord(k)[0] == phi_coord(k)[1] for k in range(Sprod)))
check("(g) Φ bijective on monomials", len({phi_coord(k)[0] for k in range(Sprod)}) == Sprod)
# convolution preserved
fa = [rng.randint(-3, 3) for _ in range(Sprod)]; ga = [rng.randint(-3, 3) for _ in range(Sprod)]
cyc = [sum(fa[i] * ga[(j - i) % Sprod] for i in range(Sprod)) for j in range(Sprod)]
fm = {}; gm = {}
for k in range(Sprod):
    fm[phi_coord(k)[0]] = fa[k]; gm[phi_coord(k)[0]] = ga[k]
okc = True
for k in range(Sprod):
    bk = phi_coord(k)[0]
    acc = 0
    for k1 in range(Sprod):
        b1 = phi_coord(k1)[0]
        b2 = tuple((bk[i] - b1[i]) % sl[i] for i in range(3))
        acc += fm[b1] * gm[b2]
    okc &= acc == cyc[k]
check("(g) Φ preserves cyclic convolution (3,5,7)", okc)

print("\nALL OK" if ok_all else "\nSOME CHECKS FAILED")
