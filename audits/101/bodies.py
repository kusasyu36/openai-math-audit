"""Numerical test of Theorem 1.1: L_K <= L_simplex(n) on explicit bodies (own code)."""
import mpmath as mp, random, itertools, math
mp.mp.dps = 30
M = mp.mpf

def Lsimp_formula(n):
    n = M(n)
    return mp.factorial(n)**(1/n) / ((n+1)**((n+1)/(2*n)) * mp.sqrt(n+2))

def simplex_moments(V):
    """V: list of n+1 vertices (mp vectors). returns vol, int x, int x x^T"""
    n = len(V)-1
    Mx = mp.matrix([[V[i][j]-V[0][j] for j in range(n)] for i in range(1, n+1)])
    vol = abs(mp.det(Mx))/mp.factorial(n)
    S = [sum(V[i][j] for i in range(n+1)) for j in range(n)]
    first = [vol*S[j]/(n+1) for j in range(n)]
    sec = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            sec[a, b] = vol/((n+1)*(n+2))*(sum(V[i][a]*V[i][b] for i in range(n+1)) + S[a]*S[b])
    return vol, first, sec

def L_from_simplices(simps, n):
    vol = M(0); fst = [M(0)]*n; sec = mp.matrix(n, n)
    for V in simps:
        v, f, s2 = simplex_moments(V)
        vol += v; fst = [fst[j]+f[j] for j in range(n)]; sec += s2
    mu = [x/vol for x in fst]
    C = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            C[a, b] = sec[a, b]/vol - mu[a]*mu[b]
    return (mp.det(C)/vol**2)**(1/(2*M(n)))

def hull_simplices(P, n):
    """brute-force convex hull facets for points in general position; cone from interior pt."""
    N = len(P)
    cen = [sum(p[j] for p in P)/N for j in range(n)]
    simps = []
    for comb in itertools.combinations(range(N), n):
        A = mp.matrix([[P[i][j] for j in range(n)] + [1] for i in comb])
        # normal: solve via null-space using cofactors
        normal = []
        for col in range(n+1):
            sub = mp.matrix([[A[r, cc] for cc in range(n+1) if cc != col] for r in range(n)])
            normal.append((-1)**col*mp.det(sub))
        val = lambda x: sum(normal[j]*x[j] for j in range(n)) + normal[n]
        sides = [val(P[i]) for i in range(N) if i not in comb]
        if all(x >= -mp.mpf(10)**-20 for x in sides) or all(x <= mp.mpf(10)**-20 for x in sides):
            simps.append([P[i] for i in comb] + [cen])
    return simps

Ls = {n: Lsimp_formula(n) for n in range(1, 13)}
print('n, L_simplex(formula), L_simplex(vertex calc), 1/e - L')
for n in range(1, 9):
    V = [[M(0)]*n] + [[M(1) if j == i else M(0) for j in range(n)] for i in range(n)]
    Lv = L_from_simplices([V], n) if n > 1 else 1/mp.sqrt(12)
    print(n, mp.nstr(Ls[n], 15), mp.nstr(Lv, 15), mp.nstr(1/mp.e - Ls[n], 6))
print('large n check: L_simplex(n) < 1/e:', all(Lsimp_formula(n) < 1/mp.e for n in list(range(1, 200))+[10**3, 10**5, 10**7]))

print('\n--- closed-form bodies, ratio L_K / L_simplex (must be <1) ---')
def Lball(n): n = M(n); return (n+2)**-0.5 * (mp.pi**(n/2)/mp.gamma(n/2+1))**(-1/n)
def Lcross(n): n = M(n); return mp.factorial(n)**(1/n)/mp.sqrt(2*(n+1)*(n+2))
def Lprod(L1, k, L2, l): return (L1**(2*k)*L2**(2*l))**(1/(2*M(k+l)))
worst = 0
for n in range(2, 13):
    rows = {'cube': 1/mp.sqrt(12), 'ball': Lball(n), 'cross': Lcross(n)}
    for k in range(1, n):
        rows[f'simplex{k}xsimplex{n-k}'] = Lprod(Ls[k], k, Ls[n-k], n-k)
        rows[f'simplex{k}xcube{n-k}'] = Lprod(Ls[k], k, 1/mp.sqrt(12), n-k)
        rows[f'simplex{k}xball{n-k}'] = Lprod(Ls[k], k, Lball(n-k) if n-k > 1 else 1/mp.sqrt(12), n-k)
    best = max(rows.items(), key=lambda kv: kv[1])
    worst = max(worst, best[1]/Ls[n])
    print(n, ' '.join(f'{k}:{mp.nstr(v/Ls[n],5)}' for k, v in list(rows.items())[:3]), '| max', best[0], mp.nstr(best[1]/Ls[n], 8))
print('max ratio over closed-form/product bodies:', mp.nstr(worst, 10))

print('\n--- profile bodies K={(x,t): x in rho(t)K0, 0<=t<=1}, K0 in R^{n-1} ---')
def L_profile(n, rho, K0):
    """K0: (vol, detcov) of centered base in dim n-1"""
    v0, dc0 = K0
    I = lambda pw: mp.quad(lambda t: rho(t)**pw, [0, 1])
    In1, In1p = I(n-1), I(n+1)
    m1 = mp.quad(lambda t: t*rho(t)**(n-1), [0, 1])/In1
    m2 = mp.quad(lambda t: t*t*rho(t)**(n-1), [0, 1])/In1
    vol = v0*In1
    det = (m2-m1**2)*dc0*(In1p/In1)**(n-1)
    return (det/vol**2)**(1/(2*M(n)))
def base(kind, k):
    if k == 1: return (M(1), M(1)/12)
    if kind == 'cube': return (M(1), M(12)**-k)
    if kind == 'ball':
        v = mp.pi**(M(k)/2)/mp.gamma(M(k)/2+1); return (v, (M(1)/(k+2))**k)
    if kind == 'simplex': L = Ls[k]; return (M(1), L**(2*k))  # vol 1 normalisation
    if kind == 'cross': L = Lcross(k); return (M(1), L**(2*k))
random.seed(1)
worstp = (0, None)
for n in range(2, 8):
    for kind in ('simplex', 'cube', 'ball', 'cross'):
        K0 = base(kind, n-1)
        fams = [('cone', lambda t: t), ('frustum.5', lambda t: 1-t/2), ('frustum.9', lambda t: 1-0.9*t),
                ('frustum.99', lambda t: 1-0.99*t), ('cyl', lambda t: M(1)), ('dcone', lambda t: min(t, 1-t)),
                ('sqrt', lambda t: mp.sqrt(t)), ('t^.9', lambda t: t**0.9)]
        for _ in range(6):
            sl = sorted([random.uniform(-3, 3) for _ in range(3)], reverse=True)
            ints = [random.uniform(0.2, 1.5) for _ in range(3)]
            fams.append(('pl', (lambda sl, ints: (lambda t: max(M(1e-12), min(ints[i]+sl[i]*t for i in range(3)))))(sl, ints)))
        for name, rho in fams:
            r = L_profile(n, rho, K0)/Ls[n]
            if r > worstp[0]: worstp = (r, (n, kind, name))
            if name == 'cone' and kind == 'simplex':
                print(f'n={n} cone over simplex ratio (should be 1): {mp.nstr(r, 15)}')
print('max ratio over profile bodies:', mp.nstr(worstp[0], 12), worstp[1])

print('\n--- random polytopes (brute-force hull), ratio L_K/L_simplex ---')
for n in (2, 3, 4):
    worst = 0
    for trial in range(40 if n < 4 else 12):
        mode = trial % 3
        if mode == 0:   # random points in cube
            P = [[M(random.random()) for _ in range(n)] for _ in range(n+3+trial % 4)]
        elif mode == 1:  # simplex + small perturbation points
            S = [[M(0)]*n] + [[M(1) if j == i else M(0) for j in range(n)] for i in range(n)]
            P = [[x + M(random.gauss(0, 0.02)) for x in v] for v in S] + [[M(random.random())*0.3+0.1 for _ in range(n)] for _ in range(2)]
            P.append([M(0.5+random.uniform(0, 0.05))]*1 + [M(0.5)] + [M(0)]*(n-2) if n >= 2 else [M(0.5)])
        else:  # gaussian points
            P = [[M(random.gauss(0, 1)) for _ in range(n)] for _ in range(n+2+trial % 5)]
        simps = hull_simplices(P, n)
        r = L_from_simplices(simps, n)/Ls[n]
        worst = max(worst, r)
    print(f'n={n}: max L_K/L_simplex over random polytopes = {mp.nstr(worst, 10)}')
