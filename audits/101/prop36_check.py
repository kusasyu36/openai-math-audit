"""Direct numerical test of Proposition 3.6 (3.27): for each eigenvalue pair (ai,aj),
LHS-RHS is a quadratic form in (G,g,Z); check it is negative semidefinite.
Independent of the paper's proof (4.14)-(4.23)."""
import mpmath as mp, math, sys
mp.mp.dps = 40
M = mp.mpf
al = 3/40; s = 1/(1-al); beta = 7/25; c = 1.5; gam = .75; eps = .15
SQ2 = mp.sqrt(2)

def raw(a):
    a = M(a); phi = mp.npdf(a); Pm = mp.erfc(a/SQ2)/2; q = phi/Pm
    t = -mp.log1p(-mp.erfc(-a/SQ2)/2) if a < 0 else -mp.log(Pm)
    d = q-a; k = t/q; b = 1-k*d; l = k-q+a*b
    f = k*(1+b); fp = b*(1+b)+k*l
    p = b+k*l-l**2/M(al); e = fp-2*b**2-2*l**2/M(al)
    return dict(a=a, logq=mp.log(q), k=k, b=b, l=l, f=f, fp=fp, p=p, e=e)

def m0(a):
    if a < -3: return 0
    if a < -2: return .19
    if a < -1: return .38
    return .68

def form(Fi, Fj):
    ai, aj = Fi['a'], Fj['a']
    if ai == aj:
        dk, db, df = Fj['b'], Fj['l'], Fj['fp']   # divided difference = derivative
    else:
        h = ai-aj
        dk = (Fi['k']-Fj['k'])/h; db = (Fi['b']-Fj['b'])/h; df = (Fi['f']-Fj['f'])/h
    out = [[0.0]*3 for _ in range(3)]
    def add(vec, w):
        for x in range(3):
            for y in range(3):
                out[x][y] += w*vec[x]*vec[y]
    for (X, dK, dB, dF) in [(Fj, dk-Fj['b'], db-Fj['l'], df-Fj['fp']),
                            (Fi, dk-Fi['b'], db-Fi['l'], df-Fi['fp'])]:
        bj, lj, ej, pj = (float(X[n]) for n in ('b','l','e','p'))
        dK, dB = float(dK), float(dB)
        psi = float(dF) - 2*bj*dK - 2*lj/al*dB
        D = (bj, 0.0, dK)
        add((D[0], -beta, D[2]), s/2)
        add((D[0], -gam, D[2]), (1+s/8)/2)
        add((0, 1, 0), gam**2*(s-1)/2)
        add((lj, 0, dB), 1/al/2)
        add((ej, c*bj, psi), 1/(4*(1-eps)*pj)/2)
    # RHS
    add((1, 0, 0), -.43); add((-.5, 1, 0), -2.72)
    # cross .26*u*v, u=G, v=g-G/2  -> .13*(u v^T + v u^T)
    u = (1, 0, 0); v = (-.5, 1, 0)
    for x in range(3):
        for y in range(3):
            out[x][y] -= .13*(u[x]*v[y]+v[x]*u[y])
    dl = abs(float(Fi['logq']-Fj['logq']))
    r = (1+m0(min(float(ai), float(aj))))/16*(1+dl**2/30+dl**4/1680)
    out[2][2] -= r*dl**2
    return out, dl

def maxeig(A):
    m = mp.matrix(A)
    E = mp.eigsy(m, eigvals_only=True)
    return max(float(x) for x in E)

step = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
lo, hi = -12.0, 20.0
pts = [lo + i*step for i in range(int(round((hi-lo)/step))+1)]
# densify near bin boundaries
for bnd in (-3, -2, -1, .5, 4):
    pts += [bnd + o for o in (-1e-3, -1e-6, 0, 1e-6, 1e-3)]
pts = sorted(set(round(x, 9) for x in pts))
FF = [raw(x) for x in pts]
worst = (-1e9, None); worst_rel = (-1e9, None)
cnt = 0
for i in range(len(FF)):
    for j in range(i+1, len(FF)):
        A, dl = form(FF[i], FF[j])
        ev = maxeig(A)
        cnt += 1
        scale = max(abs(A[x][x]) for x in range(3)) + 1e-300
        if ev > worst[0]: worst = (ev, (pts[i], pts[j], dl))
        if ev/scale > worst_rel[0]: worst_rel = (ev/scale, (pts[i], pts[j], dl))
print('pairs tested', cnt)
print('max eigenvalue of (LHS-RHS):', worst)
print('max relative eigenvalue:', worst_rel)
