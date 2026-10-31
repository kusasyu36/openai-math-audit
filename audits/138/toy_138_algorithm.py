"""
Audit 138: small-n implementation of the paper's pipeline, checked against brute force.
Everything here is written from the TeX description (setup.tex, checksum.tex,
aliases.tex, implementation.tex); nothing from the OpenAI repository is executed.

Checks:
 T1  Isolation lemma (setup.tex:32-61): equivalence + isolation frequency >= 11/12.
 T2  Deficiency lemma (setup.tex:78-138): increments in [0,1], superadditivity,
     exchange bound, Var d(V_q) <= q (exact enumeration).
 T3  Fourier identity (checksum.tex:52-66): V_rho = (1/P) sum_k J_rho(k).
 T4  Alias extraction (aliases.tex:70-189): representation + small-modulus filter +
     exact-integer exclusion.  Every saved set is a genuine alias; every visited pair
     has D != 0 and P | D (aliases.tex:217-228); coverage grows with trials.
 T5  Checksum estimator mechanics (checksum.tex:313-425): sum_k mu(k) = L exactly,
     exact unbiasedness  E[X] = v_pi over covered frequencies, with real f_*, g_* scores.
 T6  End-to-end decision (implementation.tex:11-64) vs brute force.
"""
import random, math, cmath, itertools, sys
from fractions import Fraction
random.seed(138)

def is_prime(p):
    if p < 2: return False
    for d in range(2, int(p**0.5)+1):
        if p % d == 0: return False
    return True
def rand_prime(lo, hi):
    ps = [p for p in range(lo, hi) if is_prime(p)]
    return random.choice(ps)

# ---------------------------------------------------------------- T1 isolation
def transform(a, t):
    n = len(a); N0 = 1
    while N0 < 12*n: N0 *= 2
    assert 12*n <= N0 < 24*n
    u = [random.randint(1, N0) for _ in range(n)]
    C = 1 + n*N0
    c = [C*ai + ui for ai, ui in zip(a, u)]
    targets = [C*t + s for s in range(n*N0 + 1)]
    return c, targets, N0

def subset_sums(w):
    n = len(w); out = [0]*(1 << n)
    for mask in range(1, 1 << n):
        low = mask & -mask; i = low.bit_length()-1
        out[mask] = out[mask ^ low] + w[i]
    return out

def test_isolation(trials=400):
    n = 10; bad = 0; equiv_fail = 0
    for _ in range(trials):
        a = [random.randint(1, 6) for _ in range(n)]        # many equal sums -> many solutions
        sols = subset_sums(a); t = sols[random.randrange(1, 1 << n)]
        orig_yes = t in sols
        c, targets, N0 = transform(a, t)
        cs = subset_sums(c)
        cnt = {}
        for m in range(1 << n): cnt[cs[m]] = cnt.get(cs[m], 0) + 1
        trans_yes = any(ts in cnt for ts in targets)
        if trans_yes != orig_yes: equiv_fail += 1
        if orig_yes and not any(cnt.get(ts, 0) == 1 for ts in targets): bad += 1
    print(f"T1 isolation: equivalence failures {equiv_fail}/{trials}; no-isolated-bin {bad}/{trials} "
          f"(allowed <= 1/12 = {trials/12:.1f}); #solutions per instance is large here")
    return equiv_fail == 0 and bad <= trials/12

# ---------------------------------------------------------------- T2 deficiency
def deficiency(c, V):
    if not V: return 0.0
    w = [c[i] for i in V]
    return len(V) - math.log2(len(set(subset_sums(w))))

def test_deficiency():
    n = 8; okc = True
    for _ in range(20):
        c = [random.randint(1, 5) for _ in range(n)]
        for V in itertools.chain.from_iterable(itertools.combinations(range(n), k) for k in range(n+1)):
            V = set(V); dV = deficiency(c, V)
            okc &= dV >= -1e-12
            for i in range(n):
                if i not in V:
                    inc = deficiency(c, V | {i}) - dV
                    okc &= -1e-9 <= inc <= 1 + 1e-9
            # superadditivity with a random disjoint U
            U = set(random.sample([i for i in range(n) if i not in V], random.randint(0, n-len(V))))
            okc &= deficiency(c, V | U) >= deficiency(c, V) + deficiency(c, U) - 1e-9
        # Var d(V_q) <= q exactly over all q-subsets
        for q in range(n+1):
            vals = [deficiency(c, set(S)) for S in itertools.combinations(range(n), q)]
            mu = sum(vals)/len(vals); var = sum((v-mu)**2 for v in vals)/len(vals)
            okc &= var <= q + 1e-9
        # monotone expectation in q
        Ed = [sum(deficiency(c, set(S)) for S in itertools.combinations(range(n), q))/math.comb(n, q) for q in range(n+1)]
        okc &= all(Ed[q] <= Ed[q+1] + 1e-9 for q in range(n))
    print("T2 deficiency lemma (increments, superadditivity, Var<=q, monotone E):", "OK" if okc else "FAIL")
    return okc

# ---------------------------------------------------------------- T3 Fourier identity
def eP(P, j): return cmath.exp(2j*math.pi*(j % P)/P)
def checksum_exact(b, phi, P, rho):
    n = len(b); tot = 0
    bs = subset_sums(b); ps = subset_sums(phi)
    for m in range(1 << n):
        if bs[m] % P == rho % P: tot += eP(P, ps[m])
    return tot
def J(b, phi, P, rho, k):
    val = eP(P, -k*rho)
    for bi, pi in zip(b, phi): val *= 1 + eP(P, pi + k*bi)
    return val
def test_fourier():
    n, P = 9, 37; b = [random.randrange(P) for _ in range(n)]; phi = [random.randrange(P) for _ in range(n)]
    worst = 0
    for rho in range(P):
        V = checksum_exact(b, phi, P, rho)
        Vf = sum(J(b, phi, P, rho, k) for k in range(P))/P
        worst = max(worst, abs(V-Vf))
    print(f"T3 Fourier identity: max |V_rho - (1/P)sum_k J_rho(k)| = {worst:.2e}")
    return worst < 1e-8

# ---------------------------------------------------------------- T4 alias extraction
def alias_trial(c, b, P, t_s, q, m, Pprime_range):
    """One trial of aliases.tex:92-189 for fixed (s,q). Returns (saved aliases, stats)."""
    n = len(c); T_sq = (t_s % P) + q*P
    M = sorted(random.sample(range(n), m)); Mset = set(M)
    F = {i for i in range(n) if random.random() < 0.5}
    sgn = [(-1 if i in F else 1) for i in range(n)]
    bF = [sgn[i]*b[i] for i in range(n)]; cF = [sgn[i]*c[i] for i in range(n)]
    TF = T_sq - sum(b[i] for i in F); tF = t_s - sum(c[i] for i in F)
    rest = [i for i in range(n) if i not in Mset]
    saved = set(); visited = raw = 0; viol = 0
    for k in range(m//2 + 1):
        # balance l (aliases.tex:122-135)
        l = round((n - m + math.log2(math.comb(m, m//2 - k)) - math.log2(math.comb(m, k)))/2)
        l = max(0, min(n-m, l)); L, Rr = rest[:l], rest[l:]
        lo, hi = Pprime_range
        Pp = rand_prime(lo, hi); jp = random.randrange(Pp)
        def records(side, size):
            recs = []
            for A in itertools.chain.from_iterable(itertools.combinations(side, r) for r in range(len(side)+1)):
                for X in itertools.combinations(M, size):
                    red = sum(bF[i] for i in A) + sum(bF[i] for i in X)
                    unred = sum(cF[i] for i in A) + sum(cF[i] for i in X)
                    recs.append((A, X, red, unred))
            return recs
        left = records(L, k); right = records(Rr, m//2 - k)
        raw += len(left)*len(right)
        leftF = [r for r in left if r[2] % Pp == jp]
        rightF = [r for r in right if r[2] % Pp == (TF - jp) % Pp]
        # group right by reduced weight, then by unreduced (the exclusion is an interval; here a dict)
        grp = {}
        for r in rightF: grp.setdefault(r[2], []).append(r)
        for (A, X, red, unred) in leftF:
            for (B, Y, red2, unred2) in grp.get(TF - red, []):
                if unred + unred2 == tF: continue           # exact-integer exclusion
                visited += 1
                # aliases.tex:217-228: D != 0 and P | D must hold for every visited pair
                mu = {}
                for i in A + B + X + Y: mu[i] = mu.get(i, 0) + 1
                D = sum(sgn[i]*mu[i]*c[i] for i in mu) - t_s + sum(c[i] for i in F)
                if D == 0 or D % P != 0: viol += 1
                if set(X) & set(Y): continue
                I = frozenset(set(A) | set(B) | set(X) | set(Y)) ^ frozenset(F)
                cI = sum(c[i] for i in I)
                if cI % P == t_s % P and cI != t_s: saved.add(I)
                else: viol += 1
    return saved, (raw, visited, viol)

def test_aliases(n=12, m=4, trials=25):
    a = [random.randint(1, 1 << 8) for _ in range(n)]; t = sum(random.sample(a, 5))
    c, targets, N0 = transform(a, t)
    P = rand_prime(64, 128); b = [ci % P for ci in c]
    cs = subset_sums(c)
    allok = True; totraw = totvis = 0
    results = []
    for s in random.sample(range(len(targets)), 4):
        t_s = targets[s]
        truth = {frozenset(i for i in range(n) if mask >> i & 1) for mask in range(1 << n)
                 if cs[mask] % P == t_s % P and cs[mask] != t_s}
        found = set(); cov = []
        for tr in range(trials):
            for q in range(0, sum(b)//P + 1):
                sv, (raw, vis, viol) = alias_trial(c, b, P, t_s, q, m, (2, 4))
                allok &= viol == 0 and sv <= truth
                found |= sv; totraw += raw; totvis += vis
            cov.append(len(found)/max(1, len(truth)))
        results.append((len(truth), cov[0], cov[len(cov)//2], cov[-1]))
    for r in results:
        print(f"T4 aliases: |A_s(P)|={r[0]:4d}; coverage after 1/{trials//2}/{trials} trials = {r[1]:.2f}/{r[2]:.2f}/{r[3]:.2f}")
    print(f"T4 aliases: false alias or (D==0 or P∤D) violations: {'none' if allok else 'FOUND'};"
          f" raw pairs {totraw}, visited {totvis} (ratio {totvis/totraw:.4f}; 1/P~{1/P:.4f}, P={P})")
    return allok

# ---------------------------------------------------------------- T5 estimator mechanics
def finite_filter_scores(P):
    rc = [0.851, 0.466, 0.325, 0.253, 0.224, -0.314]; sc = [1.281, -0.938, 0.731, -0.559, 0.391, 0.110]
    kappa = 1.0315
    ang = [2*math.pi*j/P for j in range(P)]
    Fv = [2**sum(cj*math.cos((j+1)*x) for j, cj in enumerate(rc)) for x in ang]
    Gv = [2**sum(cj*math.cos((j+1)*x) for j, cj in enumerate(sc)) for x in ang]
    f = [v/(sum(Fv)/P) for v in Fv]; g = [v/(sum(Gv)/P) for v in Gv]
    fk = [v**kappa for v in f]; gk = [v**kappa for v in g]
    fs = [v/(sum(fk)/P) for v in fk]; gs = [v/(sum(gk)/P) for v in gk]
    lfs = [math.log2(v) for v in fs]; lgs = [math.log2(v) for v in gs]
    # D, A, E per difference x (finite circle)
    DAE = []
    for x in range(P):
        prod = [f[(x+w) % P]*g[w] for w in range(P)]; Q = sum(prod)/P; K = [p/Q for p in prod]
        D = sum(K[w]*math.log2(K[w]) for w in range(P))/P
        A = sum(K[w]*lfs[(x+w) % P] for w in range(P))/P
        E = sum(K[w]*lgs[w] for w in range(P))/P
        DAE.append((D, A, E))
    return lfs, lgs, DAE

def test_estimator(n=8, P=31, R=600, grid=0.1):
    b = [random.randrange(P) for _ in range(n)]; phi = [random.randrange(P) for _ in range(n)]
    lfs, lgs, DAE = finite_filter_scores(P)
    B = 8                                                             # digit base, k = h*B + l
    H = (P + B - 1)//B
    # one block of size n; pattern = rounded (D,A,E) averages; alpha = 1
    def pattern(k):
        D = sum(DAE[(phi[i] + k*b[i]) % P][0] for i in range(n))/n
        A = sum(DAE[(phi[i] + k*b[i]) % P][1] for i in range(n))/n
        E = sum(DAE[(phi[i] + k*b[i]) % P][2] for i in range(n))/n
        return (math.ceil(D/grid)*grid, math.floor(A/grid)*grid, math.floor(E/grid)*grid)
    pats = {}
    for k in range(P): pats.setdefault(pattern(k), []).append(k)
    rho = random.randrange(P)
    Vexact = sum(J(b, phi, P, rho, k) for k in range(P))/P
    total_est = 0; covered_all = True; ident_ok = True; worst_unbias = 0
    for pi, ks in pats.items():
        Dp, Ap, Ep = pi
        tables = [[random.randrange(P) for _ in range(n)] for _ in range(R)]   # alpha = 1: no '*' labels
        def passes_hi(h, tb): return sum(lfs[(phi[i] + h*B*b[i] - tb[i]) % P] for i in range(n))/n >= Ap
        def passes_lo(l, tb): return sum(lgs[(-l*b[i] - tb[i]) % P] for i in range(n))/n >= Ep
        hi = {r: [h for h in range(H) if passes_hi(h, tables[r])] for r in range(R)}
        lo = {r: [l for l in range(B) if passes_lo(l, tables[r])] for r in range(R)}
        L = sum(len(hi[r])*len(lo[r]) for r in range(R))                 # eq checksum-match-count
        mu = {}
        for r in range(R):
            for h in hi[r]:
                for l in lo[r]: mu[h*B + l] = mu.get(h*B + l, 0) + 1
        ident_ok &= sum(mu.values()) == L                                # eq checksum-multiplicity consistency
        cov = [k for k in ks if mu.get(k, 0) > 0]
        covered_all &= len(cov) == len(ks)
        ncov = globals().setdefault('_ncov', [0, 0]); ncov[0] += len(cov); ncov[1] += len(ks)
        # exact expectation of the sample value X over a uniform match
        if L > 0:
            EX = sum((mu[k]/L)*(L/(P*mu[k]))*J(b, phi, P, rho, k) for k in cov)
            v_cov = sum(J(b, phi, P, rho, k) for k in cov)/P
            worst_unbias = max(worst_unbias, abs(EX - v_cov))
            total_est += EX
    print(f"T5 estimator: covered frequencies {_ncov[0]}/{_ncov[1]}")
    print(f"T5 estimator: #patterns={len(pats)}, sum_k mu(k)=L identity {'OK' if ident_ok else 'FAIL'}, "
          f"all frequencies covered: {covered_all}, max |E[X]-v_pi(covered)| = {worst_unbias:.1e}, "
          f"|sum_pi E[X_pi] - V_rho| = {abs(total_est - Vexact):.1e} (|V_rho|={abs(Vexact):.3f})")
    return ident_ok and worst_unbias < 1e-9 and (not covered_all or abs(total_est - Vexact) < 1e-9)

# ---------------------------------------------------------------- T6 end to end
def decide(a, t, m=4, trials=40, use_exact_V=True):
    n = len(a); c, targets, N0 = transform(a, t)
    P = rand_prime(64, 128); b = [ci % P for ci in c]
    phi = [random.randrange(P) for _ in range(n)]
    cs = subset_sums(c); ps = subset_sums(phi)
    # checksum (exact via Fourier; the sampling estimator's variance is only meaningful asymptotically)
    Vhat = {}
    for s, t_s in enumerate(targets):
        rho = t_s % P
        if rho not in Vhat:
            Vhat[rho] = sum(J(b, phi, P, rho, k) for k in range(P))/P
    complete = True; answer = False
    for s, t_s in enumerate(targets):
        found = set()
        for tr in range(trials):
            for q in range(0, sum(b)//P + 1):
                sv, _ = alias_trial(c, b, P, t_s, q, m, (2, 4)); found |= sv
        truth = {frozenset(i for i in range(n) if mask >> i & 1) for mask in range(1 << n)
                 if cs[mask] % P == t_s % P and cs[mask] != t_s}
        complete &= found == truth
        Bhat = Vhat[t_s % P] - sum(eP(P, sum(phi[i] for i in I)) for I in found)
        if abs(Bhat) >= 0.5: answer = True
    return answer, complete

def test_end_to_end(n=10, inst=6):
    okc = 0; comp = 0; tot = 0
    for j in range(inst):
        a = [random.randint(1, 1 << 7) for _ in range(n)]
        sums = set(subset_sums(a))
        if j % 2 == 0: t = sum(random.sample(a, 4))
        else:
            t = random.randint(1, sum(a))
            while t in sums: t = random.randint(1, sum(a))      # force a genuine NO instance
        truth = t in sums
        ans, complete = decide(a, t)
        tot += 1; okc += (ans == truth); comp += complete
        print(f"T6 instance {j}: truth={truth} algorithm={ans} aliases complete={complete}")
    print(f"T6 end-to-end: correct {okc}/{tot}, alias lists complete {comp}/{tot}")
    return okc >= 2*tot/3

if __name__ == "__main__":
    res = [test_isolation(), test_deficiency(), test_fourier(), test_aliases(), test_estimator(), test_end_to_end()]
    print("TOY SUMMARY:", sum(res), "/", len(res), "checks passed")
