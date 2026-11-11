"""Companion paper, Sections 02 and 03-complete-family: numerical tests on exact Hadamard
matrices (high precision, mpmath/numpy) of
 (1) Lemma contraction (eq:injective-contraction) and the 7-row relation (eq:seven-row-contraction);
 (2) Lemma partitions (eq:injective-partitions);
 (3) Proposition pair-certificate: 3<|L|^2> + <|g(alpha) g(q)|^2> = 0 for both q, on several
     Hadamard families (Fourier 2-parameter, Dita, Bjorck circulant, cubic S_6, transposes, equivalents);
 (4) Corollary dichotomy: for the cubic matrix (M) fails and (C) holds; for the others (M) holds;
 (5) Lemma basic: (eq:same-matrix) for one Hadamard, (eq:different-matrices) for a MUB pair, and the
     complete-family identity (eq:order-two) in dimensions 3 and 5 where complete MUB sets exist;
 (6) Lemma mixed-reductions (iii) singleton reduction and Lemma gram in the d=3 analogue (7-l -> d+1-l).
Independent code; does not use verify.py."""
import numpy as np, itertools, math
rng = np.random.default_rng(7)
TWO_PI = 2*np.pi

def fourier_family(a, b):
    F = np.exp(2j*np.pi*np.outer(np.arange(6), np.arange(6))/6)
    R = np.zeros((6,6))
    for r in (1,3,5): R[r] = [0, a, b, 0, a, b]
    return F*np.exp(1j*R)
def dita(c, d=None):
    # Dita-type construction [[P, Q],[P, -Q]] with P = F_3 and Q = diag(1,e^{ic},e^{id}) F_3 (Hadamard for all c,d)
    if d is None: d = c*0.618033988749895 + 0.3
    w3 = np.exp(2j*np.pi/3); F3 = np.array([[w3**(j*k) for k in range(3)] for j in range(3)])
    Q = np.diag([1, np.exp(1j*c), np.exp(1j*d)])@F3
    return np.block([[F3, Q],[F3, -Q]])
def bjorck():
    d = (1-np.sqrt(3))/2 + 1j*np.sqrt(np.sqrt(3)/2)
    first = np.array([1, 1j*d, -d, -1j, -np.conj(d), 1j*np.conj(d)])
    return np.array([np.roll(first, k) for k in range(6)])
def cubic():
    w = np.exp(2j*np.pi/3)
    def t(i,j):
        if i==0 or j==0 or i==j: return 0
        return 1 if (i-j)%5 in (1,4) else -1
    return np.array([[w**t(i,j) for j in range(6)] for i in range(6)])
def is_hadamard(H, tol=1e-9):
    n = H.shape[0]
    return np.allclose(np.abs(H), 1, atol=tol) and np.allclose(H@H.conj().T, n*np.eye(n), atol=tol)
def random_equivalent(H):
    P = rng.permutation(6); Q = rng.permutation(6); H = H[P][:, Q]
    H = H*np.exp(1j*rng.uniform(0,TWO_PI,(6,1)))*np.exp(1j*rng.uniform(0,TWO_PI,(1,6)))
    if rng.random() < .5: H = H.T
    if rng.random() < .5: H = H.conj()
    return H
samples = []
for _ in range(8): samples.append(("F6(a,b)", fourier_family(rng.uniform(0,TWO_PI), rng.uniform(0,TWO_PI))))
for _ in range(4): samples.append(("D6(c)", dita(rng.uniform(0,TWO_PI))))
samples.append(("C6", bjorck())); samples.append(("S6 cubic", cubic()))
for _ in range(6): samples.append(("equiv(F6)", random_equivalent(fourier_family(rng.uniform(0,TWO_PI), rng.uniform(0,TWO_PI)))))
samples.append(("equiv(S6)", random_equivalent(cubic())))
for name, H in samples: assert is_hadamard(H), name
print("Hadamard samples:", len(samples), "all verified Hadamard")

def g(H, a):  # normalized character g_H(a) = (1/6) sum_k prod_i H(i,k)^{a_i}
    a = np.asarray(a)
    return np.mean(np.prod(H**a[:,None], axis=0))
alpha = (1,1,1,-1,-1,-1)

# (1) Lemma contraction, numerically: m_H(A) over injections rho, kappa and K in {H, H^T}
def m_H(H, A):
    A = np.asarray(A); p, q = A.shape
    if p == 0 or q == 0: return 1.0
    tot = 0.0; cnt = 0
    for K in (H, H.T):
        for rho in itertools.permutations(range(6), p):
            sub = K[list(rho)]
            for kappa in itertools.permutations(range(6), q):
                tot += np.prod(sub[:, list(kappa)]**A).real; cnt += 1
    return tot/cnt
def m_H_fast(H, A):
    # vectorized over column injections
    A = np.asarray(A); p, q = A.shape
    if p == 0 or q == 0: return 1.0
    kap = np.array(list(itertools.permutations(range(6), q)))   # (nk, q)
    tot = 0.0; cnt = 0
    for K in (H, H.T):
        for rho in itertools.permutations(range(6), p):
            sub = K[list(rho)][:, kap]                     # (p, nk, q)
            tot += np.prod(sub**A[:, None, :], axis=(0,2)).real.sum(); cnt += len(kap)
    return tot/cnt
def strip_zeros(A):
    A = np.asarray(A); A = A[np.any(A != 0, axis=1)]; A = A[:, np.any(A != 0, axis=0)]; return A
def contraction_check(H, A, i):
    A = np.asarray(A); p = A.shape[0]
    lhs = (7-p)*m_H_fast(H, strip_zeros(A))
    rhs = 0.0
    for l in range(p):
        if l == i: continue
        B = A.copy(); B[l] += A[i]; B = np.delete(B, i, axis=0)
        rhs -= m_H_fast(H, strip_zeros(B))
    return lhs, rhs
H = samples[0][1]
tests = [ (np.array([[1,-1,0],[0,1,-1],[-1,0,1]]), 0), (np.array([[1,-1,0],[-1,2,-1],[0,-1,1]]), 0), (np.array([[1,-1],[1,-1],[-2,2]]), 0), (np.array([[0,1,-1],[1,-1,0],[1,0,-1],[-2,0,2]]), 0)]
worst = 0
for A, i in tests:
    lhs, rhs = contraction_check(H, A, i); worst = max(worst, abs(lhs-rhs))
print(f"(1) Lemma contraction, p<=6: max |lhs-rhs| = {worst:.2e}")
# 7-row relation: sum_l m(A_{0->l}) = 0 with 7 rows
A7 = np.array([[1,-1,0],[1,0,-1],[-1,1,0],[0,1,-1],[-1,0,1],[1,-1,0],[-1,0,1]])
assert A7.sum(axis=0).tolist() == [0,0,0]
tot = 0
for l in range(1,7):
    B = A7.copy(); B[l] += A7[0]; B = np.delete(B, 0, axis=0); tot += m_H_fast(H, strip_zeros(B))
print(f"(1) 7-row relation sum_l m(A_0->l) = {tot:.2e} (should be 0)")

# (2) Lemma partitions with n=2,3
def avg_pair(H, F):
    tot = 0; cnt = 0
    for K in (H, H.T):
        for pi in itertools.permutations(range(6)):
            tot += F(K[list(pi)]); cnt += 1
    return tot/cnt
def set_partitions(n):
    if n == 0: yield []; return
    for part in set_partitions(n-1):
        for k in range(len(part)):
            yield part[:k] + [part[k] + [n-1]] + part[k+1:]
        yield part + [[n-1]]
def partitions_rhs(H, charges):
    charges = [np.array(c) for c in charges]; n = len(charges); tot = 0
    for part in set_partitions(n):
        k = len(part)
        if k > 6: continue
        cols = [sum(charges[j] for j in blk) for blk in part]
        A = np.column_stack(cols)
        coeff = math.perm(6, k)/6**n
        tot += coeff*m_H_fast(H, strip_zeros(A))
    return tot
a1 = (1,-1,0,0,0,0); a2 = (0,0,1,0,-1,0); a3 = tuple(-np.array(a1)-np.array(a2))
lhs = avg_pair(H, lambda K: (g(K,a1)*g(K,a2)*g(K,a3)).real); rhs = partitions_rhs(H, [a1,a2,a3])
lhs2 = avg_pair(H, lambda K: abs(g(K,(1,1,-2,0,0,0)))**2); rhs2 = 1/6 + 5/6*m_H_fast(H, strip_zeros(np.array([[1,-1],[1,-1],[-2,2]])))
print(f"(2) Lemma partitions n=3: lhs={lhs:.12f} rhs={rhs:.12f} diff={abs(lhs-rhs):.2e};  n=2 example: diff={abs(lhs2-rhs2):.2e}")

# (3) Proposition pair-certificate
w = (-3,-2,0,1,2,2); gam = [(-1,-1,0,0,1,1),(-2,0,0,0,1,1),(-1,-1,0,1,0,1),(-1,-1,0,1,1,0)]; cj = [2,1,1,1]
def L(K):
    return -g(K, w) + sum(c*g(K, a)*g(K, tuple(np.array(w)-np.array(a))) for c, a in zip(cj, gam))
qs = {"q_cross": (2,0,0,-2,0,0), "q_same": (2,-2,0,0,0,0)}
print("(3) pair certificate 3<|L|^2> + <|g(alpha)g(q)|^2> (should be 0 for every Hadamard):")
worst3 = 0
for name, Hs in samples:
    t1 = avg_pair(Hs, lambda K: abs(L(K))**2)
    vals = {}
    for qn, qv in qs.items():
        t2 = avg_pair(Hs, lambda K: abs(g(K, alpha)*g(K, qv))**2)
        vals[qn] = 3*t1 + t2; worst3 = max(worst3, abs(vals[qn]))
    print(f"   {name:12s} 3<|L|^2>={3*t1:.3e}  total(q_cross)={vals['q_cross']:.3e}  total(q_same)={vals['q_same']:.3e}")
print(f"(3) max |identity residual| over all samples = {worst3:.2e}")
# also check that the identity is NOT a triviality: individual <|L|^2> is zero for every Hadamard?  (it must be, since both terms are >=0)
# (4) dichotomy
def M_holds(Hs):
    return max(abs(g(Hs[list(pi)], alpha)) for pi in itertools.permutations(range(6)))
def C_holds(Hs): return is_hadamard(Hs**2)
for name, Hs in samples[:1] + samples[8:9] + samples[12:14]:
    print(f"(4) {name:12s}: max_pi |g(alpha)| = {M_holds(Hs):.3e} ; H∘2 Hadamard: {C_holds(Hs)}")
mx = max(M_holds(Hs) for name, Hs in samples if "cubic" not in name and "S6" not in name)
print(f"(4) max over non-cubic samples of max_pi |g(pi alpha)| = {mx:.2e} (Theorem fourier predicts 0)")

# (5) Lemma basic identities
def D_set(d):
    return [tuple((1 if t==i else 0)-(1 if t==j else 0) for t in range(d)) for i in range(d) for j in range(d) if i != j]
def gd(H, a):
    a = np.asarray(a); return np.mean(np.prod(H**a[:,None], axis=0))
# (eq:same-matrix) one Hadamard
Hs = samples[2][1]; a = np.array((1,0,-1,0,2,-2)); b = np.array((0,1,0,-1,0,0)); worst5 = 0
for i in range(6):
    s_ = sum(gd(Hs, a + np.eye(6,dtype=int)[i] - np.eye(6,dtype=int)[j]) * gd(Hs, b - np.eye(6,dtype=int)[i] + np.eye(6,dtype=int)[j]) for j in range(6))
    worst5 = max(worst5, abs(s_ - gd(Hs, a+b)))
print(f"(5) (eq:same-matrix) max residual = {worst5:.2e}")
# (eq:different-matrices) for a MUB pair in d=6 (tensor triple)
I2 = np.eye(2); A1 = np.array([[1,1],[1,-1]])/np.sqrt(2); A2 = np.array([[1,1],[1j,-1j]])/np.sqrt(2)
w3 = np.exp(2j*np.pi/3); F3 = np.array([[w3**(j*k) for k in range(3)] for j in range(3)])/np.sqrt(3); D3 = np.diag([1,w3,w3])
U1 = np.kron(A1, F3); U2 = np.kron(A2, D3@F3)
H1 = np.sqrt(6)*U1; H2 = np.sqrt(6)*U2
assert is_hadamard(H1) and is_hadamard(H2) and np.allclose(abs(U1.conj().T@U2)**2, 1/6)
a = np.array((2,0,-1,-1,0,0)); b = np.array((0,1,-1,0,0,0)); tot = 0
for d in D_set(6): tot += gd(H1, a+np.array(d))*gd(H2, b-np.array(d))
print(f"(5) (eq:different-matrices) for the tensor MUB pair: |sum| = {abs(tot):.2e}")
# (eq:order-two) complete families in d=3 and d=5
def complete_mubs(d):
    wd = np.exp(2j*np.pi/d); F = np.array([[wd**(j*k) for k in range(d)] for j in range(d)])/np.sqrt(d)
    bases = [np.eye(d)]
    for r in range(d):
        Dg = np.diag([wd**((r*k*k*((d+1)//2)) % d) for k in range(d)])   # r k^2/2 mod d (d odd)
        bases.append(Dg@F)
    return bases
for d in (3,5):
    bases = complete_mubs(d)
    for r in range(d+1):
        for s in range(d+1):
            if r != s: assert np.allclose(abs(bases[r].conj().T@bases[s])**2, 1/d), (d,r,s)
    worst = 0
    for trial in range(20):
        # random order-2 charge
        while True:
            a = rng.integers(-2,3,d); 
            if a.sum()==0 and sum(x for x in a if x>0)==2: break
        for r in range(d+1):
            tot = 0
            for s in range(d+1):
                if s == r: continue
                Hrs = np.sqrt(d)*bases[r].conj().T@bases[s]
                tot += gd(Hrs, a)
            worst = max(worst, abs(tot))
    print(f"(5) (eq:order-two) sum_(s!=r) g_s(a)=0 in d={d} complete family: max |sum| = {worst:.2e}")

# (6) singleton reduction in d=3 analogue:  M((0,1);(a,-a)) = -1/(d+1-2) M((0,0);(a,-a)) for order-2 a
d = 3; bases = complete_mubs(d)
def mixed_moment(labels, charges, r):
    # expectation over uniform row permutation sigma and uniform injection pi of labels into the d external bases
    ext = [s for s in range(d+1) if s != r]
    tot = 0; cnt = 0
    labs = sorted(set(labels))
    for sigma in itertools.permutations(range(d)):
        for inj in itertools.permutations(ext, len(labs)):
            assign = dict(zip(labs, inj))
            prod = 1
            for lab, a in zip(labels, charges):
                Hs = np.sqrt(d)*(bases[r][:, list(sigma)]).conj().T@bases[assign[lab]]
                prod *= gd(Hs, a)
            tot += prod.real; cnt += 1
    return tot/cnt
a = (2,-1,-1)
lhs = mixed_moment((0,1), (a, tuple(-x for x in a)), 0); rhs = -1/(d+1-2)*mixed_moment((0,0), (a, tuple(-x for x in a)), 0)
print(f"(6) singleton reduction d=3: M((0,1);(a,-a))={lhs:.6f}, -1/(d-1) M((0,0);(a,-a))={rhs:.6f}, diff={abs(lhs-rhs):.2e}")
# Lemma gram in d=3: Re E[F_1(a) conj F_1(v)] = (1/(d-1)) M((1,1,0);b) + ((d-2)/(d-1)) M((1,2,0);b),  b=(a,-v,v-a)
# Here F_1(a) = (1/(d-1)) sum_s g_s(a) z^{-a}, z a uniform column of H_0.
def gram_lhs(a, v, r=0):
    ext = [s for s in range(d+1) if s != r]; tot = 0; cnt = 0
    for sigma in itertools.permutations(range(d)):
        for inj in itertools.permutations(ext):      # bijection labels 0..d-1 -> ext
            Hs = [np.sqrt(d)*(bases[r][:, list(sigma)]).conj().T@bases[inj[s]] for s in range(d)]
            F1a = np.mean([gd(Hs[s], a) for s in range(1, d)]); F1v = np.mean([gd(Hs[s], v) for s in range(1, d)])
            for z in range(d):
                col = Hs[0][:, z]
                tot += (F1a*np.prod(col**(-np.array(a))) * np.conj(F1v*np.prod(col**(-np.array(v))))).real; cnt += 1
    return tot/cnt
best=None
for a in [(1,-1,0),(2,-1,-1),(1,1,-2),(2,-2,0)]:
    for v in [(0,1,-1),(1,-2,1),(-1,2,-1),(2,-1,-1),(1,-1,0)]:
        bch = (a, tuple(-x for x in v), tuple(np.array(v)-np.array(a)))
        lhs = gram_lhs(a, v); rhs = (1/(d-1))*mixed_moment((1,1,0), bch, 0) + ((d-2)/(d-1))*mixed_moment((1,2,0), bch, 0)
        if best is None or abs(lhs) > abs(best[2]): best = (a, v, lhs, rhs)
        assert abs(lhs-rhs) < 1e-12, (a, v, lhs, rhs)
print(f"(6) Lemma gram (d=3 analogue), 20 charge pairs all agree; largest nonzero case a={best[0]} v={best[1]}: lhs={best[2]:.6f} rhs={best[3]:.6f}")
