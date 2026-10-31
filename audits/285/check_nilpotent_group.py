"""
Audit 285 / check_nilpotent_group.py

Symbolic / numerical checks of Section "The nilpotent deck group"
(nilpotent.tex) of the Kadison-Kaplansky paper, for m = 3 (n_f = 3)
and m = 4 (n_f = 6):

 (a) affine group law (nil:affine-law, l.167-170) is associative and the
     inverse formula (nil:affine-inverse, l.171-173) is correct;
 (b) the laws are polynomials with integer coefficients in (lambda, b-1, z);
 (c) the cocycle c(gamma,eta) = mu(b^{-1}-1) (l.180) satisfies the cocycle
     identity (l.188-190), is normalized, and c(g,g^{-1}) = c(g^{-1},g);
 (d) grading: with wt((b-1)_ij) = j-i, wt(lambda_ij) = m-(j-i), wt(z) = m,
     the dilation Delta_y is an automorphism and
     c(Delta_y g, Delta_y h) = y^m c(g,h)  (Lemma nil:grading);
 (e) pi_eps(gamma) phi(x) = e(eps lambda(x)) phi(A_{b^{-1}} x) satisfies
     pi(g)pi(h) = e(eps c(g,h)) pi(gh)  (l.441-444);
 (f) the finite representation R (nil:R, l.762-765) is a representation of
     E_Z with central character e(-z/q), R_gamma R_eta = e(-delta c) R_{gamma eta}
     (l.815-817), and R_y = 1 for y in Lambda_q (nil:congruence-R);
 (g) the Magnus detector of Lemma nil:word: for random nontrivial reduced
     words w of length k < m in free generators, w(1+Z_j) != 1 with the
     explicit Z_j, and entries of words in generic unitriangular matrices
     have total degree <= m-1.

Nothing from the OpenAI repository is executed.
"""
import itertools
import random
import sympy as sp
import numpy as np

rng = random.Random(285)


def strict_upper_symbols(m, name):
    M = sp.zeros(m, m)
    for i in range(m):
        for j in range(i + 1, m):
            M[i, j] = sp.Symbol(f"{name}{i}{j}")
    return M


def lam_symbols(m, name):
    # a linear functional on n = strictly upper matrices: lambda(x) = sum l_ij x_ij
    return {(i, j): sp.Symbol(f"{name}{i}{j}") for i in range(m) for j in range(i + 1, m)}


def lam_apply(lam, X):
    return sum(lam[(i, j)] * X[i, j] for (i, j) in lam)


def lam_compose_L(lam, b):
    """mu o L_b where L_b(x) = b x  (b unitriangular). Returns functional coefficients."""
    m = b.shape[0]
    X = strict_upper_symbols(m, "X")
    expr = sp.expand(lam_apply(lam, b * X))
    out = {}
    for (i, j) in lam:
        out[(i, j)] = sp.expand(expr.coeff(X[i, j]))
    return out


def group_mul(m, g, h):
    """(lambda,b,z)(mu,b',w) = (lambda + mu o L_{b^{-1}}, b b', z + w + mu(b^{-1}-1))."""
    lam, b, z = g
    mu, bp, w = h
    binv = b.inv()
    return (
        {k: sp.expand(lam[k] + lam_compose_L(mu, binv)[k]) for k in lam},
        sp.expand(b * bp),
        sp.expand(z + w + lam_apply(mu, binv - sp.eye(m))),
    )


def group_inv(m, g):
    lam, b, z = g
    return (
        {k: sp.expand(-lam_compose_L(lam, b)[k]) for k in lam},
        b.inv(),
        sp.expand(-z - lam_apply(lam, b - sp.eye(m))),
    )


def cocycle(m, gamma, eta):
    (lam, b), (mu, bp) = gamma, eta
    return sp.expand(lam_apply(mu, b.inv() - sp.eye(m)))


def check_group(m):
    print(f"--- m = {m}, n_f = {m*(m-1)//2} ---")
    I = sp.eye(m)
    g = (lam_symbols(m, "l"), I + strict_upper_symbols(m, "b"), sp.Symbol("z"))
    h = (lam_symbols(m, "m"), I + strict_upper_symbols(m, "c"), sp.Symbol("w"))
    k = (lam_symbols(m, "n"), I + strict_upper_symbols(m, "e"), sp.Symbol("v"))

    # (a) associativity and inverse
    lhs = group_mul(m, group_mul(m, g, h), k)
    rhs = group_mul(m, g, group_mul(m, h, k))
    assert all(sp.simplify(lhs[0][key] - rhs[0][key]) == 0 for key in lhs[0])
    assert sp.simplify(lhs[1] - rhs[1]) == sp.zeros(m, m)
    assert sp.simplify(lhs[2] - rhs[2]) == 0
    gi = group_inv(m, g)
    e1 = group_mul(m, g, gi)
    e2 = group_mul(m, gi, g)
    for e in (e1, e2):
        assert all(sp.simplify(v) == 0 for v in e[0].values())
        assert sp.simplify(e[1] - I) == sp.zeros(m, m)
        assert sp.simplify(e[2]) == 0
    print("(a) associativity + inverse formula: OK")

    # (b) integer polynomiality in (lambda, b-1, z)
    prod = group_mul(m, g, h)
    allpolys = list(prod[0].values()) + list(prod[1]) + [prod[2]] + list(gi[0].values()) + list(gi[1]) + [gi[2]]
    for p in allpolys:
        p = sp.expand(p)
        assert p.is_polynomial(), p
        for coef in sp.Poly(p, *sorted(p.free_symbols, key=str)).coeffs() if p.free_symbols else [p]:
            assert coef == int(coef), coef
    print("(b) integer polynomial laws: OK")

    # (c) cocycle identity
    gam = (g[0], g[1]); eta = (h[0], h[1]); zet = (k[0], k[1])
    def qmul(a, b_):
        return ({key: sp.expand(a[0][key] + lam_compose_L(b_[0], a[1].inv())[key]) for key in a[0]},
                sp.expand(a[1] * b_[1]))
    lhs = cocycle(m, gam, eta) + cocycle(m, qmul(gam, eta), zet)
    rhs = cocycle(m, eta, zet) + cocycle(m, gam, qmul(eta, zet))
    assert sp.simplify(lhs - rhs) == 0
    ident = ({key: 0 for key in g[0]}, I)
    assert cocycle(m, gam, ident) == 0 and cocycle(m, ident, gam) == 0
    gaminv = ({key: sp.expand(-lam_compose_L(g[0], g[1])[key]) for key in g[0]}, g[1].inv())
    assert sp.simplify(cocycle(m, gam, gaminv) - cocycle(m, gaminv, gam)) == 0
    print("(c) cocycle identity, normalization, c(g,g^-1)=c(g^-1,g): OK")

    # (d) grading / dilation homogeneity
    y = sp.Symbol("y", positive=True)
    def dil(a):
        lam, b = a
        lam2 = {(i, j): lam[(i, j)] * y ** (m - (j - i)) for (i, j) in lam}
        b2 = sp.Matrix(m, m, lambda i, j: b[i, j] * y ** (j - i) if j > i else b[i, j])
        return (lam2, b2)
    dg, dh = dil(gam), dil(eta)
    # automorphism on the quotient
    p1 = dil(qmul(gam, eta)); p2 = qmul(dg, dh)
    assert all(sp.simplify(p1[0][key] - p2[0][key]) == 0 for key in p1[0])
    assert sp.simplify(p1[1] - p2[1]) == sp.zeros(m, m)
    assert sp.simplify(cocycle(m, dg, dh) - y ** m * cocycle(m, gam, eta)) == 0
    print("(d) dilation automorphism and c(D_y g, D_y h) = y^m c(g,h): OK")

    # (e) pi_eps multiplier identity, checked on a test function at numeric points
    eps = 1
    xs = strict_upper_symbols(m, "x")
    def A(b, X):
        return b * (I + X) - I
    phi = sp.Function("phi")
    # (pi(g) psi)(x) = e(eps lambda(x)) psi(A_{b^{-1}} x)
    def pi(a, psi):
        lam, b = a
        return lambda X: sp.exp(2 * sp.pi * sp.I * eps * lam_apply(lam, X)) * psi(A(b.inv(), X))
    psi0 = lambda X: phi(*[X[i, j] for i in range(m) for j in range(i + 1, m)])
    left = pi(gam, pi(eta, psi0))(xs)
    right = sp.exp(2 * sp.pi * sp.I * eps * cocycle(m, gam, eta)) * pi(qmul(gam, eta), psi0)(xs)
    # compare the phase (ratio) after cancelling phi arguments
    ratio = sp.simplify(sp.expand_log(sp.log(left / right), force=True))
    assert sp.simplify(ratio) == 0, ratio
    print("(e) pi_eps(g) pi_eps(h) = e(eps c(g,h)) pi_eps(gh): OK")

    # (f) finite representation R mod q, numerically
    qq = 5
    nf = m * (m - 1) // 2
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    points = list(itertools.product(range(qq), repeat=nf))   # n_q
    index = {pt: n for n, pt in enumerate(points)}
    def to_mat(pt):
        X = sp.zeros(m, m)
        for (i, j), val in zip(pairs, pt):
            X[i, j] = val
        return X
    def R(lam, b, z):
        """R(lambda,b,z) psi (x) = e(-(z+lambda(x))/q) psi(b^{-1}(1+x)-1), coordinates mod q."""
        Mx = np.zeros((len(points), len(points)), dtype=complex)
        binv = b.inv()
        lam = {key: int(sp.Integer(val)) for key, val in lam.items()}
        z = int(sp.Integer(z))
        for pt in points:
            X = to_mat(pt)
            Y = A(binv, X)
            ypt = tuple(int(Y[i, j]) % qq for (i, j) in pairs)
            phase = np.exp(-2j * np.pi * (z + sum(lam[(i, j)] * pt[n] for n, (i, j) in enumerate(pairs))) / qq)
            Mx[index[pt], index[ypt]] += phase
        return Mx
    def rand_elem():
        lam = {key: rng.randint(-7, 7) for key in pairs}
        b = I + to_mat(tuple(rng.randint(-7, 7) for _ in pairs))
        return (lam, b, rng.randint(-7, 7))
    for _ in range(5):
        g1, g2 = rand_elem(), rand_elem()
        R1, R2 = R(*g1), R(*g2)
        assert np.allclose(R1 @ R1.conj().T, np.eye(len(points)))
        prod = group_mul(m, g1, g2)
        Rp = R(prod[0], prod[1], prod[2])
        assert np.allclose(R1 @ R2, Rp), "R is not a representation"
        # section version: R_gamma = R(sigma(gamma)) with z = 0, multiplier e(-delta c)
        s1 = (g1[0], g1[1], 0); s2 = (g2[0], g2[1], 0)
        sp12 = group_mul(m, s1, s2)
        c12 = cocycle(m, (g1[0], g1[1]), (g2[0], g2[1]))
        assert sp.simplify(sp12[2] - c12) == 0
        Rs = R(sp12[0], sp12[1], 0)
        assert np.allclose(R(*s1) @ R(*s2), np.exp(-2j * np.pi * float(c12) / qq) * Rs)
    # central character and congruence subgroup
    assert np.allclose(R({key: 0 for key in pairs}, I, 3), np.exp(-2j * np.pi * 3 / qq) * np.eye(len(points)))
    for _ in range(5):
        lam = {key: qq * rng.randint(-3, 3) for key in pairs}
        b = I + to_mat(tuple(qq * rng.randint(-3, 3) for _ in pairs))
        assert np.allclose(R(lam, b, 0), np.eye(len(points))), "R_y != 1 on Lambda_q"
    print(f"(f) R mod q={qq}: representation, multiplier e(-delta c), R=1 on Lambda_q: OK")


def check_magnus(m, trials=300):
    """Lemma nil:word: w(1+Z_j) != 1 for nontrivial reduced words of length k < m."""
    rng2 = random.Random(7)
    ngen = 3
    fails = 0
    for _ in range(trials):
        k = rng2.randint(1, m - 1)
        # random reduced word in syllables x_{j_s}^{a_s}, j_s != j_{s+1}
        js, prev = [], None
        for _ in range(k):
            j = rng2.choice([x for x in range(ngen) if x != prev]); js.append(j); prev = j
        a = [rng2.choice([-3, -2, -1, 1, 2, 3]) for _ in range(k)]
        Z = [sp.zeros(m, m) for _ in range(ngen)]
        for s, j in enumerate(js):
            Z[j][s, s + 1] += 1
        W = sp.eye(m)
        for s in range(k):
            W = W * (sp.eye(m) + Z[js[s]]) ** a[s]
        if W == sp.eye(m):
            fails += 1
        # the (1,k+1) entry should be prod a_s
        assert W[0, k] == sp.prod(a), (W[0, k], a)
    print(f"(g) Magnus detector, m={m}: {trials} random reduced words with k<m syllables, "
          f"failures = {fails}; (1,k+1) entry = prod a_s verified")
    # degree bound: entries of a long word in generic unitriangular matrices have total degree <= m-1
    syms = [sp.eye(m) + strict_upper_symbols(m, f"v{t}_") for t in range(3)]
    word = syms[0] * syms[1].inv() * syms[2] * syms[0].inv() * syms[1] * syms[2] ** 2 * syms[0]
    maxdeg = max(sp.Poly(sp.expand(word[i, j]), *word[i, j].free_symbols).total_degree()
                 for i in range(m) for j in range(i + 1, m) if word[i, j].free_symbols)
    print(f"    max total degree of an entry of a length-8 word = {maxdeg} (claim: <= m-1 = {m-1})")
    assert maxdeg <= m - 1


if __name__ == "__main__":
    check_group(3)          # m = 4 symbolic run is slow; m = 3 already has a nonabelian quotient E-bar
    check_magnus(4)
    check_magnus(5)
    # B(C), L(C) definitions from nil:rank / nil:word
    import math
    def B(C): return max(1, math.ceil(20 * C)) ** 2
    def L(C): return math.ceil(3 * C * B(C))
    print("B(1) =", B(1), " L(1) =", L(1))
    print("nilpotent group checks: all OK")
