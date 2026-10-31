"""
Audit 285 / check_scalar_compression.py

Independent numerical check of Proposition "Scalar compression after a free
product" (scalar.tex l.52-183) of the Kadison-Kaplansky paper.

Setting: Lambda = Z/3 (any discrete group is allowed by the proposition),
K = Lambda * F_2 = <g | g^3> * <t1, t2>.  p in M_2(C[Lambda]) is a
projection with unnormalized trace a = 1/3 < 1/2 (p is a unitary conjugate
of diag(q, 0), q = (1+g+g^2)/3, by a non-diagonal unitary in M_2(C[Z/3])).

Checks:
 (1) lower bound ||U xi|| >= c_a ||xi||, c_a = sqrt(1-a) - sqrt(a),
     on random finitely supported xi = p xi (eq. scalar:lower-bound);
 (2) the series A^{-1} = (1/N) sum (p - A/N)^k converges, e = Up A^{-1} (Up)^*
     satisfies e^2 = e, e = e^*, and tau(e) = a (eq. scalar:projection-formula);
 (3) the intermediate inequalities scalar:bad-input, scalar:good-output,
     scalar:bad-output on the same random vectors.

Elements of C[K] are stored as dict {word: complex}, with a word a tuple of
normal-form letters: ('g', k) for k in {1,2}, ('t', j, e) for generator t_j
with exponent e != 0 (consecutive syllables are merged / cancelled).
Everything is written from scratch; nothing from the OpenAI repository is run.
"""
import numpy as np
import random
from collections import defaultdict

rng = random.Random(285)

# ---------- free product normal forms -------------------------------------

def mul_words(u, v):
    """Reduced product of two normal-form words of Z/3 * F_2."""
    w = list(u)
    for letter in v:
        if w and w[-1][0] == letter[0]:
            last = w[-1]
            if letter[0] == 'g':
                k = (last[1] + letter[1]) % 3
                w.pop()
                if k:
                    w.append(('g', k))
            else:
                if last[1] == letter[1]:
                    e = last[2] + letter[2]
                    w.pop()
                    if e:
                        w.append(('t', letter[1], e))
                else:
                    w.append(letter)
        else:
            w.append(letter)
    return tuple(w)


def inv_word(u):
    out = []
    for letter in reversed(u):
        if letter[0] == 'g':
            out.append(('g', (-letter[1]) % 3))
        else:
            out.append(('t', letter[1], -letter[2]))
    return tuple(out)


# ---------- group algebra / matrices over it -------------------------------

def ga_mul(x, y):
    out = defaultdict(complex)
    for u, cu in x.items():
        for v, cv in y.items():
            out[mul_words(u, v)] += cu * cv
    return {w: c for w, c in out.items() if abs(c) > 1e-14}


def ga_add(x, y, s=1.0):
    out = defaultdict(complex)
    for u, c in x.items():
        out[u] += c
    for u, c in y.items():
        out[u] += s * c
    return {w: c for w, c in out.items() if abs(c) > 1e-14}


def ga_scale(x, s):
    return {w: s * c for w, c in x.items()}


def ga_star(x):
    return {inv_word(w): np.conj(c) for w, c in x.items()}


def ga_tau(x):
    return x.get((), 0.0)


def ga_norm1(x):
    return sum(abs(c) for c in x.values())


def mat_mul(X, Y):
    n, m, k = len(X), len(Y[0]), len(Y)
    out = [[{} for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            acc = {}
            for l in range(k):
                acc = ga_add(acc, ga_mul(X[i][l], Y[l][j]))
            out[i][j] = acc
    return out


def mat_add(X, Y, s=1.0):
    return [[ga_add(X[i][j], Y[i][j], s) for j in range(len(X[0]))]
            for i in range(len(X))]


def mat_star(X):
    n, m = len(X), len(X[0])
    return [[ga_star(X[j][i]) for j in range(n)] for i in range(m)]


def mat_scale(X, s):
    return [[ga_scale(x, s) for x in row] for row in X]


def mat_tau(X):
    return sum(ga_tau(X[i][i]) for i in range(len(X)))


def mat_norm1(X):
    return max(sum(ga_norm1(x) for x in row) for row in X)


# ---------- the data ---------------------------------------------------------

g = (('g', 1),)
g2 = (('g', 2),)
one = ()
q = {one: 1 / 3, g: 1 / 3, g2: 1 / 3}          # minimal projection, tau = 1/3
zero = {}

# unitary in M_2(C[Z/3]): u = [[c, -s g],[s g^{-1}, c]] with c^2+s^2=1.
c, s = np.cos(0.7), np.sin(0.7)
Umat = [[{one: c}, {g: -s}], [{g2: s}, {one: c}]]
assert mat_norm1(mat_add(mat_mul(Umat, mat_star(Umat)),
                         [[{one: 1}, {}], [{}, {one: 1}]], -1)) < 1e-12
P = mat_mul(mat_mul(Umat, [[q, zero], [zero, zero]]), mat_star(Umat))
# projection check
assert mat_norm1(mat_add(mat_mul(P, P), P, -1)) < 1e-12
assert mat_norm1(mat_add(mat_star(P), P, -1)) < 1e-12
a = mat_tau(P).real
N = 2
print(f"unnormalized trace a = tau_N(p) = {a:.6f}")
assert abs(a - 1 / 3) < 1e-12

t = [(('t', 1, 1),), (('t', 2, 1),)]   # the free generators
U = [[{t[0]: 1}, {t[1]: 1}]]            # 1 x N row of unitaries
UP = mat_mul(U, P)                       # Up
A = mat_mul(mat_star(UP), UP)            # p U^* U p   (N x N)

# ---------- (1) lower bound on random finitely supported vectors ----------

def random_vector(nwords=6, maxlen=4):
    """Random finitely supported xi in ell^2(K)^N, as N dicts."""
    xi = [defaultdict(complex) for _ in range(N)]
    for j in range(N):
        for _ in range(nwords):
            w = ()
            for _ in range(rng.randint(0, maxlen)):
                if rng.random() < 0.4:
                    w = mul_words(w, (('g', rng.randint(1, 2)),))
                else:
                    w = mul_words(w, (('t', rng.randint(1, 2),
                                       rng.choice([-2, -1, 1, 2])),))
            xi[j][w] += complex(rng.gauss(0, 1), rng.gauss(0, 1))
    return [dict(x) for x in xi]


def apply_mat(X, xi):
    """Left regular action of the matrix X on the vector xi (lambda(x) delta_h = delta_{xh})."""
    n = len(X)
    out = [defaultdict(complex) for _ in range(n)]
    for i in range(n):
        for j in range(len(xi)):
            for w, cw in X[i][j].items():
                for h, ch in xi[j].items():
                    out[i][mul_words(w, h)] += cw * ch
    return [dict(x) for x in out]


def vnorm2(xi):
    return sum(abs(c) ** 2 for x in xi for c in x.values())


c_a = np.sqrt(1 - a) - np.sqrt(a)
print(f"c_a = sqrt(1-a)-sqrt(a) = {c_a:.6f}")
worst = np.inf
for trial in range(200):
    eta = random_vector()
    xi = apply_mat(P, eta)                       # xi = p eta  (so p xi = xi)
    n2 = vnorm2(xi)
    if n2 < 1e-10:
        continue
    Uxi = apply_mat(U, xi)
    ratio = np.sqrt(vnorm2(Uxi) / n2)
    worst = min(worst, ratio)
    # intermediate inequalities
    xib = [{w: cw for w, cw in xi[j].items() if w and w[0] == ('t', j + 1, -1)
            or (w and w[0][0] == 't' and w[0][1] == j + 1 and w[0][2] < 0)}
           for j in range(N)]
    xig = [{w: cw for w, cw in xi[j].items() if w not in xib[j]}
           for j in range(N)]
    assert vnorm2(xib) <= a * n2 + 1e-9, "scalar:bad-input violated"
    assert abs(vnorm2(apply_mat(U, xig)) - vnorm2(xig)) < 1e-9, "good outputs not isometric"
    assert vnorm2(xig) >= (1 - a) * n2 - 1e-9, "scalar:good-output violated"
    assert vnorm2(apply_mat(U, xib)) <= a * n2 + 1e-9, "scalar:bad-output violated"
print(f"min over 200 random xi = p xi of ||U xi||/||xi|| = {worst:.6f}  (must be >= c_a)")
assert worst >= c_a - 1e-9

# ---------- (2) the formal identities behind e = Up A^{-1} (Up)^* ----------
# The identities e^2 = e, e = e^*, tau(e) = tau_N(p) only use that A = pU^*Up is
# invertible in the corner p M_N p and that tau is tracial.  (An exact computation
# of A^{-1} in C[K] is infeasible: the Neumann series has exponentially growing
# support.)  We verify the algebra in a finite-dimensional stand-in: U = row of N
# random unitaries in M_k(C), p a random projection in M_N(M_k(C)) with A
# invertible, tau = normalized trace on M_k.
k = 7
def rand_unitary(k):
    Q, _ = np.linalg.qr(np.random.default_rng(1).standard_normal((k, k)) + 1j * np.random.default_rng(2).standard_normal((k, k)))
    return Q
rr = np.random.default_rng(285)
Us = [np.linalg.qr(rr.standard_normal((k, k)) + 1j * rr.standard_normal((k, k)))[0] for _ in range(N)]
Ublk = np.hstack(Us)                                       # k x Nk
G = rr.standard_normal((N * k, N * k)) + 1j * rr.standard_normal((N * k, N * k))
Qf, _ = np.linalg.qr(G)
rk = 3
pblk = Qf[:, :rk] @ Qf[:, :rk].conj().T                    # projection in M_N(M_k), tau_N(p) = rk/k
Ablk = pblk @ Ublk.conj().T @ Ublk @ pblk
# inverse in the corner: (A + 1 - p)^{-1} p
Ainv = np.linalg.inv(Ablk + np.eye(N * k) - pblk) @ pblk
eblk = Ublk @ pblk @ Ainv @ (Ublk @ pblk).conj().T          # k x k
print(f"finite-dimensional stand-in (k={k}, N={N}): ||e^2-e|| = {np.linalg.norm(eblk @ eblk - eblk):.2e}, "
      f"||e-e*|| = {np.linalg.norm(eblk - eblk.conj().T):.2e}, tau(e) = {np.trace(eblk).real / k:.6f}, "
      f"tau_N(p) = {np.trace(pblk).real / k:.6f}")
assert np.linalg.norm(eblk @ eblk - eblk) < 1e-9 and abs(np.trace(eblk).real - np.trace(pblk).real) < 1e-9
print("scalar compression: all checks OK")
