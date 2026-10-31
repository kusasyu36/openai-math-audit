"""
Audit 285 / check_model_band.py

Finite numerical models for two mechanisms of the Kadison-Kaplansky paper.

Part A -- Lemma spec:model (sampling.tex l.65-135).
  Finite tower: layers i = 0..3, sides A_i (size |A|) and B_i (size |B|),
  random (d, D')-biregular incidences between A_i and B_k, |i-k| <= 1,
  normalized by (d D')^{-1/2}; window weights h_i with sum h_i^2 = 1 and
  two adjacent active layers; transverse projection p of rank rho on C^n.
  Model O = sum_{ik} U_ik (x) h_i h_k p  + adjoints.
  Claims checked: eigenvalue +1 with multiplicity exactly rho, eigenvalue -1
  with multiplicity rho, all other eigenvalues of modulus <= 1-chi where
  1-chi = max second singular value of the U_ik; the diagonal of the upper
  projection at a vertex of A_i is h_i^2 p /(2|A_i|) (so the sum over all
  vertices of the "identity diagonal" is tr p, not |A_i| tr p).

Part B -- the sampling mechanism of Theorem sam:realization (one voltage and
  one weight per edge approximating "incidence (x) convolution projection").
  Deck group F = Z/n (finite, so the regular representation is a direct sum of
  characters = Fourier fibers).  K = sum_g k(g) lambda(g) is the convolution
  projection onto r of the n characters.  Each edge e: a -> b gets an
  independent voltage z_e ~ rho (rho proportional to |k|) and weight
  w_e = k(z_e)/rho(z_e), so that E[w_e lambda(z_e)] = K.  The cover operator
  T = (dD')^{-1/2} sum_e w_e |a><b| (x) lambda(z_e) + adj is compared with the
  model U (x) K: we report ||T - U (x) K|| for growing d and check the band
  near +1 (count of eigenvalues in [1-c, 1+c] equals r * (number of active
  constant modes)).

Everything is self-written; nothing from the OpenAI repository is executed.
"""
import numpy as np

rng = np.random.default_rng(285)


def random_biregular(nA, nB, d, Dp):
    """Random (d, Dp)-biregular bipartite 0/1 incidence (nA x nB), nA*d == nB*Dp.
    Built from d random perfect matchings of A-slots onto B-slots; parallel edges
    are kept as multiplicities (the paper also counts slots separately)."""
    assert nA * d == nB * Dp
    M = np.zeros((nA, nB))
    perm = rng.permutation(np.repeat(np.arange(nB), Dp))   # nB*Dp = nA*d slots
    for a in range(nA):
        for b in perm[a * d:(a + 1) * d]:
            M[a, b] += 1
    assert np.all(M.sum(axis=1) == d) and np.all(M.sum(axis=0) == Dp)
    return M


def part_a():
    print("=== Part A: model band (Lemma spec:model) ===")
    nA, nB, d, Dp = 64, 32, 8, 16
    n, rho = 6, 2
    layers = 4
    U = {}
    second = 0.0
    for i in range(layers):
        for k in range(layers):
            if abs(i - k) <= 1:
                M = random_biregular(nA, nB, d, Dp) / np.sqrt(d * Dp)
                U[(i, k)] = M
                sv = np.linalg.svd(M, compute_uv=False)
                assert abs(sv[0] - 1) < 1e-9
                second = max(second, sv[1])
    chi = 1 - second
    print(f"d={d}, D'={Dp}, |A|={nA}, |B|={nB}: max second singular value = {second:.4f} (chi = {chi:.4f})")
    # window weights: active layers 1,2
    h = np.zeros(layers); h[1], h[2] = np.cos(0.4), np.sin(0.4)
    assert abs(np.sum(h ** 2) - 1) < 1e-12
    # transverse projection of rank rho on C^n
    Q, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
    p = Q[:, :rho] @ Q[:, :rho].conj().T
    # assemble O on (A_0..A_3, B_0..B_3) (x) C^n
    dimA, dimB = layers * nA, layers * nB
    O = np.zeros(((dimA + dimB) * n, (dimA + dimB) * n), dtype=complex)
    for (i, k), M in U.items():
        blk = np.kron(M, h[i] * h[k] * p)
        r0, c0 = i * nA * n, (dimA + k * nB) * n
        O[r0:r0 + nA * n, c0:c0 + nB * n] = blk
        O[c0:c0 + nB * n, r0:r0 + nA * n] = blk.conj().T
    ev, evec = np.linalg.eigh(O)
    n_plus = np.sum(np.abs(ev - 1) < 1e-8)
    n_minus = np.sum(np.abs(ev + 1) < 1e-8)
    rest = np.max(np.abs(ev[(np.abs(ev - 1) > 1e-8) & (np.abs(ev + 1) > 1e-8)]))
    print(f"multiplicity of +1: {n_plus} (claim rho = {rho}); of -1: {n_minus}; "
          f"max |other eigenvalue| = {rest:.4f} (claim <= 1-chi = {1-chi:.4f})")
    assert n_plus == rho and n_minus == rho and rest <= 1 - chi + 1e-9
    # upper projection and its diagonal
    E = evec[:, np.abs(ev - 1) < 1e-8] @ evec[:, np.abs(ev - 1) < 1e-8].conj().T
    # diagonal block at a vertex of A_1: expected h_1^2 p / (2 |A|)
    i = 1; a = 5
    r0 = (i * nA + a) * n
    blk = E[r0:r0 + n, r0:r0 + n]
    err = np.max(np.abs(blk - h[i] ** 2 * p / (2 * nA)))
    print(f"diagonal block of E at a vertex of A_1 vs h_1^2 p/(2|A|): max error {err:.2e}")
    assert err < 1e-9
    tot = sum(np.trace(E[(v) * n:(v + 1) * n, v * n:(v + 1) * n]).real for v in range(dimA + dimB))
    print(f"sum over all vertices of Tr_n E(v,v) = {tot:.6f} (= tr p = rho = {rho})")
    assert abs(tot - rho) < 1e-9
    print("Part A OK\n")


def part_b():
    print("=== Part B: one voltage + one weight per edge approximating incidence (x) convolution projection ===")
    nF = 24              # deck group Z/nF
    r = 5                # rank of the convolution projection (5 of 24 characters)
    chars = np.arange(r)  # selected characters
    g = np.arange(nF)
    # k(g) = (1/nF) sum_{c in chars} e(c g / nF): convolution projection coefficients
    k = np.array([np.sum(np.exp(2j * np.pi * chars * gg / nF)) / nF for gg in g])
    # regular representation of Z/nF: lambda(g) = shift matrix
    def lam(gg):
        return np.roll(np.eye(nF), gg, axis=0)
    Kop = sum(k[gg] * lam(gg) for gg in g)
    assert np.allclose(Kop @ Kop, Kop) and np.allclose(Kop, Kop.conj().T)
    assert abs(np.trace(Kop).real - r) < 1e-9
    rho_dist = np.abs(k) / np.sum(np.abs(k))        # importance sampling law
    print(f"deck group Z/{nF}, projection K of rank {r}, ||k||_1 = {np.sum(np.abs(k)):.3f}")

    nB_base = 16
    for d in [16, 64, 256, 1024]:
        Dp = 2 * d
        nB = nB_base * 2
        nA = nB * Dp // d                              # nA*d == nB*Dp  -> nA = 2 nB
        M = random_biregular(nA, nB, d, Dp)
        # build T with independent (voltage, weight) per edge slot
        T_blk = np.zeros((nA * nF, nB * nF), dtype=complex)
        for a in range(nA):
            for b in range(nB):
                mult = int(M[a, b])
                for _ in range(mult):
                    z = rng.choice(nF, p=rho_dist)
                    w = k[z] / rho_dist[z]
                    T_blk[a * nF:(a + 1) * nF, b * nF:(b + 1) * nF] += w * lam(z)
        T_blk /= np.sqrt(d * Dp)
        O_blk = np.kron(M / np.sqrt(d * Dp), Kop)
        err = np.linalg.norm(T_blk - O_blk, 2)
        sv2 = np.linalg.svd(M / np.sqrt(d * Dp), compute_uv=False)[1]
        # self-adjoint doubling and the band near +-1
        dim = (nA + nB) * nF
        T = np.zeros((dim, dim), dtype=complex)
        T[:nA * nF, nA * nF:] = T_blk
        T[nA * nF:, :nA * nF] = T_blk.conj().T
        ev = np.linalg.eigvalsh(T)
        band = np.sum(ev > 0.75)
        below = ev[ev <= 0.75]
        rest = np.max(below) if below.size else 0.0
        print(f"d={d:4d} D'={Dp:4d} |A|={nA:3d} |B|={nB:3d}: ||T - U(x)K|| = {err:.4f}  (sqrt(d)*err = {np.sqrt(d)*err:.2f}),"
              f" second s.v. of U = {sv2:.3f}; #eigenvalues > 0.75: {band} (model: r = {r});"
              f" largest eigenvalue <= 0.75: {rest:.3f}")
    print("Part B: the error decays like C/sqrt(d); once it is below the model gap, the upper band has exactly the model multiplicity r\n")


if __name__ == "__main__":
    part_a()
    part_b()
