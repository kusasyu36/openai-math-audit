"""Independent exact diagonalization of the periodic spin-1 Heisenberg ring.

Checks for audit 268:
  * E_0(n), gap gamma_n (next distinct level) for n = 4..NMAX, untwisted.
  * Shifted partition functions Z_n(b,g) = Tr exp[-b(H_n(g) + a n)] for
    g = 1, P (theta = pi), C (theta = 2 pi/3) at b = 21/2 and 49/4, compared
    with the table in initialization.tex (Lemma thermal-inputs).
  * Sector moments I_n, O_n, N_n (eq. sector-inversion) and their signs.
  * Hankel / log-convexity consistency of the moment sequence.
Everything is built from Kronecker products of the spin-1 matrices; nothing is
transcribed from the paper's matrix descriptions.
"""
import sys, time
import numpy as np
import scipy.sparse as sp
from fractions import Fraction

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
a = Fraction(700741, 500000)
af = float(a)

# spin-1 matrices in weight basis (m = +1, 0, -1)
Sz = np.diag([1.0, 0.0, -1.0])
Sp = np.zeros((3, 3)); Sp[0, 1] = np.sqrt(2); Sp[1, 2] = np.sqrt(2)
Sm = Sp.T
I3 = np.eye(3)

def site_op(op, j, n):
    mats = [sp.identity(3, format='csr')] * n
    mats[j] = sp.csr_matrix(op)
    out = mats[0]
    for m in mats[1:]:
        out = sp.kron(out, m, format='csr')
    return out

def hamiltonian(n, theta):
    """H_n(g) with g = exp(-i theta S^z) on the seam bond (n-1,0)."""
    H = sp.csr_matrix((3**n, 3**n), dtype=complex)
    SZ = [site_op(Sz, j, n) for j in range(n)]
    SPp = [site_op(Sp, j, n) for j in range(n)]
    SMm = [site_op(Sm, j, n) for j in range(n)]
    for j in range(n - 1):
        H = H + SZ[j] @ SZ[j + 1] + 0.5 * (SPp[j] @ SMm[j + 1] + SMm[j] @ SPp[j + 1])
    # seam: sum_alpha S^alpha_{n-1} (g S^alpha g^{-1})_0, g S^pm g^{-1} = e^{∓ i theta} S^pm
    ph = np.exp(-1j * theta)
    H = H + SZ[n - 1] @ SZ[0] + 0.5 * (ph * SPp[n - 1] @ SMm[0] + np.conj(ph) * SMm[n - 1] @ SPp[0])
    return H.tocsr(), SZ

def spectrum(n, theta):
    H, SZ = hamiltonian(n, theta)
    sztot = sum(SZ).diagonal().real
    eigs = []
    for m in range(0, n + 1):
        idx = np.where(np.abs(sztot - m) < 1e-9)[0]
        if len(idx) == 0:
            continue
        Hs = H[idx][:, idx].toarray()
        if abs(theta - np.pi) < 1e-12 or theta == 0:
            Hs = Hs.real
        w = np.linalg.eigvalsh(Hs)
        eigs.append(w)
        if m > 0:
            eigs.append(w)  # sector -m has the same spectrum
    return np.sort(np.concatenate(eigs))

def Z(ev, n, b):
    return float(np.sum(np.exp(-b * (ev + af * n))))

table_105 = {6: (8945303, 346734, 939987), 8: (3741281, 567167, 984108), 10: (2327547, 721090, 1002613)}
table_1225 = {4: (124885379, 65757), 5: (2925, 2960742), 6: (12870614, 286599), 7: (54069, 1920414),
              8: (4639283, 510140), 9: (195455, 1514846), 10: (2654823, 675419), 11: (377629, 1314544),
              12: (1900547, 787405)}

results = {}
for n in range(4, NMAX + 1):
    t0 = time.time()
    ev = {}
    for name, th in (('1', 0.0), ('P', np.pi), ('C', 2 * np.pi / 3)):
        ev[name] = spectrum(n, th)
    e = ev['1']
    E0 = e[0]
    distinct = e[e > E0 + 1e-7]
    E1 = distinct[0]
    mult0 = int(np.sum(e < E0 + 1e-7))
    print(f"n={n:2d}  dim={len(e)}  E0={E0:.10f}  E0/n={E0/n:.10f}  E1-E0={E1-E0:.10f}  "
          f"ground multiplicity={mult0}  E0+an={E0+af*n:+.6e}  ({time.time()-t0:.1f}s)")
    # lower bound check -3n/2 <= H <= n
    assert e[0] >= -1.5 * n - 1e-9 and e[-1] <= n + 1e-9
    for name in 'PC':
        assert ev[name][0] >= -1.5 * n - 1e-9 and ev[name][-1] <= n + 1e-9
    for b, tab, tol in ((10.5, table_105, 1e-5), (12.25, table_1225, 30e-6)):
        zs = {name: Z(ev[name], n, b) for name in ev}
        results[(n, b)] = zs
        line = f"   b={b:5.2f}: Z1={zs['1']:.9f} ZP={zs['P']:.9f} ZC={zs['C']:.9f}"
        if n in tab:
            c = tab[n]
            d1 = zs['1'] - c[0] * 1e-6; dP = zs['P'] - c[1] * 1e-6
            ok = abs(d1) < tol and abs(dP) < tol
            line += f" | table diff: Z1 {d1:+.2e}, ZP {dP:+.2e}"
            if len(c) == 3:
                dC = zs['C'] - c[2] * 1e-6; ok = ok and abs(dC) < tol
                line += f", ZC {dC:+.2e}"
            line += "  " + ("OK" if ok else "MISMATCH")
        print(line)
        if n % 2 == 0:
            z1, zP, zC = zs['1'], zs['P'], zs['C']
            J = (z1 + 3 * zP) / 4; N = (z1 - zP) / 4; I = (z1 + 3 * zP + 8 * zC) / 12; O = (z1 + 3 * zP - 4 * zC) / 12
            print(f"            sectors: I={I:.7f} O={O:.7f} N={N:.7f} J={J:.7f}  all>=0: {min(I,O,N)>=0}")

# moment consistency at b = 49/4 (untwisted and P-twisted)
print("\nMoment-sequence consistency at b=49/4 (own values, n=4..NMAX):")
for name in ('1',):  # only the untwisted trace is a positive moment sequence (Z_n(P) = I+2O-N)
    zs = [results[(n, 12.25)][name] for n in range(4, NMAX + 1)]
    ns = list(range(4, NMAX + 1))
    # Cauchy-Schwarz Z_{n+m}^2 <= Z_{2n} Z_{2m} for the real eigenvalue list
    for i in range(len(ns)):
        for j in range(i + 1, len(ns)):
            if ns[i] % 2 == 0 and ns[j] % 2 == 0:   # Cauchy-Schwarz Z_{(i+j)/2}^2 <= Z_i Z_j needs even i, j (eigenvalues may be negative)
                k = (ns[i] + ns[j]) // 2
                if k in ns:
                    zk = zs[ns.index(k)]
                    ok = zk**2 <= zs[i] * zs[j] * (1 + 1e-12)
                    if not ok:
                        print(f"  g={name}: Z_{k}^2 <= Z_{ns[i]} Z_{ns[j]} FAILS")
    # Hankel PSD: moments m_4..m_{NMAX}; matrix [m_{i+j}] for i,j in 2..NMAX/2
    kmax = NMAX // 2
    Hk = np.array([[zs[ns.index(i + j)] for j in range(2, kmax + 1)] for i in range(2, kmax + 1)])
    w = np.linalg.eigvalsh(Hk)
    print(f"  g={name}: Hankel [Z_(i+j)] i,j=2..{kmax} eigenvalues min={w.min():.3e} (PSD expected since Z_n = sum lambda^n)")
# same with the paper's table (n = 4..12) 
print("Hankel test on the paper's b=49/4 table (n=4..12):")
for col, name in ((0, '1'),):
    m = {n: table_1225[n][col] * 1e-6 for n in range(4, 13)}
    Hk = np.array([[m[i + j] for j in range(2, 7)] for i in range(2, 7)])
    w = np.linalg.eigvalsh(Hk)
    print(f"  g={name}: eigenvalues of 5x5 Hankel: {w}")
    # J and N lists: J = (Z + 3 Z_P)/4, N = (Z - Z_P)/4 must each be moment sequences
mJ = {n: (table_1225[n][0] + 3 * table_1225[n][1]) / 4 * 1e-6 for n in range(4, 13)}
mN = {n: (table_1225[n][0] - table_1225[n][1]) / 4 * 1e-6 for n in range(4, 13)}
for lab, m in (('J', mJ), ('N', mN)):
    Hk = np.array([[m[i + j] for j in range(2, 7)] for i in range(2, 7)])
    print(f"  sector {lab}: Hankel eigenvalues: {np.linalg.eigvalsh(Hk)}")
print("done")
