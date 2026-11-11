"""Exact diagonalization of the periodic spin-1 ring using total-S^z and
(twisted) translation symmetry, to reach n = 11, 12 (audit 268).

Z_n(b,g) = Tr exp[-b(H_n(g) + a n)] for g = 1 (theta=0) and P (theta=pi) at
b = 49/4, compared with the table at initialization.tex:40-52.  For n <= 8 the
block spectra are compared with a dense diagonalization as a regression test.
Usage: python -I check_268_ed12.py n [theta_list]
"""
import sys, time
import numpy as np
from fractions import Fraction

n = int(sys.argv[1])
thetas = [float(eval(t, {"pi": np.pi})) for t in sys.argv[2].split(',')] if len(sys.argv) > 2 else [0.0, np.pi]
a = float(Fraction(700741, 500000))
table_1225 = {4: (124885379, 65757), 5: (2925, 2960742), 6: (12870614, 286599), 7: (54069, 1920414),
              8: (4639283, 510140), 9: (195455, 1514846), 10: (2654823, 675419), 11: (377629, 1314544),
              12: (1900547, 787405)}
table_105 = {6: (8945303, 346734, 939987), 8: (3741281, 567167, 984108), 10: (2327547, 721090, 1002613)}
pow3 = 3 ** np.arange(n)[::-1]

def encode(words):  # words: (m, n) array of letters in {-1,0,1}
    return (words + 1) @ pow3

def all_words_sector(M):
    import itertools
    # generate words with sum M by iterating over all 3^n (n <= 12 -> 531441, fine)
    allw = np.array(np.unravel_index(np.arange(3**n), (3,)*n)).T - 1
    return allw[allw.sum(1) == M]

def sector_spectrum_blocks(M, theta):
    words = all_words_sector(M)
    m = len(words)
    codes = encode(words)
    pos = {int(c): i for i, c in enumerate(codes)}
    # symmetry V|x> = exp(+i theta x_{n-1}) |s(x)>, s = right rotation
    # compute orbit data: rep index, shift r, phase phi with V^r|rep> = phi |x>
    rep = -np.ones(m, dtype=int); shift = np.zeros(m, dtype=int); phase = np.zeros(m, dtype=complex); period = np.zeros(m, dtype=int)
    for i in range(m):
        if rep[i] >= 0: continue
        x = words[i].copy(); ph = 1.0 + 0j; r = 0
        orbit = []
        while True:
            j = pos[int(encode(x))]
            if r > 0 and j == i: break
            orbit.append((j, r, ph))
            ph = ph * np.exp(1j * theta * x[n-1])  # V|x> = exp(+i theta x_{n-1}) |s(x)> commutes with H_n(g) in our seam convention
            x = np.roll(x, 1)
            r += 1
        p = r
        for (j, rr, pph) in orbit:
            rep[j] = i; shift[j] = rr; phase[j] = pph; period[j] = p
    reps = np.where(rep == np.arange(m))[0]
    # H action: list of (z, y, h) nonzeros, H|y> = sum_z h |z>
    diag = np.einsum('ij,ij->i', words, np.roll(words, 1, axis=1)).astype(complex)  # sum_j x_{j-1} x_j
    rows, cols, vals = [], [], []
    for j in range(n):
        k = (j + 1) % n
        for (dj, dk, ph) in ((1, -1, np.exp(-1j*theta) if j == n-1 else 1.0), (-1, 1, np.exp(1j*theta) if j == n-1 else 1.0)):
            ok = (words[:, j] + dj <= 1) & (words[:, j] + dj >= -1) & (words[:, k] + dk <= 1) & (words[:, k] + dk >= -1)
            src = np.where(ok)[0]
            z = words[src].copy(); z[:, j] += dj; z[:, k] += dk
            zc = encode(z)
            zi = np.array([pos[int(c)] for c in zc], dtype=int)
            rows.append(zi); cols.append(src); vals.append(np.full(len(src), ph, dtype=complex))
    rows = np.concatenate(rows).astype(int); cols = np.concatenate(cols).astype(int); vals = np.concatenate(vals).astype(complex)
    # keep only columns y that are representatives (we need <x;K|H|y;K> for y in reps)
    isrep = (rep == np.arange(m))
    sel = isrep[cols]
    rows, cols, vals = rows[sel], cols[sel], vals[sel]
    repidx = -np.ones(m, dtype=int); repidx[reps] = np.arange(len(reps))
    eigs = []
    for K in range(n):
        lamK = np.exp(1j*theta*M/n) * np.exp(2j*np.pi*K/n)
        active = reps[(K % (n // period[reps])) == 0]
        if len(active) == 0: continue
        aidx = -np.ones(m, dtype=int); aidx[active] = np.arange(len(active))
        Hk = np.zeros((len(active), len(active)), dtype=complex)
        # diagonal
        Hk[np.arange(len(active)), np.arange(len(active))] = diag[active]
        # off-diagonal: <rep(z);K|H|y;K> += sqrt(p_y/p_rep(z)) h conj(phi_z) lamK^{r_z}
        ysel = aidx[cols] >= 0
        r_, c_, v_ = rows[ysel], cols[ysel], vals[ysel]
        xr = rep[r_]
        keep = aidx[xr] >= 0
        r_, c_, v_, xr = r_[keep], c_[keep], v_[keep], xr[keep]
        coef = v_ * np.conj(phase[r_]) * lamK ** shift[r_] * np.sqrt(period[c_] / period[xr])
        np.add.at(Hk, (aidx[xr], aidx[c_]), coef)
        herm = np.abs(Hk - Hk.conj().T).max()
        if herm > 1e-9:
            print(f"   WARNING: block M={M} K={K} not Hermitian: {herm:.2e}", flush=True)
        eigs.append(np.linalg.eigvalsh(Hk))
    return np.concatenate(eigs), m

for theta in thetas:
    t0 = time.time()
    allev = []
    for M in range(0, n + 1):
        ev, m = sector_spectrum_blocks(M, theta)
        assert len(ev) == m, (len(ev), m)
        allev.append(ev)
        if M > 0: allev.append(ev)
        print(f"   theta={theta:.4f} M={M} dim={m} done ({time.time()-t0:.0f}s)", flush=True)
    ev = np.sort(np.concatenate(allev))
    assert len(ev) == 3**n
    if n <= 8:
        # regression against dense diagonalization
        import scipy.sparse as sp
        Sz = np.diag([1.0, 0.0, -1.0]); Sp = np.zeros((3,3)); Sp[0,1] = Sp[1,2] = np.sqrt(2); Sm = Sp.T
        def op(o, j):
            mats = [sp.identity(3, format='csr')]*n; mats[j] = sp.csr_matrix(o)
            out = mats[0]
            for q in mats[1:]: out = sp.kron(out, q, format='csr')
            return out
        H = sp.csr_matrix((3**n, 3**n), dtype=complex)
        for j in range(n):
            k = (j+1) % n
            ph = np.exp(-1j*theta) if j == n-1 else 1.0
            H = H + op(Sz,j)@op(Sz,k) + 0.5*(ph*op(Sp,j)@op(Sm,k) + np.conj(ph)*op(Sm,j)@op(Sp,k))
        evd = np.linalg.eigvalsh(H.toarray())
        print(f"   regression vs dense: max |diff| = {np.abs(evd - ev).max():.2e}", flush=True)
    E0 = ev[0]; E1 = ev[ev > E0 + 1e-7][0]
    print(f"n={n} theta={theta:.4f}: E0={E0:.10f} E0/n={E0/n:.10f} gap={E1-E0:.10f} mult={int(np.sum(ev < E0+1e-7))}", flush=True)
    for b, tab, tol in ((12.25, table_1225, 30e-6), (10.5, table_105, 1e-5)):
        Z = float(np.sum(np.exp(-b*(ev + a*n))))
        line = f"   b={b}: Z = {Z:.9f}"
        if n in tab and abs(theta) < 1e-12:
            line += f"  table {tab[n][0]*1e-6}  diff {Z - tab[n][0]*1e-6:+.2e}  {'OK' if abs(Z-tab[n][0]*1e-6) < tol else 'MISMATCH'}"
        if n in tab and abs(theta - np.pi) < 1e-12:
            line += f"  table {tab[n][1]*1e-6}  diff {Z - tab[n][1]*1e-6:+.2e}  {'OK' if abs(Z-tab[n][1]*1e-6) < tol else 'MISMATCH'}"
        if n in tab and len(tab[n]) == 3 and abs(theta - 2*np.pi/3) < 1e-9:
            line += f"  table {tab[n][2]*1e-6}  diff {Z - tab[n][2]*1e-6:+.2e}  {'OK' if abs(Z-tab[n][2]*1e-6) < tol else 'MISMATCH'}"
        print(line, flush=True)
    print(f"   [{time.time()-t0:.0f}s]", flush=True)
