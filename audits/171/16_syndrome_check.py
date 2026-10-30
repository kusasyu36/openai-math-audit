"""Brute-force check of the syndrome / late-class construction in
16-low-modes-geometry-and-cell-states.tex, lines 49-77 (Lemma lem:late-classes-cells).

H = F_2^m with |H| = H_* a power of two, n <= H_* < 2n.
H0 = linear hyperplane {x : top bit = 0}.
IDs a_1..a_n distinct, containing all of H0 (plus n - H_*/2 elements outside H0).
L = subspace of size r, NOT contained in H0.
s(z) = XOR_{k: z_k=1} a_k.  Odd role late iff s(z) in L; class = value of s(z).
Claims checked:
 (C1) each class {odd z : s(z)=t} has size 2^{n-1}/H_*  (for every t in H)
 (C2) each even v has at most one neighbour in each late class
 (C3) p(v) = #late neighbours satisfies r/2 <= p(v) <= r
 (C4) with classes ordered alternating H0 / non-H0, every suffix of length j
      contains between floor(j/2) and j neighbours of v
 (C5) if internal axes have IDs in distinct cosets mod L, each v has
      at most one late *internal* neighbour
Also a negative control: L contained in H0 can give p(v) < r/2.
Everything is exact integer arithmetic.
"""
import itertools, random, sys

def span(gens):
    S = {0}
    for g in gens:
        S |= {x ^ g for x in S}
    return S

def run(n, m_r, seed, L_in_H0=False, verbose=False):
    rng = random.Random(seed)
    Hs = 1
    while Hs < n:
        Hs *= 2
    assert n <= Hs < 2 * n
    m = Hs.bit_length() - 1
    top = 1 << (m - 1) if m >= 1 else 0
    H = list(range(Hs))
    H0 = [x for x in H if not (x & top)]
    out = [x for x in H if x & top]
    extra = n - len(H0)
    assert 0 <= extra <= len(out)
    ids_set = H0 + rng.sample(out, extra)
    # random subspace L of dim m_r, not contained in H0 (unless L_in_H0)
    while True:
        gens = [rng.randrange(1, Hs) for _ in range(m_r)]
        L = span(gens)
        if len(L) != 2 ** m_r:
            continue
        inH0 = all(not (t & top) for t in L)
        if inH0 == L_in_H0:
            break
    r = len(L)
    # internal axes: choose up to a few IDs in H0 from distinct cosets mod L
    cosets_used = set()
    internal_ids = []
    for x in rng.sample(H0, len(H0)):
        c = min(x ^ t for t in L)
        if c not in cosets_used:
            cosets_used.add(c)
            internal_ids.append(x)
        if len(internal_ids) >= 3:
            break
    # assign IDs to coordinates in random order
    rng.shuffle(ids_set)
    a = ids_set
    internal_coords = {k for k in range(n) if a[k] in internal_ids}
    # ordering of classes: alternate in-H0 / out-H0
    Lin = sorted(t for t in L if not (t & top))
    Lout = sorted(t for t in L if t & top)
    order = []
    for x, y in itertools.zip_longest(Lin, Lout):
        if x is not None: order.append(x)
        if y is not None: order.append(y)
    # syndromes of all words
    syn = [0] * (1 << n)
    for z in range(1, 1 << n):
        low = z & -z
        k = low.bit_length() - 1
        syn[z] = syn[z ^ low] ^ a[k]
    par = [bin(z).count('1') & 1 for z in range(1 << n)]
    ok = {k: True for k in ['C1', 'C2', 'C3', 'C4', 'C5']}
    # C1
    cnt = {}
    for z in range(1 << n):
        if par[z] == 1:
            cnt[syn[z]] = cnt.get(syn[z], 0) + 1
    exp_size = 2 ** (n - 1) // Hs
    assert 2 ** (n - 1) % Hs == 0
    if not (len(cnt) == Hs and all(c == exp_size for c in cnt.values())):
        ok['C1'] = False
    minp = 10 ** 9
    for v in range(1 << n):
        if par[v] != 0:
            continue
        per_class = {}
        late_int = 0
        for k in range(n):
            t = syn[v] ^ a[k]
            if t in L:
                per_class[t] = per_class.get(t, 0) + 1
                if k in internal_coords:
                    late_int += 1
        if any(c > 1 for c in per_class.values()):
            ok['C2'] = False
        p = sum(per_class.values())
        minp = min(minp, p)
        if not (r // 2 <= p <= r):
            ok['C3'] = False
        for j in range(r + 1):
            suf = order[r - j:]
            c = sum(per_class.get(t, 0) for t in suf)
            if not (j // 2 <= c <= j):
                ok['C4'] = False
        if late_int > 1:
            ok['C5'] = False
    return dict(n=n, Hs=Hs, r=r, L_in_H0=L_in_H0, minp=minp, **ok)

if __name__ == '__main__':
    allgood = True
    for n in range(4, 15):
        for m_r in range(1, 4):
            Hs = 1
            while Hs < n: Hs *= 2
            if 2 ** m_r > Hs:
                continue
            for seed in range(3):
                res = run(n, m_r, seed)
                good = all(res[k] for k in ['C1', 'C2', 'C3', 'C4', 'C5'])
                allgood &= good
                if not good or seed == 0:
                    print(res)
    print('ALL POSITIVE CHECKS PASS:', allgood)
    # negative control: L inside H0
    print('--- negative control (L contained in H0) ---')
    for n in [5, 6, 7, 9, 10, 11, 12, 13]:
        for seed in range(3):
            try:
                res = run(n, 2, seed, L_in_H0=True)
            except AssertionError:
                continue
            print(res)
