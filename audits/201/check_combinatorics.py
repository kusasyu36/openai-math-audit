"""Audit 201: combinatorial checks for Sections 9 and 10.
(1) Directed-cycle trace test (proof of Thm 9.1): tr(S^omega) = l/|O| if omega in orbit of omega0, else 0.
(2) Prop 10.1 annihilator: with offsets s1:{1,2}, s2:{3,4} mod 11, the identifications
    c[i, s_c B] = c[j, B s_c] on every edge i<-j of color c (i-j in offsets of c)
    have equivalence classes = cyclic orbits of words (independent of vertex).
(3) No return path of 1 or 2 forward steps; reverse edges uniquely matched.
"""
import itertools
import numpy as np

def rot_orbit(w):
    return {w[k:] + w[:k] for k in range(len(w))}

# (1)
ok1 = True
for l in range(1, 8):
    words = list(itertools.product((1, 2), repeat=l))
    for w0 in words:
        S = {c: np.zeros((l, l), dtype=np.int64) for c in (1, 2)}
        for j, c in enumerate(w0):
            S[c][j, (j + 1) % l] = 1
        O = rot_orbit(w0)
        for w in words:
            M = np.eye(l, dtype=np.int64)
            for c in w:
                M = M @ S[c]
            tr = int(np.trace(M))
            exp = (l // len(O)) if w in O else 0
            if tr != exp:
                ok1 = False
print("(1) cycle trace test l<=7:", ok1)

# (2)
m = 11
off = {1: (1, 2), 2: (3, 4)}
ok2 = True
for l in range(1, 9):
    words = list(itertools.product((1, 2), repeat=l))
    parent = {}
    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent.setdefault(parent[x], parent[x])
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    for i in range(m):
        for w in words:
            find((i, w))
    for c in (1, 2):
        for r in off[c]:
            for j in range(m):
                i = (j + r) % m
                for B in itertools.product((1, 2), repeat=l - 1):
                    union((i, (c,) + B), (j, B + (c,)))
    classes = {}
    for i in range(m):
        for w in words:
            classes.setdefault(find((i, w)), set()).add((i, w))
    expected = {frozenset((i, w) for i in range(m) for w in rot_orbit(w0)) for w0 in words}
    got = {frozenset(v) for v in classes.values()}
    if got != expected:
        ok2 = False
        print("mismatch at l", l, len(got), len(expected))
print("(2) annihilator classes = cyclic orbits x all vertices, l<=8:", ok2)

# (3)
offs = [1, 2, 3, 4]
print("(3) one/two-step returns:", [s for s in offs if s % m == 0] + [a + b for a in offs for b in offs if (a + b) % m == 0])
rev = {(-r) % m: r for r in offs}
print("    reverse offsets -> unique forward:", rev, "distinct:", len(rev) == 4)
# also: which other m would work?
good = [mm for mm in range(2, 30) if all((a + b) % mm for a in offs for b in offs) and all(a % mm for a in offs) and len({(-r) % mm for r in offs} | set(o % mm for o in offs)) == 8]
print("    moduli with no 1/2-step returns and reverse!=forward:", good)
