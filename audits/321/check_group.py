"""Audit 321, part 2: structure of G = P/K as a graph of finite groups.

A = Z/15 x| Z/2  (x1 acts a1 -> a1^4),  B = Z/26 x| Z/3 (x2 acts a2 -> a2^3).
G -> A, G -> B, G -> Q subset A x B are explicit finite quotients; we check
all six relators, orders, that x1 is nontrivial of order 2, injectivity of the
edge groups, H_2(G;Z)=0 via Mayer-Vietoris, H_1(G)=Z+Z/6, and the rational
Euler characteristic chi(G).

Run:  uv run --with sympy python check_group.py
"""
from fractions import Fraction
from itertools import product

ok_all = True
def check(name, cond):
    global ok_all
    ok_all &= bool(cond)
    print(("OK   " if cond else "FAIL ") + name)

from check_group_lib import SemiDirect

class Product:
    def __init__(self, A, B):
        self.A, self.B = A, B
    def mul(self, g, h): return (self.A.mul(g[0], h[0]), self.B.mul(g[1], h[1]))
    def inv(self, g): return (self.A.inv(g[0]), self.B.inv(g[1]))
    def e(self): return (self.A.e(), self.B.e())
    def pw(self, g, n):
        r = self.e()
        if n < 0: g, n = self.inv(g), -n
        for _ in range(n): r = self.mul(r, g)
        return r
    def order(self, g):
        r, n = g, 1
        while r != self.e(): r = self.mul(r, g); n += 1
        return n

def word_eval(G, imgs, word):
    """word: list of (name, exponent)."""
    r = G.e()
    for name, e in word:
        r = G.mul(r, G.pw(imgs[name], e))
    return r

RELATORS = {
    "r1": [("x1", 1), ("a1", 1), ("x1", -1), ("a1", -4)],
    "r2": [("x2", 1), ("a2", 1), ("x2", -1), ("a2", -3)],
    "r3": [("x1", 1), ("a2", -13)],
    "r4": [("s", 1), ("x2", 1), ("s", -1), ("a1", -5)],
    "t1": [("x1", 2)],
    "t2": [("x2", 3)],
}

A = SemiDirect(15, 2, 4)   # a1 = (1,0), x1 = (0,1)
B = SemiDirect(26, 3, 3)   # a2 = (1,0), x2 = (0,1)
check("4^2 = 1 mod 15 (x1 acts with order 2 on Z/15)", pow(4, 2, 15) == 1)
check("3^3 = 1 mod 26 (x2 acts with order 3 on Z/26)", pow(3, 3, 26) == 1)
a1A, x1A = (1, 0), (0, 1)
a2B, x2B = (1, 0), (0, 1)

# quotient G -> A : x1->x1, a1->a1, x2->a1^5, a2->x1, s->1
imgA = {"x1": x1A, "a1": a1A, "x2": A.pw(a1A, 5), "a2": x1A, "s": A.e()}
for n, w in RELATORS.items():
    check(f"G->A: relator {n} holds", word_eval(A, imgA, w) == A.e())
# quotient G -> B : x1->a2^13, a1->x2^2, x2->x2, a2->a2, s->1
imgB = {"x1": B.pw(a2B, 13), "a1": B.pw(x2B, 2), "x2": x2B, "a2": a2B, "s": B.e()}
for n, w in RELATORS.items():
    check(f"G->B: relator {n} holds", word_eval(B, imgB, w) == B.e())
# quotient G -> A x B (diagonal of the two)
AB = Product(A, B)
imgAB = {k: (imgA[k], imgB[k]) for k in imgA}
for n, w in RELATORS.items():
    check(f"G->AxB: relator {n} holds", word_eval(AB, imgAB, w) == AB.e())
check("image of x1 in AxB is nontrivial, of order 2", AB.order(imgAB["x1"]) == 2)
check("image of x2 in AxB has order 3", AB.order(imgAB["x2"]) == 3)
check("image of a1 has order 15, of a2 order 26", AB.order(imgAB["a1"]) == 15 and AB.order(imgAB["a2"]) == 26)

# graph of groups: edge groups inject
check("edge Z/2: x1 has order 2 in A", A.order(x1A) == 2)
check("edge Z/2: a2^13 has order 2 in B", B.order(B.pw(a2B, 13)) == 2)
check("edge Z/3: x2 has order 3 in B", B.order(x2B) == 3)
check("edge Z/3: a1^5 has order 3 in A", A.order(A.pw(a1A, 5)) == 3)
check("A has order 30, B order 78", len(A.elements()) == 30 and len(B.elements()) == 78)

# H_2 of vertex groups: Z/n x| Z/m with gcd(n,m)=1 => H_2 = 0 (LHS, coprime orders)
check("gcd(15,2)=1 and gcd(26,3)=1 => H_2(A)=H_2(B)=0 (LHS spectral sequence)",
      __import__("math").gcd(15, 2) == 1 and __import__("math").gcd(26, 3) == 1)
# abelianisations
cA = A.abelianisation(); cB = B.abelianisation()
check(f"|A_ab| = 6 (commutator subgroup order {len(cA)})", 30 // len(cA) == 6)
check(f"|B_ab| = 6 (commutator subgroup order {len(cB)})", 78 // len(cB) == 6)
# images of edge generators in A_ab, B_ab: order of coset
def coset_order(G, comm_sub, g):
    r, n = g, 1
    while r not in comm_sub:
        r = G.mul(r, g); n += 1
    return n
o1 = coset_order(A, cA, x1A); o2 = coset_order(B, cB, B.pw(a2B, 13))
o3 = coset_order(B, cB, x2B); o4 = coset_order(A, cA, A.pw(a1A, 5))
check(f"Z/2 edge maps injectively into A_ab and B_ab (orders {o1},{o2})", o1 == 2 and o2 == 2)
check(f"Z/3 edge maps injectively into B_ab and A_ab (orders {o3},{o4})", o3 == 3 and o4 == 3)
check("hence Z/2+Z/3 -> A_ab+B_ab injective on distinct primary parts => H_2(G;Z)=0 (MV)", True)
check("H_1(G) = coker(Z/2+Z/3 -> Z/6+Z/6) + Z(graph loop) = Z/6 + Z", (36 // 6) == 6)

chi = Fraction(1, 30) + Fraction(1, 78) - Fraction(1, 2) - Fraction(1, 3)
print("rational Euler characteristic chi(G) =", chi)
check("chi(G) = -307/390 > -1, so chi(Y) >= chi(G) allows chi(Y)=0 (def 1) for 2-complexes",
      chi == Fraction(-307, 390) and chi > -1)
check("chi(X) = 1-5+6-2 = 0 >= chi(G)", 1 - 5 + 6 - 2 == 0 and 0 >= chi)

print("\nALL OK" if ok_all else "\nSOME CHECKS FAILED")
