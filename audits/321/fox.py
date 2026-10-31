"""Free-group words, Fox derivatives, and evaluation in representations.

Independent implementation (no code from the OpenAI repository is used).

Words are tuples of (generator_index, exponent in {+1,-1}).
Group-ring elements are dicts  {reduced_word: integer coefficient}.
Fox derivative convention (as in the paper, 04-construction.tex l.115-123):
    d(uv)/dg = du/dg + u dv/dg,  dg/dg = 1,  d(g^-1)/dg = -g^-1.
"""
from __future__ import annotations
from collections import defaultdict


def reduce_word(w):
    out = []
    for g, e in w:
        if out and out[-1][0] == g and out[-1][1] == -e:
            out.pop()
        else:
            out.append((g, e))
    return tuple(out)


def inv(w):
    return tuple((g, -e) for g, e in reversed(w))


def mul(*ws):
    r = ()
    for w in ws:
        r = reduce_word(r + tuple(w))
    return r


def power(w, n):
    if n >= 0:
        return mul(*([w] * n)) if n else ()
    return mul(*([inv(w)] * (-n)))


def gen(i, e=1):
    return ((i, e),)


def conj(a, b):
    """a b a^-1"""
    return mul(a, b, inv(a))


def comm(a, b):
    """[a,b] = a b a^-1 b^-1 (paper's convention, 04-construction.tex l.34)"""
    return mul(a, b, inv(a), inv(b))


def fox(word, j):
    """Fox derivative d(word)/dx_j as a group-ring element (dict word->int)."""
    res = defaultdict(int)
    prefix = ()
    for g, e in word:
        if g == j:
            if e == 1:
                res[prefix] += 1
            else:
                res[mul(prefix, gen(g, -1))] -= 1
        prefix = mul(prefix, gen(g, e))
    return {w: c for w, c in res.items() if c != 0}


def evaluate(elem, rep):
    """Evaluate a group-ring element under a representation.

    rep(word) must return something supporting + and * (e.g. sympy scalars
    or matrices).  Returns sum coeff * rep(word)."""
    total = None
    for w, c in elem.items():
        v = c * rep(w)
        total = v if total is None else total + v
    return total


def scalar_rep(images):
    """Rank-one representation from a list of generator images (scalars)."""
    def rep(w):
        v = 1
        for g, e in w:
            v = v * (images[g] if e == 1 else images[g] ** -1)
        return v
    return rep


def fox_matrix(relators, ngens, rep, zero=0):
    """Rows = relators, columns = generators, entries evaluated under rep."""
    rows = []
    for r in relators:
        row = []
        for j in range(ngens):
            d = fox(r, j)
            row.append(evaluate(d, rep) if d else zero)
        rows.append(row)
    return rows
