# -*- coding: utf-8 -*-
"""
Follow-up to stage3_c_groups (layer D), written AFTER its official run.  No criterion of layer D is changed.
In that run the control ND2 (psi without r^-1) was not rejected in 72 of the 2169 groups where my registered
liveness test ("D_0 does not divide delta'(r^-1 - 1)", i.e. the two characters differ somewhere on O) said it
was live; in those 72 groups the control gave the SAME value (relative difference 1e-16).
Question examined here: do the two characters coincide on the SUPPORT of d_sigma (the only place where psi is
evaluated)?  For every group of the same configuration list: compare psi and psi-without-r^-1 on the first
20000 support points m = lambda^4 l of the table of d_sigma.
Expected if my explanation is right: 'live on O but identical on the support' in exactly 72 groups, and those
are the groups not rejected.  This is a diagnosis of a dead control, not a test of the paper.
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage3_c_nd2_liveness.py > audits/003dual/stage3_c_nd2_liveness.out 2>/dev/null
"""
import sys, itertools
from stage2_lib import *
print(__doc__)
LINT = 18; setup = Setup((LINT, 0))
P7, P13, M5, P19 = (1, 3), (4, 3), (-5, 0), cl.split_prime(19)
CONF = [(nm, {}) for nm in "ABCDEF"] \
     + [(nm, {P7: j}) for nm in "ABCEF" for j in range(6)] \
     + [(nm, jd) for nm in "ABF" for jd in ({P13: 1}, {P13: 4}, {P13: 0}, {M5: 1}, {M5: 2}, {P19: 3}, {P19: 5})] \
     + [(nm, jd) for nm in "ABC" for jd in ({P7: 1, P13: 5}, {P7: 0, P13: 4}, {P7: 2, P13: 0}, {P7: 4, P13: 4}, {P7: 1, P13: 1})]
th = Theta(1000000)
PH = {nm: phi_hat(phi_table(nm, LINT), LINT) for nm in "ABCDEF"}
tot = {}
for nm, jd in CONF:
    primes = list(jd); must = [p for p in primes if jd[p] != 0]; opt = [p for p in primes if jd[p] == 0]
    for h0, f0 in PH[nm].items():
        if abs(f0) < 1e-13: continue
        for kopt in range(len(opt) + 1):
            for extra in itertools.combinations(opt, kopt):
                act = [p for p in primes if p in must or p in extra]
                cfg = setup.build(h0, act, [ONE] * len(act))
                liveO = not divides(cfg["D0"], mul(cfg["dl"], sub(cfg["rinv"], ONE)))
                m = th.d[cfg["sigma"]].m[:20000]
                dsup = float(np.max(np.abs(psi_char(cfg)(m) - AddChar(neg(cfg["dl"]), cfg["D0"])(m))))
                key = (nm, str({norm(p): j for p, j in jd.items()}), str(sorted(norm(p) for p in act)), cfg["case"], norm(cfg["c0"]))
                d = tot.setdefault(key, [0, 0, 0]); d[0] += 1; d[1] += liveO; d[2] += liveO and dsup < 1e-9
nO = sum(d[1] for d in tot.values()); nD = sum(d[2] for d in tot.values())
print("groups: %d ; live by the registered test (characters differ on O): %d ; of these, identical on the support of d_sigma: %d" % (sum(d[0] for d in tot.values()), nO, nD))
print("where (Psi_0, twists, active primes, case, N(c_0)) : #groups dead-on-support")
for key in sorted(tot):
    if tot[key][2]: print("   ", key, ":", tot[key][2], "of", tot[key][0])
