import sys, math, time; sys.path.insert(0, '.')
import numpy as np
import stage7_lib as L; PP=int(sys.argv[1])
def run(v0, M, n_side, shift, Cmax=40):
    nums = L.dual_indices(M); t=time.time()
    A,b,sc,info = L.build_system(nums, v0, L.sample_points(n_side, v0, shift), Cmax, per_point=PP)
    tau, sv, rank, resid = L.solve(A,b,sc)
    print(f"v0={v0} M={M} unk={len(nums)} rows={len(b)} chi={info['chi']} rep={info['rep']} skip={info['skipped']} {info['skipped_rep']} rank={rank} sv[0]={sv[0]:.3g} sv[-1]={sv[-1]:.3g} cond={sv[0]/sv[-1]:.3g} resid={resid:.2e} t={time.time()-t:.0f}s", flush=True)
    return dict(zip(nums, tau)), sv
s1, sv1 = run(0.30, 8.0, 30, (0,0))
s2, sv2 = run(0.26, 9.0, 32, (0.37,0.61))
common = [m for m in s1 if m in s2]
mu = {m: abs(L.mu_of([m])[0]) for m in common}
for Mc in (1.0, 2.0, 3.0, 4.0):
    d = [abs(s1[m]-s2[m]) for m in common if mu[m] <= Mc]
    mag = [abs(s1[m]) for m in common if mu[m] <= Mc]
    print(f"|mu|<={Mc}: n={len(d)} max|diff|={max(d):.2e} max|tau|={max(mag):.3g}")

rows = sorted(((abs(s1[m]-s2[m]), m, mu[m], abs(s1[m])) for m in common if mu[m] <= 1.5), reverse=True)[:15]
for r in rows: print("diff %.2e m=%s |mu|=%.3f |tau|=%.3g" % r)
small = sorted(((abs(s1[m]), m, mu[m]) for m in common if mu[m] <= 1.0))[:12]
for r in small: print("small |tau|=%.2e m=%s |mu|=%.3f" % r)
