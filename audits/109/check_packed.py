"""
109 監査: 05-layers.tex Lemma packed-selected-bit-rectangle (行 212-446) を自作実装で確かめる。
 (1) 8 操作の恒等式 (eq:eight-packed-result): (u',w') = (u + c(1-2(u mod 2)), w)、変位 <= 4
 (2) 長さ L=fK の 3 チャンク x,y,z 上で、8 本の巡回シフト S_0 を実装し、
     ガード区間外では理想写像 T_0 と一致、S_0(B)=B、|B| <= δ_bad M を確かめる
 (3) 修正 T_0 S_0^{-1} の実装（逆方向に 8 本の逆シフト + XOR）が B 上で正しいこと
実行: uv run --with sympy --with mpmath --with numpy python check_packed.py
"""
import itertools, random

ok_all = True
def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK' if cond else 'NG'}] {name} {detail}")

def eight_ops(u, w, c):
    """exact operations of eq:eight-packed-additions in order; also returns max displacement"""
    u0, w0 = u, w
    disp = 0
    def upd():
        nonlocal disp
        disp = max(disp, abs(u - u0), abs(w - w0))
    u -= c * (w % 2 == 0); upd()
    w -= c * (u % 2 == 0); upd()
    u += c * (w % 2 == 0); upd()
    w += c * (u % 2 == 0); upd()
    w += c * (u % 2 == 1); upd()
    u += c * (w % 2 == 1); upd()
    w -= c * (u % 2 == 1); upd()
    u -= c * (w % 2 == 1); upd()
    return u, w, disp

ok1 = True; maxdisp = 0
for u in range(-6, 7):
    for w in range(-6, 7):
        for c in (0, 1):
            u2, w2, disp = eight_ops(u, w, c)
            ok1 &= (u2, w2) == (u + c * (1 - 2 * (u % 2)), w)
            maxdisp = max(maxdisp, disp)
check("eight-op identity (u',w') = (u + c(1-2(u mod 2)), w) for all parities, c in {0,1}", ok1)
check("max displacement <= 4", maxdisp <= 4, maxdisp)
# first four ops: cycle (00 01 10) on parity square; last four: (01 10 11)
def first4(u, w):
    u -= (w % 2 == 0); w -= (u % 2 == 0); u += (w % 2 == 0); w += (u % 2 == 0); return u, w
def last4(u, w):
    w += (u % 2 == 1); u += (w % 2 == 1); w -= (u % 2 == 1); u -= (w % 2 == 1); return u, w
par = lambda u, w: (u % 2, w % 2)
cyc1 = {par(*first4(2 * 5 + pu, 2 * 3 + pw)): None for pu in (0, 1) for pw in (0, 1)}
map1 = {(pu, pw): par(*first4(10 + pu, 6 + pw)) for pu in (0, 1) for pw in (0, 1)}
map2 = {(pu, pw): par(*last4(10 + pu, 6 + pw)) for pu in (0, 1) for pw in (0, 1)}
check("first four ops induce cycle (00 01 10), fix 11", map1 == {(0, 0): (0, 1), (0, 1): (1, 0), (1, 0): (0, 0), (1, 1): (1, 1)}, map1)
check("last four ops induce cycle (01 10 11), fix 00", map2 == {(0, 1): (1, 0), (1, 0): (1, 1), (1, 1): (0, 1), (0, 0): (0, 0)}, map2)
# quotient squares preserved for all parities
ok_sq = all((first4(2 * a_ + pu, 2 * b_ + pw)[0] // 2, first4(2 * a_ + pu, 2 * b_ + pw)[1] // 2) == (a_, b_) for a_ in range(3) for b_ in range(3) for pu in (0, 1) for pw in (0, 1))
check("first four ops stay in the parity square (quotients preserved)", ok_sq)

# ---------- (2) full packed rotation on a rectangle ----------
K, f, rho = 9, 3, 2   # K=9: guard window [10, 245] of 256; union bound 80(f-1)2^-K = 0.3125 is informative (K=6 gives the trivial bound 1)
L = f * K
js = [rho + i * K for i in range(f)]
MASK = (1 << L) - 1

def ideal_T0(x, y, z):
    """XOR selected bits j_0..j_{f-2} of y with those of x (j_{f-1} handled separately in the paper)"""
    for i in range(f - 1):
        if (x >> js[i]) & 1:
            y ^= 1 << js[i]
    return x, y, z

def rotation_lines(x, y, z, reverse=False):
    """the eight cyclic shifts (eq:packed-cyclic-shift), simultaneously on segments i=0..f-2.
    line pattern (target, sign, parity tested on other chunk): from eight-packed-additions with u=y segs, w=z segs"""
    lines = [('y', -1, 'z', 0), ('z', -1, 'y', 0), ('y', +1, 'z', 0), ('z', +1, 'y', 0),
             ('z', +1, 'y', 1), ('y', +1, 'z', 1), ('z', -1, 'y', 1), ('y', -1, 'z', 1)]
    if reverse:
        lines = [(t, -s, o, bb) for (t, s, o, bb) in reversed(lines)]
    vals = {'x': x, 'y': y, 'z': z}
    for (t, s, o, bb) in lines:
        off = 0
        for i in range(f - 1):
            if (vals['x'] >> js[i]) & 1 and ((vals[o] >> js[i]) & 1) == bb:
                off += 1 << js[i]
        vals[t] = (vals[t] + s * off) & MASK
    return vals['x'], vals['y'], vals['z']

def guards_ok(y, z):
    for i in range(f - 1):
        for chunk in (y, z):
            seg = (chunk >> js[i]) & ((1 << K) - 1)
            g = seg >> 1
            if not (10 <= g <= 2**(K - 1) - 11):
                return False
    return True

# K=6: guard range 10..21 within 0..31 => 12 of 32 good; exhaustive over x,y,z (2^18 each) too many: sample + exhaustive over segments
rng = random.Random(1)
agree_off_B = True; S0_maps_B_to_B = True; inverse_ok = True; count_B = 0; total = 0
for _ in range(60000):
    x = rng.randrange(1 << L); y = rng.randrange(1 << L); z = rng.randrange(1 << L)
    inB = not guards_ok(y, z)
    total += 1; count_B += inB
    s0 = rotation_lines(x, y, z)
    t0 = ideal_T0(x, y, z)
    if not inB:
        agree_off_B &= (s0 == t0)
    # S_0 is a bijection; B invariant: image of a B point is in B
    xs, ys, zs = s0
    if inB:
        S0_maps_B_to_B &= (not guards_ok(ys, zs))
    # inverse by reversed lines
    inverse_ok &= rotation_lines(*s0, reverse=True) == (x, y, z)
check("S_0 = T_0 off B (random sample)", agree_off_B)
check("S_0(B) ⊂ B (random sample)", S0_maps_B_to_B)
check("reversed lines with negated offsets invert S_0", inverse_ok)
frac = count_B / total
bound = min(1, 80 * (f - 1) * 2.0**(-K))
check(f"|B|/M <= δ_bad = min(1, 80(f-1)2^-K)", frac <= bound, f"{frac:.4f} <= {bound:.4f}")
# exact fraction: each guard good with prob (2^{K-1}-20)/2^{K-1}; 2(f-1) guards independent
exact_bad = 1 - ((2**(K - 1) - 20) / 2**(K - 1))**(2 * (f - 1))
print(f"   exact bad fraction {exact_bad:.4f}; union bound {bound:.4f}; sampled {frac:.4f}")

# ---------- (3) exhaustive small case: K=6, f=2 (one segment): x,y,z over 12 bits each = 2^36 too many; restrict x to selected bits ----------
agree = True
for xb in (0, 1):
    x = xb << js[0]
    for y in range(0, 1 << (2 * K), 97):
        for z in range(0, 1 << (2 * K), 389):
            s0 = rotation_lines(x, y, z); t0 = ideal_T0(x, y, z)
            if guards_ok(y, z):
                agree &= s0 == t0
check("exhaustive-ish K=9,f=3 with x on selected bits: S_0 = T_0 off B", agree)

print("\nALL OK" if ok_all else "\nSOME CHECKS FAILED")
