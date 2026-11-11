"""Independent parser for the saved-run Record files of OpenAI preprint 266.

Reads only data (binary records, centers.txt, commands.json, results/*.txt).
Does NOT execute any OpenAI code. Checks the claims of execution.tex
(Table execute:counts, record format, radius replacement, batch partition,
bad-file subsequence property, fallback input order, final modulo-4 partition,
Root constant ranges from refinement.tex:609-621 / rounding.tex:992-997).
"""
import struct, sys, json, os, math, hashlib
ROOT = "data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/verification/saved-run/execution"
D = os.path.join(ROOT, "data")

def read_records(path):
    with open(path, "rb") as f:
        buf = f.read()
    (n,) = struct.unpack_from("<Q", buf, 0)
    off = 8
    recs = []
    for _ in range(n):
        (k,) = struct.unpack_from("<i", buf, off); off += 4
        p = struct.unpack_from("<36d", buf, off); off += 36 * 8
        T = []
        for _a in range(k):
            T.append(struct.unpack_from("<36d", buf, off)); off += 36 * 8
        C, n0, idv, ort, hor, mnorm, froot, r = struct.unpack_from("<8d", buf, off); off += 64
        recs.append(dict(k=k, p=p, T=T, C=C, n0=n0, id=idv, ort=ort, hor=hor, mnorm=mnorm, froot=froot, r=r))
    assert off == len(buf), (path, off, len(buf))
    return recs

def rec_key(rec):
    # identity of a record ignoring radius (used for subsequence / batch checks)
    return (rec["k"], rec["p"], tuple(rec["T"]), rec["C"], rec["n0"], rec["id"], rec["ort"], rec["hor"], rec["mnorm"], rec["froot"])

def validate(recs, name, rmax):
    bad = 0
    for rec in recs:
        ok = (0 <= rec["k"] <= 10 and all(math.isfinite(x) and abs(x) <= 4 for x in rec["p"])
              and all(abs(t) <= 1.01 for row in rec["T"] for t in row)
              and 0 <= rec["C"] < 20 and 0 <= rec["n0"] < 1e-6 and 0 <= rec["froot"] < 1e-6
              and 0 <= rec["id"] < 1e-4 and 0 <= rec["ort"] < 1e-5 and 0 <= rec["hor"] < 1e-5
              and 0 <= rec["mnorm"] < 2 and 0 <= rec["r"] <= rmax + 1e-15)
        if not ok:
            bad += 1
    rs = [rec["r"] for rec in recs]
    ks = [rec["k"] for rec in recs]
    print(f"{name}: {len(recs)} records; range-check failures={bad}; r in [{min(rs):.6g},{max(rs):.6g}] (cap {rmax}); k histogram={ {k: ks.count(k) for k in sorted(set(ks))} }")
    return bad

def is_ordered_subsequence(sub, seq):
    it = iter(seq)
    return all(any(s == x for x in it) for s in sub)

total_bad = 0
coarse = []
for i in range(4):
    recs = read_records(os.path.join(D, f"coarse{i}.bin"))
    total_bad += validate(recs, f"coarse{i}", 0.125)
    coarse.append(recs)
    assert os.path.getsize(os.path.join(D, f"bad{i}.txt")) == 0, f"bad{i}.txt not empty"
print("A outputs total:", sum(map(len, coarse)), "(paper: 764+774+757+788 =", 764 + 774 + 757 + 788, ")")

mid = read_records(os.path.join(D, "mid.bin")); total_bad += validate(mid, "mid", 0.07)
fine = read_records(os.path.join(D, "fine.bin")); total_bad += validate(fine, "fine", 0.016)
fb1 = read_records(os.path.join(D, "fb1.bin")); total_bad += validate(fb1, "fb1", 0.011)
for f in ["midbad.bin", "finebad.bin", "fb1bad.bin", "fb1cbad0.bin", "fb1cbad1.bin", "fb1cbad2.bin", "fb1cbad3.bin"]:
    assert os.path.getsize(os.path.join(D, f)) == 8 and read_records(os.path.join(D, f)) == [], f
print("zero-count bad files OK (8 bytes each): midbad finebad fb1bad fb1cbad0-3")

# batches: f{j}.bin concatenation == fine.bin with radius replaced by .016
B = os.path.join(D, "batches")
nb = 0
while os.path.exists(os.path.join(B, f"f{nb}.bin")):
    nb += 1
print("number of batches:", nb, "(paper: 88; 87*256+107 =", 87 * 256 + 107, ")")
cat = []
sizes = []
for j in range(nb):
    recs = read_records(os.path.join(B, f"f{j}.bin"))
    sizes.append(len(recs))
    cat.extend(recs)
print("batch sizes:", sorted(set(sizes)), "last:", sizes[-1], "sum:", sum(sizes))
assert len(cat) == len(fine)
mismatch = sum(1 for a, b in zip(cat, fine) if rec_key(a) != rec_key(b))
radius_repl = all(abs(a["r"] - 0.016) < 1e-15 for a in cat)
enlarged = all(a["r"] >= b["r"] - 1e-15 for a, b in zip(cat, fine))
print(f"batch concat == fine (ignoring radius): mismatches={mismatch}; all batch radii == .016: {radius_repl}; replacement enlarges: {enlarged}")

# bad files per batch: ordered subsequence of its batch
total_badrecs = 0
nonempty = []
for j in range(nb):
    bad = read_records(os.path.join(B, f"bad{j}.bin"))
    fj = read_records(os.path.join(B, f"f{j}.bin"))
    if bad:
        assert is_ordered_subsequence([rec_key(x) for x in bad], [rec_key(x) for x in fj]), j
        nonempty.append((f"bad{j}.bin", len(bad)))
        total_badrecs += len(bad)
print("initial C unresolved records:", total_badrecs, "(paper: 163); nonempty bad files:", len(nonempty), "(paper: 15)")
print("nonempty bad files (C-locale sort):", sorted(n for n, _ in nonempty))

# commands.json: fallback input order and final partition
cmd = json.load(open(os.path.join(ROOT, "commands.json")))["children"]
fb_inputs = cmd["C.Bfallback"]["inputs"]
print("fallback inputs match sorted nonempty bad files:", fb_inputs == ["batches/" + n for n in sorted(n for n, _ in nonempty)])
fb_in_recs = []
for name in fb_inputs:
    fb_in_recs.extend(read_records(os.path.join(D, name)))
print("fallback parents:", len(fb_in_recs), "(paper: 163); fallback outputs:", len(fb1), "(paper: 1127)")
# results: C.final{i} done counts
def parse_result_done(path):
    last = [l for l in open(path).read().splitlines() if l.startswith("result")][-1]
    return last
tot = 0
for i in range(4):
    line = parse_result_done(os.path.join(ROOT, "results", f"C.final{i}.txt"))
    done = int(line.split("done")[1].split()[0]); tot += done
    print(f"final C{i}:", line.strip())
print("final C total done:", tot, "== fb1 count:", tot == len(fb1))
# C batch results: exits vs bad counts
rc = {}
for j in range(nb):
    c = cmd[f"C.batch{j}"]
    line = parse_result_done(os.path.join(ROOT, "results", f"C.batch{j}.txt"))
    badn = int(line.split("bad")[1].split()[0])
    done = int(line.split("done")[1].split()[0])
    assert done == sizes[j], (j, done, sizes[j])
    assert (c["returncode"] == 1) == (badn > 0), j
    bfile = len(read_records(os.path.join(B, f"bad{j}.bin")))
    assert bfile == badn, (j, bfile, badn)
    cats = line.split("(")[1].rstrip(")").split()
    assert sum(int(x.lstrip("graphcapoverclusterr")) for x in cats) == badn
    rc[c["returncode"]] = rc.get(c["returncode"], 0) + 1
print("C batch exits:", rc, "(paper: 73 exits 0, 15 exits 1); all done==batch size, bad==file count, categories sum OK")
# cap categories
caps = 0
for j in range(nb):
    line = parse_result_done(os.path.join(ROOT, "results", f"C.batch{j}.txt"))
    caps += int(line.split("cap")[1].split()[0])
print("graph-cap failures:", caps, "(paper: all 163 graph caps)")

# centers.txt  (format per covout.cpp:141: x bins b[0..4], then s, then y tail bins y2..y5)
with open(os.path.join(D, "centers.txt")) as f:
    header = f.readline().split()
    rows = [list(map(int, l.split())) for l in f if l.strip()]
print("centers header:", header, "rows:", len(rows), "(paper: 249,878)")
N = 112
ok = all(len(r) == 10 and 0 <= r[5] <= 4 and all(0 <= b <= 111 for b in r[:5] + r[6:]) for r in rows)
print("all rows: 10 ints, s in 0..4, bins in 0..111:", ok)
xs = set(tuple(r[:5]) for r in rows)
print("distinct five-bin x tuples among centers:", len(xs), "(program printed 'xs 22829' = x lists passing filters, not all need appear)")
print("x bins nondecreasing in every row:", all(r[0] <= r[1] <= r[2] <= r[3] <= r[4] for r in rows))
print("orientation filter a<=l-d+1 (seed.tex:187-189) for all rows:", all(r[0] <= r[4] - r[3] + 1 for r in rows))
def first_minimizer(b):
    s = 0
    for j in range(5):
        if abs(2 * b[j] + 1 - N) < abs(2 * b[s] + 1 - N): s = j
    return s
print("s == first minimizer of |2b_s+1-N| (seed.tex:192-193) for all rows:", all(r[5] == first_minimizer(r[:5]) for r in rows))
print("y tail bins nondecreasing:", all(r[6] <= r[7] <= r[8] <= r[9] for r in rows))
# midpoint gap filters (seed.tex:231-233) using the mk gap formula (seed.tex:244-250)
def midgap(bins):
    q = sorted(bins)
    g = max((q[0] + .5) / N, 1 - (q[4] + .5) / N)
    for j in range(4): g = max(g, (q[j + 1] - q[j]) / N)
    return g
gx = all(midgap(r[:5]) <= 1 - (r[4] - 1) / N + 1e-9 for r in rows)
gy = all(midgap([r[r[5]]] + r[6:]) >= 1 - (r[4] + 2) / N - 1e-9 for r in rows)
print("gap filters: g_x,mid <= 1-(l-1)/N:", gx, "; g_y,mid >= 1-(l+2)/N:", gy)
print("TOTAL range-check failures:", total_bad)
