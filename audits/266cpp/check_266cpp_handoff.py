"""Check the recorded C-stage invocations against the saved data (own script; reads data only).
- all 88 initial C argv identical in parameters, part=0, parts=1, each reads exactly its own f<j>.bin
- each C.batch result 'done' count == record count of its batch file; bad count == count in bad<j>.bin
- final C: parts 0..3 of 4 over fb1.bin; done counts sum to fb1 count and match the modulo partition
- A: parts 0..3 of 4 over centers.txt; 'done' counts match the modulo partition of the center count
- B: argv thresholds as in execution.tex; 'result N bad 0' matches output counts"""
import json, os, re, struct, sys
E = ("data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-"
     "September-24-2026/verification/saved-run/execution")
ch = json.load(open(os.path.join(E, "commands.json")))["children"]
ok = True
def chk(name, cond, extra=""):
    global ok
    ok &= bool(cond)
    print(f"[{'OK' if cond else 'FAIL'}] {name} {extra}")
def count(path):
    with open(os.path.join(E, "data", path), "rb") as f:
        return struct.unpack("<Q", f.read(8))[0]
def last_result(name):
    lines = open(os.path.join(E, "results", name + ".txt")).read().strip().split("\n")
    return [l for l in lines if l.startswith("result") or l.startswith("total")][-1]
print("children:", len(ch))
batches = sorted((k for k in ch if k.startswith("C.batch")), key=lambda s: int(s[7:]))
chk("88 initial C invocations", len(batches) == 88 and [int(b[7:]) for b in batches] == list(range(88)))
params = {tuple(ch[b]["argv"][:10]) for b in batches}
chk("all 88 argv[0:10] identical", len(params) == 1, str(params))
chk("  == crun 12 4 1 .042 1 1 5000000 0 1", params == {("./build/crun", "12", "4", "1", ".042", "1", "1", "5000000", "0", "1")})
io_ok = all(ch[f"C.batch{j}"]["argv"][10:] == [f"batches/bad{j}.bin", f"batches/f{j}.bin"] for j in range(88))
chk("each batch j: output batches/bad<j>.bin, single input batches/f<j>.bin", io_ok)
tot = totbad = 0; per_ok = True; rc_ok = True; nonempty = []
for j in range(88):
    n = count(f"batches/f{j}.bin"); nb = count(f"batches/bad{j}.bin")
    m = re.search(r"result done(\d+) bad(\d+) \(graph(\d+) cap(\d+) over(\d+) cluster(\d+)\)", last_result(f"C.batch{j}"))
    d, b, g, c, o, cl = map(int, m.groups())
    per_ok &= (d == n and b == nb and g + c + o + cl == b)
    rc_ok &= (ch[f"C.batch{j}"]["returncode"] == (1 if nb else 0)) and ch[f"C.batch{j}"]["completed"]
    tot += n; totbad += nb
    if nb: nonempty.append(f"batches/bad{j}.bin")
chk("each batch: logged done == records in f<j>.bin, logged bad == records in bad<j>.bin", per_ok)
chk("each batch: return code 1 iff bad>0, completed", rc_ok)
fd = open(os.path.join(E, "data", "batches", "fdone.txt")).read().split()
chk("fdone.txt == (88, total batch records) and total == fine.bin count", fd == ["88", str(tot)] and tot == count("fine.bin"), f"total={tot} bad={totbad}")
fb = ch["C.Bfallback"]["argv"]
chk("fallback argv = brun .011 .0038 fb1.bin fb1bad.bin - + all nonempty bad files in C-locale order",
    fb[:6] == ["./build/brun", ".011", ".0038", "fb1.bin", "fb1bad.bin", "-"] and fb[6:] == sorted(nonempty), f"{len(fb)-6} files")
m = re.search(r"result (\d+) bad (\d+) .*parent(\d+)", last_result("C.Bfallback"))
chk("fallback: parents == initial C bad total, output == fb1.bin count, bad 0",
    int(m.group(3)) == totbad and int(m.group(1)) == count("fb1.bin") and int(m.group(2)) == 0 and count("fb1bad.bin") == 0, m.group(0)[:60])
nfb = count("fb1.bin"); s = 0; fin_ok = True
for i in range(4):
    a = ch[f"C.final{i}"]["argv"]
    fin_ok &= a == ["./build/crun", "12", "4", "1", ".042", "1", "1", "5000000", str(i), "4", f"fb1cbad{i}.bin", "fb1.bin"]
    m = re.search(r"result done(\d+) bad(\d+)", last_result(f"C.final{i}"))
    expect = len(range(i, nfb, 4))
    fin_ok &= int(m.group(1)) == expect and int(m.group(2)) == 0 and count(f"fb1cbad{i}.bin") == 0 and ch[f"C.final{i}"]["returncode"] == 0
    s += int(m.group(1))
chk("final C: argv parts i/4 over fb1.bin, done == modulo-class size, bad 0, rc 0; sum == fb1 count", fin_ok and s == nfb, f"sum={s}")
with open(os.path.join(E, "data", "centers.txt")) as f:
    ncent = int(f.readline().split()[1])
a_ok = True; sa = 0; outs = 0
for i in range(4):
    a = ch[f"A{i}"]["argv"]
    a_ok &= a == ["./build/arun", "centers.txt", str(i), "4", f"coarse{i}.bin", f"bad{i}.txt"]
    m = re.search(r"total done(\d+) .*global(\d+) failedBalls(\d+)", last_result(f"A{i}"))
    a_ok &= int(m.group(1)) == len(range(i, ncent, 4)) and int(m.group(3)) == 0 and int(m.group(2)) == count(f"coarse{i}.bin")
    a_ok &= ch[f"A{i}"]["returncode"] == 0 and os.path.getsize(os.path.join(E, "data", f"bad{i}.txt")) == 0
    sa += int(m.group(1)); outs += count(f"coarse{i}.bin")
chk("A: argv parts i/4, done == modulo-class size of centers, failedBalls 0, global == coarse count, bad txt empty", a_ok and sa == ncent, f"centers={ncent} coarse total={outs}")
bm = ch["Bmid"]["argv"]; bf = ch["Bfine"]["argv"]
chk("Bmid argv", bm == ["./build/brun", ".07", ".013", "mid.bin", "midbad.bin", "-", "coarse0.bin", "coarse1.bin", "coarse2.bin", "coarse3.bin"])
chk("Bfine argv", bf == ["./build/brun", ".016", ".005", "fine.bin", "finebad.bin", "batches/f", "mid.bin"])
for nm, outf, badf, nin in (("Bmid", "mid.bin", "midbad.bin", outs), ("Bfine", "fine.bin", "finebad.bin", count("mid.bin"))):
    m = re.search(r"result (\d+) bad (\d+) .*parent(\d+)", last_result(nm))
    chk(f"{nm}: parents == input count, result == {outf} count, bad 0, rc 0",
        int(m.group(3)) == nin and int(m.group(1)) == count(outf) and int(m.group(2)) == 0 and count(badf) == 0 and ch[nm]["returncode"] == 0, m.group(0)[:50])
chk("covrun argv", ch["covrun"]["argv"] == ["./build/covrun", "centers.txt"] and ch["covrun"]["returncode"] == 0)
print("\nALL OK" if ok else "\nSOME CHECK FAILED")
