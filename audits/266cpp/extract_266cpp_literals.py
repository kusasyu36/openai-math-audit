"""Extract numeric literals and math-library call sites from the distributed C++ sources.
Reads the sources as text only (no compilation, no execution of the OpenAI code)."""
import re, sys, os, collections
BASE = ("data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-"
        "September-24-2026/verification/saved-run/compiled-build/source")
FILES = ["seed.cpp", "covout.cpp", "ub.cpp", "patchx.cpp", "atlas.cpp", "stageB.cpp", "stageC.cpp", "abi-fp-probe.cpp"]
# floating literals: have a '.', or exponent, or hex-float
FLOAT_RE = re.compile(r"(?<![\w.])(0x[0-9a-fA-F.]+p[+-]?\d+|\d+\.\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?|\d+[eE][+-]?\d+)(?![\w.])")
INT_RE = re.compile(r"(?<![\w.])(\d+)(?:u|UL|ull|ULL|UL)?(?![\w.])")
FUNCS = ["pow", "sqrt", "nearbyint", "floor", "ceil", "atan2", "sin", "cos", "tan", "exp", "log", "fma", "hypot",
         "copysign", "isfinite", "fesetround", "fegetround", "atof", "atoi", "strtoull", "abs", "min", "max",
         "__builtin_ctz", "__builtin_ctzll", "__builtin_popcountll", "_mm_getcsr", "cbrt", "asin", "acos", "atan", "round", "lround", "trunc", "fmod", "remainder", "ldexp", "frexp"]
func_re = re.compile(r"(?<![\w.])(" + "|".join(re.escape(f) for f in FUNCS) + r")\s*\(")
kw_re = re.compile(r"long\s+double|float\b|__float128|_Float|-ffast|#pragma|volatile|reinterpret_cast|union\b")
tot = 0
flo = collections.defaultdict(list)
ints = collections.defaultdict(list)
calls = collections.defaultdict(list)
kws = []
for fn in FILES:
    with open(os.path.join(BASE, fn)) as f:
        lines = f.read().split("\n")
    if lines and lines[-1] == "":
        lines = lines[:-1]
    tot += len(lines)
    for ln, line in enumerate(lines, 1):
        code = line
        if code.lstrip().startswith("#include"):
            continue
        # strip string literals
        code_ns = re.sub(r'"(?:\\.|[^"\\])*"', '""', code)
        spans = []
        for m in FLOAT_RE.finditer(code_ns):
            flo[m.group(1)].append(f"{fn}:{ln}")
            spans.append(m.span())
        for m in INT_RE.finditer(code_ns):
            if any(a <= m.start() < b for a, b in spans):
                continue
            ints[m.group(1)].append(f"{fn}:{ln}")
        for m in func_re.finditer(code_ns):
            calls[m.group(1)].append(f"{fn}:{ln}")
        if kw_re.search(code_ns):
            kws.append(f"{fn}:{ln}: {line.strip()[:110]}")
print("total source lines", tot)
def key(s):
    try:
        return float.fromhex(s) if s.startswith("0x") else float(s)
    except Exception:
        return 0.0
print("\n== floating literals: %d distinct, %d occurrences ==" % (len(flo), sum(len(v) for v in flo.values())))
for k in sorted(flo, key=key):
    print(f"{k:>24}  x{len(flo[k]):<3} {' '.join(flo[k])}")
print("\n== integer literals >= 7 (smaller ones are loop bounds/indices): ==")
for k in sorted(ints, key=lambda s: int(s)):
    if int(k) >= 7:
        print(f"{k:>12}  x{len(ints[k]):<3} {' '.join(ints[k][:40])}")
print("\n== small integer literal counts (0..6) ==")
print({k: len(v) for k, v in sorted(ints.items(), key=lambda kv: int(kv[0])) if int(k) < 7})
print("\n== math/library call sites ==")
for f in FUNCS:
    if f in calls and f not in ("abs", "min", "max"):
        print(f"{f:>22} x{len(calls[f]):<3} {' '.join(calls[f])}")
for f in ("abs", "min", "max"):
    print(f"{f:>22} x{len(calls[f])}")
print("\n== keywords (long double / float / pragma / volatile / casts) ==")
for k in kws:
    print(k)
