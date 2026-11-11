# -*- coding: utf-8 -*-
"""Builds and caches the theta coefficient table (kubota.ThetaData) used by stages 2-4.
No check is made here.  Cache: audits/003dual/cache/thetadata_<Nmax>.pkl
Run: UV_CACHE_DIR=$TMPDIR uv run --no-project --with numpy --with scipy --with mpmath python -B -u audits/003dual/stage2_buildcache.py <Nmax>
"""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "029reflect"))
import kubota
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
t0 = time.time()
td = kubota.ThetaData(N, os.path.join(HERE, "cache"))
print("ThetaData(%d): %d entries, %d primes, |g|-1 max %.1e, cube %.1e, %.0fs" % (N, len(td.A), td.checks.get("nprimes", -1), td.checks["abs"], td.checks["cube"], time.time() - t0))
