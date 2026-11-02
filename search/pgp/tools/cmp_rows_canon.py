#!/usr/bin/env python3
"""cmp_rows_canon.py w h CANON_OUTDIR ORACLE_ROWS_FILE : rotations of the canonical-key rows must give exactly the oracle rows"""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pgplib import *
w, h = int(sys.argv[1]), int(sys.argv[2])
types = parse_types(w, h)
m = types["P"].m
eng = set()
rows = read_rows(sys.argv[3], types)
for (t, xi, xo, yi, yo) in rows:
    for r in range(m):
        for e1 in (0, 1):
            for e2 in (0, 1):
                eng.add((s2([xi] + [xo[(j + r + e1 * (m // 2)) % m] for j in range(m)]), s2([yi] + [yo[(j + r + e2 * (m // 2)) % m] for j in range(m)])))
ora = set(tuple(ln.split()) for ln in open(sys.argv[4]))
print(f"canonical rows {len(rows)} -> closure under rotations and reversals {len(eng)}  oracle {len(ora)}  equal={eng == ora}  only_engine={len(eng - ora)} only_oracle={len(ora - eng)}")
sys.exit(0 if eng == ora else 1)
