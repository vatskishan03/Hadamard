#!/usr/bin/env python3
"""cmp_rows.py w h ENGINE_OUTDIR ORACLE_ROWS_FILE : engine rows (all-keys mode) vs oracle rows, as sets of (X, Y) strings"""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pgplib import *
w, h = int(sys.argv[1]), int(sys.argv[2])
types = parse_types(w, h)
P = types["P"]
eng = set()
for (t, xi, xo, yi, yo) in read_rows(sys.argv[3], types):
    eng.add((s2([xi] + xo), s2([yi] + yo)))
ora = set()
for ln in open(sys.argv[4]):
    a, b = ln.split()
    ora.add((a, b))
print(f"engine rows {len(eng)}  oracle rows {len(ora)}  equal={eng == ora}  only_engine={len(eng - ora)} only_oracle={len(ora - eng)}")
sys.exit(0 if eng == ora else 1)
