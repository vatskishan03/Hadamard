#!/usr/bin/env python3
"""count_pairs_oracle.py ORACLE_ROWS m : independent Python support join on the oracle's complete row list: number of ordered
pairs (r1, r2) of sign-normalised rows with complementary supports; raw (a,b) count = 16 x that."""
import sys
from collections import Counter
m = int(sys.argv[2])
rows = [tuple(l.split()) for l in open(sys.argv[1])]
def sup(s): return (s[0] != '0', tuple(c != '0' for c in s[1:]))
cnt = Counter((sup(x), sup(y)) for x, y in rows)
def comp(s): return (not s[0], tuple(not c for c in s[1]))
tot = 0
for x, y in rows:
    tot += cnt[(comp(sup(x)), comp(sup(y)))]
print(f"rows {len(rows)}; ordered (r1,r2) pairs with complementary supports: {tot}; raw ordered (a,b) = {16*tot}")
