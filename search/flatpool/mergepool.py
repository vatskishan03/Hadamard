#!/usr/bin/env python3
"""mergepool.py DIR [MAXN] -> DIR/pool.bin = unique canonical necklaces of all finished chunks (c*.done), optionally first MAXN chunks' worth"""
import sys, os, numpy as np
d = sys.argv[1]
arrs = []
for fn in sorted(os.listdir(d), key=lambda s: (len(s), s)):
    if fn.endswith('.done'): arrs.append(np.fromfile(os.path.join(d, fn[:-5] + '.bin'), dtype=np.uint64))
u = np.unique(np.concatenate(arrs)); u.tofile(os.path.join(d, 'pool.bin')); print(d, len(u))
