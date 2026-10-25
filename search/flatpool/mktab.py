#!/usr/bin/env python3
"""mktab.py v ulen [gen] -> tab_<v>.txt for flat.c"""
import sys, numpy as np
import gsflat
v, ulen = int(sys.argv[1]), int(sys.argv[2])
gen = int(sys.argv[3]) if len(sys.argv) > 3 else None
S = gsflat.build(v, ulen, gen)
with open(f"tab_{v}.txt", "w") as f:
    f.write(f"{v} {S['m']} {S['C']} {ulen} {S['pm']}\n")
    f.write(" ".join(repr(float(x)) for x in S['eta'].real) + "\n")
    f.write(" ".join(repr(float(x)) for x in S['eta'].imag) + "\n")
print("wrote", f"tab_{v}.txt", "m", S['m'], "C", S['C'], "rho", S['rho'])
