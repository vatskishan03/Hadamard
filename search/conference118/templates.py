#!/usr/bin/env python3
"""Write templates.txt: for each cubic template graph (K4, the prism, K33, the cube Q3) the solvable classes of
(edge signs S, diagonal D) that ./tmpl lists (see tmpl.c). Format per template: a line `T name n m count`, a line with
the edge list, then one line `S D` per class."""
import subprocess

TEMPLATES = {"K4": (4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
             "prism": (6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]),
             "K33": (6, [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]),
             "Q3": (8, [(0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (1, 3), (4, 6), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7)])}

with open("templates.txt", "w") as f:
    for name, (n, E) in TEMPLATES.items():
        out = subprocess.run(["./tmpl", str(n)] + [str(x) for e in E for x in e], capture_output=True, text=True,
                             check=True).stdout.splitlines()
        cls = [l.split() for l in out[1:] if l.startswith("S=")]
        f.write(f"T {name} {n} {len(E)} {len(cls)}\n" + " ".join(f"{a} {b}" for a, b in E) + "\n")
        f.writelines(f"{s[2:]} {d[2:]}\n" for s, d in cls)
