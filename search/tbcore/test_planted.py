"""End-to-end control on real data: a known solution planted into small random lists must be found by the production code path.

usage: python3 test_planted.py [SOLUTION.json] [--out DIR] [--size L]
Takes a recorded solution (default data/SOLUTION_0_329_743_226805_104054.json; its tuple index is stored in the file), builds the four
lists exactly as s1_prod.py does -- the exhaustive A and B lists of its tuple (zA, zB), C and D lists of L random fibre samples
(s1_prepare.sample_C, fixed seed) with the true C and D appended at the end -- turns them into class vectors (classvec), writes
the binary input of join4 and runs it. Passes iff
  * the recorded (A, B) are found in the exhaustive lists and the planted (a, b, c, d) is among the solutions printed by join4;
  * every solution printed satisfies sum of the four PAF class vectors == target (zs.ZS.paf_classes, not the hash);
  * the planted solution is also found with 1 and with 32 key chunks (the setting of the production run).
"""
import argparse, json, os, subprocess, sys
import numpy as np
import tbio
from s1_prepare import *

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("solution", nargs="?", default=os.path.join(tbio.DATA, "SOLUTION_0_329_743_226805_104054.json"))
ap.add_argument("--out", default=None, help="output directory (default: current directory)")
ap.add_argument("--size", type=int, default=20000)
args = ap.parse_args()
out = tbio.outdir(args.out)
sol = json.load(open(args.solution))
XA0, XB0, XC0, XD0 = [np.array(sol[k], dtype=np.int8) for k in "ABCD"]
f0, f1 = enum_A(); XA = build(f0, f1); z5, q0, q1 = zq_compress(XA)
mQ = (q0 == 0) & (q1 == 3)
def zclass(X):
    z, a, b = zq_compress(X[None, :]); return tuple(int(v) for v in z[0]), int(a[0]), int(b[0])
(zA, qa0, qa1), (zB, qb0, qb1), (zC, qc0, qc1), (zD, qd0, qd1) = [zclass(x) for x in (XA0, XB0, XC0, XD0)]
print("Z29-level classes of the four sequences:", zA, zB, zC, zD, " Q-compressions:", (qa0, qa1), (qb0, qb1), (qc0, qc1), (qd0, qd1))
assert (qa0, qa1) == (0, 3) and (qb0, qb1) == (0, 3) and (qc0, qc1) == (-1, 1) and (qd0, qd1) == (-1, 1)
tuples = None
try:
    tuples = tbio.load_tuples(out)
    print("tuple index recorded in the file:", sol["tuple"], " matches the classes:", tuples[sol["tuple"]][:4] == (zA, zB, zC, zD))
except FileNotFoundError:
    print("(no tuples.json in --out: tuple index not cross-checked)")
LA = XA[mQ & (z5 == np.array(zA)).all(1)]; LB = XA[mQ & (z5 == np.array(zB)).all(1)]
ia = [i for i in range(len(LA)) if np.array_equal(LA[i], XA0)]; ib = [i for i in range(len(LB)) if np.array_equal(LB[i], XB0)]
assert len(ia) == 1 and len(ib) == 1, "recorded A or B not in the exhaustive list"
rng = np.random.default_rng(12345)
L = args.size
fc0, fc1, _ = sample_C(zC, (-1, 1), L, rng); XC = build(fc0, fc1)
fd0, fd1, _ = sample_C(zD, (-1, 1), L, rng); XD = build(fd0, fd1)
XC = np.concatenate([XC, XC0[None, :]]); XD = np.concatenate([XD, XD0[None, :]])
planted = (ia[0], ib[0], L, L)
VA, VB, VC, VD = [classvec(x) for x in (LA, LB, XC, XD)]
print("lists |A| = %d, |B| = %d, |C| = |D| = %d; planted indices %s" % (len(LA), len(LB), len(XC), planted))
exe = tbio.build_join4(out)
fn = os.path.join(out, "planted_tmp.bin")
tbio.write_binary_lists(fn, VA, VB, VC, VD, Tvec)
ok = True
for lch in (0, 5):
    r = subprocess.run([exe, fn, str(lch)], capture_output=True, text=True, check=True)
    sols = [tuple(int(x) for x in l.split()[1:5]) for l in r.stdout.splitlines() if l.startswith("SOL")]
    good = True
    for a, b, c, d in sols:
        tot = sum(st.paf_classes(X[None, :])[0] for X in (LA[a], LB[b], XC[c], XD[d]))
        good &= bool((tot == st.tvec).all())
    found = planted in sols
    print("  join4 with %d key chunks: %d solutions, planted found: %s, all verified against the target: %s" % (1 << lch, len(sols), found, good))
    ok &= found and good
os.remove(fn)
print("test_planted:", "ALL OK" if ok else "FAILURES")
sys.exit(0 if ok else 1)
