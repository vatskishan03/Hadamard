#!/usr/bin/env python3
"""driver.py -- Python driver around src/pgp (C).

  driver.py plan   W H                      exact numbers of candidates stage 1 will enumerate (about 7 s at (97,3))
  driver.py run    W H WORKDIR [--threads T] [--no-cert]    stage 1 (resumable) -> stage 2 -> classes -> certificate
  driver.py status W H WORKDIR              progress from WORKDIR/done.log (needs the plan; computed once and cached)

W prime = 1 (mod 4), H odd, H | W - 1, (W - 1)/H <= 32.  Types are always swap-reduced (the swap a <-> b is a symmetry).  Everything is
resumable: a bucket (type, key) is recorded in WORKDIR/done.log only after its rows file is written; a rerun skips recorded buckets."""
import argparse, os, re, subprocess, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
PGP = os.path.join(HERE, "..", "src", "pgp")
def sh(cmd, **kw): return subprocess.run(cmd, check=True, **kw)
def ensure_binary():
    if not os.path.exists(PGP):
        sh(["clang", "-O3", "-march=native", "-pthread", "-o", PGP, os.path.join(HERE, "..", "src", "pgp.c"), "-lm"])
def plan(w, h):
    r = sh([PGP, "plan", str(w), str(h), "--swap-reduce"], capture_output=True, text=True)
    return r.stdout
def parse_plan(txt):
    m = re.search(r"TOTAL leaves ([0-9.e+]+).*?buckets ([0-9]+)", txt)
    return float(m.group(1)), int(m.group(2))
ap = argparse.ArgumentParser()
ap.add_argument("cmd", choices=["plan", "run", "status"]); ap.add_argument("w", type=int); ap.add_argument("h", type=int); ap.add_argument("workdir", nargs="?")
ap.add_argument("--threads", type=int, default=os.cpu_count()); ap.add_argument("--no-cert", action="store_true")
a = ap.parse_args()
ensure_binary()
if a.cmd == "plan":
    print(plan(a.w, a.h))
elif a.cmd == "status":
    cache = os.path.join(a.workdir, "plan.txt")
    if not os.path.exists(cache):
        os.makedirs(a.workdir, exist_ok=True); open(cache, "w").write(plan(a.w, a.h))
    tl, tb = parse_plan(open(cache).read())
    nb = nl = nr = 0; cpu = 0.0
    for ln in open(os.path.join(a.workdir, "done.log")):
        f = ln.split(); nb += 1; nr += int(f[2]); nl += int(f[3]); cpu += float(f[6]) if len(f) > 6 else 0
    print(f"buckets {nb}/{tb} ({100*nb/tb:.1f}%)  leaves {nl:.4g}/{tl:.4g} ({100*nl/tl:.1f}%)  rows so far {nr}  thread-CPU {cpu/3600:.2f} h"
          + (f"  remaining at the measured average rate: {(tl-nl)/(nl/cpu)/3600/a.threads:.2f} h on {a.threads} threads" if cpu > 0 and nl > 0 else ""))
else:
    w, h, d = a.w, a.h, a.workdir
    os.makedirs(d, exist_ok=True)
    t0 = time.time()
    sh([PGP, "rows", str(w), str(h), d, "--swap-reduce", "--threads", str(a.threads)])
    with open(os.path.join(d, "solutions.txt"), "w") as f:
        sh([PGP, "join", str(w), str(h), d, "--swap-reduce"], stdout=f)
    cl = os.path.join(d, "classes.json")
    with open(cl, "w") as f:
        sh(["python3.11", os.path.join(HERE, "extract_pairs.py"), str(w), str(h), d, "--swap-reduce"], stdout=f)
    classes = json.load(open(cl))["classes"]
    print(f"stage 1 + 2 finished in {time.time()-t0:.0f} s: {len(classes)} equivalence class(es) of periodic Golay pairs of length {2*w} with |H| = {h}")
    q = 4 * w + 1
    if classes and not a.no_cert and all(q % p for p in range(2, int(q ** .5) + 1)):
        sh(["python3.11", os.path.join(HERE, "make_cert.py"), "--q", str(q), "--pair", cl, "--class", "0", "--out", os.path.join(d, "cert")])
