#!/usr/bin/env python3
"""Run the low-energy k-tree (lowe_sz.py) on the candidates of lowe_scan.py, best first, one at a time.

Score = 10^logEin / (1 + (pool_C + pool_D) / 1e9). Skips orders already certified in ../complex-hadamard, candidates
already run (recorded in WORK/lowe/dispatch_done.txt) and those with logEin < MINEIN. Every verified hit is
certified with ssxx_certify.py (which writes ../complex-hadamard/C<2N>.json).

    python3 search/lowe_dispatch.py SCAN.jsonl LOGFILE [--minein -0.7] [--L 3e7] [--nsig 3e5] [--threads 3]
"""
import argparse, json, os, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import workdir as W


def solved():
    return {json.loads(p.read_text())["N"] for p in (ROOT / "complex-hadamard").glob("C*.json")}


def key(r):
    return f"{r['N']}|{'b' if r['bordered'] else 'c'}|{','.join(map(str, r['gens']))}|{','.join(map(str, r['rows']))}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan"); ap.add_argument("log")
    ap.add_argument("--minein", type=float, default=-0.7)
    ap.add_argument("--L", type=float, default=3e7)
    ap.add_argument("--nsig", type=float, default=3e5)
    ap.add_argument("--threads", type=int, default=3)
    a = ap.parse_args()
    donef = W.WORK / "lowe" / "dispatch_done.txt"; donef.parent.mkdir(exist_ok=True)
    hits = W.WORK / "logs" / "ssxx_hits.jsonl"
    while True:
        done = set(donef.read_text().split()) if donef.exists() else set()
        rows = []
        if Path(a.scan).exists():
            for l in open(a.scan):
                try:
                    r = json.loads(l)
                except Exception:
                    continue
                rows.append(r)
        sol = solved()
        skipf = W.WORK / "lowe" / "dispatch_skipN.txt"            # N handled elsewhere (e.g. by a queue)
        skip = {int(x) for x in skipf.read_text().split()} if skipf.exists() else set()
        cands = [r for r in rows if r["N"] not in sol and r["N"] not in skip and key(r) not in done
                 and r["logEin"] >= a.minein]
        if not cands:
            time.sleep(60); continue
        cands.sort(key=lambda r: -(10 ** r["logEin"]) / (1 + (r["sizes"][2] + r["sizes"][3]) / 1e9))
        r = cands[0]
        with open(donef, "a") as f:
            f.write(key(r) + "\n")
        cmd = [sys.executable, "-u", str(HERE / "lowe_sz.py"), str(r["N"]), *map(str, r["gens"]),
               "--rows", ",".join(map(str, r["rows"])), "--L", str(a.L), "--nsig", str(a.nsig),
               "--threads", str(a.threads), "--log", str(hits)] + (["--bordered"] if r["bordered"] else [])
        with open(a.log, "a") as f:
            f.write(f"--- {time.strftime('%H:%M')} dispatch {key(r)} logE {r['logE']} logEin {r['logEin']}\n"); f.flush()
            p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        txt = open(a.log).read().rsplit("--- ", 1)[-1]
        if f"verified: complex Hadamard matrix of order {2 * r['N']}" in txt:
            c = subprocess.run([sys.executable, str(HERE / "ssxx_certify.py"), str(hits), str(r["N"])],
                               capture_output=True, text=True)
            with open(a.log, "a") as f:
                f.write(f"FOUND N = {r['N']}: certified: {c.stdout.strip()[-200:]} {c.stderr.strip()[-200:]}\n")


if __name__ == "__main__":
    main()
