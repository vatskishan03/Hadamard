#!/usr/bin/env python3
"""Queue runner for lowe_sz.py: NSLOTS slots, each claims the first pending ticket of the queue and runs
`lowe_sz.py TICKET --auto --budget LEFT` until the runner's budget is spent or the queue is empty. Every call stops
by itself (exit code 75) before LEFT runs out, and a later runner continues the same ticket from its saved lists and
join state, so a run never outlives its runner and resumes after an interruption.

    python3 search/runq.py BUDGET_SECONDS [NSLOTS] [--queue FILE]

The queue (default $HADAMARD_WORK/lowe_queue.txt) holds one ticket per line, the arguments of lowe_sz.py:

    341 2 --rows 1,21,9,29 --nsig 3e5 --threads 1 --kt4-procs 2       pending
    run 341 2 --rows ...                                               claimed by a live runner
    done 341 2 --rows ...                                              finished (or a matrix was found)
    fail(RC) 341 2 --rows ...                                          lowe_sz.py exited with RC
    # comment

A runner releases every claim when it starts (those of a killed runner), so use one runner per queue. Logs go to
$HADAMARD_WORK/logs/runq_<slot>.log and hits to $HADAMARD_WORK/logs/ssxx_hits.jsonl; a verified hit is certified at
once with ssxx_certify.py.
"""
import argparse, fcntl, os, subprocess, sys, threading, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import workdir as W

REPO = Path(__file__).resolve().parent
STATES = ("run", "done")


class Queue:
    def __init__(self, path):
        self.path = Path(path)

    def locked(self, fn):
        with open(str(self.path) + ".lock", "w") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            lines = self.path.read_text().splitlines()
            res = fn(lines)
            tmp = Path(str(self.path) + ".tmp")
            tmp.write_text("\n".join(lines) + "\n")
            os.replace(tmp, self.path)
            return res


def pending(l):
    t = l.strip()
    return bool(t) and not t.startswith("#") and t.split()[0] not in STATES and not t.startswith("fail(")


def reset(lines):
    for i, l in enumerate(lines):
        if l.startswith("run "):
            lines[i] = l[4:]


def claim(lines):
    for i, l in enumerate(lines):
        if pending(l):
            lines[i] = "run " + l
            return l
    return None


def release(q, ticket, rc):
    def f(lines):
        for i, l in enumerate(lines):
            if l == "run " + ticket:
                lines[i] = ticket if rc == 75 else ("done " + ticket if rc == 0 else f"fail({rc}) " + ticket)
                return
    q.locked(f)


def worker(q, slot, end, logs, hits):
    log = logs / f"runq_{slot}.log"
    env = dict(os.environ, HADAMARD_WORK=str(W.WORK))
    while end - time.time() > 900:
        t = q.locked(claim)
        if t is None:
            return
        left = int(end - time.time() - 120)
        with open(log, "a") as f:
            f.write(f"--- {time.strftime('%H:%M')} {t} (budget {left}s)\n"); f.flush()
            start = f.tell()
            rc = subprocess.run([sys.executable, "-u", "lowe_sz.py", *t.split(), "--auto", "--budget", str(left),
                                 "--log", str(hits)], cwd=REPO, stdout=f, stderr=subprocess.STDOUT, env=env).returncode
        N = int(t.split()[0])
        with open(log) as f:
            f.seek(start); txt = f.read()
        if f"verified: complex Hadamard matrix of order {2 * N}" in txt:
            with open(log, "a") as f:
                subprocess.run([sys.executable, "ssxx_certify.py", str(hits), str(N)], cwd=REPO, stdout=f,
                               stderr=subprocess.STDOUT, env=env)
                f.write(f"FOUND N = {N}\n")
            rc = 0
        release(q, t, rc)
        if rc == 75:
            return


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("budget", type=float); ap.add_argument("nslots", type=int, nargs="?", default=2)
    ap.add_argument("--queue", default=str(W.WORK / "lowe_queue.txt"))
    a = ap.parse_args()
    q = Queue(a.queue)
    logs = W.WORK / "logs"; logs.mkdir(exist_ok=True)
    hits = logs / "ssxx_hits.jsonl"
    end = time.time() + a.budget
    q.locked(reset)
    ths = [threading.Thread(target=worker, args=(q, k + 1, end, logs, hits)) for k in range(a.nslots)]
    for th in ths:
        th.start()
    for th in ths:
        th.join()
    lines = q.path.read_text().splitlines()
    print(time.strftime("%H:%M"), "runner done;", sum(pending(l) for l in lines), "tickets pending;",
          sum(l.startswith("done ") for l in lines), "done;", sum(l.startswith("fail(") for l in lines), "failed")
    for l in lines:
        if l.strip() and not l.startswith("#"):
            print("  ", l)


if __name__ == "__main__":
    main()
