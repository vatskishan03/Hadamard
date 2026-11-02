#!/usr/bin/env python3
"""timed.py LABEL cmd args... : run a command, append 'label user=.. sys=.. real=.. rc=..' to logs/cpu_ledger.txt (rusage of the child)."""
import os, sys, time, subprocess
label = sys.argv[1]; cmd = sys.argv[2:]
t0 = time.time()
p = subprocess.Popen(cmd)
_, status, ru = os.wait4(p.pid, 0)
real = time.time() - t0
rc = os.waitstatus_to_exitcode(status)
line = f"{label} user={ru.ru_utime:.2f} sys={ru.ru_stime:.2f} real={real:.2f} rc={rc}"
ledger = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs", "cpu_ledger.txt")
with open(ledger, "a") as f: f.write(line + "\n")
print("[cpu] " + line, file=sys.stderr)
sys.exit(rc if rc >= 0 else 1)
