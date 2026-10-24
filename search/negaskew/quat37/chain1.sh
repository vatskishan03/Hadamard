#!/bin/bash
# sequential chain of the quaternary runs behind the nu = 7 and nu = 9 verdicts of NOTES section 90 (one core, resume-safe: run37.py skips the shapes already in its jsonl; 943 is skipped if runs/p943_b.out has a STATS line).
# Needs runs/lists37_cheap2.jsonl and runs/lists_t5_nu7.jsonl (lister.py output; data/runs/ holds the originals).  About 1.5 core-hours, dominated by 943.
cd "$(dirname "$0")"
export QP=${QP:-./qp37q}
python3 -u run37.py runs/lists37_cheap2.jsonl runs/run_cheap.jsonl 2400 7,9 usym=1 tsym=1 sym=1 lcap=40000000 >> runs/run_cheap.log 2>&1
if ! grep -q "^STATS" runs/p943_b.out 2>/dev/null; then ./go.sh p943_b shapes/p943.txt l2c=1 first=1 maxprint=2 usym=1 tsym=1 sym=1 lcap=50000000; fi
python3 -u run37.py runs/lists_t5_nu7.jsonl runs/run_t5_nu7.jsonl 3000 7 usym=1 tsym=1 sym=1 lcap=40000000 >> runs/run_t5_nu7.log 2>&1
