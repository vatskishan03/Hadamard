#!/bin/bash
# usage: supervise.sh NAME DEADLINE_EPOCH 'K' q 'hunits' SLOTFILE TOPN ORDER ENGINE   (env EX_* passed through)
# restarts exh38.py (resumable) until 'ALL SLOTS FINISHED' appears in its log or the deadline has passed; pid of this supervisor in runs/NAME.sup.pid, logs runs/NAME.log and runs/NAME.sup.log
cd "$(dirname "$0")/.."
mkdir -p runs
NAME=$1; DL=$2; shift 2
echo $$ > runs/$NAME.sup.pid
ulimit -v 3500000
export OMP_NUM_THREADS=1
while true; do
  now=$(date +%s); left=$(( DL - now )); [ $left -le 30 ] && { echo "$(date -u +%T) deadline" >> runs/$NAME.sup.log; break; }
  echo "$(date -u +%T) start exh38 $NAME secs=$left" >> runs/$NAME.sup.log
  nice -n 5 python3 -u exh38.py "$NAME" "$left" "$@" >> runs/$NAME.log 2>&1
  rc=$?
  echo "$(date -u +%T) exh38 exited rc=$rc; last log line: $(tail -n 1 runs/$NAME.log | cut -c1-200)" >> runs/$NAME.sup.log
  tail -n 1 runs/$NAME.log | grep -q "ALL SLOTS FINISHED" && break
  [ -f runs/$NAME.STOP ] && break
  sleep 2
done
