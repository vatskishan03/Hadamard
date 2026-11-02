#!/bin/bash
# Resumable wrapper for cp_lowe_dd.py (the dedup k-tree) as a scheduler "cmd" job.
#   dd_run.sh WORKROOT <cp_lowe_dd.py arguments...>
# Every call writes round_NN.log in WORKROOT, rebuilds done_sigmas.txt from the "kt4: sigma n/N (a,b,c,d,e)" lines of the earlier
# rounds and passes it as --skip-sigmas (same mechanism as ops/sched.py for the ktree kind).  HADAMARD_WORK = WORKROOT, so the hits
# go to WORKROOT/skew_hits.jsonl.  Exit code of cp_lowe_dd.py is returned (0 = all chosen split values searched or a hit found).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
W=$1; shift
mkdir -p "$W"; cd "$W" || exit 3
export HADAMARD_WORK="$W"
# list files on tmpfs like the production work/tmp (-> /dev/shm/hadamard-kt); removed on exit
SHM=/dev/shm/hadamard-dd-$(basename "$W"); rm -rf "$SHM"; mkdir -p "$SHM"; rm -rf "$W/tmp"; ln -sfn "$SHM" "$W/tmp"
trap 'rm -rf "$SHM"' EXIT
trap 'exit 143' TERM INT
n=$(ls round_*.log 2>/dev/null | wc -l); n=$((n+1)); log=$(printf 'round_%02d.log' $n)
python3 - <<'PY'
import re, glob
pat = re.compile(r'kt4: sigma \d+/\d+ \(([-0-9,]+)\)')
seen = []
for f in sorted(glob.glob('round_*.log')):
    for ln in open(f, errors='replace'):
        m = pat.match(ln)
        if m: seen.append(m.group(1))
seen = list(dict.fromkeys(seen))
open('done_sigmas.txt', 'w').write(''.join('(' + s + ')\n' for s in seen))
print('done split values:', len(seen))
PY
SKIP=""; [ -s done_sigmas.txt ] && SKIP="--skip-sigmas $W/done_sigmas.txt"
echo "$(date) round $n: python3 cp_lowe_dd.py $* $SKIP" >> cmd.txt
python3 "$HERE/cp_lowe_dd.py" "$@" $SKIP > "$log" 2>&1
rc=$?
tail -n 3 "$log"
exit $rc
