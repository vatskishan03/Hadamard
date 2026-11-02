#!/bin/sh
# sp_v247_H9_g68x87: v=247 H=<68,87> |H|=9 O_A=18 O_B=35 C=17 |Gamma|=48 patterns=12  NS=3
# exact three-list propus search (propus239/p3, read-only); resumable; hits (reduced orbit masks, one line per Gamma-orbit rep) in hits/<id>_<tag>.txt
# after a hit:  python3 $P3/p3cert.py --v 247 --gens 68,87 --hits <hits file> --out DIR --matrix ; then repo verify/verify.py on DIR
P3=$(cd "$(dirname "$0")/.." && pwd)
W=${W:-$PWD/work}
ID=sp_v247_H9_g68x87
THREADS=${THREADS:-8}
cd "$P3" || exit 2
mkdir -p "$W/ckpt" "$W/hits" "$W/probs" "$W/done"
run() {   # rowstag rows ns [shard i/n]
  RT=$1; RW=$2; NSV=$3; SH=${4:-}
  if [ -n "$SH" ]; then ST="_s$(echo "$SH" | tr / o)"; SA="-shard $SH"; else ST=""; SA=""; fi
  K="${ID}_${RT}${ST}"
  [ -e "$W/done/$K" ] && return 0
  [ -s "$W/probs/${ID}_$RT.txt" ] || python3 -B "$P3/jobs/mkprob.py" 247 68,87 --rows "$RW" --outdir "$W/probs" --tag "$ID" > /dev/null || exit 3
  R=""; [ -s "$W/ckpt/$K.ck" ] && R="-resume"
  ./p3 "$W/probs/${ID}_$RT.txt" -t "$THREADS" -ns "$NSV" $SA -ckpt "$W/ckpt/$K.ck" $R -hits "$W/hits/$K.txt" || exit 4
  touch "$W/done/$K"
}
run r5_1_31 5,1,31 3
run r5_13_25 5,13,25 3
run r7_7_29 7,7,29 3
run r7_17_19 7,17,19 3
run r17_13_19 17,13,19 3
run r17_17_11 17,17,11 3
run r19_1_25 19,1,25 3
run r19_7_23 19,7,23 3
run r19_13_17 19,13,17 3
run r19_17_7 19,17,7 3
run r29_7_7 29,7,7 3
run r31_1_5 31,1,5 3
rm -f "$W/probs/${ID}_"*.txt
exit 0
