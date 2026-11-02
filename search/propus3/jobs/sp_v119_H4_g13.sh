#!/bin/sh
# sp_v119_H4_g13: v=119 H=<13> |H|=4 O_A=20 O_B=32 C=19 |Gamma|=48 patterns=8  NS=3
# exact three-list propus search (propus239/p3, read-only); resumable; hits (reduced orbit masks, one line per Gamma-orbit rep) in hits/<id>_<tag>.txt
# after a hit:  python3 $P3/p3cert.py --v 119 --gens 13 --hits <hits file> --out DIR --matrix ; then repo verify/verify.py on DIR
P3=$(cd "$(dirname "$0")/.." && pwd)
W=${W:-$PWD/work}
ID=sp_v119_H4_g13
THREADS=${THREADS:-4}
cd "$P3" || exit 2
mkdir -p "$W/ckpt" "$W/hits" "$W/probs" "$W/done"
run() {   # rowstag rows ns [shard i/n]
  RT=$1; RW=$2; NSV=$3; SH=${4:-}
  if [ -n "$SH" ]; then ST="_s$(echo "$SH" | tr / o)"; SA="-shard $SH"; else ST=""; SA=""; fi
  K="${ID}_${RT}${ST}"
  [ -e "$W/done/$K" ] && return 0
  [ -s "$W/probs/${ID}_$RT.txt" ] || python3 -B "$P3/jobs/mkprob.py" 119 13 --rows "$RW" --outdir "$W/probs" --tag "$ID" > /dev/null || exit 3
  R=""; [ -s "$W/ckpt/$K.ck" ] && R="-resume"
  ./p3 "$W/probs/${ID}_$RT.txt" -t "$THREADS" -ns "$NSV" $SA -ckpt "$W/ckpt/$K.ck" $R -hits "$W/hits/$K.txt" || exit 4
  touch "$W/done/$K"
}
run r1_15_5 1,15,5 3
run r3_11_15 3,11,15 3
run r5_9_17 5,9,17 3
run r5_15_1 5,15,1 3
run r13_3_17 13,3,17 3
run r15_11_3 15,11,3 3
run r17_3_13 17,3,13 3
run r17_9_5 17,9,5 3
rm -f "$W/probs/${ID}_"*.txt
exit 0
