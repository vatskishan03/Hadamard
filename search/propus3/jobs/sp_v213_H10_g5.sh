#!/bin/sh
# sp_v213_H10_g5: v=213 H=<5> |H|=10 O_A=16 O_B=30 C=15 |Gamma|=28 patterns=8  NS=3
# exact three-list propus search (propus239/p3, read-only); resumable; hits (reduced orbit masks, one line per Gamma-orbit rep) in hits/<id>_<tag>.txt
# after a hit:  python3 $P3/p3cert.py --v 213 --gens 5 --hits <hits file> --out DIR --matrix ; then repo verify/verify.py on DIR
P3=$(cd "$(dirname "$0")/.." && pwd)
W=${W:-$PWD/work}
ID=sp_v213_H10_g5
THREADS=${THREADS:-4}
cd "$P3" || exit 2
mkdir -p "$W/ckpt" "$W/hits" "$W/probs" "$W/done"
run() {   # rowstag rows ns [shard i/n]
  RT=$1; RW=$2; NSV=$3; SH=${4:-}
  if [ -n "$SH" ]; then ST="_s$(echo "$SH" | tr / o)"; SA="-shard $SH"; else ST=""; SA=""; fi
  K="${ID}_${RT}${ST}"
  [ -e "$W/done/$K" ] && return 0
  [ -s "$W/probs/${ID}_$RT.txt" ] || python3 -B "$P3/jobs/mkprob.py" 213 5 --rows "$RW" --outdir "$W/probs" --tag "$ID" > /dev/null || exit 3
  R=""; [ -s "$W/ckpt/$K.ck" ] && R="-resume"
  ./p3 "$W/probs/${ID}_$RT.txt" -t "$THREADS" -ns "$NSV" $SA -ckpt "$W/ckpt/$K.ck" $R -hits "$W/hits/$K.txt" || exit 4
  touch "$W/done/$K"
}
run r7_19_9 7,19,9 3
run r9_11_23 9,11,23 3
run r9_19_7 9,19,7 3
run r11_1_27 11,1,27 3
run r11_19_3 11,19,3 3
run r13_11_21 13,11,21 3
run r27_1_11 27,1,11 3
run r29_1_3 29,1,3 3
rm -f "$W/probs/${ID}_"*.txt
exit 0
