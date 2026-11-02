#!/bin/sh
# sp_v267_H11_g4: v=267 H=<4> |H|=11 O_A=14 O_B=27 C=13 |Gamma|=288 patterns=2  NS=3
# exact three-list propus search (propus239/p3, read-only); resumable; hits (reduced orbit masks, one line per Gamma-orbit rep) in hits/<id>_<tag>.txt
# after a hit:  python3 $P3/p3cert.py --v 267 --gens 4 --hits <hits file> --out DIR --matrix ; then repo verify/verify.py on DIR
P3=$(cd "$(dirname "$0")/.." && pwd)
W=${W:-$PWD/work}
ID=sp_v267_H11_g4
THREADS=${THREADS:-4}
cd "$P3" || exit 2
mkdir -p "$W/ckpt" "$W/hits" "$W/probs" "$W/done"
run() {   # rowstag rows ns [shard i/n]
  RT=$1; RW=$2; NSV=$3; SH=${4:-}
  if [ -n "$SH" ]; then ST="_s$(echo "$SH" | tr / o)"; SA="-shard $SH"; else ST=""; SA=""; fi
  K="${ID}_${RT}${ST}"
  [ -e "$W/done/$K" ] && return 0
  [ -s "$W/probs/${ID}_$RT.txt" ] || python3 -B "$P3/jobs/mkprob.py" 267 4 --rows "$RW" --outdir "$W/probs" --tag "$ID" > /dev/null || exit 3
  R=""; [ -s "$W/ckpt/$K.ck" ] && R="-resume"
  ./p3 "$W/probs/${ID}_$RT.txt" -t "$THREADS" -ns "$NSV" $SA -ckpt "$W/ckpt/$K.ck" $R -hits "$W/hits/$K.txt" || exit 4
  touch "$W/done/$K"
}
run r1_23_3 1,23,3 3
run r3_23_1 3,23,1 3
rm -f "$W/probs/${ID}_"*.txt
exit 0
