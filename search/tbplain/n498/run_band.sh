#!/bin/bash
# Restartable band run of band.c (see its header).  Resumes by skipping the s0 buckets already listed in OUT/NAME.prog.
# usage: run_band.sh OUT NAME CLO:CHI DLO:DHI [S0LO S0HI]        e.g.  run_band.sh out b16 16:16 16:16 -3 -3
cd "$(dirname "$0")" || exit 1
OUT=$1; NAME=$2; CB=$3; DB=$4; LO=${5:--1000}; HI=${6:-1000}
mkdir -p "$OUT"
ulimit -v 3000000; export OMP_NUM_THREADS=1
echo "$(date -u) start $NAME C=$CB D=$DB s0=[$LO,$HI]" >> "$OUT/$NAME.runlog"
CBAND=$CB DBAND=$DB MAXB=20000000 exec nice -n 12 ./band 41 "$OUT/$NAME" "$LO" "$HI" >> "$OUT/$NAME.stderr" 2>&1
