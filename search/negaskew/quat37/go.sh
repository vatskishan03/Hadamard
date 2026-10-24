#!/bin/bash
# usage: go.sh TAG shapefile qp37q-args...   (single job, logs in runs/TAG.{out,err}; run detached with setsid nohup for jobs longer than 10 minutes)
TAG=$1; SH=$2; shift 2
cd "$(dirname "$0")"
mkdir -p runs
ulimit -v 3500000; export OMP_NUM_THREADS=1
nice -n 5 ${QP:-./qp37q} $SH "$@" > runs/$TAG.out 2> runs/$TAG.err
echo "EXIT $? $(date -u +%T)" >> runs/$TAG.out
