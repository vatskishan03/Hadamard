#!/bin/sh
# Exhaustive search of the propus structure v = 239, H = <10> (order 7), rows (15, 1, 27)  ->  symmetric Hadamard order 956.
#   sh run239.sh            run (or continue after an interruption: the checkpoint is picked up automatically)
#   THREADS=32 NS=3 SHARD=0/1 sh run239.sh
# Needs only a C compiler for the search (prob239.txt is shipped; python3 + numpy only to regenerate it or to certify a hit).
# Memory: the D list is ~7.2e8 entries x 16 bytes = ~11.5 GB (plus ~0.5 GB of counters); use SHARD=i/n (n processes, or n
# machines) to divide it: each shard re-enumerates D but keeps only its share of the sigma groups.
set -e
cd "$(dirname "$0")"
THREADS=${THREADS:-$(getconf _NPROCESSORS_ONLN)}
NS=${NS:-3}
SHARD=${SHARD:-0/1}
TAG=$(echo "$SHARD" | tr '/' '_')
mkdir -p work
[ -x ./p3 ] || cc -O3 -march=native -ffast-math -pthread p3.c -o p3 -lm
RESUME=""
[ -s "work/ck239_$TAG.txt" ] && RESUME="-resume"
echo "threads=$THREADS ns=$NS shard=$SHARD $RESUME"
./p3 prob239.txt -t "$THREADS" -ns "$NS" -shard "$SHARD" -ckpt "work/ck239_$TAG.txt" $RESUME -hits "work/hits239_$TAG.txt" 2>&1 | tee -a "work/run239_$TAG.log"
echo "finished. Hits (orbit masks, one line per Gamma-orbit representative) are in work/hits239_$TAG.txt"
echo "If it is non-empty:  python3 p3cert.py --v 239 --gens 10 --hits work/hits239_$TAG.txt --out work/cert --matrix"
echo "                     python3 /path/to/Hadamard/verify/verify.py \$PWD/work/cert/symmetric/956"
