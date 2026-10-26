#!/bin/bash
# Reproduce the n = 955 core and the Hadamard matrices of order 3820 from the start and compare with the recorded data.
#   usage: bash reproduce.sh [OUT] [--border]
#   OUT: working directory (default ./out); about 5 minutes on one core, about 1 GB; --border also reruns the border search
#   (1.5 minutes) and compares data/border_k8_pattern2.json. Exit status 0 iff every comparison below holds.
set -e
OUT=out; BORDER=0
for a in "$@"; do case "$a" in --border) BORDER=1;; *) OUT=$a;; esac; done
HERE=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$OUT"
HITA=SOLUTION_0_329_743_226805_104054.json      # the two hits of iteration 0 (tuple 60), indices (a, b, c, d) in the lists
HITB=SOLUTION_0_746_137_184467_29180.json
sha() { sha256sum "$1" | cut -d' ' -f1; }
want() { [ "$1" = "$2" ] && echo "  == $3" || { echo "  MISMATCH: $3 ($1 != $2)"; exit 1; }; }

nice -n 10 python3 "$HERE/s1_tuples.py" --out "$OUT"                    # 5 s:    3316 tuples
nice -n 10 python3 "$HERE/s1_tstats.py" --out "$OUT"                    # 2 s:    list sizes per tuple
nice -n 10 python3 "$HERE/s1_prod.py" 400000 4 101 5 900 --out "$OUT"   # 4 min:  iteration 0 = tuple 60, join, two hits
want "$(sha "$OUT/$HITA")" "$(sha "$HERE/data/$HITA")" "regenerated $HITA is byte-identical to data/"
want "$(sha "$OUT/$HITB")" "$(sha "$HERE/data/$HITB")" "regenerated $HITB is byte-identical to data/"
python3 "$HERE/check_core.py" "$OUT/$HITA" "$OUT/$HITB"                 # independent exact check of the sequences
build() {   # build NAME HIT K478_SHA H_INT8_SHA H_CSV_SHA
  mkdir -p "$OUT/$1"
  python3 "$HERE/assemble_K.py" "$OUT/$2" --out "$OUT/$1"
  want "$(sha "$OUT/$1/K478.npy")" "$3" "K478.npy sha256 ($1)"
  python3 "$HERE/cyk955.py" --out "$OUT/$1" | tee "$OUT/$1/cyk955.log"
  want "$(grep 'int8 rowmajor' "$OUT/$1/cyk955.log" | awk '{print $NF}')" "$4" "sha256 of the int8 matrix H ($1)"
  want "$(grep 'canonical CSV' "$OUT/$1/cyk955.log" | awk '{print $NF}')" "$5" "sha256 of the canonical CSV of H ($1)"
}
build hitB "$HITB" fdf17efee93da270ad2b7095de8c1bca0fb41148fe2e57ebdb352618a408bfa0 \
  66f69fadacffb8d838b4355d3524e5b0941770697a276afd76aa97908e650c75 145a49849d9b6b01e9ec0cba0ba0486f89e859fd2bbc926a663337f9f9587f49
build hitA "$HITA" dae04363eae2d2914cc9224061f8e9219d8bf305dfd020a1ee2ac40a669b9374 \
  c9b9cc08f2e750b3830217f0327c77dd8bde3f8a800ad959e8b78e48341e470b 79c29912373c20336c82e8d51e71f6fd165af674971b975367b8b20579b1482e
if [ "$BORDER" = 1 ]; then
  nice -n 10 python3 "$HERE/border8.py" --out "$OUT/border"
  want "$(sha "$OUT/border/border_k8_pattern2.json")" "$(sha "$HERE/data/border_k8_pattern2.json")" "regenerated border_k8_pattern2.json is byte-identical to data/"
fi
echo "reproduce.sh: ALL OK"
