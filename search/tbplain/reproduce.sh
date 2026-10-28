#!/bin/bash
# Reproduction of README.md section 5: build, small tests, the hit of QH(358), the hit of QH(498), each ending in the comparison of the
# regenerated certificate with complex-hadamard/C<n>.json and in verify/verify_complex.py run on the regenerated certificate.
#   bash reproduce.sh [OUT] [tests|358|358alt|498|sym|all]    OUT defaults to ./out, the stage to all; REPO = repository root (default ../..)
#   358alt compares its hit with the one of the 358 stage (OUT/p358); the other stages are independent of each other
# One core, nice, 3 GB of address space.  Timings: README section 5.
set -e
HERE=$(cd "$(dirname "$0")" && pwd); OUT=$(mkdir -p "${1:-$HERE/out}" && cd "${1:-$HERE/out}" && pwd); WHAT=${2:-all}
REPO=$(cd "${REPO:-$HERE/../..}" && pwd)
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RUN="nice -n 12"
T0=$SECONDS; step() { echo; echo "== $* [$((SECONDS - T0)) s]"; }

# the build commands of README.md section 5 (a binary newer than its source is kept)
build358() { (cd "$HERE/n358"
  [ pjoin -nt pjoin.c ] || gcc -O2 -march=native -o pjoin pjoin.c
  [ tower -nt tower.c ] || gcc -O2 -march=native -o tower tower.c); }
build498() { (cd "$HERE/n498"
  [ band -nt band.c ] || gcc -O2 -march=native -o band band.c -lm
  [ mlib.so -nt mlib.c ] || gcc -O2 -march=native -shared -fPIC -o mlib.so mlib.c -lm); }

# verify_complex.py is run on the regenerated certificate: a copy of the repository's verifier sits next to it (the script reads
# ../complex-hadamard/C<n>.json relative to itself); then the certificate is compared with the one in the repository.
check_cert() {
  local n=$1 d=$OUT/c$1
  mkdir -p "$d/verify"; cp "$REPO/verify/verify_complex.py" "$d/verify/"
  step "verify_complex.py on the regenerated certificate C$n.json"
  python3 "$d/verify/verify_complex.py" $n
  python3 - "$d/complex-hadamard/C$n.json" "$REPO/complex-hadamard/C$n.json" <<'PY'
import json, sys
a, b = (json.load(open(f)) for f in sys.argv[1:3])
assert a["sha256"] == b["sha256"], "matrix digest differs from the repository"
assert a == b, "certificate differs from the repository"
print("regenerated C%d.json equals the repository certificate: sha256 %s" % (a["order"], a["sha256"]))
PY
  cmp "$d/complex-hadamard/C$n.json" "$REPO/complex-hadamard/C$n.json" && echo "C$n.json is byte-identical to complex-hadamard/C$n.json"
}

tests() {
  step "build"; build358; build498
  cd "$HERE/n358"
  step "n358/structure.py"; $RUN python3 structure.py
  step "n358/test_pjoin.py (pjoin against numpy brute force)"; $RUN python3 test_pjoin.py
  step "n358/small_analogue.py (N = 51: exhaustive cores, border, QH(102))"; $RUN python3 small_analogue.py --out "$OUT/c102"
  mkdir -p "$OUT/c102/verify"; cp "$REPO/verify/verify_complex.py" "$OUT/c102/verify/"; python3 "$OUT/c102/verify/verify_complex.py" 102
  step "n358/test_tower.py (planted control of tower.c at real size)"; $RUN python3 test_tower.py 3e4 0
  cd "$HERE/n498"
  step "n498/mcomp.py (M-part quadratic forms against the exact PAF)"; $RUN python3 mcomp.py
  step "n498/test_mlib.py (mlib.c against numpy)"; $RUN python3 test_mlib.py
  step "n498/small_analogue.py (G = Z2 x Z3: brute-force cores, border, QH(18))"; $RUN python3 small_analogue.py --out "$OUT/c18"
  mkdir -p "$OUT/c18/verify"; cp "$REPO/verify/verify_complex.py" "$OUT/c18/verify/"; python3 "$OUT/c18/verify/verify_complex.py" 18
  step "n498/band PLANTM=1 (planted M completions of real P-tuples, band 20:99, s0 in [-12, 12])"
  mkdir -p "$OUT/plant498"; rm -f "$OUT"/plant498/plant.*
  PLANTM=1 CBAND=20:99 DBAND=20:99 MAXB=20000000 $RUN ./band 41 "$OUT/plant498/plant" -12 12 2> "$OUT/plant498/plant.stderr"; grep PLANT "$OUT/plant498/plant.prog"
}

r358() {
  step "build"; build358
  cd "$HERE/n358"
  step "QH(358): R-tuple rank 1 (seed 2), (R,w)-tuple rank 1291 only"
  rm -rf "$OUT/p358"
  $RUN python3 produce.py 18,0,14,14 1 1e5 3000 "$OUT/p358" 2 --only 1291
  step "QH(358): core -> border -> K"
  $RUN python3 "$HERE/assemble.py" 358 "$OUT/p358/hits/hit_r1_t1291_0.json" --out "$OUT/c358" --text
  check_cert 358
}

# the same hit selected by explicit indices (robust against a different NumPy sort order): R-tuple index 19593, class indices of stage 1
r358alt() {
  step "build"; build358
  cd "$HERE/n358"
  step "QH(358): the same hit by --ridx 19593 --quad 3889,3410,63496,63837"
  rm -rf "$OUT/p358alt"
  $RUN python3 produce.py 18,0,14,14 1 1e5 3000 "$OUT/p358alt" 2 --ridx 19593 --quad 3889,3410,63496,63837
  python3 - "$OUT/p358/hits/hit_r1_t1291_0.json" "$OUT/p358alt/hits/hit_r1_t1291_0.json" <<'PY'
import json, sys
a, b = (json.load(open(f)) for f in sys.argv[1:3])
assert a["seqs"] == b["seqs"] and a["quad"] == b["quad"], "the two selections give different hits"
print("the hit of --only 1291 and the hit of --ridx/--quad are the same core")
PY
}

# the 32-fold symmetry reduction of band.c on band 20:99 (s0 in [-12, 12]): all 672 P-tuples without reduction against the 21 canonical ones
sym498() {
  step "build"; build498
  cd "$HERE/n498"; mkdir -p "$OUT/sym498"; rm -f "$OUT"/sym498/*
  step "n498/band without symmetry reduction (DUMPT=1 NOCANA NOCANB NOREFL NODEDUP)"
  DUMPT=1 NOCANA=1 NOCANB=1 NOREFL=1 NODEDUP=1 CBAND=20:99 DBAND=20:99 MAXB=20000000 $RUN ./band 41 "$OUT/sym498/full" -12 12 2> "$OUT/sym498/full.stderr"
  step "n498/band with the reductions (DUMPT=1)"
  DUMPT=1 CBAND=20:99 DBAND=20:99 MAXB=20000000 $RUN ./band 41 "$OUT/sym498/keep" -12 12 2> "$OUT/sym498/keep.stderr"
  step "n498/check_sym.py"
  $RUN python3 check_sym.py "$OUT/sym498/full.sol" "$OUT/sym498/keep.sol"
}

r498() {
  step "build"; build498
  cd "$HERE/n498"
  step "QH(498): band n_C = n_D = 16, bucket s0 = -3"
  mkdir -p "$OUT/p498"; rm -f "$OUT"/p498/b16.*
  ./run_band.sh "$OUT/p498" b16 16:16 16:16 -3 -3
  tail -3 "$OUT/p498/b16.stderr"; head -c 600 "$OUT/p498/b16.prog"; echo; test -s "$OUT/p498/b16.sol"
  step "QH(498): SOL lines -> sequences over G, exact check"
  $RUN python3 verify_sol.py "$OUT/p498/b16.sol" --out "$OUT/p498/verified_solutions.json" | tail -2
  step "QH(498): core -> border -> K"
  $RUN python3 "$HERE/assemble.py" 498 "$OUT/p498/verified_solutions.json" --index 0 --out "$OUT/c498" --text
  check_cert 498
}

case $WHAT in
  tests) tests;;
  358) r358;;
  358alt) r358alt;;
  498) r498;;
  sym) sym498;;
  all) tests; r358; r358alt; r498; sym498;;
  *) echo "usage: reproduce.sh [OUT] [tests|358|358alt|498|sym|all]"; exit 2;;
esac
step "reproduce.sh $WHAT: ALL OK"
