#!/bin/bash
# validate_all.sh : the validation suite behind the README.  Writes logs/validation_log.txt.  About 100 CPU-seconds, dominated by (73,3).
set -u
cd "$(dirname "$0")/.."
LOG=logs/validation_log.txt
REPO=${HADAMARD_REPO:-$(cd ../.. && pwd)}
exec > >(tee "$LOG") 2>&1
echo "# validation log  $(date -u +%Y-%m-%dT%H:%MZ)  host: $(uname -sm)  cc: $(clang --version | head -1)"
echo "# repository checkout read from: $REPO (read-only); git HEAD $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null)"
mkdir -p work out
cc="clang -O3 -march=native -pthread -o"
$cc src/pgp src/pgp.c -lm && clang -O3 -o src/oracle src/oracle.c && echo "compiled src/pgp src/oracle"
T="python3 tools/timed.py"
echo; echo "== 1. self-tests: tables, parity keys, rotation covariance, Fourier tables, support index (against direct computation)"
for wh in "5 1" "13 3" "13 1" "37 3" "41 5" "61 5" "73 9" "113 7" "101 5" "97 3"; do ./src/pgp selftest $wh; done
echo; echo "== 2. rows (stage 1) against the brute-force oracle (all ternary H-invariant pairs, PAF from the definition on Z_w)"
for wh in "5 1" "41 5" "13 1" "37 3"; do set -- $wh; w=$1; h=$2
  ./src/oracle rows $w $h work/or_${w}_${h}.txt > /dev/null
  rm -rf work/e_${w}_${h} work/c_${w}_${h}
  ./src/pgp rows $w $h work/e_${w}_${h} --all-keys --no-rev --quiet > /dev/null
  echo -n "(w,h)=($w,$h) all keys, no reversal reduction: "; python3 tools/cmp_rows.py $w $h work/e_${w}_${h} work/or_${w}_${h}.txt
  ./src/pgp rows $w $h work/c_${w}_${h} --quiet > /dev/null
  echo -n "(w,h)=($w,$h) canonical keys + reversal reduction, closed under rotations/reversals: "; python3 tools/cmp_rows_canon.py $w $h work/c_${w}_${h} work/or_${w}_${h}.txt
done
echo; echo "== 3. Golay pairs (stage 1 + stage 2) against the brute-force oracle over all H-invariant +-1 pairs on Z_2w (ordered pairs)"
for wh in "5 1" "13 3" "41 5" "29 7" "73 9" "37 9"; do set -- $wh; w=$1; h=$2
  [ -f work/op_${w}_${h}.txt ] || ./src/oracle pgp $w $h work/op_${w}_${h}.txt > /dev/null
  for sr in "" "--swap-reduce"; do
    python3 tools/solve.py $w $h work/s_${w}_${h} $sr --oracle work/op_${w}_${h}.txt | grep -v "stage1 done\|^    \|class repr" | tr '\n' ' '; echo
  done
done
echo; echo "== 4. stage 2 against an independent Python support join on the oracle's complete row list, (13,1), m = 12"
python3 tools/count_pairs_oracle.py work/or_13_1.txt 12
./src/pgp join 13 1 work/e_13_1 2>/dev/null | python3 -c "
import sys
s=set(tuple(l.split()[:12]) for l in sys.stdin); print('engine (all keys, no reversal reduction): distinct (A,B,C,D) quadruples =',len(s))"
echo; echo "== 5. positive controls with known / expected results"
echo "-- (37,3): the repository's note says 35328 pairs for w = 37, H = <10> of order 3"
python3 tools/solve.py 37 3 work/s_37_3 --swap-reduce | grep -v "stage1 done\|^    \|class repr"
echo "-- the Dokovic-Kotsireas pair of length 74 (symmetric/596/certificate.json) among the images of the found solutions"
python3 tools/check_known.py 37 3 work/s_37_3 $REPO/symmetric/596/certificate.json --swap-reduce
echo "-- (113,7): Dokovic-Kotsireas found length 226 with |H| = 7"
python3 tools/solve.py 113 7 work/s_113_7 --swap-reduce | grep -v "stage1 done\|^    \|class repr"
echo "-- (61,3): length 122"
python3 tools/solve.py 61 3 work/s_61_3 --swap-reduce | grep -v "stage1 done\|^    \|class repr"
echo; echo "== 6. threads, checkpoint and resume: (61,3) in one run, then in two runs (first 100 buckets, then resume) with 4 threads"
rm -rf work/r_61_3 && ./src/pgp rows 61 3 work/r_61_3 --swap-reduce --threads 4 --quiet --limit 100 > /dev/null
echo "after the first run: $(wc -l < work/r_61_3/done.log) buckets done"
./src/pgp rows 61 3 work/r_61_3 --swap-reduce --threads 4 --quiet > /dev/null
echo "after the resume: $(wc -l < work/r_61_3/done.log) buckets done (single run: $(wc -l < work/s_61_3/done.log))"
python3 tools/solve.py 61 3 work/r_61_3 --swap-reduce --keep | grep "classes\|join matches"
echo; echo "== 7. the control (73,3): length 146 over Z_2 x Z_73, |H| = 3 (exhaustive; 4 threads)"
rm -rf work/s_73_3
$T rows_73_3_final ./src/pgp rows 73 3 work/s_73_3 --swap-reduce --threads 4 --quiet
python3 tools/solve.py 73 3 work/s_73_3 --swap-reduce --keep | grep -v "^    \|class repr\|stage1 done"
python3.11 tools/extract_pairs.py 73 3 work/s_73_3 --swap-reduce > out/pgp146_h3_classes.json
python3.11 -c "
import json; d=json.load(open('out/pgp146_h3_classes.json'))
print('PGP(146) classes:',len(d['classes']),' exact integer PAF check:',all(c['exact_ok'] for c in d['classes']),' numpy FFT check:',all(c['fft_ok'] for c in d['classes']))"
echo; echo "== 8. rows found at m = 32 (sample run on real (97,3) buckets) are exact rows"
python3 - <<'PY'
import sys; sys.path.insert(0,'tools')
from pgplib import *
types=parse_types(97,3,True); P=types["P"]
rows=read_rows('work/t_97_3b',types); bad=0
for (t,xi,xo,yi,yo) in rows:
    X=P.expand(xi,xo); Y=P.expand(yi,yo); px=paf_cyc(X); py=paf_cyc(Y)
    bad += not (px[0]+py[0]==97 and all(px[s]+py[s]==0 for s in range(1,97)))
print("m=32 sample rows:",len(rows),"all exact rows (PAF_X+PAF_Y = 97 delta on Z_97):",bad==0)
PY
echo; echo "== 9. from a found pair to the symmetric Hadamard matrix and a certificate accepted by verify/verify.py (builder step)"
python3.11 tools/extract_pairs.py 37 3 work/s_37_3 --swap-reduce > out/pgp74_h3_classes.json
rm -rf out/cert_test_149 out/cert_test_101 out/cert_own_149 out/cert_own_293
echo "-- q = 149 from the repository's pair: must reproduce the repository's certificate"
python3.11 tools/make_cert.py --q 149 --from-cert $REPO/symmetric/596/certificate.json --out out/cert_test_149 --compare $REPO/symmetric/596/certificate.json | cut -c1-260
echo "-- q = 101 (group Z_2 x Z_5^2) from the repository's pair"
python3.11 tools/make_cert.py --q 101 --from-cert $REPO/symmetric/404/certificate.json --out out/cert_test_101 --compare $REPO/symmetric/404/certificate.json | cut -c1-260
echo "-- q = 149 from a pair found by this engine (class 5 of the 27 classes at (37,3))"
python3.11 tools/make_cert.py --q 149 --pair out/pgp74_h3_classes.json --class 5 --out out/cert_own_149 | cut -c1-260
echo "-- q = 293 (order 1172) from the new PGP(146) pair found by this engine (a test artefact: v = 293 is not census-undecided)"
python3.11 tools/make_cert.py --q 293 --pair out/pgp146_h3_classes.json --class 0 --out out/cert_own_293 | cut -c1-260
echo; echo "repository working tree changes made by this suite: $(git -C "$REPO" status --short | wc -l) (must be 0)"
echo; echo "== CPU ledger total: $(python3 - <<'PY'
tot=0
for ln in open('logs/cpu_ledger.txt'):
    if 'user=' in ln:
        d=dict(x.split('=') for x in ln.split()[1:]); tot+=float(d['user'])+float(d['sys'])
print(round(tot,1),'CPU-s (timed commands only)')
PY
)"
