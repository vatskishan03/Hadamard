#!/bin/bash
# validate_extra.sh : exhaustive class counts at the multiplier structures of the Dokovic-Kotsireas pairs (arXiv:1310.5773, Sect. 3), appended to logs/validation_log.txt.
# DK list: 74 (H of order 3: 2 solutions), 82 (order 5: 2), 122 (order 5: 1), 164 (order 5: 3), 202 (order 5: 1), 226 (order 7: 2) -- "non-equivalent" in their sense.
cd "$(dirname "$0")/.."
{
echo; echo "== 10. exhaustive class counts (this engine, all H-invariant pairs with the stated |H|), next to the numbers of solutions listed by Dokovic-Kotsireas"
for wh in "37 3 74(DK:2)" "41 5 82(DK:2)" "61 5 122(DK:1)" "101 5 202(DK:1)" "113 7 226(DK:2)" "61 3 122" "73 3 146" "109 9 218" "181 9 362"; do
  set -- $wh
  r=$(python3 tools/solve.py $1 $2 work/sc_$1_$2 --swap-reduce | grep "equivalence classes" | sed -E 's/.*found representatives: ([0-9]+) .*/\1/')
  echo "(w,h)=($1,$2) length $3 : $r equivalence class(es)"
done
} | tee -a logs/validation_log.txt
