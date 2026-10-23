#!/bin/bash
# Battery of the validation of the fast level 2 (skf2) against sl2 (val_skf.py): 22 slots of seven shapes
# (Z17 x F7 with H_K = 1: 7 slots, Z19 x F7 <7>, Z13 x F7 <3>, Z13 x F3 <3>: 4 slots each, Z7 x F11 <2>, Z43 x F11 <6>, Z35 x F7 <11>: 1 slot each).
# Every RESULT line must end with 'count/status mismatches 0, solution-set mismatches 0' (two slots print NOINST: no level-1 instance in the sample).
# Run from any directory (needs libgl1.so, sl2_*, skf2_* built in this directory); 7 to 30 minutes on one core.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[17]" 7 "[]" "[2, 1, 4, 6]" "[2, 0, 2, 5]" 20000 20000 20000 300 4) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[17]" 7 "[]" "[2, 0, 4, 7]" "[0, 0, 3, 6]" 20000 20000 20000 300 356) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[17]" 7 "[]" "[4, 2, 3, 5]" "[2, 1, 3, 5]" 20000 20000 20000 300 502) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[17]" 7 "[]" "[0, 1, 1, 11]" "[0, 1, 1, 7]" 20000 20000 20000 300 97) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[17]" 7 "[]" "[2, 0, 5, 7]" "[2, 0, 3, 6]" 20000 20000 20000 300 503) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[17]" 7 "[]" "[0, 3, 3, 10]" "[0, 3, 3, 9]" 20000 20000 20000 300 711) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[17]" 7 "[]" "[8, 0, 3, 6]" "[8, 0, 3, 6]" 20000 20000 20000 300 815) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[19]" 7 "[7]" "[0, 1, 3, 7]" "[0, 0, 3, 0]" 3000 10000 10000 400 267) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[19]" 7 "[7]" "[6, 0, 3, 10]" "[6, 0, 3, 10]" 3000 10000 10000 400 947) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[19]" 7 "[7]" "[6, 3, 4, 4]" "[4, 3, 4, 4]" 3000 10000 10000 400 797) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[19]" 7 "[7]" "[6, 0, 3, 7]" "[6, 0, 0, 7]" 3000 10000 10000 400 286) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[13]" 7 "[3]" "[0, 0, 3, 10]" "[0, 0, 3, 10]" 3000 10000 10000 400 467) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[13]" 7 "[3]" "[0, 1, 3, 4]" "[0, 0, 3, 0]" 3000 10000 10000 400 924) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[13]" 7 "[3]" "[0, 1, 1, 7]" "[0, 1, 1, 3]" 3000 10000 10000 400 494) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[13]" 7 "[3]" "[0, 1, 4, 4]" "[0, 1, 0, 4]" 3000 10000 10000 400 116) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[13]" 3 "[3]" "[12, 4, 9, 13]" "[12, 1, 3, 9]" 3000 10000 10000 300 629) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[13]" 3 "[3]" "[6, 3, 6, 10]" "[6, 3, 3, 0]" 3000 10000 10000 300 697) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[13]" 3 "[3]" "[6, 4, 7, 13]" "[6, 1, 7, 3]" 3000 10000 10000 300 892) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[13]" 3 "[3]" "[12, 4, 6, 9]" "[12, 0, 3, 3]" 3000 10000 10000 300 482) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[7]" 11 "[2]" "[0, 1, 1, 1]" "[0, 0, 1, 1]" 3000 10000 10000 300 847) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[43]" 11 "[6]" "[6,3,4,7]" "[6,0,4,7]" 3000 30000 30000 1500 3) 2>&1 | tail -3
(ulimit -v 3500000; nice -n 6 python3 val_skf.py "[35]" 7 "[11]" "[2,3,3,10]" "[0,0,0,1]" 256 922 5036 400 4) 2>&1 | tail -3
