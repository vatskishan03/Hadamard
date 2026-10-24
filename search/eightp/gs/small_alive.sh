#!/bin/bash
# exact-pattern production for the tau-alive shapes with <= 11 orbits (EXHAUSTED verdicts when all lists are small)
cd "$(dirname "$0")"; mkdir -p runs; export OMP_NUM_THREADS=1; ulimit -v 3500000
[ -e runs/sa_1198_0.done ] || (GS_PCACHE=runs/sa_1198_0_pc.pkl GS_EXACT_CAP=3e4 nice -n 5 timeout 420 python3 prod_tau.py spence '[3,3,7]' 19 '[[[2,1,1],1],[[1,2,1],1],[[1,1,3],1]]' 360 7 runs/sa_1198_0 2000 200000 > runs/sa_1198_0.log 2>&1; touch runs/sa_1198_0.done)
[ -e runs/sa_1198_1.done ] || (GS_PCACHE=runs/sa_1198_1_pc.pkl GS_EXACT_CAP=3e4 nice -n 5 timeout 420 python3 prod_tau.py spence '[9,19]' 7 '[[[2,1],1],[[1,4],1]]' 360 8 runs/sa_1198_1 2000 200000 > runs/sa_1198_1.log 2>&1; touch runs/sa_1198_1.done)
[ -e runs/sa_1198_2.done ] || (GS_PCACHE=runs/sa_1198_2_pc.pkl GS_EXACT_CAP=3e4 nice -n 5 timeout 420 python3 prod_tau.py spence '[9,19]' 7 '[[[2,1],-1],[[1,4],1]]' 360 9 runs/sa_1198_2 2000 200000 > runs/sa_1198_2.log 2>&1; touch runs/sa_1198_2.done)
[ -e runs/sa_1198_3.done ] || (GS_PCACHE=runs/sa_1198_3_pc.pkl GS_EXACT_CAP=3e4 nice -n 5 timeout 420 python3 prod_tau.py spence '[7,19]' 9 '[[[2,1],1],[[1,4],1]]' 360 10 runs/sa_1198_3 2000 200000 > runs/sa_1198_3.log 2>&1; touch runs/sa_1198_3.done)
[ -e runs/sa_1198_4.done ] || (GS_PCACHE=runs/sa_1198_4_pc.pkl GS_EXACT_CAP=3e4 nice -n 5 timeout 420 python3 prod_tau.py spence '[7,9]' 19 '[[[2,1],1],[[1,2],1]]' 360 11 runs/sa_1198_4 2000 200000 > runs/sa_1198_4.log 2>&1; touch runs/sa_1198_4.done)
[ -e runs/sa_1678_5.done ] || (GS_PCACHE=runs/sa_1678_5_pc.pkl GS_EXACT_CAP=3e4 nice -n 5 timeout 420 python3 prod_tau.py spence '[13,43]' 3 '[[[4,1],1],[[1,9],1]]' 360 12 runs/sa_1678_5 2000 200000 > runs/sa_1678_5.log 2>&1; touch runs/sa_1678_5.done)
[ -e runs/sa_1198_6.done ] || (GS_PCACHE=runs/sa_1198_6_pc.pkl GS_EXACT_CAP=3e4 nice -n 5 timeout 420 python3 prod_tau.py spence '[7,9]' 19 '[[[3,1],1],[[1,4],1]]' 360 13 runs/sa_1198_6 2000 200000 > runs/sa_1198_6.log 2>&1; touch runs/sa_1198_6.done)
[ -e runs/sa_1198_7.done ] || (GS_PCACHE=runs/sa_1198_7_pc.pkl GS_EXACT_CAP=3e4 nice -n 5 timeout 420 python3 prod_tau.py spence '[7,9]' 19 '[[[3,1],1],[[1,8],1]]' 360 14 runs/sa_1198_7 2000 200000 > runs/sa_1198_7.log 2>&1; touch runs/sa_1198_7.done)
echo finished > runs/small_alive.finished
