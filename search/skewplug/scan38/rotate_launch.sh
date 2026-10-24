#!/bin/bash
# usage: rotate_launch.sh TAG DEADLINE_EPOCH 'PASSES' IDXLIST   (one rotation worker; pid in ../runs/rotate_TAG.pid, log ../runs/rotate_TAG.log)
cd "$(dirname "$0")"; mkdir -p ../runs; echo $$ > ../runs/rotate_$1.pid
exec python3 -u rotate.py "$2" "$3" "$4" >> ../runs/rotate_$1.log 2>&1
