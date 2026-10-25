#!/bin/sh
# Reproduce search/data/conference118_5fe566972f21790b.txt in a scratch directory (about 10 s).
set -e
SRC=$(cd "$(dirname "$0")" && pwd)
W=${1:-$(mktemp -d)}
mkdir -p "$W" && cd "$W"
cp "$SRC"/tmpl.c "$SRC"/inst.cpp "$SRC"/xc.cpp "$SRC"/templates.py "$SRC"/build.py .
cc -O2 -o tmpl tmpl.c && c++ -O2 -o inst inst.cpp && c++ -O2 -o xc xc.cpp
python3 templates.py
./inst
./xc 1 > xc.out
python3 build.py xc.out
cmp conference118_5fe566972f21790b.txt "$SRC"/../data/conference118_5fe566972f21790b.txt && echo "identical to search/data/conference118_5fe566972f21790b.txt"
