# Paley plug-ins in the plain Szöllősi array

The sampler that found the quaternary Hadamard matrix of order 418 ([../../complex-hadamard/C418.json](../../complex-hadamard/C418.json);
[../NOTES.md](../NOTES.md) §85). It searches the plain (not zero-diagonal) Szöllősi array over G = K × F_q, K cyclic,
with A, B, C, D = ±1, A and B symmetric and Σ PAF = 4|G| δ, in the Paley plug-in class: each block is
s ⊗ I + t ⊗ Q + e ⊗ J with Q the Paley matrix of F_q and a few defects (positions of K where the block is constant on
F_q*). The K-level system and the defect identity e₊(q + 1) + e₋(q − 3) = 4|K| are as in §81–85.

| File | Contents |
|---|---|
| `l1.c`, `l1.py` | level 1: exact A, B lists against random C, D configurations of a given defect weight, class-sorted hash join on the trivial-F_q-character features |
| `pl2.c`, `l2io.py` | level 2: exact completion of the signs t off the defects (row-sum filter, spectral filter, meet in the middle); any odd q |
| `samp_l1.py`, `prod_worker.sh` | one level-1 and level-2 round per process; the restartable production loop |
| `certify.py` | a K-level solution to the four sequences on G, the exact checks, and the certificate C<2|G|>.json |
| `pk.py`, `probe.py`, `est.py`, `exact2.py`, `samp.py`, `idx.py`, `buildq.py`, `mitm.py`, `cmitm.c` | shapes, the K-level equations, exact small-case solvers and estimators |

Build (in this directory):

```
gcc -O3 -march=native -shared -fPIC -o libl1.so l1.c -lm
gcc -O3 -march=native -shared -fPIC -o libcmitm.so cmitm.c
gcc -O3 -march=native -o pl2_64 pl2.c -lm
gcc -O3 -march=native -DPS=32 -o pl2_32 pl2.c -lm
mkdir -p runs tmp
```

The run that found order 418 (two workers, about 0.6 CPU-hours to the first hit):

```
./prod_worker.sh p418w1 '[19]' 11 '[]' "1:1:2:3@1.2e7,1:1:2:3@1.2e7,1:1:1:4@1e7,1:1:3:3@1.2e7" 12000000 $(date +%s) p418
```

A hit is written to `runs/` as a pickle of K-level solutions; `python3 certify.py <pickle> 0 '[19]' 11 '[]'` writes
`out/C418.json`.
