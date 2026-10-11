# Search tools

These are the programs that found the Goethals–Seidel orders in this repository.
[METHOD.md](../METHOD.md) gives the background.

| Program | Role |
|---|---|
| `triage.py` | ranks multiplier structures for a target n by predicted matching work |
| `generate.py` | generates candidate sequences by batched local search (JAX; runs on TPU, GPU or CPU) |
| `dedupe.c` | removes duplicate candidates |
| `match.c` | extended k-tree matcher; prints complementary quadruples |
| `certify.py` | writes `certificate.json`, `sequences.txt` and the matrix for a match |

## Build

```
pip install numpy jax
cc -O2 search/dedupe.c -o dedupe
cc -O3 -march=native -fopenmp search/match.c -o match -lm
```

## Example: order 4492

**1. Pick a structure.**

```
python3 search/triage.py 1123
```

The best row is H = ⟨60⟩ (order 33, 35 orbits, 17 classes) with row sums
(1, 1, 1, 67).

**2. Generate candidates for each row sum.**

```
python3 search/generate.py 1123 60 --rowsum 1  --count 160000000 --out r1.bin  --greedy 3 --polish 1 --width 4 --rounds 0
python3 search/generate.py 1123 60 --rowsum 67 --count 160000000 --out r67.bin --greedy 3 --polish 1 --width 4 --rounds 0
```

Each command reports the number of classes; for 1123 it is 17.

**3. Remove duplicates.**

```
./dedupe 17 r1.bin u1.bin
./dedupe 17 r67.bin u67.bin
```

**4. Match.** Roles A and D form the left pairs. B supplies both middle
sequences, so the repeated row sum goes in the second file. The target is
−(p mod 4).

```
./match -s 6 -t -3 u1.bin u1.bin u67.bin 17 > solutions.txt
```

**5. Write and check the certificate.**

```
python3 search/certify.py 1123 60 "$(head -1 solutions.txt)" 4492
```

## Notes

- **Candidate depth.** Small structures (pools of up to about 10¹⁰ sequences) need
  only light descent, as in the example. For orders like 2092, use the default
  settings, which add iterated local search. These run at about 1.3·10⁵ candidates
  per second on one TPU v6e chip for p = 523.
- **Memory.** `match` holds about 20 bytes per list element plus the pair lists for
  one value of τ. Lists of 10⁹ elements fit in 172 GB.
- **Choosing `-s`.** Choose the clamp `-s` so that the expected number of pairs per
  τ (printed at start-up) is close to the list size.
