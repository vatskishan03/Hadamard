# Skew Paley plug-ins in the Goethals–Seidel array

The sampler that found the skew-Hadamard matrix of order 476 ([../../skew/476](../../skew/476);
[../NOTES.md](../NOTES.md) §87). It searches skew Goethals–Seidel quadruples over G = K × F_q, q ≡ 3 (mod 4), invariant
under H = {1} × (squares of F_q*), the Paley products of [../NOTES.md](../NOTES.md) §15. Each block is
s ⊗ I + t ⊗ Q + e ⊗ J with Q the Paley matrix of F_q; the first block is of skew type (s and e odd under k ↦ −k off 0,
s(0) = 1, e(0) = 0, t even), B, C, D are free. §15 solved the case of at most five w-positions (defects); this code
handles many (nine to seventeen for 476).

| File | Contents |
|---|---|
| `sk.py` | the shape: orbits, the skew block, the K-level equations (I), (J), (Q), exact lists of the skew block, checks |
| `sk_pat.py` | allocations allowed by the defect identity (q + 1)e₊ + (q − 3)e₋ = 4\|K\|, their splits into the two defect types, exact row-sum patterns |
| `gl1.c`, `gl1.py` | level 1: the exact list of the skew block against sampled lists of B, C, D, class-restricted hash join on the trivial-F_q-character equation (U) |
| `gl2.c`, `sl2io.py`, `gl2io.py` | level 2: exact completion of the signs t off the defects, with the skew block's own t-orbits |
| `sk_samp.py`, `sk_worker.py` | one level-1 and level-2 round; the production loop (seeded rounds, resumable, stops at the first hit) |
| `sk_cert.py` | a K-level solution to four sequences over Z_m × Z_q, checked with the builders of `verify/verify.py`, written to `out/` |
| `cyclic.py` | the same certificate over the cyclic group Z_mq (the format of [../../skew/820](../../skew/820)) |
| `bf.py`, `sk_exh.py`, `validate_skew.py` | G-level and K-level brute force, the exhaustive two-level pipeline, the validation battery |
| `indep_check.py` | an independent numpy check of a certificate over Z₁₇ × Z₇ |
| `gsk.py`, `pk.py`, `est.py`, `idx.py`, `sk_filter.py` | shared shape code (from the real Goethals–Seidel and Szöllősi plug-in searches), an estimator, a level-1 prefilter |

Build (in this directory):

```
gcc -O3 -march=haswell -shared -fPIC -o libgl1.so gl1.c -lm
for ps in 32 64 128; do gcc -O3 -march=haswell -DPS=$ps -o sl2_$ps gl2.c -lm; done
mkdir -p runs tmp
```

The run that found order 476 (one core, about 40 s; the rounds are seeded, so it is deterministic):

```
python3 -u sk_worker.py w119a 1 900 "[17]" 7 "[[2,2,2,3]]" 2000 20000 20000 30 1.0
python3 cyclic.py out/skew_476_w119a.json out/skew_476_cyclic.json
```

The arguments are the worker name, the seed base, a time limit in seconds, K, q, the defect allocations (or `E9,10,11`
for every allocation of those total weights), the list sizes of B, C, D, the number of level-1 instances passed to
level 2 per round, and the fraction of them kept. The first command writes `out/skew_476_w119a.json` (group form over
Z₁₇ × Z₇); the second gives the certificate of [../../skew/476](../../skew/476) without its `sha256` field.

`python3 validate_skew.py` re-runs the validation battery (31 shapes; G-level brute force, K-level brute force and the
two-level pipeline agree).
