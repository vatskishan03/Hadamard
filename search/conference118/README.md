# A symmetric conference matrix of order 118 by an exact cover

The programs behind the conference matrix `search/data/conference118_5fe566972f21790b.txt` ([../NOTES.md](../NOTES.md)
§95), which gives [symmetric/940](../../symmetric/940) and [hadamard/11292](../../hadamard/11292). A symmetric
conference matrix C of order 118 (zero diagonal, ±1 elsewhere, C² = 117I) is the same as a conference graph
SRG(117, 58, 28, 29). Its existence was open in the literature we read; Brouwer's table now cites J. Y. Z. Ho (email,
October 2026) for one, which is not public. Ours was found independently.

## The ansatz

Points: ℤ₁₃ × AG(2, 3) and one point ∞. C = [[0, eᵀ], [e, Q]], with Q of order 117 made of 13 × 13 blocks of order 9:

- diagonal blocks J − I − 2R, R the graph of the two directions (1, 1), (1, 2) of AG(2, 3) (degree 4);
- the block (i, j), i ≠ j, equal to −χ(i − j)(J − 2M), χ the Legendre symbol mod 13 and M the incidence matrix
  {(x, y) : ⟨a, x⟩ + ⟨b, y⟩ = c} of a hyperplane of AG(2, 3)², with directions a, b and constant c.

On the characters ψ of ℤ₃² the condition Q² = 117I − J splits. The trivial character gives −3 times Paley's core of
order 13, whose square is 117I − 9J, as needed. A nontrivial ψ belongs to one of the four directions of AG(2, 3); on
the two characters of a direction, the diagonal block acts as 3D with D = +1 for the directions (0, 1), (1, 0) and
D = −1 for (1, 1), (1, 2), and the block (i, j) maps the direction a at i to the direction b at j with modulus 6 and a
cube-root phase. Call (i, direction) a flag: 52 flags, each of degree 3 (twelve blocks at i over four flags), so the
blocks define a cubic graph on the flags, and Q² = 117I − J becomes, on each of its components, W² = 13I for
W = D + 2Φ over ℤ[ω], Φ the signed phases of the edges. Each point carries two flags with D = +1 and two with D = −1.

The search therefore has two levels:

1. `tmpl.c`: for each small cubic template (K₄, the prism, K₃,₃, the cube) the classes (edge signs S, diagonal D) for
   which phases with W² = 13I exist: 48, 64, 64 and 256 (`templates.py` writes them to `templates.txt`).
2. `inst.cpp`: every placement of a template class on points of ℤ₁₃ whose edge signs are the Paley signs χ(p_v − p_u),
   with distinct point pairs and each (point, D) used at most twice (585 for K₄, 17160 prisms, 624 K₃,₃, 196560
   cubes). `xc.cpp`: an exact cover of the 78 pairs of points of ℤ₁₃ by such instances, each (point, D) used at most
   twice, hence exactly twice, since a point lies on twelve pairs and a flag covers three (the pair with fewest
   candidates first).

With K₄ components alone the first cover is found after 13,589 nodes (0.05 s). Thirteen K₄ that cover the 78 pairs of
ℤ₁₃ once are the lines of PG(2, 3). `build.py` then chooses, for each component, the phases (a tree gauge and
exhaustion over 3³), assigns the four flag classes at each point to the four directions, builds the blocks, checks
C = Cᵀ, zero diagonal and C² = 117I exactly, and normalises C (∞ first, first row and column +).

## Reproduction

```
search/conference118/run.sh [WORKDIR]
```

compiles the three programs, runs the pipeline in WORKDIR (a new temporary directory by default) and compares the result
with `search/data/conference118_5fe566972f21790b.txt` (identical; about 10 s). There is no randomness. The signed
automorphism group of the matrix is {±I} (checked in §95).

| File | Contents |
|---|---|
| `tmpl.c` | the template classes: for a cubic graph, all (S, D, phases) with W² = 13I over ℤ[ω] |
| `templates.py` | runs `tmpl` on K₄, the prism, K₃,₃ and the cube; writes `templates.txt` |
| `inst.cpp` | the instances on ℤ₁₃ with Paley signs; writes `insts.bin` (6.9 MB) |
| `xc.cpp` | the exact cover; `./xc 1` uses K₄ only (bitmask: 1 K₄, 2 prism, 4 K₃,₃, 8 cube) |
| `build.py` | the conference matrix from the first cover; checks and writes `conference118_<sha256 prefix>.txt` |
| `run.sh` | the whole pipeline |
