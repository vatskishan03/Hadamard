"""Write skew/<order>/README.md from the certificate.

    python3 search/skew_readme.py ORDER M "search sentence" [COMMIT "command" "randomness and caps" "result"]

M is the exponent listed for v = ORDER/4 in Table 5 of Cati-Pasechnik (smallest m with a known skew-Hadamard
matrix of order 2^m v); M = 0 for v > 999, where the Djokovic-Seberry census is cited instead."""
import json, sys
from pathlib import Path
from sympy import factorint
ROOT = Path(__file__).resolve().parent.parent
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def main():
    order, m, search = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    repro = sys.argv[4:8] if len(sys.argv) >= 8 else None
    d = ROOT / "skew" / str(order)
    c = json.loads((d / "certificate.json").read_text())
    v = c["v"]; gens = c.get("multipliers", [c.get("multiplier")])
    f = factorint(v)
    fac = " · ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))
    H = {1}
    while True:
        g2 = H | {h * g % v for h in H for g in gens}
        if g2 == H: break
        H = g2
    H = sorted(H)
    Hs = "⟨" + ", ".join(map(str, gens)) + "⟩"
    rs = c["rowsums"]
    sq = " + ".join(f"{r}²" for r in rs)
    mu = c.get("spin_multiplier")
    spin_par = ""
    if mu:
        spin_par = (f"**Spin multiplier.** μ = {mu} has order 3 modulo {v} (μ³ ∈ H). It fixes x₁, and the other blocks satisfy "
                    f"x₃(t) = x₂({mu}t) and x₄(t) = x₂({mu}²t).\n\n")
    if m > 0:
        status_block = (f"To our knowledge this is the first skew-Hadamard matrix of order {order}.\n\n"
                        f"- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, "
                        f"lists {v}({m}). So the smallest known skew order of the form 2ᵐ · {v} was {2 ** m * v}, and {order} was unknown.\n"
                        f"- A literature search in October 2026 found no later construction.\n")
        if m > 3:
            status_block += ("\nThe doubling construction (Seberry) then also gives skew-Hadamard matrices of orders "
                             + ", ".join(str(2 ** k * v) for k in range(3, m)) + ".\n")
    else:
        status_block = (f"To our knowledge this is the first skew-Hadamard matrix of order {order}.\n\n"
                        f"- Cati–Pasechnik Table 5 stops at v = 999. The Đoković–Seberry census of skew orders 4v for odd v < 2500 "
                        f"(May 2025) lists {v} among its 737 undecided values.\n"
                        f"- A literature search in October 2026 found no later construction.\n")
    vsub = str(v).translate(SUB)
    text = f"""# Skew-Hadamard matrix of order {order}

{order} = 4 · {v}{"" if len(f) == 1 and list(f.values())[0] == 1 else f", with {v} = {fac}"}{" (prime)" if len(f) == 1 and list(f.values())[0] == 1 else ""}.

## Construction

A skew-Hadamard matrix is a Hadamard matrix H with H + Hᵀ = 2I.

**Difference family.** Four ±1 sequences x₁, …, x₄ on ℤ{vsub} satisfy two conditions:

- each is constant on the orbits of the multiplier subgroup
  H = {Hs} ≤ ℤ{vsub}*, of order {len(H)}{":" if len(H) <= 40 else ""}{" {" + ", ".join(map(str, H)) + "}" if len(H) <= 40 else ""};
- their periodic autocorrelations sum to zero at every nonzero shift.

`certificate.json` lists, for each sequence, one representative r of each orbit rH on which the sequence is −1.

**Skew first block.** The first sequence is of skew type: x₁(0) = 1 and x₁(−t) = −x₁(t) for t ≠ 0. This is possible because −1 ∉ H, so negation pairs the orbits.

{spin_par}**Row sums.** The row sums are {", ".join(map(str, rs))}, and {sq} = {order}.

**The matrix.** Let A, B, C, D be the circulant matrices with these first rows (`sequences.txt`) and R the back-diagonal identity. The Goethals–Seidel array

```
    [  A     BR     CR     DR  ]
H = [ -BR    A     DᵀR   -CᵀR ]
    [ -CR   -DᵀR    A     BᵀR ]
    [ -DR    CᵀR  -BᵀR    A   ]
```

is then a Hadamard matrix of order {order}. Because A + Aᵀ = 2I, it satisfies H + Hᵀ = 2I (Wallis–Whiteman).

## Search

{search}

{{repro_section}}## Status

{status_block}
## Verification

```
python3 verify/verify.py skew/{order}
gzip -dc skew/{order}/H{order}.csv.gz | ./check_matrix
```

SHA-256 of `H{order}.csv`:
`{c["sha256"]}`
"""
    repro_section = ""
    if repro:
        commit, cmd, caps, result = repro
        repro_section = (f"## Reproduction\n\n| | |\n|---|---|\n| Code | commit `{commit}` |\n| Command | `{cmd}` |\n"
                         f"| Randomness and caps | {caps} |\n| Result of that run | {result} |\n\n"
                         "The tools write hit logs, temporary files and compiled helpers to `$HADAMARD_WORK` (default `<tmp>/hadamard-work`). "
                         "The certificate does not depend on the search: `verify/verify.py` rebuilds and checks the matrix from `certificate.json` alone.\n\n")
    text = text.replace("{repro_section}", repro_section)
    (d / "README.md").write_text(text)
    print(f"wrote {d / 'README.md'}")


if __name__ == "__main__":
    main()
