"""Write skew/<order>/README.md from the certificate.

    python3 attack/skew_readme.py ORDER M "search sentence"

M is the exponent listed for v = ORDER/4 in Table 5 of Cati-Pasechnik (smallest m with a known skew-Hadamard
matrix of order 2^m v)."""
import json, sys
from pathlib import Path
from sympy import factorint
ROOT = Path(__file__).resolve().parent.parent
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def main():
    order, m, search = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
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

**Row sums.** The row sums are {", ".join(map(str, rs))}, and {sq} = {order}.

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

## Status

To our knowledge this is the first skew-Hadamard matrix of order {order}.

- Cati and Pasechnik, *A database of constructions of Hadamard matrices*, arXiv:2411.18897v2 (Aug. 2025), Table 5, lists {v}({m}). So the smallest known skew order of the form 2ᵐ · {v} was {2 ** m * v}, and {order} was unknown.
- A literature search in October 2026 found no later construction.

The doubling construction (Seberry) then also gives skew-Hadamard matrices of orders {", ".join(str(2 ** k * v) for k in range(3, m))}{"" if m > 3 else ""}.

## Verification

```
python3 verify/verify.py skew/{order}
gzip -dc skew/{order}/H{order}.csv.gz | ./check_matrix
```

SHA-256 of `H{order}.csv`:
`{c["sha256"]}`
"""
    if m <= 3:
        text = text.replace("\nThe doubling construction (Seberry) then also gives skew-Hadamard matrices of orders .\n", "")
    (d / "README.md").write_text(text)
    print(f"wrote {d / 'README.md'}")


if __name__ == "__main__":
    main()
