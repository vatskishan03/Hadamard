#!/usr/bin/env python3
"""Print the table of complex-hadamard/README.md from the certificates C<n>.json (order, N, factorisation,
multiplier group, row sums).

    python3 search/ssxx_table.py
"""
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "complex-hadamard"


def factor(n):
    out, p = [], 2
    while p * p <= n:
        while n % p == 0:
            out.append(p); n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def fmt_factor(n):
    f = factor(n)
    if len(f) == 1:
        return str(n)
    parts = []
    for p in sorted(set(f)):
        k = f.count(p)
        parts.append(f"{p}" + ("" if k == 1 else "^" + str(k)))
    return f"{n} = " + " · ".join(parts)


def main():
    rows = []
    for f in sorted(ROOT.glob("C*.json"), key=lambda p: int(p.stem[1:])):
        c = json.loads(f.read_text())
        N = c["N"]
        src = c.get("source")
        H = f"⟨{', '.join(map(str, c['multipliers']))}⟩" + (f", order {c['subgroup_order']}" if c.get("subgroup_order")
                                                           else "")
        r = c["rowsums"]
        rows.append(f"| {c['order']} | {fmt_factor(N)} | {src or H} | {r[0]}, {r[1]}; {r[2]}, {r[3]} | [{f.name}]({f.name}) |")
    print("| order 2N | N | multiplier group / source | row sums (a, b; c, d) | file |")
    print("|---:|---|---|---|---|")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
