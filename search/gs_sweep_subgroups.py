from math import gcd
def subgroups(n):
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    def closure(gens):
        S, todo = {1}, [1]
        while todo:
            x = todo.pop()
            for g in gens:
                y = x * g % n
                if y not in S: S.add(y); todo.append(y)
        return frozenset(S)
    cyc = {}
    for u in units: cyc.setdefault(closure([u]), u)
    subs = {H: (g,) for H, g in cyc.items()}
    items = list(cyc.items())
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            subs.setdefault(closure([items[i][1], items[j][1]]), (items[i][1], items[j][1]))
    return subs

