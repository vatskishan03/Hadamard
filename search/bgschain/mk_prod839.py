#!/usr/bin/env python3
"""Write the n = 839 production configs configs/n839_P{1,2,3}.json (G = Z4^3 x Z13, k = 8; the three live patterns; cap 2000 on the last level, 40 children per node in random
order, exact mod-8 screen)."""
import json, sys, os
import mkcfg2
HERE = os.path.dirname(os.path.abspath(__file__))
PATS = {'P1': ([2, 4, 4, 4], [4, -2, -2, 2]), 'P2': ([2, 4, 4, 4], [-4, -2, 2, 2]), 'P3': ([0, 4, 4, 4], [2, -2, 2, 4])}
def make(nm, cap=2000, child_cap=40, bdepth=2, tower='Z4,2Z4,Odd'):
    a, d = PATS[nm]; A = [a[i] + d[i] for i in range(4)]
    cfg = mkcfg2.build([13], [[3]], A, a, tower=tower, split='', bdepth=bdepth, name=f'n839_{nm}')
    cfg['normalizer'] = json.load(open(os.path.join(HERE, 'configs', 'normalizer_839.json')))
    cfg['cap'] = [0, 0, cap]; cfg['child_cap'] = [0, 0, child_cap]
    cfg['maxfibre'] = 4_000_000; cfg['maxstates'] = 80_000_000; cfg['dp_cache_states'] = 3_000_000
    cfg['pattern'] = dict(a=a, d=d, A=A, rowsums=[8 * a[i] + d[i] for i in range(4)])
    cfg['screen_mods'] = [[2, 3]]; cfg['score'] = True; cfg['leaflog'] = True; cfg['child_order'] = 'random'
    return cfg
if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'configs')      # output directory (default: the shipped configs/)
    for nm in PATS:
        cfg = make(nm); json.dump(cfg, open(os.path.join(out, f'n839_{nm}.json'), 'w'), indent=1); print(nm, cfg['pattern'])
