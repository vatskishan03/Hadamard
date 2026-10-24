"""cmp_logs.py -- compares logs of plain_spence3.py / plain_nc3.py with recorded ones, ignoring the time field 'secs'.
usage: cmp_logs.py [--subset] RECORDED NEW.log [NEW.log ...]      RECORDED = a directory (the log of the same name is taken) or a file
Default: the two logs have the same lines in the same order (records compared as dictionaries).  --subset: every record of NEW occurs in RECORDED (a shorter run of the same sampler).
Prints one line per log and 'k of n identical'."""
import sys, os, ast, re
args = [a for a in sys.argv[1:] if a != '--subset']; subset = '--subset' in sys.argv
rec, news = args[0], args[1:]
def load(f):
    out = []
    for l in open(f):
        l = l.rstrip('\n')
        if l.startswith('{'):
            d = ast.literal_eval(l); d.pop('secs', None); out.append(str(sorted(d.items())))
        else: out.append(re.sub(r", 'secs': [0-9.]+", '', l))
    return out
ok = 0
for f in news:
    r = rec if os.path.isfile(rec) else os.path.join(rec, os.path.basename(f))
    a, b = load(f), load(r)
    if subset: same = all(x in set(b) for x in a if x.startswith('['))
    else: same = (a == b)
    ok += same
    print('%-24s %d lines (recorded %d): %s' % (os.path.basename(f), len(a), len(b), 'IDENTICAL' if same else 'DIFFERENT'))
print('%d of %d identical' % (ok, len(news)))
