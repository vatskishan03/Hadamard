"""control of the scan (in search/skewplug: python3 scan38/ctl539.py): 539 = Z7 x Z7 x F11 with the scalars <(2,2)> of order 3 (the shape of skew/2156), all slots of the top energy window; output in est/.  New in this port (the control of the notes, with the K-level shape given as a matrix group)."""
import sys, json, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import scan38
rows = json.load(open(os.path.join(HERE, 'data', 'ctl539_rows.json')))
info = scan38.run_shape(0, rows[0], 600, 5, 60)
print(json.dumps(info))
