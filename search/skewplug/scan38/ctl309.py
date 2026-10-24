"""control of the scan (in search/skewplug: python3 scan38/ctl309.py); output in est_ctl/ or est/"""
import sys, json
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import scan38
rows=json.load(open(os.path.join(HERE, 'data', 'ctl309_rows.json')))
info=scan38.run_shape(0, rows[0], 120, 5, 60, wmode='lo')
print(json.dumps(info))
