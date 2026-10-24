"""control of the scan (in search/skewplug: python3 scan38/ctl119.py); output in est_ctl/ or est/"""
import sys, json
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import scan38
rows=json.load(open(os.path.join(HERE, 'data', 'ctl119_rows.json')))
for i,r in enumerate(rows):
    info=scan38.run_shape(i, r, 150, 5, 60, outdir='est_ctl')
    print(json.dumps({k:info.get(k) for k in ('name','nu','Emin','Emax','pool','processed','status')}), flush=True)
