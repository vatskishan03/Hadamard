"""The Table 5 targets of the octic class: entries v(m) with m >= 5 (Table 5 of Cati and Pasechnik, ../t5.json) for which skew/ has no matrix of
order 2^k v with k <= 4 (a skew order 16v is the smallest the octic ring gives, so m <= 4 is nothing to improve).  Of the 74 entries with m >= 5,
13 are already improved by the repository (191, 277 -> 3; 377, 403, 463, 601, 933 -> 2; 579 -> 3; 589, 611, 817, 897, 899 -> 4): 61 targets remain.
usage: python3 targets.py      writes runs/octic_targets.json (rows (v, m, k) with k the best exponent known or null) and prints the lists"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
T5 = {int(a): b for a, b in json.load(open(os.path.join(HERE, '..', 't5.json'))).items()}
ents = sorted(T5.items())
# already-known skew orders in the repository (skew/ dirs): v odd part, k
have={}
for d in os.listdir(os.path.join(HERE, '..', '..', '..', 'skew')):
    if d.isdigit():
        N=int(d); k=0
        while N%2==0: N//=2; k+=1
        have[N]=min(have.get(N,99),k)
# plus hadamard/3732 skew
have[933]=min(have.get(933,99),2)
tg=[]
for v,m in ents:
    if m<5: continue
    k=have.get(v,None)
    if k is not None and k<=4: continue
    tg.append((v,m,k))
print(len(tg),'targets (m>=5 and not already improved to <=4 by skew/):')
print([v for v,m,k in tg])
print('already improved (m>=5 in Table 5, repo has skew order 2^k v, k<=4):',[(v,m,have[v]) for v,m in ents if m>=5 and v in have and have[v]<=4])
os.makedirs('runs',exist_ok=True)
json.dump(tg,open('runs/octic_targets.json','w'))
