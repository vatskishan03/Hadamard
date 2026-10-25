import itertools, numpy as np, pickle, sys
p=29; n=8
qr={x*x%p for x in range(1,p)}
chi=lambda x: 0 if x%p==0 else (1 if x%p in qr else -1)
hits=pickle.load(open('k8hits_29.pkl','rb'))
sq=sorted(qr)
allsets=set()
for h in hits:
    for a in sq:
        for b in range(p):
            allsets.add(tuple(sorted((a*x+b)%p for x in h)))
print("distinct 8-subsets with the template spectrum:",len(allsets))
cls=pickle.load(open('k8_classes.pkl','rb')); (target,)=cls
def canon_key(pts,plus):
    minus=[x for x in pts if x not in plus]
    best=None
    for p1 in itertools.permutations(plus):
        for p2 in itertools.permutations(minus):
            perm=list(p1)+list(p2)
            key=tuple(chi(perm[b]-perm[a])*chi(perm[c]-perm[b])*chi(perm[c]-perm[a]) for a,b,c in itertools.combinations(range(n),3))
            if best is None or key<best: best=key
    return best
# quick: triple-product function must match; check labelings on orbit representatives only (one per affine orbit)
reps={}
for s in allsets:
    can=min(tuple(sorted((a*x+b)%p for x in s)) for a in sq for b in range(p))
    reps.setdefault(can,[]).append(s)
print("affine orbits:",len(reps))
good={}
for can in reps:
    labs=[]
    for plus in itertools.combinations(can,4):
        if canon_key(can,list(plus))==target: labs.append(plus)
    good[can]=labs
    print(can, "valid D-labelings:",len(labs))
pickle.dump((reps,good),open('k8place_29.pkl','wb'))
