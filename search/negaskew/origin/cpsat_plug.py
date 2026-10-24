"""CP-SAT model of the origin-defect plug-in problem at the K-orbit level (independent of qp.c).
 Unknowns: for each block X (A skew, B,C,D free) and orbit o >= 1: alpha[X][o] (x0 exponent), beta[X][o] (x1 exponent; x2 = x1 + 2, E-type).
 Origin data (S, tau, eps) enumerated outside from the cost pattern.  Equations (exact, lags = orbit representatives l_j, j with j < neg[j]):
   (d)  sum_X [P_s + q P_t](l_j) = 0 ;  (Q) sum_X [CC_ts - CC_st](l_j) = 0 ;  (J) sum_X Phi_X(l_j) = sum_X P_t(l_j)
 usage: python3 cpsat_plug.py m "[hK]" q A B C D [tlimit] [maxorigin]     (A,B,C,D = origin costs)
"""
import sys, os, json, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from ortools.sat.python import cp_model
from qshape import kshape
U = [1, 1j, -1, -1j]
def origin_options(cost, q, skew):
    """list of (x0e, x1e, x2e) exponent triples of the origin with the given cost (free blocks: x0e = 0)"""
    out = []
    for x0e in ([0, 3] if skew else [0]):
        for x1e in range(4):
            for x2e in range(4):
                x0, x1, x2 = U[x0e], U[x1e], U[x2e]
                if skew and abs(x2 - 1j * x1.conjugate()) > 1e-9: continue
                e = (x1 + x2) / 2
                c = (q - 1) * abs(e) ** 2 + 2 * (x0 * e.conjugate()).real
                if abs(c - cost) < 1e-9: out.append((x0e, x1e, x2e))
    return out
def build_and_solve(m, hK, q, origins, tlimit=60, seed=1, sym=True, log=False, fixA=None, fix=None):
    sh = kshape([m], hK); nu = sh['nu']; neg = sh['neg']; N = sh['N']
    model = cp_model.CpModel()
    # origin constants (complex, exact in halves)
    S0 = []; T0 = []; E0 = []
    for (x0e, x1e, x2e) in origins:
        x0, x1, x2 = U[x0e], U[x1e], U[x2e]
        e = (x1 + x2) / 2; t = (x1 - x2) / 2
        S0.append(x0 - e); T0.append(t); E0.append(e)
    alpha = [[None] * nu for _ in range(4)]; beta = [[None] * nu for _ in range(4)]
    oh_a = [[None] * nu for _ in range(4)]; oh_b = [[None] * nu for _ in range(4)]
    def onehot(v, name):
        ys = [model.NewBoolVar(f'{name}_{e}') for e in range(4)]
        model.Add(sum(ys) == 1)
        model.Add(v == sum(e * ys[e] for e in range(4)))
        return ys
    for X in range(4):
        for o in range(1, nu):
            alpha[X][o] = model.NewIntVar(0, 3, f'a{X}_{o}'); beta[X][o] = model.NewIntVar(0, 3, f'b{X}_{o}')
            oh_a[X][o] = onehot(alpha[X][o], f'A{X}_{o}'); oh_b[X][o] = onehot(beta[X][o], f'B{X}_{o}')
    # skew block A: x0(-k) = i conj x0(k): alpha' = 1 - alpha ;  x1(-k) such that x2(-k) = i conj x1(k): E-type x1(-k) = -x2(-k)... derive: x2 = x1+2 => x1(-k) = x2(-k)+2 = (i conj x1(k)) + 2 -> beta' = 1 - beta + 2 = 3 - beta
    for o in range(1, nu):
        p = neg[o]
        if o < p:
            for e in range(4):
                model.Add(oh_a[0][p][(1 - e) % 4] == oh_a[0][o][e]); model.Add(oh_b[0][p][(3 - e) % 4] == oh_b[0][o][e])
    # helper: product indicator of exponents  (u_a conj u_b) -> exponent (ea - eb) mod 4 ; returns list of 4 literals-expressions via auxiliary bools
    cache = {}
    def prod(ohx, ohy, key):
        # returns [expr_e] e=0..3: indicator that exponent(x) - exponent(y) == e mod 4
        if key in cache: return cache[key]
        zs = {}
        res = [[] for _ in range(4)]
        for ea in range(4):
            for eb in range(4):
                z = model.NewBoolVar(f'z_{key}_{ea}{eb}')
                model.AddBoolOr([ohx[ea].Not(), ohy[eb].Not(), z]); model.AddImplication(z, ohx[ea]); model.AddImplication(z, ohy[eb])
                res[(ea - eb) % 4].append(z)
        out = [sum(r) for r in res]
        cache[key] = out
        return out
    # scale: all coefficients doubled to stay integral: 2*S etc are Gaussian integers.
    def cx2(z): return (int(round(2 * z.real)), int(round(2 * z.imag)))
    # a complex linear expression = (re_expr, im_expr), integer-scaled by 2
    def add_c(acc, coef, ind_by_exp, sign_exp_shift=0):
        """acc += coef * i^(exp) summed over exps with indicators ind_by_exp[exp] (exp already includes shift)"""
        for e in range(4):
            u = U[e]
            cr, ci = coef.real, coef.imag
            zr = (coef * u)
            acc[0].append((int(round(2 * zr.real)), ind_by_exp[e])); acc[1].append((int(round(2 * zr.imag)), ind_by_exp[e]))
    def lin(terms): return sum(c * x for c, x in terms)
    # block values by orbit: s_o: ('var', X, o) for o>=1, ('const', S0[X]) for o=0
    def val_terms(X, o, kind):
        """one-hot representation of the unit value at orbit o (kind 's' or 't'): returns ('c', const) or ('v', oh)"""
        if o == 0: return ('c', S0[X] if kind == 's' else T0[X])
        return ('v', oh_a[X][o] if kind == 's' else oh_b[X][o])
    def bil(X, kx, a, ky, b, acc, weight):
        """acc += weight * x_a * conj(y_b), x of kind kx at orbit a, y of kind ky at orbit b"""
        vx = val_terms(X, a, kx); vy = val_terms(X, b, ky)
        if vx[0] == 'c' and vy[0] == 'c': 
            c = weight * vx[1] * vy[1].conjugate(); acc[0].append((int(round(2 * c.real)), 1)); acc[1].append((int(round(2 * c.imag)), 1)); return
        if vx[0] == 'c':   # const * conj(u_b): exponent -eb
            ind = [vy[1][(-e) % 4] for e in range(4)]
            for e in range(4):
                c = weight * vx[1] * U[e]; acc[0].append((int(round(2 * c.real)), ind[e])); acc[1].append((int(round(2 * c.imag)), ind[e]))
            return
        if vy[0] == 'c':   # u_a * conj(const)
            for e in range(4):
                c = weight * vy[1].conjugate() * U[e]; acc[0].append((int(round(2 * c.real)), vx[1][e])); acc[1].append((int(round(2 * c.imag)), vx[1][e]))
            return
        key = (X, kx, a, ky, b)
        ind = prod(vx[1], vy[1], key)
        for e in range(4):
            c = weight * U[e]; acc[0].append((int(round(2 * c.real)), ind[e])); acc[1].append((int(round(2 * c.imag)), ind[e]))
    reps = [j for j in range(1, nu) if j < neg[j]]
    for j in reps:
        Nj = N[j]
        accD = ([], []); accQ = ([], []); accJ = ([], [])
        for X in range(4):
            for a in range(nu):
                for b in range(nu):
                    n = Nj[a][b]
                    if not n: continue
                    bil(X, 's', a, 's', b, accD, n)
                    bil(X, 't', a, 't', b, accD, n * q)
                    bil(X, 't', a, 's', b, accQ, n); bil(X, 's', a, 't', b, accQ, -n)
                    bil(X, 't', a, 't', b, accJ, -n)     # -P_t on the J side
            # Phi_X(l) = eps conj(u(l)) + conj(eps) u(-l):  u(l)= s at orbit(l)=j ; u(-l) at neg[j]
            e = E0[X]
            if abs(e) > 1e-9:
                vj = val_terms(X, j, 's'); vn = val_terms(X, neg[j], 's')
                for ee in range(4):
                    c = e * U[(-ee) % 4]; accJ[0].append((int(round(2 * c.real)), vj[1][ee])); accJ[1].append((int(round(2 * c.imag)), vj[1][ee]))
                    c = e.conjugate() * U[ee]; accJ[0].append((int(round(2 * c.real)), vn[1][ee])); accJ[1].append((int(round(2 * c.imag)), vn[1][ee]))
        for acc in (accD, accQ, accJ):
            for part in (0, 1):
                # separate constants (indicator == 1 int) from vars
                const = sum(c for c, x in acc[part] if isinstance(x, int)); expr = sum(c * x for c, x in acc[part] if not isinstance(x, int))
                model.Add(expr + const == 0)
    # symmetry breaking: order of B,C,D by alpha vector (only if origins equal)
    if sym:
        for X, Y in ((1, 2), (2, 3)):
            if origins[X] == origins[Y]:
                # lexicographic on alpha then beta  (simple weighted sum comparison)
                wX = sum((4 ** (o - 1)) * alpha[X][o] for o in range(1, nu)); wY = sum((4 ** (o - 1)) * alpha[Y][o] for o in range(1, nu))
                model.Add(wX <= wY)
    if fix is not None:
        for X in range(4):
            for o in range(1, nu):
                model.Add(alpha[X][o] == fix[X][0][o])
                if fix[X][1] is not None: model.Add(beta[X][o] == fix[X][1][o])
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = tlimit; solver.parameters.num_search_workers = 1; solver.parameters.random_seed = seed
    solver.parameters.log_search_progress = log
    t0 = time.time(); st = solver.Solve(model); dt = time.time() - t0
    res = dict(status=solver.StatusName(st), time=round(dt, 2))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = []
        for X in range(4):
            x0 = [origins[X][0]] + [solver.Value(alpha[X][o]) for o in range(1, nu)]
            b = [None] + [solver.Value(beta[X][o]) for o in range(1, nu)]
            x1 = [origins[X][1]] + b[1:]; x2 = [origins[X][2]] + [(v + 2) % 4 for v in b[1:]]
            sol.append(''.join(map(str, x0)) + '|' + ''.join(map(str, x1)) + '|' + ''.join(map(str, x2)))
        res['sol'] = 'SOL ' + ' '.join(sol)
    return res
if __name__ == '__main__':
    m = int(sys.argv[1]); hK = json.loads(sys.argv[2]); q = int(sys.argv[3]); costs = list(map(int, sys.argv[4:8]))
    tl = float(sys.argv[8]) if len(sys.argv) > 8 else 60; maxo = int(sys.argv[9]) if len(sys.argv) > 9 else 10**9
    opts = [origin_options(costs[0], q, True)] + [origin_options(c, q, False) if c != 0 else [(0, t, (t + 2) % 4) for t in range(4)] for c in costs[1:]]
    print('origin options per block', [len(o) for o in opts], flush=True)
    n = 0
    for org in itertools.product(*opts):
        n += 1
        if n > maxo: break
        r = build_and_solve(m, hK, q, org, tl)
        print(org, r, flush=True)
        if 'sol' in r: break
