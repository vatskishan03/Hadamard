"""Run lg_pt on one block of a prime structure at the Gaussian-model threshold for a target ball size.
  gen_run.py v e cfg k Lball seconds R Tmin Tmax nprop threads [seed] [maxdistinct]"""
import sys, os, time, json, resource, subprocess, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, lg_lib as LG, cyclo_product as CP
from scipy import stats, optimize
WORK = os.environ['HADAMARD_WORK']

def tau_for(p, r, nc):
    f = (lambda x: stats.ncx2.logcdf(x, r, nc) - math.log(p)) if nc > 1e-9 else (lambda x: stats.chi2.logcdf(x, r) - math.log(p))
    return optimize.brentq(f, 1e-3, 10 * r + 10 * nc)

def cpu_children():
    r = resource.getrusage(resource.RUSAGE_CHILDREN); return r.ru_utime + r.ru_stime

if __name__ == '__main__':
    v, e, cfg, k = [int(x) for x in sys.argv[1:5]]
    Lb = float(sys.argv[5]); secs = float(sys.argv[6]); R = int(sys.argv[7]); Tmin = float(sys.argv[8]); Tmax = float(sys.argv[9]); nprop = int(sys.argv[10]); th = int(sys.argv[11])
    seed = int(sys.argv[12]) if len(sys.argv) > 12 else 1
    mode_arg = sys.argv[14] if len(sys.argv) > 14 else '0'; lam_arg = sys.argv[15] if len(sys.argv) > 15 else '0'
    S, y, N = LG.make_structure(v, e, cfg)
    t0 = time.time(); mom = LG.moments(S, N, y, nsamp=int(os.environ.get('NSAMP', '300000'))); print('moments %.1fs' % (time.time() - t0), flush=True)
    ck, W, Z = LG.energy_setup(mom, k)
    delta = W @ (mom[k]['mu'] - ck); nc = float(delta @ delta)
    size = 2.0 ** (e // 2) if k == 0 else math.comb(e, (e + mom[k]['s']) // 2)
    p = Lb / size
    tau = tau_for(p, Z, nc)
    print('block %d size %.3g ball target %.3g p=%.3g rank %d nc %.2f tau %.3f' % (k, size, Lb, p, Z, nc, tau), flush=True)
    tag = '%d_%d_%d_k%d' % (v, e, cfg, k)
    inp = WORK + '/g_%s.in' % tag; outp = WORK + '/g_%s.out' % tag
    LG.export_block(inp, S, N, mom, k, ck=ck, W=W)
    maxd = int(sys.argv[13]) if len(sys.argv) > 13 else int(Lb)
    c0 = cpu_children(); t0 = time.time()
    r = subprocess.run(['./lg_pt', inp, outp, repr(tau), repr(tau), str(secs), '1000000000', str(th), str(seed), str(R), str(Tmin), str(Tmax), str(nprop), str(maxd), WORK + '/g_%s.stat' % tag, mode_arg, lam_arg], capture_output=True, text=True)
    print(r.stdout.strip()); print('\n'.join(r.stderr.strip().splitlines()[-R - 3:]), flush=True)
    print('cpu %.1f wall %.1f' % (cpu_children() - c0, time.time() - t0))
    print(open(WORK + '/g_%s.stat' % tag).read().replace('\n', ' | ')[:600])
