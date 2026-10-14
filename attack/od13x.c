/* od13x.c: exhaustive search for skew-type OD(4t; 1, 3, 4t - 4) in Goethals-Seidel form (see od13.c, od13b.c)
 * that are invariant under a multiplier group H <= Z_t^* (all four blocks fixed by g -> hg).
 *
 * H must not move any s to -s (W0 is skew), i.e. -1 is not in H modulo any divisor d > 1 of t on which H acts.  The
 * pairs {s, -s} split into H-orbits; each block is determined by its values on one pair per orbit.  W0: one sign per
 * pair orbit.  (W1, W2, W3): one of the 20 patterns with S1 + S2 + S3 = 0 per pair orbit (od13b.c).  Row sums:
 * sum over pairs of Si = (a, b, -a - b) with a^2 + ab + b^2 = (t - 1)/2, up to permutations and a global sign, so all
 * targets are enumerated.  The search is a meet in the middle: the PAF vectors (on orbit representatives) of all W0 go
 * into a hash table; every pattern assignment with the right S sums is looked up.
 *
 * Usage: od13x t h1 [h2 ...]   (generators of H)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int t, h, ng, gens[8];
static int Hel[512], nH;
static int porb[256];         /* pair orbit id of pair p (1..h) */
static int norb, orep[256], osize[256];
static int pat[20][6], psym[20][3];
static int maxdepth;

static int pairof(int g) { g %= t; if (g < 0) g += t; return g <= h ? g : t - g; }

/* a sequence on Z_t determined by per-orbit values (vals[o][0] at rep, vals[o][1] at -rep) */
static void expand(int vals[][2], int *W) {
    memset(W, 0, sizeof(int) * t);
    for (int o = 0; o < norb; o++) {
        int r = orep[o];
        for (int k = 0; k < nH; k++) {
            int g = (int)((long)Hel[k] * r % t);
            W[g] = vals[o][0];
            W[(t - g) % t] = vals[o][1];
        }
    }
}

static int nreps;   /* number of distinct PAF constraints: one per pair orbit */
static void paf_reps(const int *W, long *out) {
    for (int o = 0; o < norb; o++) {
        int s = orep[o];
        long a = 0;
        for (int g = 0; g < t; g++) a += W[g] * W[(g + s) % t];
        out[o] = a;
    }
}

/* hash table of W0 PAF vectors */
#define HB 22
static long *hkeys; static int *hval; static int hused;
static uint64_t hashv(const long *v) {
    uint64_t x = 1469598103934665603ULL;
    for (int i = 0; i < norb; i++) { x ^= (uint64_t)(v[i] + 100000); x *= 1099511628211ULL; }
    return x;
}
static long *w0paf; static int nw0;

static int pa[256];      /* pattern per orbit */
static long found = 0, visited = 0;
static int tgt[3];
static int W1[3][1024];

static void report(int w0idx) {
    int vals[256][2], W[1024];
    printf("SOL t=%d\n", t);
    for (int o = 0; o < norb; o++) { int sg = (w0idx >> o) & 1 ? -1 : 1; vals[o][0] = sg; vals[o][1] = -sg; }
    expand(vals, W);
    for (int g = 0; g < t; g++) printf("%d%c", g == 0 ? 1 : 3 * W[g], g == t - 1 ? '\n' : ' ');
    for (int i = 0; i < 3; i++) {
        for (int o = 0; o < norb; o++) { vals[o][0] = pat[pa[o]][2 * i]; vals[o][1] = pat[pa[o]][2 * i + 1]; }
        expand(vals, W);
        for (int g = 0; g < t; g++) printf("%d%c", g == 0 ? 2 : 3 * W[g], g == t - 1 ? '\n' : ' ');
    }
    fflush(stdout);
}

static void leaf(void) {
    visited++;
    int vals[256][2], W[1024];
    long v[256], tot[256];
    for (int o = 0; o < norb; o++) tot[o] = 0;
    for (int i = 0; i < 3; i++) {
        for (int o = 0; o < norb; o++) { vals[o][0] = pat[pa[o]][2 * i]; vals[o][1] = pat[pa[o]][2 * i + 1]; }
        expand(vals, W);
        paf_reps(W, v);
        for (int o = 0; o < norb; o++) tot[o] += v[o];
    }
    for (int o = 0; o < norb; o++) tot[o] = -tot[o];
    uint64_t hv = hashv(tot) & ((1u << HB) - 1);
    while (hval[hv] >= 0) {
        if (!memcmp(&hkeys[hv * norb], tot, sizeof(long) * norb)) {
            /* all W0 with this vector: scan (few) */
            for (int k = 0; k < nw0; k++)
                if (!memcmp(&w0paf[(long)k * norb], tot, sizeof(long) * norb)) { found++; report(k); }
            return;
        }
        hv = (hv + 1) & ((1u << HB) - 1);
    }
}


/* ---- incremental search: PAF contributions of pattern choices on pairs of orbits ---- */
static int ocount;                     /* number of orbits (= norb) */
static int opos[256][1024], onpos[256]; /* positions of orbit o: first half h*rep (value index 0), then -h*rep (1) */
static int oidx[256][1024];
static long *T2;                       /* T2[((o1*norb + o2)*20 + c1)*20 + c2][k] */
static long cur[64][256];              /* running PAF vectors by depth */
static int ssum[64][3];
static int order_[256];

static inline long *T2at(int o1, int o2, int c1, int c2) {
    return &T2[((((long)o1 * norb + o2) * 20 + c1) * 20 + c2) * norb];
}

static void build_tables(void) {
    if (T2) return;
    for (int o = 0; o < norb; o++) {
        int r = orep[o], n = 0;
        int mark[1024]; memset(mark, 0, sizeof(int) * t);
        for (int k = 0; k < nH; k++) {
            int g = (int)((long)Hel[k] * r % t);
            if (!mark[g]) { mark[g] = 1; opos[o][n] = g; oidx[o][n] = 0; n++; }
            int ng_ = (t - g) % t;
            if (!mark[ng_]) { mark[ng_] = 1; opos[o][n] = ng_; oidx[o][n] = 1; n++; }
        }
        onpos[o] = n;
    }
    T2 = calloc((size_t)norb * norb * 400 * norb, sizeof(long));
    static int pos2orb[1024], pos2idx[1024];
    for (int g = 0; g < t; g++) pos2orb[g] = -1;
    for (int o = 0; o < norb; o++) for (int j = 0; j < onpos[o]; j++) { pos2orb[opos[o][j]] = o; pos2idx[opos[o][j]] = oidx[o][j]; }
    for (int o1 = 0; o1 < norb; o1++) for (int o2 = 0; o2 < norb; o2++)
        for (int c1 = 0; c1 < 20; c1++) for (int c2 = 0; c2 < 20; c2++) {
            long *T = T2at(o1, o2, c1, c2);
            for (int k = 0; k < norb; k++) {
                int sft = orep[k];
                long a = 0;
                /* sum over g in orbit o1, g + s in orbit o2 */
                for (int j = 0; j < onpos[o1]; j++) {
                    int g = opos[o1][j], g2 = (g + sft) % t;
                    if (pos2orb[g2] != o2) continue;
                    for (int i = 0; i < 3; i++)
                        a += (long)pat[c1][2 * i + oidx[o1][j]] * pat[c2][2 * i + pos2idx[g2]];
                }
                T[k] = a;
            }
        }
}

static void leaf2(int depth) {
    visited++;
    long tot[256];
    for (int o = 0; o < norb; o++) tot[o] = -cur[depth][o];
    uint64_t hv = hashv(tot) & ((1u << HB) - 1);
    while (hval[hv] >= 0) {
        if (!memcmp(&hkeys[hv * norb], tot, sizeof(long) * norb)) {
            for (int k = 0; k < nw0; k++)
                if (!memcmp(&w0paf[(long)k * norb], tot, sizeof(long) * norb)) { found++; report(k); }
            return;
        }
        hv = (hv + 1) & ((1u << HB) - 1);
    }
}

static void rec2(int d) {
    if (d == norb) {
        if (ssum[d][0] == tgt[0] && ssum[d][1] == tgt[1] && ssum[d][2] == tgt[2]) leaf2(d);
        return;
    }
    int rem = 0; for (int k = d; k < norb; k++) rem += osize[k];
    for (int i = 0; i < 3; i++) if (abs(tgt[i] - ssum[d][i]) > rem) return;
    for (int c = 0; c < 20; c++) {
        pa[d] = c;
        for (int i = 0; i < 3; i++) ssum[d + 1][i] = ssum[d][i] + osize[d] * psym[c][i];
        long *nx = cur[d + 1];
        for (int k = 0; k < norb; k++) nx[k] = cur[d][k];
        long *Td = T2at(d, d, c, c);
        for (int k = 0; k < norb; k++) nx[k] += Td[k];
        for (int j = 0; j < d; j++) {
            long *A = T2at(j, d, pa[j], c), *B = T2at(d, j, c, pa[j]);
            for (int k = 0; k < norb; k++) nx[k] += A[k] + B[k];
        }
        rec2(d + 1);
    }
}

static void rec(int o, int s0, int s1, int s2) {
    if (o == norb) {
        if (s0 == tgt[0] && s1 == tgt[1] && s2 == tgt[2]) leaf();
        return;
    }
    /* remaining capacity check: each orbit contributes at most osize in absolute value per coordinate */
    int rem = 0; for (int k = o; k < norb; k++) rem += osize[k];
    if (abs(tgt[0] - s0) > rem || abs(tgt[1] - s1) > rem || abs(tgt[2] - s2) > rem) return;
    for (int c = 0; c < 20; c++) {
        pa[o] = c;
        rec(o + 1, s0 + osize[o] * psym[c][0], s1 + osize[o] * psym[c][1], s2 + osize[o] * psym[c][2]);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: od13x t h1 [h2 ...]\n"); return 1; }
    t = atoi(argv[1]); h = (t - 1) / 2;
    ng = argc - 2; for (int i = 0; i < ng; i++) gens[i] = atoi(argv[2 + i]);
    /* group closure */
    nH = 1; Hel[0] = 1;
    for (int changed = 1; changed;) {
        changed = 0;
        for (int i = 0; i < nH; i++) for (int j = 0; j < ng; j++) {
            int x = (int)((long)Hel[i] * gens[j] % t), seen = 0;
            for (int k = 0; k < nH; k++) if (Hel[k] == x) seen = 1;
            if (!seen) { Hel[nH++] = x; changed = 1; }
        }
    }
    /* pair orbits; reject if some h maps a pair element s to -s */
    for (int p = 1; p <= h; p++) porb[p] = -1;
    norb = 0;
    for (int p = 1; p <= h; p++) {
        if (porb[p] >= 0) continue;
        int sz = 0;
        for (int k = 0; k < nH; k++) {
            int g = (int)((long)Hel[k] * p % t);
            if (g == t - p) { fprintf(stderr, "H maps %d to -%d: not allowed\n", p, p); return 1; }
            int q = pairof(g);
            if (porb[q] < 0) { porb[q] = norb; sz++; }
        }
        orep[norb] = p; osize[norb] = sz; norb++;
    }
    /* consistency: every element g of the orbit of rep has W(g) = W(rep); check the stabiliser acts trivially */
    fprintf(stderr, "t=%d |H|=%d pair orbits %d:", t, nH, norb);
    for (int o = 0; o < norb; o++) fprintf(stderr, " %d(x%d)", orep[o], osize[o]);
    fprintf(stderr, "\n");
    int npat = 0;
    for (int m = 0; m < 64; m++) {
        int v[6], sy = 0;
        for (int k = 0; k < 6; k++) v[k] = (m >> k) & 1 ? 1 : -1;
        for (int i = 0; i < 3; i++) sy += v[2 * i] + v[2 * i + 1];
        if (sy) continue;
        memcpy(pat[npat], v, sizeof v);
        for (int i = 0; i < 3; i++) psym[npat][i] = (v[2 * i] + v[2 * i + 1]) / 2;
        npat++;
    }
    /* W0 table */
    nw0 = 1 << norb;
    if (norb > 24) { fprintf(stderr, "too many orbits\n"); return 1; }
    w0paf = malloc(sizeof(long) * (size_t)nw0 * norb);
    hkeys = malloc(sizeof(long) * (size_t)(1u << HB) * norb);
    hval = malloc(sizeof(int) * (1u << HB));
    for (unsigned i = 0; i < (1u << HB); i++) hval[i] = -1;
    for (int k = 0; k < nw0; k++) {
        int vals[256][2], W[1024];
        for (int o = 0; o < norb; o++) { int sg = (k >> o) & 1 ? -1 : 1; vals[o][0] = sg; vals[o][1] = -sg; }
        expand(vals, W);
        paf_reps(W, &w0paf[(long)k * norb]);
        uint64_t hv = hashv(&w0paf[(long)k * norb]) & ((1u << HB) - 1);
        int dup = 0;
        while (hval[hv] >= 0) {
            if (!memcmp(&hkeys[hv * norb], &w0paf[(long)k * norb], sizeof(long) * norb)) { dup = 1; break; }
            hv = (hv + 1) & ((1u << HB) - 1);
        }
        if (!dup) { memcpy(&hkeys[hv * norb], &w0paf[(long)k * norb], sizeof(long) * norb); hval[hv] = k; hused++; }
    }
    fprintf(stderr, "W0: %d sequences, %d distinct PAF vectors\n", nw0, hused);
    /* targets: (a, b, -a-b) with a^2 + ab + b^2 = h.  Permuting W1, W2, W3 and negating all three preserve the
     * conditions, and act transitively on the 12 sign/order variants of one solution (a, b), so one target per
     * unordered {a, b, -a-b} up to sign suffices. */
    int ntg = 0, done[64][2]; int nd = 0;
    for (int a = -64; a <= 64; a++) for (int b = -64; b <= 64; b++) {
        if (a * a + a * b + b * b != h) continue;
        int v[3] = {abs(a), abs(b), abs(a + b)}, x;
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2 - i; j++) if (v[j] > v[j + 1]) { x = v[j]; v[j] = v[j + 1]; v[j + 1] = x; }
        int seen = 0;
        for (int k = 0; k < nd; k++) if (done[k][0] == v[0] && done[k][1] == v[1]) seen = 1;
        if (seen) continue;
        done[nd][0] = v[0]; done[nd][1] = v[1]; nd++;
        tgt[0] = a; tgt[1] = b; tgt[2] = -a - b;
        ntg++;
        long f0 = found, v0 = visited;
        build_tables();
        rec2(0);
        fprintf(stderr, "target (%d, %d, %d): %ld leaves, %ld solutions\n", a, b, -a - b, visited - v0, found - f0);
    }
    fprintf(stderr, "targets %d, total leaves %ld, solutions %ld\n", ntg, visited, found);
    return 0;
}
