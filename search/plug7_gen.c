// plug7_gen.c (search/NOTES.md section 82): fast level-1 generator (E0 instances) for the plug-in class Z_p x F_7, H = squares, q = 3 mod 4, unbordered.
// Per x0-quad (A,B symmetric with x0(0)=0; C,D free; sampled with sum PAF_x0(k) = 0 mod 3 for all k, necessary):
//   * exact scalar class conditions from the trivial characters:
//       R0: sum_X (sigma_X + 6 eps_X)^2 = 4m-2        (E0 at the trivial K-character; sigma_X = sum x0_X, eps_X = sum e_X)
//       R1: sum_X (sigma_X - eps_X)^2 + 7 (tau_C^2 + tau_D^2) = 4m-2 solvable (level-2 row sums; tau_A = tau_B = 0)
//   * defect-pattern lists (u = x0 + 6e) tagged by (eps, np); sublist joins only for feasible (eps, np) class combinations,
//     exact hash join on (PAF_u(1..h), |u|^2) = (0..0, 4m-2).
// Output: instances in the X0 format of the sampler of search/NOTES.md section 78 (program not included), with Gaussian log-weight.  usage: plug7_gen p a b c d nquads seed logWthr outfile
// The run that found the core of order 322 (hadamard/2572): GEN_QMIN=10 GEN_QMAX=10 GEN_EFFMIN=3.5e3 GEN_TLIMIT=150
//   plug7_gen 23 1 1 4 4 100000000 SEED -1e9 out.g; then the instances (lines "I ... | X0 | X0 | X0 | X0") are written
//   as "p n" followed by the four X0 rows of each instance (out.l) and completed by plug7_l2 out.l 20.
// env: GEN_TLIMIT (seconds), GEN_STAT (print per-quad statistics)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MAXP 31
#define MAXH 15
static int p, h, MT;
static uint64_t rs = 88172645463325252ULL;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static int X0[4][MAXP];
static double binom_[40][40]; static double WTHR = 0; static int QMIN = 0; static double EFFMIN = 0; static int QMAX = 1000000; static int CDMAX = 1000; static int ABMIN = 0;
static double wcalc(int nvC, int nvD, int Q) { double w = 0; for (int tc = -nvC; tc <= nvC; tc += 2) { int rem = Q - tc * tc; if (rem < 0) continue; for (int td = -nvD; td <= nvD; td += 2) if (td * td == rem) w += binom_[nvC][(nvC + tc) / 2] * binom_[nvD][(nvD + td) / 2]; } return w; }
static uint64_t coef[MAXH + 1];      // hash coefficients for keys 0..h (h = norm coordinate)

typedef struct { uint64_t h, pat; uint8_t tag; } Ent;     // tag = (eps+16) | np<<6 ... see below
typedef struct { Ent *e; long n, cap; } List;

// ---- per-block state for DFS
static int lagk[MAXP][MAXP];            // reduced lag between positions
static uint64_t hv[4][MAXP][2];         // hash increment of a defect (pos, sign index 0:+1 1:-1), excluding couplings
static uint64_t hcpl[MAXH + 1];         // 36*coef_k
static int ptype[4][MAXP][2];           // 1 if p-type (e == x0)
static int dpos[16], dsg[16], dnd; static uint64_t curH; static int cur_eps, cur_np;
static List *curL; static int needE[200]; static int nneedE; static int stepU;
static inline int eps_reachable(int eps, int r) { for (int i = 0; i < nneedE; i++) { int d = needE[i] - eps; if (d % stepU) continue; d /= stepU; if (d < 0) d = -d; if (d <= r && ((d - r) & 1) == 0) return 1; } return 0; }

static void push_leaf(void) {
    if (curL->n == curL->cap) { curL->cap = curL->cap ? curL->cap * 2 : 4096; curL->e = realloc(curL->e, sizeof(Ent) * curL->cap); if (!curL->e) { fprintf(stderr, "oom\n"); exit(3); } }
    Ent *x = &curL->e[curL->n++];
    x->h = curH; uint64_t pt = 0; for (int i = 0; i < dnd && i < 10; i++) pt |= (uint64_t)((dpos[i] << 1) | (dsg[i] < 0)) << (6 * i);
    x->pat = pt; x->tag = (uint8_t)(((cur_eps + 16) & 31) | (cur_np << 5));    // eps in [-16,15] -> 5 bits, np <= 7 -> 3 bits
}
static inline void add_d(int bi, int j, int sgn) {
    int si = sgn > 0 ? 0 : 1;
    uint64_t dh = hv[bi][j][si];
    for (int i = 0; i < dnd; i++) { int k = lagk[j][dpos[i]]; dh += (sgn * dsg[i] > 0) ? hcpl[k] : (0 - hcpl[k]); }
    curH += dh; dpos[dnd] = j; dsg[dnd] = sgn; dnd++; cur_eps += sgn; cur_np += ptype[bi][j][si];
}
static inline void rem_d(int bi, int j, int sgn) {
    dnd--; int si = sgn > 0 ? 0 : 1;
    uint64_t dh = hv[bi][j][si];
    for (int i = 0; i < dnd; i++) { int k = lagk[j][dpos[i]]; dh += (sgn * dsg[i] > 0) ? hcpl[k] : (0 - hcpl[k]); }
    curH -= dh; cur_eps -= sgn; cur_np -= ptype[bi][j][si];
}
static void rec_cd(int bi, int start, int depth, int cnt) {
    if (depth == cnt) { push_leaf(); return; }
    for (int j = start; j <= p - (cnt - depth); j++)
        for (int sgn = 1; sgn >= (depth == 0 ? 1 : -1); sgn -= 2) {
            if (!eps_reachable(cur_eps + sgn, cnt - depth - 1)) continue;
            add_d(bi, j, sgn); rec_cd(bi, j + 1, depth + 1, cnt); rem_d(bi, j, sgn); }
}
static void rec_ab(int bi, int start, int depth, int cnt) {
    if (depth == cnt) { push_leaf(); return; }
    for (int j = start; j <= h - (cnt - depth - 1); j++)
        for (int sgn = 1; sgn >= -1; sgn -= 2) {
            if (!eps_reachable(cur_eps + 2 * sgn, cnt - depth - 1)) continue;
            add_d(bi, j, sgn); add_d(bi, p - j, sgn); rec_ab(bi, j + 1, depth + 1, cnt); rem_d(bi, p - j, sgn); rem_d(bi, j, sgn); }
}
static uint64_t base_hash(int bi) {
    const int *x0 = X0[bi]; uint64_t H = 0; int nn = 0;
    for (int k = 1; k <= h; k++) { int a = 0; for (int j = 0; j < p; j++) a += x0[j] * x0[(j + k) % p]; H += coef[k - 1] * (uint64_t)(int64_t)a; }
    for (int j = 0; j < p; j++) nn += x0[j] * x0[j];
    return H + coef[h] * (uint64_t)(int64_t)nn;
}
static void prep_block(int bi) {
    const int *x0 = X0[bi];
    for (int j = 0; j < p; j++) for (int si = 0; si < 2; si++) {
        int e = si == 0 ? 1 : -1; uint64_t H = 0;
        for (int k = 1; k <= h; k++) H += coef[k - 1] * (uint64_t)(int64_t)(6 * e * (x0[(j + k) % p] + x0[(j - k + p) % p]));
        H += coef[h] * (uint64_t)(int64_t)(12 * x0[j] * e + 36);
        hv[bi][j][si] = H; ptype[bi][j][si] = (x0[j] * e > 0);
    }
}
static void build_list(int bi, int cnt, List *L) {
    stepU = bi < 2 ? 2 : 1;
    prep_block(bi); L->n = 0; curL = L; dnd = 0; cur_eps = 0; cur_np = 0; curH = base_hash(bi);
    if (bi < 2) { add_d(bi, 0, 1); rec_ab(bi, 1, 0, cnt); rem_d(bi, 0, 1); }
    else rec_cd(bi, 0, 0, cnt);
}
// recompute exact key (h+1 ints) of a pattern
static void key_of(int bi, uint64_t pat, int ndef, int *key, int *e) {
    for (int j = 0; j < p; j++) e[j] = 0;
    for (int i = 0; i < ndef; i++) { int f = (pat >> (6 * i)) & 63; e[f >> 1] = (f & 1) ? -1 : 1; }
    int u[MAXP]; for (int j = 0; j < p; j++) u[j] = X0[bi][j] + 6 * e[j];
    for (int k = 1; k <= h; k++) { int a = 0; for (int j = 0; j < p; j++) a += u[j] * u[(j + k) % p]; key[k - 1] = a; }
    int nn = 0; for (int j = 0; j < p; j++) nn += u[j] * u[j]; key[h] = nn;
}

// ---- weight (as gen.c)
static double logweight(int e[4][MAXP]) {
    int s[4][MAXP]; int Ps[MAXH + 1]; memset(Ps, 0, sizeof(Ps));
    for (int b = 0; b < 4; b++) for (int j = 0; j < p; j++) s[b][j] = X0[b][j] - e[b][j];
    for (int b = 0; b < 4; b++) for (int k = 1; k <= h; k++) { int a = 0; for (int j = 0; j < p; j++) a += s[b][j] * s[b][(j + k) % p]; Ps[k] += a; }
    double lw = 0; int nbits = 0; int inv2 = (p + 1) / 2;
    for (int k = 1; k <= h; k++) {
        if (Ps[k] % 7 != 0) return -1e9;
        double R = -(double)Ps[k] / 7.0, mu = 0, v = 0;
        for (int b = 0; b < 4; b++) {
            int nk = 0;
            for (int j = 0; j < p; j++) if (!e[b][j] && !e[b][(j + k) % p] && !(b < 2 && (j == 0 || (j + k) % p == 0))) nk++;
            if (b < 2) {
                int j0 = (int)(((long)(p - k) * inv2) % p);
                int fixed = (!e[b][j0] && !e[b][(j0 + k) % p] && j0 != 0 && (j0 + k) % p != 0);
                if (fixed) { mu += -1; v += 2.0 * (nk - 1); } else v += 2.0 * nk;
            } else v += nk;
        }
        if (v <= 0) { if (fabs(R - mu) > 0.5) return -1e9; v = 1; }
        lw += -0.5 * (R - mu) * (R - mu) / v - 0.5 * log(v);
        double vl = 0;
        for (int b = 0; b < 4; b++) {
            if (b < 2) { for (int j = 1; j <= h; j++) if (!e[b][j]) { int w = s[b][(j - k + p) % p] - s[b][(j + k) % p]; vl += 4.0 * w * w; } }
            else for (int j = 0; j < p; j++) if (!e[b][j]) { int w = s[b][(j - k + p) % p] - s[b][(j + k) % p]; vl += (double)w * w; }
        }
        if (vl < 4) vl = 4;
        lw += -0.5 * log(vl);
    }
    for (int b = 0; b < 4; b++) { int nn = 0; for (int j = (b < 2 ? 1 : 0); j <= (b < 2 ? h : p - 1); j++) if (!e[b][j]) nn++; nbits += nn; }
    return lw + nbits * log(2.0);
}

// ---- sublist join: hash (A x C) pairs, probe (B x D)
typedef struct { uint64_t v; uint32_t i, j; } TE;
static TE *tabm = NULL; static size_t tabm_n = 0;
static uint64_t *sh[4] = {0, 0, 0, 0}; static long shn[4] = {0, 0, 0, 0};
static long join_sub(const List *LA, const uint32_t *ia, int na, const List *LB, const uint32_t *ib, int nb, const List *LC, const uint32_t *ic, int nc,
                     const List *LD, const uint32_t *id, int nd, uint64_t HT, int nA_, int nB_, int nC_, int nD_, uint32_t *out, long maxout, const int *target) {
    const List *Ls[4] = {LA, LB, LC, LD}; const uint32_t *ix[4] = {ia, ib, ic, id}; int nn[4] = {na, nb, nc, nd};
    for (int b = 0; b < 4; b++) { if (shn[b] < nn[b]) { free(sh[b]); shn[b] = nn[b] + 1024; sh[b] = malloc(8 * shn[b]); } for (int i = 0; i < nn[b]; i++) sh[b][i] = Ls[b]->e[ix[b][i]].h; }
    double pa = (double)na * nc, pb = (double)nb * nd;
    long found = 0; int swap = pa > pb;     // table side = smaller
    const uint64_t *H1 = swap ? sh[1] : sh[0], *H2 = swap ? sh[3] : sh[2], *Q1 = swap ? sh[0] : sh[1], *Q2 = swap ? sh[2] : sh[3];
    int n1 = swap ? nb : na, n2 = swap ? nd : nc, m1 = swap ? na : nb, m2 = swap ? nc : nd;
    size_t need = 1; while (need < (size_t)(2 * (swap ? pb : pa) + 2)) need <<= 1;
    if (tabm_n < need) { free(tabm); tabm_n = need; tabm = malloc(sizeof(TE) * tabm_n); }
    uint64_t cm = need - 1;
    for (size_t q = 0; q < need; q++) tabm[q].i = 0xffffffffu;
    for (int x = 0; x < n1; x++) { uint64_t hx = H1[x]; for (int y = 0; y < n2; y++) { uint64_t v = hx + H2[y]; uint64_t pos = (v * 0x9E3779B97F4A7C15ULL) >> 20 & cm; while (tabm[pos].i != 0xffffffffu) pos = (pos + 1) & cm; tabm[pos].v = v; tabm[pos].i = x; tabm[pos].j = y; } }
    for (int x = 0; x < m1; x++) {
        uint64_t hx = HT - Q1[x];
        for (int y = 0; y < m2; y++) {
            uint64_t v = hx - Q2[y]; uint64_t pos = (v * 0x9E3779B97F4A7C15ULL) >> 20 & cm;
            while (tabm[pos].i != 0xffffffffu) {
                if (tabm[pos].v == v) {
                    uint32_t q1 = tabm[pos].i, q2 = tabm[pos].j;
                    uint32_t a, b, c, d;
                    if (swap) { b = ib[q1]; d = id[q2]; a = ia[x]; c = ic[y]; } else { a = ia[q1]; c = ic[q2]; b = ib[x]; d = id[y]; }
                    int ka[MAXH + 1], kb[MAXH + 1], kc[MAXH + 1], kd[MAXH + 1], ee[MAXP];
                    key_of(0, LA->e[a].pat, nA_, ka, ee); key_of(1, LB->e[b].pat, nB_, kb, ee); key_of(2, LC->e[c].pat, nC_, kc, ee); key_of(3, LD->e[d].pat, nD_, kd, ee);
                    int ok = 1; for (int k = 0; k <= h; k++) if (ka[k] + kb[k] + kc[k] + kd[k] != target[k]) { ok = 0; break; }
                    if (ok) { if (found < maxout) { out[4 * found] = a; out[4 * found + 1] = b; out[4 * found + 2] = c; out[4 * found + 3] = d; } found++; }
                }
                pos = (pos + 1) & cm;
            }
        }
    }
    return found;
}

int main(int argc, char **argv) {
    if (argc < 10) { fprintf(stderr, "usage: plug7_gen p a b c d nquads seed logWthr outfile\n"); return 1; }
    p = atoi(argv[1]); h = (p - 1) / 2; int ca = atoi(argv[2]), cb = atoi(argv[3]), cc = atoi(argv[4]), cdd = atoi(argv[5]);
    long nq = atol(argv[6]); rs ^= (uint64_t)atol(argv[7]) * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) rnd();
    double thr = atof(argv[8]); FILE *fo = fopen(argv[9], "w");
    double tlim = getenv("GEN_TLIMIT") ? atof(getenv("GEN_TLIMIT")) : 1e30;
    int stat = getenv("GEN_STAT") != NULL; FILE *jlog = getenv("GEN_JLOG") ? fopen(getenv("GEN_JLOG"), "w") : NULL; FILE *fq = getenv("GEN2_QUADS") ? fopen(getenv("GEN2_QUADS"), "r") : NULL;
    MT = 28 * p - 2; WTHR = getenv("GEN_WTHR") ? atof(getenv("GEN_WTHR")) : 0; QMIN = getenv("GEN_QMIN") ? atoi(getenv("GEN_QMIN")) : 0; EFFMIN = getenv("GEN_EFFMIN") ? atof(getenv("GEN_EFFMIN")) : 0; QMAX = getenv("GEN_QMAX") ? atoi(getenv("GEN_QMAX")) : 1000000; CDMAX = getenv("GEN_CDMAX") ? atoi(getenv("GEN_CDMAX")) : 1000; ABMIN = getenv("GEN_ABMIN") ? atoi(getenv("GEN_ABMIN")) : 0; for (int i = 0; i < 40; i++) { binom_[i][0] = 1; for (int j = 1; j <= i; j++) binom_[i][j] = binom_[i - 1][j - 1] + (j <= i - 1 ? binom_[i - 1][j] : 0); }
    int nd_[4] = {1 + 2 * ca, 1 + 2 * cb, cc, cdd};
    int Npos = 2 * ca + 2 * cb + cc + cdd; int NP = (p - 3) - Npos;            // number of p-type positions (excluding z)
    if (NP < 0 || NP > Npos) { fprintf(stderr, "allocation infeasible (N=%d)\n", Npos); return 1; }
    for (int k = 0; k <= h; k++) coef[k] = (rnd() >> 3) | 1;
    for (int k = 1; k <= h; k++) hcpl[k] = 36 * coef[k - 1];
    for (int j = 0; j < p; j++) for (int j2 = 0; j2 < p; j2++) { int d = j - j2; if (d < 0) d += p; lagk[j][j2] = d == 0 ? 0 : (d <= h ? d : p - d); }
    // hcpl indexed by lag k (1..h) -> uses coef[k-1]
    // AB pairs by mod-3 bucket
    int nsym = 1 << h; int *Asym = malloc(sizeof(int) * nsym * p);
    for (int r = 0; r < nsym; r++) { for (int j = 1; j <= h; j++) { int v = ((r >> (j - 1)) & 1) ? -1 : 1; Asym[r * p + j] = v; Asym[r * p + p - j] = v; } Asym[r * p] = 0; }
    int pw3[MAXH + 1]; pw3[0] = 1; for (int i = 1; i <= h; i++) pw3[i] = pw3[i - 1] * 3;
    int nb = pw3[h];
    unsigned char *pafA = malloc((size_t)nsym * h);
    for (int r = 0; r < nsym; r++) for (int k = 1; k <= h; k++) { int a = 0; for (int j = 0; j < p; j++) a += Asym[r * p + j] * Asym[r * p + (j + k) % p]; pafA[(size_t)r * h + k - 1] = (unsigned char)(((a % 3) + 3) % 3); }
    int *bstart = calloc(nb + 1, sizeof(int)); uint32_t *bpairs;
    for (int r1 = 0; r1 < nsym; r1++) for (int r2 = 0; r2 < nsym; r2++) { int key = 0; for (int k = 0; k < h; k++) key += pw3[k] * ((pafA[(size_t)r1 * h + k] + pafA[(size_t)r2 * h + k]) % 3); bstart[key + 1]++; }
    for (int i = 0; i < nb; i++) bstart[i + 1] += bstart[i];
    bpairs = malloc(sizeof(uint32_t) * (size_t)nsym * nsym); { int *fill = malloc(sizeof(int) * nb); memcpy(fill, bstart, sizeof(int) * nb);
      for (int r1 = 0; r1 < nsym; r1++) for (int r2 = 0; r2 < nsym; r2++) { int key = 0; for (int k = 0; k < h; k++) key += pw3[k] * ((pafA[(size_t)r1 * h + k] + pafA[(size_t)r2 * h + k]) % 3); bpairs[fill[key]++] = ((uint32_t)r1 << 16) | (uint32_t)r2; } free(fill); }
    int target[MAXH + 1]; for (int k = 0; k < h; k++) target[k] = 0; target[h] = MT;
    uint64_t HT = 0; for (int k = 0; k <= h; k++) HT += coef[k] * (uint64_t)(int64_t)target[k];
    // eps ranges
    int emaxA[2] = {3, -1};
    long total_inst = 0, kept = 0, quads_proc = 0, quads_tried = 0, quads_feas = 0; double t0 = now(); long hist[80]; memset(hist, 0, sizeof(hist));
    double tbuild = 0, tjoin = 0; long njoins = 0; double sumpairs = 0;
    uint32_t *sol = malloc(sizeof(uint32_t) * 4 * 4000000L);
    List L[4]; memset(L, 0, sizeof(L));
    uint32_t *cls[4][256]; int ccnt[4][256];
    uint32_t *clsbuf[4]; for (int b = 0; b < 4; b++) clsbuf[b] = NULL;
    while ((fq || quads_proc < nq) && now() - t0 < tlim) {
        int chosen = -1; int feasE[100][4]; int nfe = 0; int sC = 0, sD = 0;
        if (fq) {
            int ok = 1; for (int bb = 0; bb < 4 && ok; bb++) for (int j = 0; j < p; j++) if (fscanf(fq, "%d", &X0[bb][j]) != 1) { ok = 0; break; }
            if (!ok) break;
            quads_tried++; int sg4[4]; for (int bb = 0; bb < 4; bb++) { sg4[bb] = 0; for (int j = 0; j < p; j++) sg4[bb] += X0[bb][j]; }
            nfe = 0;
            for (int sa = -ca; sa <= ca; sa += 2) for (int sb = -cb; sb <= cb; sb += 2) for (int eC = -cc; eC <= cc; eC += 2) for (int eD = -cdd; eD <= cdd; eD += 2) {
                int e4[4] = {1 + 2 * sa, 1 + 2 * sb, eC, eD};
                int r0 = 0, r1 = 0; for (int i = 0; i < 4; i++) { int a = sg4[i] + 6 * e4[i]; r0 += a * a; int b2 = sg4[i] - e4[i]; r1 += b2 * b2; }
                if (r0 != MT) continue;
                if (abs(sg4[2] - e4[2]) > CDMAX || abs(sg4[3] - e4[3]) > CDMAX) continue;
                { int ba = abs(sg4[0] - e4[0]), bb = abs(sg4[1] - e4[1]); if ((ba > bb ? ba : bb) < ABMIN) continue; }
                int q7 = MT - r1; if (q7 < 0 || q7 % 7) continue; int Q = q7 / 7; int ok2 = 0; int nvC = p - cc, nvD = p - cdd;
                double wq = wcalc(nvC, nvD, Q); ok2 = wq > 0 && wq >= WTHR && Q >= QMIN && Q <= QMAX;
                if (!ok2) continue;
                if (nfe < 100) { for (int i = 0; i < 4; i++) feasE[nfe][i] = e4[i]; nfe++; }
            }
            if (!nfe) continue; chosen = 0;
        } else {
        // sample quad
        for (int j = 0; j < p; j++) { X0[2][j] = (rnd() & 1) ? -1 : 1; X0[3][j] = (rnd() & 1) ? -1 : 1; }
        int key = 0;
        for (int k = 1; k <= h; k++) { int a = 0; for (int j = 0; j < p; j++) a += X0[2][j] * X0[2][(j + k) % p] + X0[3][j] * X0[3][(j + k) % p]; int rm = (((-a) % 3) + 3) % 3; key += pw3[k - 1] * rm; }
        for (int j = 0; j < p; j++) { sC += X0[2][j]; sD += X0[3][j]; }
        quads_tried++;
        int lo = bstart[key], hi = bstart[key + 1]; if (hi == lo) continue;
        int cand[64]; int nc_ = 0;
        for (int q = lo; q < hi && nc_ < 64; q++) cand[nc_++] = q;
        int start = (int)(rnd() % nc_);
        for (int t = 0; t < nc_ && chosen < 0; t++) {
            uint32_t pr = bpairs[cand[(start + t) % nc_]]; int rA = pr >> 16, rB = pr & 0xffff;
            int sA = 0, sB = 0; for (int j = 0; j < p; j++) { sA += Asym[rA * p + j]; sB += Asym[rB * p + j]; }
            int sg4[4] = {sA, sB, sC, sD}; nfe = 0;
            for (int sa = -ca; sa <= ca; sa += 2) for (int sb = -cb; sb <= cb; sb += 2) for (int eC = -cc; eC <= cc; eC += 2) for (int eD = -cdd; eD <= cdd; eD += 2) {
                int e4[4] = {1 + 2 * sa, 1 + 2 * sb, eC, eD};
                int r0 = 0, r1 = 0; for (int i = 0; i < 4; i++) { int a = sg4[i] + 6 * e4[i]; r0 += a * a; int b2 = sg4[i] - e4[i]; r1 += b2 * b2; }
                if (r0 != MT) continue;
                if (abs(sg4[2] - e4[2]) > CDMAX || abs(sg4[3] - e4[3]) > CDMAX) continue;
                { int ba = abs(sg4[0] - e4[0]), bb = abs(sg4[1] - e4[1]); if ((ba > bb ? ba : bb) < ABMIN) continue; }
                int q7 = MT - r1; if (q7 < 0 || q7 % 7) continue; int Q = q7 / 7; int ok = 0;
                int nvC = p - cc, nvD = p - cdd;
                double wq = wcalc(nvC, nvD, Q); ok = wq > 0 && wq >= WTHR && Q >= QMIN && Q <= QMAX;
                if (!ok) continue;
                if (nfe < 100) { for (int i = 0; i < 4; i++) feasE[nfe][i] = e4[i]; nfe++; }
            }
            if (nfe) { chosen = cand[(start + t) % nc_]; uint32_t pr2 = bpairs[chosen]; for (int j = 0; j < p; j++) { X0[0][j] = Asym[(pr2 >> 16) * p + j]; X0[1][j] = Asym[(pr2 & 0xffff) * p + j]; } }
        }
        if (chosen < 0) continue;
        }
        quads_feas++; quads_proc++;
        double tb0 = now();
        int cnts[4] = {ca, cb, cc, cdd};
        for (int b = 0; b < 4; b++) {
            nneedE = 0; for (int fe = 0; fe < nfe; fe++) { int v = feasE[fe][b], dup = 0; for (int i = 0; i < nneedE; i++) if (needE[i] == v) dup = 1; if (!dup) needE[nneedE++] = v; }
            build_list(b, cnts[b], &L[b]);
            // class buckets
            for (int t = 0; t < 256; t++) ccnt[b][t] = 0;
            for (long i = 0; i < L[b].n; i++) ccnt[b][L[b].e[i].tag]++;
            clsbuf[b] = realloc(clsbuf[b], sizeof(uint32_t) * (L[b].n + 1));
            long off = 0; for (int t = 0; t < 256; t++) { cls[b][t] = clsbuf[b] + off; off += ccnt[b][t]; }
            int fillc[256]; memset(fillc, 0, sizeof(fillc));
            for (long i = 0; i < L[b].n; i++) { int t = L[b].e[i].tag; cls[b][t][fillc[t]++] = (uint32_t)i; }
        }
        tbuild += now() - tb0;
        double tj0 = now();
        int E[4][MAXP];
        for (int fe = 0; fe < nfe; fe++) {
            // eps tags: A: eps = 1+2*s -> tag eps value stored as cur_eps (includes z defect sign +1 and pair signs *2)
            int et[4]; for (int i = 0; i < 4; i++) et[i] = (feasE[fe][i] + 16) & 31;
            // np combos
            int npmax[4] = {2 * ca, 2 * cb, cc, cdd};
            for (int na = 0; na <= npmax[0]; na += 2) for (int nbb = 0; nbb <= npmax[1]; nbb += 2) for (int nc2 = 0; nc2 <= npmax[2]; nc2++) {
                int nd2 = NP - na - nbb - nc2; if (nd2 < 0 || nd2 > npmax[3]) continue;
                int tg[4] = {et[0] | (na << 5), et[1] | (nbb << 5), et[2] | (nc2 << 5), et[3] | (nd2 << 5)};
                // A,B np counts: A's np includes z-defect? z has x0 = 0: type none (ptype uses x0*e>0 -> 0): fine
                int sz[4]; for (int i = 0; i < 4; i++) sz[i] = ccnt[i][tg[i]];
                if (!sz[0] || !sz[1] || !sz[2] || !sz[3]) continue;
                if (EFFMIN > 0 && ((double)sz[0] * sz[1] * sz[2] * sz[3]) / ((double)sz[0] * sz[2] + (double)sz[1] * sz[3]) < EFFMIN) continue;
                njoins++; sumpairs += (double)sz[0] * sz[2] + (double)sz[1] * sz[3];
                double jt0 = jlog ? now() : 0;
                long ns = join_sub(&L[0], cls[0][tg[0]], sz[0], &L[1], cls[1][tg[1]], sz[1], &L[2], cls[2][tg[2]], sz[2], &L[3], cls[3][tg[3]], sz[3], HT, nd_[0], nd_[1], nd_[2], nd_[3], sol, 4000000L, target);
                total_inst += ns;
                if (jlog) fprintf(jlog, "%d %d %d %d %ld %.3g\n", sz[0], sz[1], sz[2], sz[3], ns, (now() - jt0) * 1e6);
                for (long i = 0; i < ns && i < 4000000L; i++) {
                    for (int b = 0; b < 4; b++) { int dummy[MAXH + 1]; key_of(b, L[b].e[sol[4 * i + b]].pat, nd_[b], dummy, E[b]); }
                    double lw = logweight(E);
                    int hb = (int)floor(lw + 40); if (hb < 0) hb = 0; if (hb > 79) hb = 79; hist[hb]++;
                    if (lw >= thr) {
                        kept++;
                        fprintf(fo, "I %.4f", lw);
                        for (int b = 0; b < 4; b++) { fprintf(fo, " |"); for (int j = 0; j < p; j++) fprintf(fo, " %d", X0[b][j] + 6 * E[b][j]); }
                        fprintf(fo, "\n");
                    }
                }
            }
        }
        tjoin += now() - tj0;
        if (stat && quads_proc % 200 == 0) fprintf(stderr, "[%.1fs] tried %ld feas %ld inst %ld kept %ld tbuild %.2f tjoin %.2f joins %ld pairs/join %.3g\n", now() - t0, quads_tried, quads_feas, total_inst, kept, tbuild, tjoin, njoins, njoins ? sumpairs / njoins : 0.0);
    }
    fprintf(fo, "# DONE p=%d alloc=%d,%d,%d,%d quads_processed=%ld tried=%ld feasible=%ld instances=%ld kept=%ld time=%.2f tbuild=%.2f tjoin=%.2f joins=%ld\n", p, ca, cb, cc, cdd, quads_proc, quads_tried, quads_feas, total_inst, kept, now() - t0, tbuild, tjoin, njoins);
    fprintf(fo, "# HIST (lnW+40 bins):"); for (int i = 0; i < 80; i++) if (hist[i]) fprintf(fo, " %d:%ld", i - 40, hist[i]); fprintf(fo, "\n");
    fclose(fo);
    return 0;
}
