/*
 * hcore.c -- C kernels of the compression-chain engine (hchain.py).
 *
 *   dp_run / dp_fetch : dynamic program over the H-orbits of G.  A state is the vector of partial sums of
 *                       all "coordinates" (the quotient orbits of every compression level involved, with
 *                       exact integer weights).  Coordinates flagged `fixed` carry a target value and are
 *                       pruned exactly (reachability by the remaining weights, equality when closed);
 *                       the other (output) coordinates are listed with multiplicities.
 *   join4             : exact four-list join  V0[a] + V1[b] + V2[c] + V3[d] = T  on integer vectors, after
 *                       join4x.c of the repository (exact partition on up to NP classes, 64-bit linear hash
 *                       in small open-addressing tables, exact re-check of every hash match), as a library.
 *
 *   cc -O3 -march=native -shared -fPIC hcore.c -o libhcore.so
 *
 * Nothing is sampled: every routine is exhaustive; truncation is reported by return codes.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/mman.h>

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

/* ------------------------------------------------------------------------------------------------ DP */
static int16_t *g_k = NULL; static double *g_c = NULL; static int64_t g_n = 0; static int g_nout = 0;
static int64_t g_peak = 0;

static inline uint64_t hash_key(const int16_t *k, int Dp) {
    uint64_t h = 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < Dp; i += 4) {
        uint64_t w; memcpy(&w, k + i, 8);
        h = (h ^ w) * 0xff51afd7ed558ccdULL; h ^= h >> 32;
    }
    return h;
}

/* Returns the number of final states (>= 0), -1 when the state count exceeded max_states, -2 on allocation failure.
 * order[s]      : orbit processed at step s (a permutation of 0..O-1)
 * nco[o], co[o*Lmax+t], cw[o*Lmax+t] : the coordinates orbit o contributes to (index, weight): +-weight
 * fixed[j], tgt[j]                     : coordinate j is pruned against tgt[j] when fixed[j] != 0
 * outidx[0..nout-1]                    : the (non-fixed) coordinates reported
 * peak (out)                           : largest intermediate state count */
int64_t dp_run(int O, int D, int Lmax, const int32_t *order, const int32_t *nco, const int32_t *co, const int32_t *cw,
               const int32_t *fixed, const int32_t *tgt, int nout, const int32_t *outidx, int64_t max_states, int64_t *peak) {
    int Dp = (D + 3) & ~3;
    int32_t *rem = calloc((size_t)O * D, 4);          /* rem[s*D+j] = total weight of coordinate j left after step s */
    for (int s = O - 1; s >= 0; s--) {
        if (s < O - 1) { memcpy(rem + (size_t)s * D, rem + (size_t)(s + 1) * D, 4 * D);
            int o = order[s + 1]; for (int t = 0; t < nco[o]; t++) rem[(size_t)s * D + co[o * Lmax + t]] += abs(cw[o * Lmax + t]); }
    }
    int64_t n = 1, cap = 1;
    int16_t *ck = calloc((size_t)Dp, 2); double *cc = malloc(8); cc[0] = 1.0;
    int16_t tmp[256];
    int64_t pk = 1;
    for (int s = 0; s < O; s++) {
        int o = order[s]; int nt = nco[o];
        int64_t ncap = 2 * n; int64_t tsz = 16; while (tsz < 2 * ncap) tsz <<= 1;
        int16_t *nk = malloc((size_t)ncap * Dp * 2); double *nc = malloc((size_t)ncap * 8); int32_t *tab = malloc((size_t)tsz * 4);
        if (!nk || !nc || !tab) { free(nk); free(nc); free(tab); free(ck); free(cc); free(rem); return -2; }
        memset(tab, 0xFF, (size_t)tsz * 4);
        int64_t nn = 0; uint64_t mask = (uint64_t)tsz - 1;
        const int32_t *remS = rem + (size_t)s * D;
        for (int64_t i = 0; i < n; i++) {
            for (int sg = -1; sg <= 1; sg += 2) {
                memcpy(tmp, ck + (size_t)i * Dp, (size_t)Dp * 2);
                int ok = 1;
                for (int t = 0; t < nt; t++) {
                    int j = co[o * Lmax + t]; int v = tmp[j] + sg * cw[o * Lmax + t]; tmp[j] = (int16_t)v;
                    if (fixed[j]) { int d = tgt[j] - v; if (d < 0) d = -d; if (d > remS[j]) { ok = 0; break; } }
                }
                if (!ok) continue;
                uint64_t h = hash_key(tmp, Dp) & mask;
                for (;;) {
                    int32_t e = tab[h];
                    if (e < 0) { memcpy(nk + (size_t)nn * Dp, tmp, (size_t)Dp * 2); nc[nn] = cc[i]; tab[h] = (int32_t)nn; nn++; break; }
                    if (memcmp(nk + (size_t)e * Dp, tmp, (size_t)Dp * 2) == 0) { nc[e] += cc[i]; break; }
                    h = (h + 1) & mask;
                }
            }
        }
        free(tab); free(ck); free(cc);
        ck = nk; cc = nc; n = nn; cap = ncap;
        if (n > pk) pk = n;
        if (n > max_states) { free(ck); free(cc); free(rem); if (peak) *peak = pk; return -1; }
    }
    free(rem);
    free(g_k); free(g_c);
    g_nout = nout; g_n = n;
    g_k = malloc((size_t)(n ? n : 1) * (nout ? nout : 1) * 2); g_c = malloc((size_t)(n ? n : 1) * 8);
    for (int64_t i = 0; i < n; i++) {
        for (int q = 0; q < nout; q++) g_k[(size_t)i * nout + q] = ck[(size_t)i * Dp + outidx[q]];
        g_c[i] = cc[i];
    }
    free(ck); free(cc);
    if (peak) *peak = pk;
    (void)cap;
    return n;
}

void dp_fetch(int16_t *k, double *c) {
    if (g_n) { memcpy(k, g_k, (size_t)g_n * (g_nout ? g_nout : 1) * 2); memcpy(c, g_c, (size_t)g_n * 8); }
    free(g_k); free(g_c); g_k = NULL; g_c = NULL; g_n = 0;
}

/* ------------------------------------------------------------------------------------------------ join4 */
#define NP 3
#define MAXC 256
#define KEYLIM (1 << 22)
#define EMPTY 0xFFFFFFFFu
#define GOLD 0x9E3779B97F4A7C15ULL
#define PB 16

typedef struct { uint64_t hv; uint32_t ia, ib; } Slot;
typedef struct { int64_t n; int *v; int64_t *start; uint32_t *idx; uint64_t *hs; } Cells;

static int J_C, J_P;
static const int32_t *J_V[4], *J_T;
static int64_t J_nn[4];
static uint64_t J_R[MAXC], *J_h[4], J_hT;
static int J_pc[NP], J_W[NP], J_KW[NP], J_Tp[NP];
static int64_t J_l[MAXC], J_g[MAXC];

static int lexcmp(const int *a, const int *b) {
    for (int j = 0; j < J_P; j++) if (a[j] != b[j]) return a[j] < b[j] ? -1 : 1;
    return 0;
}
typedef struct { int v[NP]; uint32_t idx; } Rec;
static int cmprec(const void *x, const void *y) {
    const Rec *a = x, *b = y;
    int c = lexcmp(a->v, b->v); if (c) return c;
    return a->idx < b->idx ? -1 : (a->idx > b->idx);
}
static Cells make_cells(int k) {
    int64_t n = J_nn[k];
    Rec *r = malloc(sizeof(Rec) * (n + 1));
    for (int64_t i = 0; i < n; i++) {
        r[i].idx = (uint32_t)i;
        for (int j = 0; j < NP; j++) r[i].v[j] = j < J_P ? (int)((J_V[k][i * J_C + J_pc[j]] - J_l[J_pc[j]]) / J_g[J_pc[j]]) : 0;
    }
    qsort(r, n, sizeof(Rec), cmprec);
    Cells c; c.n = 0; c.v = malloc(sizeof(int) * NP * (n + 1)); c.start = malloc(8 * (n + 2)); c.idx = malloc(4 * (n + 1)); c.hs = malloc(8 * (n + 1));
    for (int64_t i = 0; i < n; i++) {
        c.idx[i] = r[i].idx; c.hs[i] = J_h[k][r[i].idx];
        if (c.n == 0 || lexcmp(c.v + (c.n - 1) * NP, r[i].v) != 0) { memcpy(c.v + c.n * NP, r[i].v, sizeof(int) * NP); c.start[c.n] = i; c.n++; }
    }
    c.start[c.n] = n; free(r);
    return c;
}
static void free_cells(Cells *c) { free(c->v); free(c->start); free(c->idx); free(c->hs); }
static int64_t lower_c(const Cells *c, const int *key) {
    int64_t a = 0, b = c->n;
    while (a < b) { int64_t m = (a + b) / 2; if (lexcmp(c->v + m * NP, key) < 0) a = m + 1; else b = m; }
    return a;
}
static int64_t upper_c(const Cells *c, const int *key) {
    int64_t a = 0, b = c->n;
    while (a < b) { int64_t m = (a + b) / 2; if (lexcmp(c->v + m * NP, key) <= 0) a = m + 1; else b = m; }
    return a;
}
static int64_t enc(const int *s) {
    int64_t id = 0;
    for (int j = 0; j < J_P; j++) { int u = s[j]; if (u < 0 || u >= J_KW[j]) return -1; id = id * J_KW[j] + u; }
    return id;
}
static void dec(int64_t id, int *s) {
    for (int j = J_P - 1; j >= 0; j--) { s[j] = (int)(id % J_KW[j]); id /= J_KW[j]; }
    for (int j = J_P; j < NP; j++) s[j] = 0;
}
static int64_t gcd64(int64_t a, int64_t b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { int64_t t = a % b; a = b; b = t; } return a; }

/* All (a,b,c,d) with V0[a]+V1[b]+V2[c]+V3[d] = T (exact).  Returns the exact number of solutions found
 * (stopping once stop_at > 0 and that many were found), writes the first maxout as 4 uint32 each into out.
 * budget = build pairs per chunk.  Returns -2 on out of memory.
 * stats[0] = build pairs, stats[1] = probe pairs available, stats[2] = probe pairs touched, stats[3] = partition classes */
int64_t join4(const int32_t *V0, int64_t n0, const int32_t *V1, int64_t n1, const int32_t *V2, int64_t n2,
              const int32_t *V3, int64_t n3, int C, const int32_t *T, uint32_t *out, int64_t maxout,
              int64_t stop_at, int64_t budget, int verbose, double *stats) {
    if (C > MAXC) return -3;
    J_C = C; J_T = T; J_V[0] = V0; J_V[1] = V1; J_V[2] = V2; J_V[3] = V3;
    J_nn[0] = n0; J_nn[1] = n1; J_nn[2] = n2; J_nn[3] = n3;
    if (n0 == 0 || n1 == 0 || n2 == 0 || n3 == 0) return 0;
    /* exact lattice normalisation of every class: v -> (v - l) / g */
    for (int c = 0; c < C; c++) {
        int64_t lo = INT64_MAX, hi = INT64_MIN, hs_ = 0; (void)hs_;
        int64_t slo = 0, shi = 0;
        for (int k = 0; k < 4; k++) {
            int64_t mn = INT64_MAX, mx = INT64_MIN;
            for (int64_t i = 0; i < J_nn[k]; i++) { int64_t x = J_V[k][i * C + c]; if (x < mn) mn = x; if (x > mx) mx = x; }
            slo += mn; shi += mx; if (mn < lo) lo = mn; if (mx > hi) hi = mx;
        }
        if (T[c] < slo || T[c] > shi) return 0;
        int64_t g = 0;
        for (int k = 0; k < 4; k++) for (int64_t i = 0; i < J_nn[k]; i++) { g = gcd64(g, J_V[k][i * C + c] - lo); if (g == 1) break; }
        if (g == 0) { if (4 * lo != T[c]) return 0; g = 1; }
        if ((((int64_t)T[c] - 4 * lo) % g) != 0) return 0;
        J_l[c] = lo; J_g[c] = g;
    }
    int map[4] = {0, 1, 2, 3};
    if ((double)J_nn[0] * J_nn[1] > (double)J_nn[2] * J_nn[3]) {       /* build side = smaller pair product */
        const int32_t *tv; int64_t tn;
        tv = J_V[0]; J_V[0] = J_V[2]; J_V[2] = tv; tn = J_nn[0]; J_nn[0] = J_nn[2]; J_nn[2] = tn;
        tv = J_V[1]; J_V[1] = J_V[3]; J_V[3] = tv; tn = J_nn[1]; J_nn[1] = J_nn[3]; J_nn[3] = tn;
        map[0] = 2; map[1] = 3; map[2] = 0; map[3] = 1;
    }
    uint64_t s = 0x9E3779B97F4A7C15ULL;
    for (int c = 0; c < C; c++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; J_R[c] = s | 1; }
    J_hT = 0; for (int c = 0; c < C; c++) J_hT += J_R[c] * (uint64_t)(int64_t)T[c];
    for (int k = 0; k < 4; k++) {
        J_h[k] = malloc(8 * (J_nn[k] + 1));
        for (int64_t i = 0; i < J_nn[k]; i++) {
            uint64_t x = 0; const int32_t *row = J_V[k] + i * C;
            for (int c = 0; c < C; c++) x += J_R[c] * (uint64_t)(int64_t)row[c];
            J_h[k][i] = x;
        }
    }
    /* partition classes: greedy by joint entropy of the normalised values (on a subsample); small joins are not partitioned */
    int64_t Wc[MAXC];
    for (int c = 0; c < C; c++) {
        int64_t mx = 0;
        for (int k = 0; k < 4; k++) for (int64_t i = 0; i < J_nn[k]; i++) { int64_t q = (J_V[k][i * C + c] - J_l[c]) / J_g[c]; if (q > mx) mx = q; }
        Wc[c] = mx + 1;
    }
    int used[MAXC]; memset(used, 0, sizeof(used));
    int Pmax = ((double)J_nn[0] * J_nn[1] + (double)J_nn[2] * J_nn[3] < 2e6) ? 0 : NP;
    {
        int64_t ns = 0, *sub[4]; int64_t cnt[4];
        uint64_t *code = malloc(8 * 4100 * 4); uint64_t *tmp = malloc(8 * 4100 * 4);
        int64_t sidx[4][1025];
        for (int k = 0; k < 4; k++) {
            int64_t st = J_nn[k] / 1024; if (st < 1) st = 1; cnt[k] = 0;
            for (int64_t i = 0; i < J_nn[k] && cnt[k] < 1024; i += st) sidx[k][cnt[k]++] = i;
            ns += cnt[k]; sub[k] = sidx[k];
        }
        for (int64_t i = 0; i < ns; i++) code[i] = 0;
        J_P = 0; double space = 1, Hprev = 0;
        while (J_P < Pmax) {
            int best = -1; double bestH = -1;
            for (int c = 0; c < C; c++) {
                if (used[c] || Wc[c] < 2 || space * (2.0 * Wc[c] - 1) > (double)KEYLIM) continue;
                int64_t m = 0;
                for (int k = 0; k < 4; k++) for (int64_t t = 0; t < cnt[k]; t++) {
                    int64_t i = sub[k][t]; int64_t q = (J_V[k][i * C + c] - J_l[c]) / J_g[c];
                    tmp[m] = code[m] * (uint64_t)Wc[c] + (uint64_t)q; m++;
                }
                /* entropy of tmp (sorted) */
                {   /* shell sort */
                    int64_t gap = m / 2;
                    while (gap > 0) { for (int64_t i = gap; i < m; i++) { uint64_t x = tmp[i]; int64_t j = i; while (j >= gap && tmp[j - gap] > x) { tmp[j] = tmp[j - gap]; j -= gap; } tmp[j] = x; } gap /= 2; }
                }
                double H = 0; int64_t run = 1;
                for (int64_t i = 1; i <= m; i++) {
                    if (i < m && tmp[i] == tmp[i - 1]) { run++; continue; }
                    double p = (double)run / m; H -= p * log(p); run = 1;
                }
                if (H > bestH + 1e-12) { bestH = H; best = c; }
            }
            if (best < 0 || bestH < Hprev + 0.15) break;
            /* commit class `best` */
            int64_t m = 0;
            for (int k = 0; k < 4; k++) for (int64_t t = 0; t < cnt[k]; t++) {
                int64_t i = sub[k][t]; int64_t q = (J_V[k][i * C + best] - J_l[best]) / J_g[best];
                code[m] = code[m] * (uint64_t)Wc[best] + (uint64_t)q; m++;
            }
            used[best] = 1; J_pc[J_P] = best; J_W[J_P] = (int)Wc[best]; J_KW[J_P] = 2 * (int)Wc[best] - 1; space *= J_KW[J_P]; J_P++; Hprev = bestH;
        }
        free(code); free(tmp);
    }
    Cells cl[4]; int64_t KS = 1;
    {   /* cell-pair work must stay a small fraction of the pair work: drop the last classes otherwise */
        double pairs = (double)J_nn[0] * J_nn[1] + (double)J_nn[2] * J_nn[3];
        for (;;) {
            for (int j = 0; j < NP; j++) J_Tp[j] = j < J_P ? (int)(((int64_t)T[J_pc[j]] - 4 * J_l[J_pc[j]]) / J_g[J_pc[j]]) : 0;
            KS = 1; for (int j = 0; j < J_P; j++) KS *= J_KW[j];
            for (int k = 0; k < 4; k++) cl[k] = make_cells(k);
            double work = (double)cl[0].n * cl[1].n + (double)cl[2].n * cl[3].n;
            if (J_P == 0 || work <= pairs / 48.0) break;
            for (int k = 0; k < 4; k++) free_cells(&cl[k]);
            J_P--;
        }
    }
    double t0 = now();
    if (verbose) fprintf(stderr, "join4: n = %ld %ld %ld %ld, C = %d, partition on %d classes (key space %ld), cells %ld %ld %ld %ld\n",
            (long)J_nn[0], (long)J_nn[1], (long)J_nn[2], (long)J_nn[3], C, J_P, (long)KS, (long)cl[0].n, (long)cl[1].n, (long)cl[2].n, (long)cl[3].n);
    int64_t *nl = calloc(KS, 8);
    for (int64_t x = 0; x < cl[0].n; x++) {
        int64_t cx = cl[0].start[x + 1] - cl[0].start[x];
        for (int64_t y = 0; y < cl[1].n; y++) {
            int sm[NP]; for (int j = 0; j < NP; j++) sm[j] = cl[0].v[x * NP + j] + cl[1].v[y * NP + j];
            nl[enc(sm)] += cx * (cl[1].start[y + 1] - cl[1].start[y]);
        }
    }
    {   /* keys that no probe pair can complete are not built */
        unsigned char *has = calloc(KS, 1);
        for (int64_t x = 0; x < cl[2].n; x++)
            for (int64_t y = 0; y < cl[3].n; y++) {
                int sm[NP]; for (int j = 0; j < NP; j++) sm[j] = J_Tp[j] - cl[2].v[x * NP + j] - cl[3].v[y * NP + j];
                int64_t id = enc(sm); if (id >= 0) has[id] = 1;
            }
        for (int64_t id = 0; id < KS; id++) if (!has[id]) nl[id] = 0;
        free(has);
    }
    double totl = 0; for (int64_t id = 0; id < KS; id++) totl += nl[id];
    if (stats) { stats[0] = totl; stats[1] = (double)J_nn[2] * J_nn[3]; stats[3] = J_P; }
    int64_t total = 0, printed = 0; double probed = 0;
    int64_t id0 = 0;
    Slot *tab = NULL; int64_t cap = 0; uint64_t *bm = NULL; int64_t bmcap = 0;
    int rc = 0; double lastp = now(), T_build = 0, T_probe = 0;
    const uint64_t GOLD2 = 0xC2B2AE3D27D4EB4FULL;
    while (id0 < KS) {
        while (id0 < KS && nl[id0] == 0) id0++;
        if (id0 >= KS) break;
        int64_t id1 = id0, acc = 0;
        while (id1 < KS && (acc == 0 || acc + nl[id1] <= budget)) { acc += nl[id1]; id1++; }
        int sb = 4; while ((1L << sb) < 2 * acc) sb++;
        int bb = 6; while ((1L << bb) < 16 * acc) bb++;
        int64_t slots = 1L << sb, bmw = (1L << bb) >> 6;
        if (slots > cap) {
            free(tab); cap = slots; tab = NULL;
            size_t bytes = (sizeof(Slot) * cap + (2u << 20) - 1) & ~(size_t)((2u << 20) - 1);
            if (posix_memalign((void **)&tab, 2u << 20, bytes) != 0) { tab = NULL; rc = -2; break; }
            madvise(tab, bytes, MADV_HUGEPAGE);
        }
        if (bmw > bmcap) { free(bm); bmcap = bmw; bm = malloc(8 * bmcap); if (!bm) { rc = -2; break; } }
        double tb0 = now();
        memset(tab, 0xFF, sizeof(Slot) * slots); memset(bm, 0, 8 * bmw);
        uint64_t smask = (uint64_t)slots - 1;
        int klo[NP], khi[NP]; dec(id0, klo); dec(id1 - 1, khi);
        for (int64_t x = 0; x < cl[0].n; x++) {
            int a[NP], b[NP];
            for (int j = 0; j < NP; j++) { a[j] = klo[j] - cl[0].v[x * NP + j]; b[j] = khi[j] - cl[0].v[x * NP + j]; }
            int64_t y0 = lower_c(&cl[1], a), y1 = upper_c(&cl[1], b);
            for (int64_t y = y0; y < y1; y++) {
                int sm[NP]; for (int j = 0; j < NP; j++) sm[j] = cl[0].v[x * NP + j] + cl[1].v[y * NP + j];
                int64_t id = enc(sm); if (id < id0 || id >= id1 || nl[id] == 0) continue;
                int64_t qs = cl[1].start[y], qe = cl[1].start[y + 1];
                for (int64_t p = cl[0].start[x]; p < cl[0].start[x + 1]; p++) {
                    uint32_t ia = cl[0].idx[p]; uint64_t hA = cl[0].hs[p];
                    for (int64_t q = qs; q < qe; q += PB) {
                        int nb = (qe - q < PB) ? (int)(qe - q) : PB; uint64_t hv[PB], sl[PB];
                        for (int t = 0; t < nb; t++) { hv[t] = hA + cl[1].hs[q + t]; sl[t] = (hv[t] * GOLD) >> (64 - sb); __builtin_prefetch(&tab[sl[t]], 1); }
                        for (int t = 0; t < nb; t++) {
                            uint64_t z = sl[t];
                            while (tab[z].ia != EMPTY) z = (z + 1) & smask;
                            tab[z].hv = hv[t]; tab[z].ia = ia; tab[z].ib = cl[1].idx[q + t];
                            uint64_t bi = (hv[t] * GOLD2) >> (64 - bb); bm[bi >> 6] |= 1ULL << (bi & 63);
                        }
                    }
                }
            }
        }
        double tb1 = now(); T_build += tb1 - tb0;
        double pr = 0;
        int stop = 0;
        for (int64_t x = 0; x < cl[2].n && !stop; x++) {
            int a[NP], b[NP];
            for (int j = 0; j < NP; j++) { a[j] = J_Tp[j] - cl[2].v[x * NP + j] - khi[j]; b[j] = J_Tp[j] - cl[2].v[x * NP + j] - klo[j]; }
            int64_t y0 = lower_c(&cl[3], a), y1 = upper_c(&cl[3], b);
            for (int64_t y = y0; y < y1 && !stop; y++) {
                int sm[NP]; for (int j = 0; j < NP; j++) sm[j] = J_Tp[j] - cl[2].v[x * NP + j] - cl[3].v[y * NP + j];
                int64_t id = enc(sm); if (id < id0 || id >= id1 || nl[id] == 0) continue;
                int64_t qs = cl[3].start[y], qe = cl[3].start[y + 1];
                pr += (double)(cl[2].start[x + 1] - cl[2].start[x]) * (double)(qe - qs);
                for (int64_t p = cl[2].start[x]; p < cl[2].start[x + 1] && !stop; p++) {
                    uint32_t ic = cl[2].idx[p]; uint64_t t = J_hT - cl[2].hs[p];
                    const uint64_t *hq = cl[3].hs;
                    for (int64_t q = qs; q < qe; q++) {
                        uint64_t want = t - hq[q];
                        uint64_t bi = (want * GOLD2) >> (64 - bb);
                        if (!((bm[bi >> 6] >> (bi & 63)) & 1)) continue;
                        uint64_t z = (want * GOLD) >> (64 - sb);
                        while (tab[z].ia != EMPTY) {
                            if (tab[z].hv == want) {
                                uint32_t ia = tab[z].ia, ib = tab[z].ib, ic_ = ic, id_ = cl[3].idx[q]; int ok = 1;
                                for (int c = 0; c < C && ok; c++)
                                    ok = J_V[0][(int64_t)ia * C + c] + J_V[1][(int64_t)ib * C + c] + J_V[2][(int64_t)ic_ * C + c] + J_V[3][(int64_t)id_ * C + c] == T[c];
                                if (ok) {
                                    total++;
                                    if (printed < maxout) {
                                        uint32_t o[4]; o[map[0]] = ia; o[map[1]] = ib; o[map[2]] = ic_; o[map[3]] = id_;
                                        memcpy(out + 4 * printed, o, 16); printed++;
                                    }
                                    if (stop_at > 0 && total >= stop_at) { stop = 1; break; }
                                }
                            }
                            z = (z + 1) & smask;
                        }
                        if (stop) break;
                    }
                }
            }
        }
        probed += pr; T_probe += now() - tb1;
        if (verbose && now() - lastp > 5.0) { lastp = now(); fprintf(stderr, "join4: keys up to %ld/%ld, build %.3g, probe pairs %.4g, solutions %ld [%.1fs]\n",
                (long)id1, (long)KS, (double)acc, probed, (long)total, now() - t0); }
        id0 = id1;
        if (stop) break;
    }
    if (verbose) fprintf(stderr, "join4: done, probe pairs %.4g, solutions %ld [%.1fs] build %.2fs probe %.2fs\n", probed, (long)total, now() - t0, T_build, T_probe);
    if (stats) stats[2] = probed;
    free(tab); free(bm); free(nl);
    for (int k = 0; k < 4; k++) { free_cells(&cl[k]); free(J_h[k]); }
    if (rc) return rc;
    return total;
}

/* ------------------------------------------------------------------------------------------------ item DP */
/* DP over "items", each with several options (delta vectors over D coordinates) and multiplicities.  Same pruning
 * as dp_run: fixed coordinates must end at tgt, pruned by the minimum / maximum the remaining items can still add.
 * Returns the number of final states (keys over outidx), -1 on state explosion, -2 on allocation failure. */
int64_t dp_items(int nitems, int D, const int64_t *optstart, const int16_t *delta, const double *mult,
                 const int32_t *fixed, const int32_t *tgt, int nout, const int32_t *outidx, int64_t max_states, int64_t *peak) {
    int Dp = (D + 3) & ~3;
    int32_t *rhi = calloc((size_t)(nitems + 1) * D, 4), *rlo = calloc((size_t)(nitems + 1) * D, 4);   /* sums over items s+1.. */
    for (int s = nitems - 1; s >= 0; s--) {
        if (s < nitems - 1) {
            int it = s + 1;
            for (int j = 0; j < D; j++) {
                int mx = INT32_MIN, mn = INT32_MAX;
                for (int64_t o = optstart[it]; o < optstart[it + 1]; o++) { int d = delta[o * D + j]; if (d > mx) mx = d; if (d < mn) mn = d; }
                rhi[(size_t)s * D + j] = rhi[(size_t)(s + 1) * D + j] + mx;
                rlo[(size_t)s * D + j] = rlo[(size_t)(s + 1) * D + j] + mn;
            }
        }
    }
    int64_t n = 1;
    int16_t *ck = calloc((size_t)Dp, 2); double *cc = malloc(8); cc[0] = 1.0;
    int16_t tmp[256]; int64_t pk = 1;
    for (int s = 0; s < nitems; s++) {
        int64_t no = optstart[s + 1] - optstart[s];
        int64_t ncap = n * no; if (ncap < 1) ncap = 1;
        /* the state count can grow up to n*no: bound the table by max_states */
        if (ncap > max_states * 2) ncap = max_states * 2;
        int64_t tsz = 16; while (tsz < 2 * ncap) tsz <<= 1;
        int16_t *nk = malloc((size_t)ncap * Dp * 2); double *nc = malloc((size_t)ncap * 8); int32_t *tab = malloc((size_t)tsz * 4);
        if (!nk || !nc || !tab) { free(nk); free(nc); free(tab); free(ck); free(cc); free(rhi); free(rlo); return -2; }
        memset(tab, 0xFF, (size_t)tsz * 4);
        int64_t nn = 0; uint64_t mask = (uint64_t)tsz - 1;
        const int32_t *hi = rhi + (size_t)s * D, *lo = rlo + (size_t)s * D;
        for (int64_t i = 0; i < n; i++) {
            for (int64_t o = optstart[s]; o < optstart[s + 1]; o++) {
                memcpy(tmp, ck + (size_t)i * Dp, (size_t)Dp * 2);
                int ok = 1; const int16_t *dl = delta + o * D;
                for (int j = 0; j < D; j++) {
                    int v = tmp[j] + dl[j]; tmp[j] = (int16_t)v;
                    if (fixed[j]) { int d = tgt[j] - v; if (d > hi[j] || d < lo[j]) { ok = 0; break; } }
                }
                if (!ok) continue;
                uint64_t h = hash_key(tmp, Dp) & mask;
                double w = cc[i] * mult[o];
                for (;;) {
                    int32_t e = tab[h];
                    if (e < 0) {
                        if (nn >= ncap) { free(nk); free(nc); free(tab); free(ck); free(cc); free(rhi); free(rlo); if (peak) *peak = nn; return -1; }
                        memcpy(nk + (size_t)nn * Dp, tmp, (size_t)Dp * 2); nc[nn] = w; tab[h] = (int32_t)nn; nn++; break;
                    }
                    if (memcmp(nk + (size_t)e * Dp, tmp, (size_t)Dp * 2) == 0) { nc[e] += w; break; }
                    h = (h + 1) & mask;
                }
            }
        }
        free(tab); free(ck); free(cc);
        ck = nk; cc = nc; n = nn;
        if (n > pk) pk = n;
        if (n > max_states) { free(ck); free(cc); free(rhi); free(rlo); if (peak) *peak = pk; return -1; }
    }
    free(rhi); free(rlo);
    free(g_k); free(g_c);
    g_nout = nout; g_n = n;
    g_k = malloc((size_t)(n ? n : 1) * (nout ? nout : 1) * 2); g_c = malloc((size_t)(n ? n : 1) * 8);
    for (int64_t i = 0; i < n; i++) {
        for (int q = 0; q < nout; q++) g_k[(size_t)i * nout + q] = ck[(size_t)i * Dp + outidx[q]];
        g_c[i] = cc[i];
    }
    free(ck); free(cc);
    if (peak) *peak = pk;
    return n;
}
