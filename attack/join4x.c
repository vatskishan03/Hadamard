/*
 * Exhaustive four-list join: every (a, b, c, d) with VA[a] + VB[b] + VC[c] + VD[d] = T on all C classes.
 *
 * Pairs (a, b) are the build side and pairs (c, d) the probe side; the sides are swapped when the build
 * side would be the larger one. Both sides are partitioned exactly on the values of up to three classes
 * (those with the largest spread). The build pairs of each partition key go into a small open-addressing
 * table keyed by a 64-bit linear hash of all classes, and every probe pair is looked up in the table of
 * its complementary key, so the tables stay in cache. Partition keys are processed in lexicographic
 * chunks that bound the memory of the build side; for a chunk, the probe cells that can complete it form
 * a lexicographic interval (the order is translation invariant), found by binary search. Hash matches
 * are re-checked exactly. Nothing is sampled: TOTAL is the exact number of solutions.
 *
 * input: int32 nA nB nC nD C; int16 T[C]; int16 VA[nA][C] VB[nB][C] VC[nC][C] VD[nD][C]   (as find4b)
 *   cc -O3 -march=native -fopenmp join4x.c -o join4x ;  ./join4x input.bin [K] [budget]
 * Prints "SOL a b c d" for the first K solutions (default 20) and "TOTAL n"; progress on stderr.
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NP 3
static int C, P;                        /* classes, partition classes used */
static int16_t *T, *V[4];
static long nn[4];
static uint64_t R[1024], *h[4], hT;
static int pc[NP], lo[NP], hi[NP];      /* partition classes and their value ranges over all lists */
static int KW[NP];                      /* ranges of pair sums: 2W - 1 */
static int Tp[NP];

typedef struct { long n; int *v; long *start; uint32_t *idx; } Cells;    /* cells sorted lexicographically */

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static int lexcmp(const int *a, const int *b) {
    for (int j = 0; j < P; j++) if (a[j] != b[j]) return a[j] < b[j] ? -1 : 1;
    return 0;
}

static int16_t pv(int k, long i, int j) { return V[k][i * C + pc[j]]; }

static int *cellv_tmp; static int cmp_tmp(const void *x, const void *y) {
    return lexcmp(cellv_tmp + (long)(*(const uint32_t *)x) * NP, cellv_tmp + (long)(*(const uint32_t *)y) * NP);
}

static Cells cells(int k) {
    long n = nn[k];
    uint32_t *ord = malloc(4 * (n + 1));
    cellv_tmp = malloc(sizeof(int) * NP * (n + 1));
    for (long i = 0; i < n; i++) { ord[i] = i; for (int j = 0; j < NP; j++) cellv_tmp[i * NP + j] = j < P ? pv(k, i, j) : 0; }
    qsort(ord, n, 4, cmp_tmp);
    Cells c = {0, malloc(sizeof(int) * NP * (n + 1)), malloc(8 * (n + 2)), ord};
    for (long i = 0; i < n; i++) {
        int *vi = cellv_tmp + (long)ord[i] * NP;
        if (c.n == 0 || lexcmp(c.v + (c.n - 1) * NP, vi) != 0) {
            memcpy(c.v + c.n * NP, vi, sizeof(int) * NP); c.start[c.n] = i; c.n++;
        }
    }
    c.start[c.n] = n;
    free(cellv_tmp);
    return c;
}

/* first cell index with v >= key (lex) */
static long lower(const Cells *c, const int *key) {
    long a = 0, b = c->n;
    while (a < b) { long m = (a + b) / 2; if (lexcmp(c->v + m * NP, key) < 0) a = m + 1; else b = m; }
    return a;
}
static long upper(const Cells *c, const int *key) {
    long a = 0, b = c->n;
    while (a < b) { long m = (a + b) / 2; if (lexcmp(c->v + m * NP, key) <= 0) a = m + 1; else b = m; }
    return a;
}

static long enc(const int *s) {           /* key triple (pair sum) -> id, or -1 when out of range */
    long id = 0;
    for (int j = 0; j < P; j++) {
        int u = s[j] - 2 * lo[j];
        if (u < 0 || u >= KW[j]) return -1;
        id = id * KW[j] + u;
    }
    return id;
}
static void dec(long id, int *s) {
    for (int j = P - 1; j >= 0; j--) { s[j] = (int)(id % KW[j]) + 2 * lo[j]; id /= KW[j]; }
    for (int j = P; j < NP; j++) s[j] = 0;
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"); int hdr[5];
    if (!f || fread(hdr, 4, 5, f) != 5) return 1;
    for (int k = 0; k < 4; k++) nn[k] = hdr[k];
    C = hdr[4];
    long K = argc > 2 ? atol(argv[2]) : 20;
    long budget = argc > 3 ? atol(argv[3]) : 8000000;      /* build pairs per chunk */
    T = malloc(2 * C); if (fread(T, 2, C, f) != (size_t)C) return 1;
    for (int k = 0; k < 4; k++) {
        V[k] = malloc(2 * nn[k] * C + 2);
        if (fread(V[k], 2, nn[k] * C, f) != (size_t)(nn[k] * C)) return 1;
    }
    fclose(f);
    for (int k = 0; k < 4; k++) if (nn[k] == 0) { printf("TOTAL 0\n"); return 0; }
    int map[4] = {0, 1, 2, 3};
    if ((double)nn[0] * nn[1] > (double)nn[2] * nn[3]) {          /* build side = smaller pair product */
        int16_t *tv; long tn;
        tv = V[0]; V[0] = V[2]; V[2] = tv; tn = nn[0]; nn[0] = nn[2]; nn[2] = tn;
        tv = V[1]; V[1] = V[3]; V[3] = tv; tn = nn[1]; nn[1] = nn[3]; nn[3] = tn;
        map[0] = 2; map[1] = 3; map[2] = 0; map[3] = 1;
    }
    uint64_t s = 0x9E3779B97F4A7C15ULL;
    for (int c = 0; c < C; c++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; R[c] = s | 1; }
    hT = 0; for (int c = 0; c < C; c++) hT += R[c] * (uint64_t)(int64_t)T[c];
    for (int k = 0; k < 4; k++) {
        h[k] = malloc(8 * (nn[k] + 1));
        for (long i = 0; i < nn[k]; i++) {
            uint64_t x = 0; for (int c = 0; c < C; c++) x += R[c] * (uint64_t)(int64_t)V[k][i * C + c];
            h[k][i] = x;
        }
    }
    /* partition classes: largest pooled variance first, while the key space stays below 2^24 */
    double var[1024]; int used[1024] = {0};
    for (int c = 0; c < C; c++) {
        var[c] = 0;
        for (int k = 0; k < 4; k++) {
            double m = 0, m2 = 0;
            for (long i = 0; i < nn[k]; i++) { double x = V[k][i * C + c]; m += x; m2 += x * x; }
            m /= nn[k]; var[c] += m2 / nn[k] - m * m;
        }
    }
    P = 0; double space = 1;
    while (P < NP && P < C) {
        int best = -1;
        for (int c = 0; c < C; c++) if (!used[c] && (best < 0 || var[c] > var[best])) best = c;
        int l = 1 << 30, u = -(1 << 30);
        for (int k = 0; k < 4; k++) for (long i = 0; i < nn[k]; i++) {
            int x = V[k][i * C + best]; if (x < l) l = x; if (x > u) u = x;
        }
        if (space * (2.0 * (u - l + 1) - 1) > (double)(1 << 24)) break;
        used[best] = 1; pc[P] = best; lo[P] = l; hi[P] = u; KW[P] = 2 * (u - l + 1) - 1; space *= KW[P]; P++;
    }
    for (int j = 0; j < NP; j++) Tp[j] = j < P ? T[pc[j]] : 0;
    long KS = (long)space;
    Cells cl[4]; for (int k = 0; k < 4; k++) cl[k] = cells(k);
    fprintf(stderr, "join4x: n = %ld %ld %ld %ld, C = %d, partition on %d classes (key space %ld), cells %ld %ld %ld %ld\n",
            nn[0], nn[1], nn[2], nn[3], C, P, KS, cl[0].n, cl[1].n, cl[2].n, cl[3].n);
    double t0 = now();
    /* build-pair counts per key */
    long *nl = calloc(KS, 8);
    for (long x = 0; x < cl[0].n; x++) {
        long cx = cl[0].start[x + 1] - cl[0].start[x];
        for (long y = 0; y < cl[1].n; y++) {
            int sm[NP]; for (int j = 0; j < NP; j++) sm[j] = cl[0].v[x * NP + j] + cl[1].v[y * NP + j];
            nl[enc(sm)] += cx * (cl[1].start[y + 1] - cl[1].start[y]);
        }
    }
    double totl = 0; for (long id = 0; id < KS; id++) totl += nl[id];
    double totr = (double)nn[2] * nn[3];
    fprintf(stderr, "join4x: build pairs %.4g, probe pairs at most %.4g [%.1fs]\n", totl, totr, now() - t0);
    long *gofs = malloc(8 * KS); unsigned char *gbits = malloc(KS);
    for (long id = 0; id < KS; id++) gofs[id] = -1;
    long total = 0, printed = 0; double probed = 0;
    long id0 = 0, chunk = 0;
    uint64_t *sh = NULL; uint32_t *sa = NULL, *sb = NULL; long cap = 0;
    while (id0 < KS) {
        while (id0 < KS && nl[id0] == 0) id0++;
        if (id0 >= KS) break;
        long id1 = id0, acc = 0, slots = 0;                     /* chunk = keys id0 .. id1-1 */
        while (id1 < KS && (acc == 0 || acc + nl[id1] <= budget)) {
            if (nl[id1]) { int b = 1; while ((1L << b) < 2 * nl[id1]) b++; gofs[id1] = slots; gbits[id1] = b; slots += 1L << b; acc += nl[id1]; }
            id1++;
        }
        if (slots > cap) {
            free(sh); free(sa); free(sb); cap = slots;
            sh = malloc(8 * cap); sa = malloc(4 * cap); sb = malloc(4 * cap);
            if (!sh || !sa || !sb) { fprintf(stderr, "join4x: out of memory (%ld slots)\n", cap); return 2; }
        }
        memset(sa, 0xFF, 4 * slots);
        int klo[NP], khi[NP]; dec(id0, klo); dec(id1 - 1, khi);
        /* insert the build pairs of this chunk (threads claim slots with a compare-and-swap) */
        #pragma omp parallel for schedule(dynamic, 1)
        for (long x = 0; x < cl[0].n; x++) {
            int a[NP], b[NP];
            for (int j = 0; j < NP; j++) { a[j] = klo[j] - cl[0].v[x * NP + j]; b[j] = khi[j] - cl[0].v[x * NP + j]; }
            long y0 = lower(&cl[1], a), y1 = upper(&cl[1], b);
            for (long y = y0; y < y1; y++) {
                int sm[NP]; for (int j = 0; j < NP; j++) sm[j] = cl[0].v[x * NP + j] + cl[1].v[y * NP + j];
                long id = enc(sm); if (id < id0 || id >= id1) continue;
                long base = gofs[id]; int bits = gbits[id]; uint64_t mask = (1UL << bits) - 1;
                for (long p = cl[0].start[x]; p < cl[0].start[x + 1]; p++) {
                    uint32_t ia = cl[0].idx[p];
                    for (long q = cl[1].start[y]; q < cl[1].start[y + 1]; q++) {
                        uint32_t ib = cl[1].idx[q];
                        uint64_t hv = h[0][ia] + h[1][ib];
                        uint64_t sl = (hv * 0x9E3779B97F4A7C15ULL) >> (64 - bits);
                        for (;;) {
                            uint32_t empty = 0xFFFFFFFFu;
                            if (__atomic_compare_exchange_n(&sa[base + sl], &empty, ia, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                                sh[base + sl] = hv; sb[base + sl] = ib; break;
                            }
                            sl = (sl + 1) & mask;
                        }
                    }
                }
            }
        }
        /* probe: cells x of list 2, cells y of list 3 with key T - x - y in [klo, khi] */
        double pr = 0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:total, pr)
        for (long x = 0; x < cl[2].n; x++) {
            int a[NP], b[NP];
            for (int j = 0; j < NP; j++) { a[j] = Tp[j] - cl[2].v[x * NP + j] - khi[j]; b[j] = Tp[j] - cl[2].v[x * NP + j] - klo[j]; }
            long y0 = lower(&cl[3], a), y1 = upper(&cl[3], b);
            for (long y = y0; y < y1; y++) {
                int sm[NP]; for (int j = 0; j < NP; j++) sm[j] = Tp[j] - cl[2].v[x * NP + j] - cl[3].v[y * NP + j];
                long id = enc(sm); if (id < id0 || id >= id1 || nl[id] == 0) continue;
                long base = gofs[id]; int bits = gbits[id]; uint64_t mask = (1UL << bits) - 1;
                pr += (double)(cl[2].start[x + 1] - cl[2].start[x]) * (cl[3].start[y + 1] - cl[3].start[y]);
                for (long p = cl[2].start[x]; p < cl[2].start[x + 1]; p++) {
                    uint32_t ic = cl[2].idx[p]; uint64_t t = hT - h[2][ic];
                    for (long q = cl[3].start[y]; q < cl[3].start[y + 1]; q++) {
                        uint32_t id_ = cl[3].idx[q];
                        uint64_t want = t - h[3][id_];
                        uint64_t sl = (want * 0x9E3779B97F4A7C15ULL) >> (64 - bits);
                        while (sa[base + sl] != 0xFFFFFFFFu) {
                            if (sh[base + sl] == want) {
                                uint32_t ia = sa[base + sl], ib = sb[base + sl]; int ok = 1;
                                for (int c = 0; c < C && ok; c++)
                                    ok = V[0][(long)ia * C + c] + V[1][(long)ib * C + c] + V[2][(long)ic * C + c] + V[3][(long)id_ * C + c] == T[c];
                                if (ok) {
                                    total++;
                                    #pragma omp critical
                                    {
                                        if (printed < K) {
                                            long o[4]; o[map[0]] = ia; o[map[1]] = ib; o[map[2]] = ic; o[map[3]] = id_;
                                            printf("SOL %ld %ld %ld %ld\n", o[0], o[1], o[2], o[3]); fflush(stdout); printed++;
                                        }
                                    }
                                }
                            }
                            sl = (sl + 1) & mask;
                        }
                    }
                }
            }
        }
        probed += pr;
        for (long id = id0; id < id1; id++) gofs[id] = -1;
        chunk++;
        if (chunk % 16 == 1 || id1 >= KS)
            fprintf(stderr, "join4x: chunk %ld, keys up to %ld/%ld, build %.3g, probe pairs %.4g, solutions %ld [%.1fs]\n",
                    chunk, id1, KS, (double)acc, probed, total, now() - t0);
        id0 = id1;
    }
    fprintf(stderr, "join4x: done, probe pairs %.4g, solutions %ld [%.1fs]\n", probed, total, now() - t0);
    printf("TOTAL %ld\n", total);
    return 0;
}
