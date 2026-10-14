/*
 * Spin-type skew search by a walk with neighbour lookups (the method of METHOD.md for order 2284).
 *
 * All skew first blocks A (2^NP of them) are enumerated and their vectors -P_A (on the C classes of
 * <K,-1>) are stored in a hash table. A tabu walk moves B by swaps of two equal-size orbits of opposite
 * sign (keeping the row sum) and minimises sum_c w_c (T_B(c) - m)^2, m the mean of -P_A. At every step
 * the T vector of *every* swap neighbour is looked up in the table; a hit is a solution.
 *
 * input (int32): OB C NA NP ; w[C] ; sizeB[OB] ; MB[C][OB][OB] ; SA[C][NA][NA] ; pairs[NP][2] ; B0[OB]
 *   (same file as spin_ls.c; MB and SA are doubled symmetric forms)
 *   cc -O3 -march=native -fopenmp spin_walk.c -o spin_walk ;  ./spin_walk in.bin seconds seed
 * prints "SOL A <NA signs> B <OB signs>".
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int OB, C, NA, NP;
static int *w, *szB, *MB, *SA, *pr, *B0;
#define MBc(c, i, j) MB[((long)(c) * OB + (i)) * OB + (j)]
#define SAc(c, i, j) SA[((long)(c) * NA + (i)) * NA + (j)]
static int *rd(FILE *f, long n) { int *p = malloc(n * 4); if (fread(p, 4, n, f) != (size_t)n) { fprintf(stderr, "short read\n"); exit(1); } return p; }

/* hash table of -P_A vectors */
static long HS; static uint64_t *hkey; static int32_t *hval; static long *PAv;   /* PAv[a*C + c] */
static uint64_t R1[64];
static inline uint64_t hvec(const long *x) { uint64_t h = 0; for (int c = 0; c < C; c++) h += R1[c] * (uint64_t)x[c]; return h | 1; }
static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; return x; }

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb");
    int h[4]; if (fread(h, 4, 4, f) != 4) return 1;
    OB = h[0]; C = h[1]; NA = h[2]; NP = h[3];
    w = rd(f, C); szB = rd(f, OB); MB = rd(f, (long)C * OB * OB); SA = rd(f, (long)C * NA * NA);
    pr = rd(f, 2L * NP); B0 = rd(f, OB);
    fclose(f);
    double seconds = atof(argv[2]); uint64_t seed = (uint64_t)atol(argv[3]);
    for (int c = 0; c < C; c++) R1[c] = mix(0x1234567ULL + 7919ULL * c) | 1;
    long nA = 1L << NP;
    PAv = malloc(nA * C * sizeof(long));
    double mean[64] = {0};
    int *A = malloc(NA * 4);
    for (long a = 0; a < nA; a++) {
        for (int i = 0; i < NA; i++) A[i] = 0;
        A[0] = 1;
        for (int k = 0; k < NP; k++) { int s = (a >> k) & 1 ? -1 : 1; A[pr[2 * k]] = s; A[pr[2 * k + 1]] = -s; }
        for (int c = 0; c < C; c++) {
            long t = 0;
            for (int i = 0; i < NA; i++) { long g = 0; for (int j = 0; j < NA; j++) g += (long)SAc(c, i, j) * A[j]; t += g * A[i]; }
            PAv[a * C + c] = -(t / 2);
            mean[c] += -(t / 2.0) / nA;
        }
    }
    HS = 1; while (HS < 4 * nA) HS <<= 1;
    hkey = calloc(HS, 8); hval = malloc(HS * 4);
    for (long a = 0; a < nA; a++) {
        uint64_t k = hvec(PAv + a * C); long p = mix(k) & (HS - 1);
        while (hkey[p] && !(hkey[p] == k && !memcmp(PAv + (long)hval[p] * C, PAv + a * C, C * sizeof(long)))) p = (p + 1) & (HS - 1);
        if (!hkey[p]) { hkey[p] = k; hval[p] = (int32_t)a; }
    }
    fprintf(stderr, "%ld skew blocks hashed; walking\n", nA);
    volatile int done = 0;
    long total_lookups = 0;
#pragma omp parallel reduction(+:total_lookups)
    {
        uint64_t rs = 0x9E3779B97F4A7C15ULL * (seed + 1 + 104729 * omp_get_thread_num());
#define RND() (rs ^= rs << 13, rs ^= rs >> 7, rs ^= rs << 17, rs)
        int *B = malloc(OB * 4), *pl = malloc(OB * 4), *mi = malloc(OB * 4);
        long *g = malloc((long)C * OB * 8), *T = malloc(C * 8), *Tn = malloc(C * 8);
        long *tabu = calloc(OB, 8);
        double t0 = omp_get_wtime();
        long it = 0;
        while (!done && omp_get_wtime() - t0 < seconds) {
            /* restart */
            memcpy(B, B0, OB * 4);
            for (int i = OB - 1; i > 0; i--) { int j = 1 + (int)(RND() % (uint64_t)i); if (szB[i] == szB[j]) { int t = B[i]; B[i] = B[j]; B[j] = t; } }
            for (int c = 0; c < C; c++) {
                long t = 0;
                for (int i = 0; i < OB; i++) { long s = 0; for (int j = 0; j < OB; j++) s += (long)MBc(c, i, j) * B[j]; g[(long)c * OB + i] = s; t += s * B[i]; }
                T[c] = t / 2;
            }
            memset(tabu, 0, OB * 8);
            for (int step = 0; step < 4000 && !done; step++, it++) {
                int np = 0, nm = 0;
                for (int i = 1; i < OB; i++) { if (B[i] > 0) pl[np++] = i; else mi[nm++] = i; }
                double best = 1e300; int bi = -1, bj = -1, ties = 0;
                for (int a = 0; a < np; a++) for (int b = 0; b < nm; b++) {
                    int i = pl[a], j = mi[b];
                    if (szB[i] != szB[j]) continue;
                    double e = 0;
                    for (int c = 0; c < C; c++) {
                        long d = (-2 * g[(long)c * OB + i] + 2 * MBc(c, i, i)) + (2 * g[(long)c * OB + j] + 2 * MBc(c, j, j)) - 4 * MBc(c, i, j);
                        Tn[c] = T[c] + d / 2;
                        double y = Tn[c] - mean[c]; e += w[c] * y * y;
                    }
                    total_lookups++;
                    uint64_t k = hvec(Tn); long p = mix(k) & (HS - 1);
                    while (hkey[p]) {
                        if (hkey[p] == k && !memcmp(PAv + (long)hval[p] * C, Tn, C * sizeof(long))) {
#pragma omp critical
                            {
                                if (!done) {
                                    done = 1;
                                    long aidx = hval[p];
                                    B[i] = -B[i]; B[j] = -B[j];
                                    printf("SOL A 1");
                                    int *Af = calloc(NA, 4); Af[0] = 1;
                                    for (int q = 0; q < NP; q++) { int s = (aidx >> q) & 1 ? -1 : 1; Af[pr[2 * q]] = s; Af[pr[2 * q + 1]] = -s; }
                                    for (int q = 1; q < NA; q++) printf(" %d", Af[q]);
                                    printf(" B"); for (int q = 0; q < OB; q++) printf(" %d", B[q]);
                                    printf("\n"); fflush(stdout);
                                }
                            }
                            break;
                        }
                        p = (p + 1) & (HS - 1);
                    }
                    if (tabu[i] > it || tabu[j] > it) continue;
                    if (e < best) { best = e; bi = i; bj = j; ties = 1; }
                    else if (e == best && (RND() % ++ties) == 0) { bi = i; bj = j; }
                }
                if (bi < 0) break;
                int xs[2] = {bi, bj};
                for (int u = 0; u < 2; u++) {
                    int x = xs[u];
                    for (int c = 0; c < C; c++) for (int r = 0; r < OB; r++) g[(long)c * OB + r] -= 2L * MBc(c, r, x) * B[x];
                    B[x] = -B[x];
                }
                for (int c = 0; c < C; c++) { long t = 0; for (int r = 0; r < OB; r++) t += g[(long)c * OB + r] * B[r]; T[c] = t / 2; }
                tabu[bi] = it + 4 + RND() % 6; tabu[bj] = it + 4 + RND() % 6;
            }
        }
    }
    fprintf(stderr, "%s; %ld neighbour lookups\n", done ? "solution" : "no solution", total_lookups);
    return 0;
}
