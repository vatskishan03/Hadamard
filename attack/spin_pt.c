/*
 * Parallel tempering for spin-type skew Goethals-Seidel families (the problem of spin_ls.c, same input file).
 *
 * State per replica: B (signs on the OB H-orbits, fixed row sum) and A (signs on the NP negation pairs of K-orbits).
 * Energy sum_c w_c R_c^2 with R_c = P_A(c) + T_B(c) (spin_ls.c). Moves: swap a + and a - B-orbit of equal size, or
 * flip an A pair; each proposal costs O(C), an accepted one O(C * OB). Metropolis at temperature T_k, replicas at
 * geometric temperatures between Tmin and Tmax, exchange of neighbours after every sweep (OB proposals per replica).
 *
 *   cc -O3 -march=native spin_pt.c -o spin_pt -lm
 *   ./spin_pt in.bin seconds seed [replicas Tmin Tmax]
 * prints "SOL A <NA signs> B <OB signs>" on success.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int OB, C, NA, NP;
static int *w, *szB, *MB, *SA, *pr, *B0;
static int NB0 = 1;
#define MBc(c, i, j) MB[((long)(c) * OB + (i)) * OB + (j)]
#define SAc(c, i, j) SA[((long)(c) * NA + (i)) * NA + (j)]
static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static int *rd(FILE *f, long n) { int *p = malloc(n * sizeof(int)); if (fread(p, sizeof(int), n, f) != (size_t)n) { fprintf(stderr, "short read\n"); exit(1); } return p; }

typedef struct { int *B, *A; long *gB, *gA, *R; long E; } Rep;

static void init_rep(Rep *r) {
    memcpy(r->B, B0 + (long)(rnd() % NB0) * OB, OB * 4);
    for (int i = 1; i < OB; i++) { int j = 1 + rnd() % (OB - 1); if (szB[i] == szB[j]) { int t = r->B[i]; r->B[i] = r->B[j]; r->B[j] = t; } }
    for (int i = 0; i < NA; i++) r->A[i] = 0;
    r->A[0] = 1;
    for (int k = 0; k < NP; k++) { int s = (rnd() & 1) ? 1 : -1; r->A[pr[2 * k]] = s; r->A[pr[2 * k + 1]] = -s; }
    r->E = 0;
    for (int c = 0; c < C; c++) {
        long t = 0, u = 0;
        for (int i = 0; i < OB; i++) { long g = 0; for (int j = 0; j < OB; j++) g += MBc(c, i, j) * r->B[j]; r->gB[(long)c * OB + i] = g; t += g * r->B[i]; }
        for (int i = 0; i < NA; i++) { long g = 0; for (int j = 0; j < NA; j++) g += SAc(c, i, j) * r->A[j]; r->gA[(long)c * NA + i] = g; u += g * r->A[i]; }
        r->R[c] = t / 2 + u / 2;
        r->E += (long)w[c] * r->R[c] * r->R[c];
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: spin_pt in.bin seconds seed [replicas Tmin Tmax]\n"); return 2; }
    FILE *f = fopen(argv[1], "rb");
    int hdr[4]; if (fread(hdr, 4, 4, f) != 4) return 1;
    OB = hdr[0]; C = hdr[1]; NA = hdr[2]; NP = hdr[3];
    w = rd(f, C); szB = rd(f, OB); MB = rd(f, (long)C * OB * OB); SA = rd(f, (long)C * NA * NA);
    pr = rd(f, 2L * NP); B0 = rd(f, OB);
    int nb0x; if (fread(&nb0x, 4, 1, f) == 1 && nb0x > 0) { free(B0); NB0 = nb0x; B0 = rd(f, (long)NB0 * OB); }
    fclose(f);
    double secs = atof(argv[2]); rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[3]) * 0x2545F4914F6CDD1DULL;
    int NR = argc > 4 ? atoi(argv[4]) : 24; double Tmin = argc > 5 ? atof(argv[5]) : 20, Tmax = argc > 6 ? atof(argv[6]) : 20000;
    double *T = malloc(NR * sizeof(double));
    for (int k = 0; k < NR; k++) T[k] = Tmin * pow(Tmax / Tmin, NR > 1 ? (double)k / (NR - 1) : 0);
    Rep *rp = malloc(NR * sizeof(Rep)); int *perm = malloc(NR * sizeof(int));    /* perm[k]: replica at temperature k */
    for (int k = 0; k < NR; k++) {
        rp[k].B = malloc(OB * 4); rp[k].A = malloc(NA * 4); rp[k].gB = malloc((long)C * OB * 8); rp[k].gA = malloc((long)C * NA * 8);
        rp[k].R = malloc(C * 8); init_rep(&rp[k]); perm[k] = k;
    }
    long *d = malloc(C * 8);
    int *plus = malloc(OB * 4), *minus = malloc(OB * 4);
    clock_t t0 = clock(); long sweeps = 0, best = -1, acc = 0, prop = 0, xacc = 0;
    long *xprop = calloc(NR, 8), *xok = calloc(NR, 8);
    while ((double)(clock() - t0) / CLOCKS_PER_SEC < secs) {
        for (int k = 0; k < NR; k++) {
            Rep *r = &rp[perm[k]]; double beta = 1.0 / T[k];
            for (int it = 0; it < OB; it++) {
                prop++;
                if (NP > 0 && (rnd() % 8) == 0) {                       /* flip an A pair */
                    int kk = rnd() % NP, o = pr[2 * kk], q = pr[2 * kk + 1];
                    long dE = 0;
                    for (int c = 0; c < C; c++) {
                        long dd = (-2L * r->A[o] * r->gA[(long)c * NA + o] + 2L * SAc(c, o, o) - 2L * r->A[q] * r->gA[(long)c * NA + q]
                                   + 2L * SAc(c, q, q) + 4L * r->A[o] * r->A[q] * SAc(c, o, q));
                        d[c] = dd; long x = r->R[c] + dd; dE += (long)w[c] * (x * x - r->R[c] * r->R[c]);
                    }
                    if (dE <= 0 || urand() < exp(-beta * dE)) {
                        int xs[2] = {o, q};
                        for (int u = 0; u < 2; u++) { int x = xs[u]; for (int c = 0; c < C; c++) for (int i = 0; i < NA; i++) r->gA[(long)c * NA + i] -= 2L * SAc(c, i, x) * r->A[x]; r->A[x] = -r->A[x]; }
                        for (int c = 0; c < C; c++) r->R[c] += d[c];
                        r->E += dE; acc++;
                    }
                } else {                                                  /* swap a + and a - orbit of equal size */
                    int i = 1 + rnd() % (OB - 1), j = 1 + rnd() % (OB - 1);
                    if (r->B[i] == r->B[j] || szB[i] != szB[j]) continue;
                    if (r->B[i] < 0) { int t = i; i = j; j = t; }           /* B_i = +1, B_j = -1 */
                    long dE = 0;
                    for (int c = 0; c < C; c++) {
                        long dd = ((-2L * r->gB[(long)c * OB + i] + 2L * MBc(c, i, i)) + (2L * r->gB[(long)c * OB + j] + 2L * MBc(c, j, j)) - 4L * MBc(c, i, j));
                        d[c] = dd; long x = r->R[c] + dd; dE += (long)w[c] * (x * x - r->R[c] * r->R[c]);
                    }
                    if (dE <= 0 || urand() < exp(-beta * dE)) {
                        int xs[2] = {i, j};
                        for (int u = 0; u < 2; u++) { int x = xs[u]; for (int c = 0; c < C; c++) for (int y = 0; y < OB; y++) r->gB[(long)c * OB + y] -= 2L * MBc(c, y, x) * r->B[x]; r->B[x] = -r->B[x]; }
                        for (int c = 0; c < C; c++) r->R[c] += d[c];
                        r->E += dE; acc++;
                    }
                }
                if (r->E == 0) {
                    printf("SOL A"); for (int a = 0; a < NA; a++) printf(" %d", r->A[a]);
                    printf(" B"); for (int b = 0; b < OB; b++) printf(" %d", r->B[b]);
                    printf("\n"); fflush(stdout);
                    fprintf(stderr, "solution after %ld sweeps, %.1f s\n", sweeps, (double)(clock() - t0) / CLOCKS_PER_SEC);
                    return 0;
                }
                if (best < 0 || r->E < best) best = r->E;
            }
        }
        for (int k = (sweeps & 1); k + 1 < NR; k += 2) {                  /* replica exchange */
            Rep *a = &rp[perm[k]], *b = &rp[perm[k + 1]];
            double x = (1.0 / T[k] - 1.0 / T[k + 1]) * (double)(a->E - b->E);
            xprop[k]++;
            if (x >= 0 || urand() < exp(x)) { int t = perm[k]; perm[k] = perm[k + 1]; perm[k + 1] = t; xacc++; xok[k]++; }
        }
        sweeps++;
        if (sweeps % 20000 == 0) {
            fprintf(stderr, "sweeps %ld best %ld cold %ld acc %.3f [%.0f s] x:", sweeps, best, rp[perm[0]].E, (double)acc / prop,
                    (double)(clock() - t0) / CLOCKS_PER_SEC);
            for (int k = 0; k + 1 < NR; k++) fprintf(stderr, " %.2f", xprop[k] ? (double)xok[k] / xprop[k] : 0.0);
            fprintf(stderr, "\n");
        }
    }
    fprintf(stderr, "no solution, best %ld, sweeps %ld\n", best, sweeps);
    return 0;
}
