/*
 * Tabu search for spin-type skew Goethals-Seidel families (see spin_skew.py).
 *
 * State: B (signs on the OB H-orbits, fixed row sum) and A (signs on NP negation pairs of K-orbits).
 * Residual on K-class c:  R_c = P_A(c) + T_B(c),  P_A = A' SA_c A / 2,  T_B = B' MB_c B / 2
 * (SA_c, MB_c are the symmetrised doubled quadratic forms, integers). Objective sum_c w_c R_c^2.
 * Moves: swap a + and a - B-orbit of equal size (keeps the row sum); flip one A pair.
 *
 * input (int32): OB C NA NP ; w[C] ; sizeB[OB] ; MB[C][OB][OB] ; SA[C][NA][NA] ; pairs[NP][2] ;
 *                B0[OB] (a start with the right row sum) [; NB0 ; B0s[NB0][OB] (optional: more starts)]
 * Each restart takes one of the starts at random and shuffles it within orbits of equal size. With a single start
 * the numbers of + signs per orbit size stay those of B0 for the whole run, which for composite v (orbits of several
 * sizes) confines the search to one row-sum split; the extra starts (uniform samples of all B with the row sum)
 * remove that restriction.
 *   cc -O3 -march=native spin_ls.c -o spin_ls ;  ./spin_ls in.bin restarts steps seed
 * prints "SOL A <NA signs> B <OB signs>" on success.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int OB, C, NA, NP;
static int *w, *szB, *MB, *SA, *pr, *B0;
#define MBc(c, i, j) MB[((long)(c) * OB + (i)) * OB + (j)]
#define SAc(c, i, j) SA[((long)(c) * NA + (i)) * NA + (j)]
static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

static int *rd(FILE *f, long n) { int *p = malloc(n * sizeof(int)); if (fread(p, sizeof(int), n, f) != (size_t)n) { fprintf(stderr, "short read\n"); exit(1); } return p; }

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb");
    int hdr[4]; if (fread(hdr, 4, 4, f) != 4) return 1;
    OB = hdr[0]; C = hdr[1]; NA = hdr[2]; NP = hdr[3];
    w = rd(f, C); szB = rd(f, OB); MB = rd(f, (long)C * OB * OB); SA = rd(f, (long)C * NA * NA);
    pr = rd(f, 2L * NP); B0 = rd(f, OB);
    int NB0 = 1, nb0x;
    if (fread(&nb0x, 4, 1, f) == 1 && nb0x > 0) { free(B0); NB0 = nb0x; B0 = rd(f, (long)NB0 * OB); }
    fclose(f);
    long restarts = atol(argv[2]), steps = atol(argv[3]); rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[4]) * 0x2545F4914F6CDD1DULL;
    int *B = malloc(OB * 4), *A = malloc(NA * 4);
    long *gB = malloc((long)C * OB * 8), *gA = malloc((long)C * NA * 8), *R = malloc(C * 8), *PA = malloc(C * 8), *TB = malloc(C * 8);
    long *tabuB = calloc(OB, 8), *tabuA = calloc(NP, 8);
    int *plus = malloc(OB * 4), *minus = malloc(OB * 4);
    long best_all = -1;
    for (long r = 0; r < restarts; r++) {
        /* random B with the same row sum: shuffle B0 within equal-size orbits */
        memcpy(B, B0 + (long)(rnd() % NB0) * OB, OB * 4);
        for (int i = 1; i < OB; i++) {
            int j = 1 + rnd() % (OB - 1);
            if (szB[i] == szB[j]) { int t = B[i]; B[i] = B[j]; B[j] = t; }
        }
        for (int i = 0; i < NA; i++) A[i] = 0;
        A[0] = 1;
        for (int k = 0; k < NP; k++) { int s = (rnd() & 1) ? 1 : -1; A[pr[2 * k]] = s; A[pr[2 * k + 1]] = -s; }
        for (int c = 0; c < C; c++) {
            long t = 0;
            for (int i = 0; i < OB; i++) { long g = 0; for (int j = 0; j < OB; j++) g += MBc(c, i, j) * B[j]; gB[(long)c * OB + i] = g; t += g * B[i]; }
            TB[c] = t / 2;
            t = 0;
            for (int i = 0; i < NA; i++) { long g = 0; for (int j = 0; j < NA; j++) g += SAc(c, i, j) * A[j]; gA[(long)c * NA + i] = g; t += g * A[i]; }
            PA[c] = t / 2;
        }
        memset(tabuB, 0, OB * 8); memset(tabuA, 0, NP * 8);
        long best = -1;
        for (long it = 1; it <= steps; it++) {
            long f0 = 0;
            for (int c = 0; c < C; c++) { R[c] = PA[c] + TB[c]; f0 += (long)w[c] * R[c] * R[c]; }
            if (best < 0 || f0 < best) best = f0;
            if (f0 == 0) {
                printf("SOL A"); for (int i = 0; i < NA; i++) printf(" %d", A[i]);
                printf(" B"); for (int i = 0; i < OB; i++) printf(" %d", B[i]);
                printf("\n"); fflush(stdout);
                fprintf(stderr, "restart %ld step %ld: solution\n", r, it);
                return 0;
            }
            int np = 0, nm = 0;
            for (int i = 1; i < OB; i++) { if (B[i] > 0) plus[np++] = i; else minus[nm++] = i; }
            long bf = -1; int bi = -1, bj = -1, ties = 0;
            for (int a = 0; a < np; a++) {
                int i = plus[a];
                for (int b = 0; b < nm; b++) {
                    int j = minus[b];
                    if (szB[i] != szB[j]) continue;
                    if (tabuB[i] > it || tabuB[j] > it) continue;
                    long fv = 0;
                    for (int c = 0; c < C; c++) {
                        /* B_i = +1, B_j = -1 before the swap */
                        long d = (-2 * gB[(long)c * OB + i] + 2 * MBc(c, i, i)) + (2 * gB[(long)c * OB + j] + 2 * MBc(c, j, j)) - 4 * MBc(c, i, j);
                        long x = R[c] + d;
                        fv += (long)w[c] * x * x;
                    }
                    if (bf < 0 || fv < bf) { bf = fv; bi = i; bj = j; ties = 1; }
                    else if (fv == bf && (rnd() % ++ties) == 0) { bi = i; bj = j; }
                }
            }
            long af = -1; int ak = -1;
            for (int k = 0; k < NP; k++) {
                if (tabuA[k] > it) continue;
                int o = pr[2 * k], q = pr[2 * k + 1];
                long fv = 0;
                for (int c = 0; c < C; c++) {
                    /* flip o and q: delta of A'SA A/2 */
                    long d = -2 * A[o] * gA[(long)c * NA + o] + 2 * SAc(c, o, o) - 2 * A[q] * gA[(long)c * NA + q] + 2 * SAc(c, q, q)
                             + 4 * A[o] * A[q] * SAc(c, o, q);
                    long x = R[c] + d;
                    fv += (long)w[c] * x * x;
                }
                if (af < 0 || fv < af) { af = fv; ak = k; }
            }
            if (bi < 0 && ak < 0) continue;
            int takeB = (ak < 0) || (bi >= 0 && bf <= af);
            if (rnd() % 100 < 3) takeB = (rnd() & 1) && bi >= 0;   /* occasional random choice of move type */
            if (takeB && bi >= 0) {
                int xs[2] = {bi, bj};
                for (int u = 0; u < 2; u++) {
                    int x = xs[u];
                    for (int c = 0; c < C; c++) for (int i = 0; i < OB; i++) gB[(long)c * OB + i] -= 2L * MBc(c, i, x) * B[x];
                    B[x] = -B[x];
                }
                for (int c = 0; c < C; c++) { long t = 0; for (int i = 0; i < OB; i++) t += gB[(long)c * OB + i] * B[i]; TB[c] = t / 2; }
                tabuB[bi] = it + 3 + rnd() % 5; tabuB[bj] = it + 3 + rnd() % 5;
            } else if (ak >= 0) {
                int xs[2] = {pr[2 * ak], pr[2 * ak + 1]};
                for (int u = 0; u < 2; u++) {
                    int x = xs[u];
                    for (int c = 0; c < C; c++) for (int i = 0; i < NA; i++) gA[(long)c * NA + i] -= 2L * SAc(c, i, x) * A[x];
                    A[x] = -A[x];
                }
                for (int c = 0; c < C; c++) { long t = 0; for (int i = 0; i < NA; i++) t += gA[(long)c * NA + i] * A[i]; PA[c] = t / 2; }
                tabuA[ak] = it + 2 + rnd() % 3;
            }
        }
        if (best_all < 0 || best < best_all) best_all = best;
        if (r % 1000 == 0) fprintf(stderr, "restart %ld best %ld (overall %ld)\n", r, best, best_all);
    }
    fprintf(stderr, "no solution, best %ld\n", best_all);
    return 0;
}
