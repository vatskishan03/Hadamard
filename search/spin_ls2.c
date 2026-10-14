/*
 * spin_ls2: the tabu search of spin_ls.c with the B-swap scan vectorised over the partner orbit j.
 * Same input format, same moves and acceptance rule (best non-tabu swap or pair flip, 3% random move type).
 * Swap scores are evaluated in float (exact for the small residuals near a solution; the exact objective is
 * recomputed in integers after every move).
 *   cc -O3 -march=native -ffast-math spin_ls2.c -o spin_ls2 ; ./spin_ls2 in.bin restarts steps seed
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
    long restarts = atol(argv[2]), steps = atol(argv[3]);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[4]) * 0x2545F4914F6CDD1DULL;
    int OBp = (OB + 7) & ~7;
    /* float copies: Mf[c][i][j] (padded), diag */
    float *Mf = aligned_alloc(64, sizeof(float) * (size_t)C * OB * OBp);
    for (int c = 0; c < C; c++) for (int i = 0; i < OB; i++) for (int j = 0; j < OBp; j++)
        Mf[((size_t)c * OB + i) * OBp + j] = j < OB ? (float)MBc(c, i, j) : 0.f;
    float *wf = malloc(sizeof(float) * C); for (int c = 0; c < C; c++) wf[c] = (float)w[c];
    int *B = malloc(OB * 4), *A = malloc(NA * 4);
    long *gB = malloc((long)C * OB * 8), *gA = malloc((long)C * NA * 8), *R = malloc(C * 8), *PA = malloc(C * 8), *TB = malloc(C * 8);
    float *E = aligned_alloc(64, sizeof(float) * (size_t)C * OBp);       /* E[c][j] = 2 g_cj + 2 M_cjj (j currently -1) */
    float *fv = aligned_alloc(64, sizeof(float) * OBp);
    long *tabuB = calloc(OB, 8), *tabuA = calloc(NP, 8);
    int *plus = malloc(OB * 4);
    unsigned char *okj = malloc(OBp);
    long best_all = -1;
    for (long r = 0; r < restarts; r++) {
        memcpy(B, B0 + (long)(rnd() % NB0) * OB, OB * 4);
        for (int i = 1; i < OB; i++) { int j = 1 + rnd() % (OB - 1); if (szB[i] == szB[j]) { int t = B[i]; B[i] = B[j]; B[j] = t; } }
        for (int i = 0; i < NA; i++) A[i] = 0;
        A[0] = 1;
        for (int k = 0; k < NP; k++) { int s = (rnd() & 1) ? 1 : -1; A[pr[2 * k]] = s; A[pr[2 * k + 1]] = -s; }
        for (int c = 0; c < C; c++) {
            long t = 0;
            for (int i = 0; i < OB; i++) { long g = 0; for (int j = 0; j < OB; j++) g += MBc(c, i, j) * B[j]; gB[(long)c * OB + i] = g; t += g * B[i]; }
            TB[c] = t / 2; t = 0;
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
            int np = 0;
            for (int j = 0; j < OBp; j++) okj[j] = 0;
            for (int i = 1; i < OB; i++) { if (B[i] > 0) plus[np++] = i; else if (tabuB[i] <= it) okj[i] = 1; }
            for (int c = 0; c < C; c++) {
                float *Ec = E + (size_t)c * OBp;
                for (int j = 0; j < OB; j++) Ec[j] = (float)(2 * gB[(long)c * OB + j] + 2 * MBc(c, j, j));
                for (int j = OB; j < OBp; j++) Ec[j] = 0.f;
            }
            float bf = INFINITY; int bi = -1, bj = -1, ties = 0;
            for (int a = 0; a < np; a++) {
                int i = plus[a];
                if (tabuB[i] > it) continue;
                for (int j = 0; j < OBp; j++) fv[j] = 0.f;
                for (int c = 0; c < C; c++) {
                    float base = (float)R[c] + (float)(-2 * gB[(long)c * OB + i] + 2 * MBc(c, i, i));
                    const float *Ec = E + (size_t)c * OBp, *Mi = Mf + ((size_t)c * OB + i) * OBp;
                    float wc = wf[c];
                    for (int j = 0; j < OBp; j++) { float x = base + Ec[j] - 4.f * Mi[j]; fv[j] += wc * x * x; }
                }
                for (int j = 1; j < OB; j++) {
                    if (!okj[j] || szB[j] != szB[i]) continue;
                    float v = fv[j];
                    if (v < bf) { bf = v; bi = i; bj = j; ties = 1; }
                    else if (v == bf && (rnd() % ++ties) == 0) { bi = i; bj = j; }
                }
            }
            long af = -1; int ak = -1;
            for (int k = 0; k < NP; k++) {
                if (tabuA[k] > it) continue;
                int o = pr[2 * k], q = pr[2 * k + 1];
                long fvA = 0;
                for (int c = 0; c < C; c++) {
                    long d = -2 * A[o] * gA[(long)c * NA + o] + 2 * SAc(c, o, o) - 2 * A[q] * gA[(long)c * NA + q] + 2 * SAc(c, q, q)
                             + 4 * A[o] * A[q] * SAc(c, o, q);
                    long x = R[c] + d; fvA += (long)w[c] * x * x;
                }
                if (af < 0 || fvA < af) { af = fvA; ak = k; }
            }
            if (bi < 0 && ak < 0) continue;
            int takeB = (ak < 0) || (bi >= 0 && bf <= (float)af);
            if (rnd() % 100 < 3) takeB = (rnd() & 1) && bi >= 0;
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
