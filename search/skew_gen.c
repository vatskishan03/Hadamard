/*
 * Low-energy candidates for one Goethals-Seidel role on a multiplier structure.
 *
 * x: signs on the O orbits with row sum r (orbit sizes sz). P_c(x) = x' S_c x / 2 (S_c symmetrised and
 * doubled, integer) is the periodic autocorrelation on class c. Each candidate is a local minimum of
 *     E(x) = sum_c w_c (P_c(x) - m)^2,     m = (r^2 - v) / (v - 1)   (the mean of P over shifts),
 * reached by steepest descent over swaps of two equal-size orbits of opposite sign from a random start.
 * Concentrating every role near its mean raises the chance that four candidates sum to zero.
 *
 * input (int32): O C v r ; sz[O] ; w[C] ; S[C][O][O] ; x0[O] (a vector with row sum r)
 *   cc -O3 -march=native -fopenmp skew_gen.c -o skew_gen
 *   ./skew_gen in.bin N maxsteps seed out.signs      -> N records of O int8 signs
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int O, C, V, Rr;
static int *sz, *w, *S, *x0;
#define Sc(c, i, j) S[((long)(c) * O + (i)) * O + (j)]

static int *rd(FILE *f, long n) { int *p = malloc(n * 4); if (fread(p, 4, n, f) != (size_t)n) { fprintf(stderr, "short read\n"); exit(1); } return p; }

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb");
    int h[4]; if (fread(h, 4, 4, f) != 4) return 1;
    O = h[0]; C = h[1]; V = h[2]; Rr = h[3];
    sz = rd(f, O); w = rd(f, C); S = rd(f, (long)C * O * O); x0 = rd(f, O);
    fclose(f);
    long N = atol(argv[2]); int maxsteps = atoi(argv[3]); uint64_t seed = (uint64_t)atol(argv[4]);
    double m = ((double)Rr * Rr - V) / (V - 1);
    int8_t *out = malloc((size_t)N * O);
#pragma omp parallel
    {
        uint64_t rs = 0x9E3779B97F4A7C15ULL * (seed + 1 + 7919 * omp_get_thread_num());
        int *x = malloc(O * 4); long *g = malloc((long)C * O * 8); double *R = malloc(C * 8);
        int *pl = malloc(O * 4), *mi = malloc(O * 4);
#pragma omp for schedule(dynamic, 64)
        for (long n = 0; n < N; n++) {
            memcpy(x, x0, O * 4);
            for (int i = O - 1; i > 0; i--) {     /* shuffle signs within equal-size orbits */
                rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
                int j = (int)(rs % (uint64_t)(i + 1));
                if (sz[i] == sz[j]) { int t = x[i]; x[i] = x[j]; x[j] = t; }
            }
            for (int c = 0; c < C; c++) {
                long t = 0;
                for (int i = 0; i < O; i++) { long s = 0; for (int j = 0; j < O; j++) s += (long)Sc(c, i, j) * x[j]; g[(long)c * O + i] = s; t += s * x[i]; }
                R[c] = t / 2 - m;
            }
            for (int step = 0; step < maxsteps; step++) {
                int np = 0, nm = 0;
                for (int i = 0; i < O; i++) { if (x[i] > 0) pl[np++] = i; else mi[nm++] = i; }
                double e0 = 0; for (int c = 0; c < C; c++) e0 += w[c] * R[c] * R[c];
                double best = e0; int bi = -1, bj = -1;
                for (int a = 0; a < np; a++) for (int b = 0; b < nm; b++) {
                    int i = pl[a], j = mi[b];
                    if (sz[i] != sz[j]) continue;
                    double e = 0;
                    for (int c = 0; c < C; c++) {
                        long d = (-2 * g[(long)c * O + i] + 2 * Sc(c, i, i)) + (2 * g[(long)c * O + j] + 2 * Sc(c, j, j)) - 4 * Sc(c, i, j);
                        double y = R[c] + d / 2.0;
                        e += w[c] * y * y;
                        if (e >= best) break;
                    }
                    if (e < best) { best = e; bi = i; bj = j; }
                }
                if (bi < 0) break;
                int xs[2] = {bi, bj};
                for (int u = 0; u < 2; u++) {
                    int k = xs[u];
                    for (int c = 0; c < C; c++) for (int i = 0; i < O; i++) g[(long)c * O + i] -= 2L * Sc(c, i, k) * x[k];
                    x[k] = -x[k];
                }
                for (int c = 0; c < C; c++) { long t = 0; for (int i = 0; i < O; i++) t += g[(long)c * O + i] * x[i]; R[c] = t / 2 - m; }
            }
            for (int i = 0; i < O; i++) out[(size_t)n * O + i] = (int8_t)x[i];
        }
        free(x); free(g); free(R); free(pl); free(mi);
    }
    FILE *fo = fopen(argv[5], "wb"); fwrite(out, 1, (size_t)N * O, fo); fclose(fo);
    return 0;
}
