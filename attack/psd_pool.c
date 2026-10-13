// All +-1 orbit vectors with a given row sum whose power spectrum stays below a bound (meet in the middle).
//
// A sequence constant on the orbits of a multiplier group has a Fourier transform X^(w) = sum_o y_o S_o(w), constant
// (up to conjugation) on the classes of frequencies; S_o(w) = sum_{x in orbit o} zeta^(w x) for one representative
// frequency per class. Two-block problems PAF_A + PAF_B = T (s != 0) need |A^(w)|^2 + |B^(w)|^2 = 2v - T for w != 0,
// so each block satisfies |X^(w)|^2 <= 2v - T (the PSD test); propus-type problems give bounds 4v and 2v.
//
// The orbits are split into two halves; every sign vector of each half is listed with its row-sum contribution and its
// transform; for each pair of halves with the right total row sum the full transform is tested class by class.
// Input on stdin (text): O C r bound, then O orbit sizes, then O*C pairs (re im) of S_o(w_c) (orbit-major).
// Output: the 64-bit masks (bit o set = orbit o has sign -1) of the surviving vectors, as raw uint64 to the file
// argv[1]; "DONE n" on stdout.
//   cc -O3 -march=native -fopenmp psd_pool.c -o psd_pool -lm
//   ./psd_pool out.bin < spec.txt
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef struct { int sum; uint32_t mask; } Half;
static int cmp(const void *a, const void *b) { const Half *x = a, *y = b; return (x->sum > y->sum) - (x->sum < y->sum); }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: psd_pool out.bin < spec\n"); return 2; }
    FILE *fo = fopen(argv[1], "wb"); if (!fo) return 1;
    int O, C, r; double bound;
    if (scanf("%d %d %d %lf", &O, &C, &r, &bound) != 4) return 1;
    if (O > 62) { fprintf(stderr, "too many orbits\n"); return 1; }
    int *sz = malloc(O * sizeof(int));
    for (int o = 0; o < O; o++) if (scanf("%d", &sz[o]) != 1) return 1;
    double *re = malloc((size_t)O * C * sizeof(double)), *im = malloc((size_t)O * C * sizeof(double));
    for (int o = 0; o < O; o++) for (int c = 0; c < C; c++) if (scanf("%lf %lf", &re[o * C + c], &im[o * C + c]) != 2) return 1;
    int h1 = O / 2, h2 = O - h1;
    long n1 = 1L << h1, n2 = 1L << h2;
    // half transforms: T[mask][c] (re, im), all signs +1 for bit 0
    double *T1 = malloc((size_t)n1 * C * 2 * sizeof(double)), *T2 = malloc((size_t)n2 * C * 2 * sizeof(double));
    Half *H1 = malloc(n1 * sizeof(Half)), *H2 = malloc(n2 * sizeof(Half));
    for (int half = 0; half < 2; half++) {
        int lo = half ? h1 : 0, m = half ? h2 : h1; long n = 1L << m;
        double *T = half ? T2 : T1; Half *H = half ? H2 : H1;
        #pragma omp parallel for schedule(static)
        for (long mask = 0; mask < n; mask++) {
            int s = 0; double *t = T + mask * C * 2;
            for (int c = 0; c < C; c++) { t[2 * c] = 0; t[2 * c + 1] = 0; }
            for (int j = 0; j < m; j++) {
                int o = lo + j, sg = (mask >> j & 1) ? -1 : 1;
                s += sg * sz[o];
                for (int c = 0; c < C; c++) { t[2 * c] += sg * re[o * C + c]; t[2 * c + 1] += sg * im[o * C + c]; }
            }
            H[mask].sum = s; H[mask].mask = (uint32_t)mask;
        }
        qsort(H, n, sizeof(Half), cmp);
    }
    long total = 0, j0 = 0;
    // walk H1 ascending and H2 for the complementary sum
    long i = 0;
    while (i < n1) {
        long i2 = i; while (i2 < n1 && H1[i2].sum == H1[i].sum) i2++;
        int need = r - H1[i].sum;
        // binary search in H2 for need
        long lo = 0, hi = n2;
        while (lo < hi) { long md = (lo + hi) / 2; if (H2[md].sum < need) lo = md + 1; else hi = md; }
        long k0 = lo; long k1 = k0; while (k1 < n2 && H2[k1].sum == need) k1++;
        if (k1 > k0) {
            #pragma omp parallel for schedule(dynamic, 64) reduction(+:total)
            for (long a = i; a < i2; a++) {
                const double *ta = T1 + (size_t)H1[a].mask * C * 2;
                for (long b = k0; b < k1; b++) {
                    const double *tb = T2 + (size_t)H2[b].mask * C * 2;
                    int ok = 1;
                    for (int c = 0; c < C; c++) {
                        double x = ta[2 * c] + tb[2 * c], y = ta[2 * c + 1] + tb[2 * c + 1];
                        if (x * x + y * y > bound + 1e-6) { ok = 0; break; }
                    }
                    if (ok) {
                        uint64_t full = (uint64_t)H1[a].mask | ((uint64_t)H2[b].mask << h1);
                        #pragma omp critical
                        fwrite(&full, 8, 1, fo);
                        total++;
                    }
                }
            }
        }
        i = i2;
    }
    (void)j0;
    fclose(fo);
    printf("DONE %ld\n", total);
    return 0;
}
