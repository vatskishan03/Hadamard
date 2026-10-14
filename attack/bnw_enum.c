/*
 * Enumerate every +-1 vector y on O orbits (Gray code) with a given row sum and write its class values
 *     v_c(y) = sum_{i,j} S_c[i][j] y_i y_j        (c = 0 .. C-1)
 * as int16 rows, for the bordered near-Williamson search (attack/bnw.py).
 *
 * Flipping y_j changes v_c by -2 y_j g_c[j], where g_c[j] = sum_{i != j} (S_c[i][j] + S_c[j][i]) y_i; the g_c are
 * updated in O(C * O) per flip. y_0 can be fixed to +1 (fix0 = 1), which removes the global sign when the row
 * sum is 0.
 *
 * input  (binary): int32 O, C, rowsum, fix0; int32 sizes[O]; int32 S[C][O][O]
 * output (binary): int16 rows of C values, then the sign vector of each row as a uint64 bitmask (bit i set = -1)
 *   cc -O3 -march=native bnw_enum.c -o bnw_enum ; ./bnw_enum in.bin vals.bin masks.bin [cap]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: bnw_enum in.bin vals.bin masks.bin [cap]\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    int32_t hdr[4];
    if (fread(hdr, 4, 4, f) != 4) return 2;
    int O = hdr[0], C = hdr[1], r = hdr[2], fix0 = hdr[3];
    if (O > 62) { fprintf(stderr, "too many orbits\n"); return 2; }
    int32_t *sz = malloc(4 * O), *S = malloc(4L * C * O * O);
    if (fread(sz, 4, O, f) != (size_t)O || fread(S, 4, (size_t)C * O * O, f) != (size_t)C * O * O) return 2;
    fclose(f);
    long cap = argc > 4 ? atol(argv[4]) : -1;
    int *Q = malloc(sizeof(int) * C * O * O);           /* Q[c][i][j] = S[c][i][j] + S[c][j][i] (i != j) */
    for (int c = 0; c < C; c++)
        for (int i = 0; i < O; i++)
            for (int j = 0; j < O; j++)
                Q[(c * O + i) * O + j] = i == j ? 0 : S[(c * O + i) * O + j] + S[(c * O + j) * O + i];
    int y[64]; long v[1024]; long g[1024 * 64];
    for (int i = 0; i < O; i++) y[i] = 1;
    long rs = 0;
    for (int i = 0; i < O; i++) rs += sz[i];
    for (int c = 0; c < C; c++) {
        long s = 0;
        for (int i = 0; i < O; i++) for (int j = 0; j < O; j++) s += S[(c * O + i) * O + j];
        v[c] = s;
        for (int j = 0; j < O; j++) { long t = 0; for (int i = 0; i < O; i++) t += Q[(c * O + i) * O + j]; g[c * O + j] = t; }
    }
    FILE *fv = fopen(argv[2], "wb"), *fm = fopen(argv[3], "wb");
    int16_t *row = malloc(2 * C);
    uint64_t mask = 0; long nout = 0;
    int first = fix0 ? 1 : 0;                              /* orbits that may flip */
    int nfree = O - first;
    uint64_t total = 1ULL << nfree;
    for (uint64_t k = 0; k < total; k++) {
        if (k) {
            int b = __builtin_ctzll(k) + first;            /* Gray code: flip orbit b */
            int yo = y[b];
            for (int c = 0; c < C; c++) {
                v[c] += -2L * yo * g[c * O + b];
                int *q = Q + (c * O + b) * O; long *gc = g + c * O;
                for (int j = 0; j < O; j++) gc[j] += -2L * yo * q[j];
            }
            y[b] = -yo; rs += -2L * yo * sz[b]; mask ^= 1ULL << b;
        }
        if (rs == r) {
            for (int c = 0; c < C; c++) row[c] = (int16_t)v[c];
            fwrite(row, 2, C, fv); fwrite(&mask, 8, 1, fm);
            if (++nout == cap) break;
        }
    }
    fclose(fv); fclose(fm);
    printf("%ld\n", nout);
    return 0;
}
