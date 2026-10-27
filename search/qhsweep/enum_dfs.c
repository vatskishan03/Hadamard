/*
 * Pruned depth-first enumeration of every +-1 vector y on O orbits with a given row sum, writing its class values
 *     v_c(y) = sum_{i,j} S_c[i][j] y_i y_j       (c = 0 .. C-1)
 * as int16 rows plus the sign mask (bit i set = -1). Same input and output format as search/bnw_enum.c, but the
 * work is proportional to the number of vectors with that row sum, not to 2^O: the orbits are taken by decreasing
 * size and a branch is cut as soon as the row sum can no longer be reached (subset-sum reachability table).
 * The class values are updated incrementally in O(C * O) per decision.
 *
 * fix0 = 1 fixes y_0 = +1 (the global sign when the row sum is 0). With nbits > 0 the first nbits DFS decisions are
 * fixed to the bits of j (bit d set = -1 for the d-th orbit in the order above): the 2^nbits slices j = 0 ..
 * 2^nbits - 1 partition the list, so one list too big for memory can be enumerated and joined slice by slice.
 *
 * input  (binary): int32 O, C, rowsum, fix0; int32 sizes[O]; int32 S[C][O][O]
 *   cc -O3 -march=native enum_dfs.c -o enum_dfs ; ./enum_dfs in.bin vals.bin masks.bin [cap [nbits j]]
 * output: vals.bin (int16 rows of C values), masks.bin (uint64 sign masks, the same order); prints the number of rows.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int O, C, r, fix0;
static int32_t *sz, *Q, *D;           /* Q[c][i][j], D[c][i] */
static int perm[64];
static unsigned char *reach[65];      /* reach[d][s + off]: sum s reachable from perm[d..O-1] */
static int off;
static FILE *fv, *fm;
static long nout = 0, cap = -1;
static int nbits = 0; static long pidx = 0;       /* slice j: the first nbits DFS decisions are the bits of j */
static int32_t *gst;                  /* gst[d][c*O + p] */
static int32_t vst[65][1024];
static int16_t *obuf; static uint64_t *mbuf; static long nbuf = 0;
#define BUF 65536

static void flush_buf(void) {
    fwrite(obuf, 2 * (size_t)C, nbuf, fv); fwrite(mbuf, 8, nbuf, fm); nbuf = 0;
}

static void dfs(int d, long rs, uint64_t mask) {
    if (cap >= 0 && nout >= cap) return;
    if (d == O) {
        int16_t *row = obuf + (size_t)nbuf * C;
        for (int c = 0; c < C; c++) row[c] = (int16_t)vst[d][c];
        mbuf[nbuf++] = mask; nout++;
        if (nbuf == BUF) flush_buf();
        return;
    }
    int p = perm[d];
    int32_t *g = gst + (size_t)d * C * O;
    for (int s = 1; s >= -1; s -= 2) {
        if (fix0 && p == 0 && s < 0) continue;
        if (d < nbits && ((s < 0) != (((pidx >> d) & 1) != 0))) continue;
        long ns = rs + (long)s * sz[p];
        long need = r - ns;
        if (need < -off || need > off || !reach[d + 1][need + off]) continue;
        for (int c = 0; c < C; c++) vst[d + 1][c] = vst[d][c] + D[c * O + p] + s * g[c * O + p];
        int32_t *g2 = gst + (size_t)(d + 1) * C * O;
        for (int c = 0; c < C; c++) {
            const int32_t *q = Q + ((size_t)c * O + p) * O; const int32_t *gc = g + (size_t)c * O; int32_t *g2c = g2 + (size_t)c * O;
            for (int e = d + 1; e < O; e++) { int j = perm[e]; g2c[j] = gc[j] + s * q[j]; }
        }
        dfs(d + 1, ns, s < 0 ? (mask | (1ULL << p)) : mask);
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: enum_dfs in.bin vals.bin masks.bin [cap [nbits j]]\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    int32_t hdr[4];
    if (!f || fread(hdr, 4, 4, f) != 4) return 2;
    O = hdr[0]; C = hdr[1]; r = hdr[2]; fix0 = hdr[3];
    if (O > 62 || C > 1000) { fprintf(stderr, "too large\n"); return 2; }
    sz = malloc(4 * O); int32_t *S = malloc(4L * C * O * O);
    if (fread(sz, 4, O, f) != (size_t)O || fread(S, 4, (size_t)C * O * O, f) != (size_t)C * O * O) return 2;
    fclose(f);
    cap = argc > 4 ? atol(argv[4]) : -1;
    if (argc > 6) { nbits = atoi(argv[5]); pidx = atol(argv[6]); }
    Q = malloc(4L * C * O * O); D = malloc(4L * C * O);
    for (int c = 0; c < C; c++) for (int i = 0; i < O; i++) {
        D[c * O + i] = S[((size_t)c * O + i) * O + i];
        for (int j = 0; j < O; j++) Q[((size_t)c * O + i) * O + j] = i == j ? 0 : S[((size_t)c * O + i) * O + j] + S[((size_t)c * O + j) * O + i];
    }
    /* order: decreasing size (stable) */
    for (int i = 0; i < O; i++) perm[i] = i;
    for (int i = 1; i < O; i++) { int x = perm[i], j = i - 1; while (j >= 0 && sz[perm[j]] < sz[x]) { perm[j + 1] = perm[j]; j--; } perm[j + 1] = x; }
    long tot = 0; for (int i = 0; i < O; i++) tot += sz[i];
    off = (int)tot;
    for (int d = 0; d <= O; d++) { reach[d] = calloc(2 * off + 1, 1); }
    reach[O][off] = 1;
    for (int d = O - 1; d >= 0; d--) {
        int s = sz[perm[d]];
        for (int x = 0; x <= 2 * off; x++) if (reach[d + 1][x]) {
            if (x + s <= 2 * off) reach[d][x + s] = 1;
            if (x - s >= 0) reach[d][x - s] = 1;
        }
    }
    gst = calloc((size_t)(O + 1) * C * O, 4);
    /* initial: v = 0, g = 0 */
    memset(vst[0], 0, sizeof(vst[0]));
    obuf = malloc(2 * (size_t)C * BUF); mbuf = malloc(8 * (size_t)BUF);
    fv = fopen(argv[2], "wb"); fm = fopen(argv[3], "wb");
    if (!fv || !fm) return 3;
    if (r + off >= 0 && r + off <= 2 * off && reach[0][r + off]) dfs(0, 0, 0);
    flush_buf();
    fclose(fv); fclose(fm);
    printf("%ld\n", nout);
    return 0;
}
