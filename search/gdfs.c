/* gdfs.c: generic depth-first list builder for orbit-structured sequences (search/gsgroup.py writes the tables).
 *
 * usage: gdfs tables.bin out.lst
 *
 * tables.bin (int32 unless noted):
 *   mode (0: +-1 values with a row-sum target, 1: ternary values with a weight bound), nv, npr, target, scale,
 *   vsize[nv], Dg[nv][npr], Cx[nv][nv][npr] (only i < j used), Lv[nv][npr]
 * A variable j with value v contributes v^2 Dg[j][k] + v sum_{i<j} val_i Cx[j][i][k] to the PAF at representative k,
 * v Lv[j][k] to the cross vector X at k, and |v| vsize[j] to the weight (mode 1) or v vsize[j] to the row sum (mode 0;
 * vsize is then the signed row-sum contribution, 0 for a skew pair).
 *
 * out.lst: int32 n, int32 D, int16 vec[n][D], uint64 masks[n][M]
 *   mode 0: D = 2 npr, vec = (PAF, X), M = 1 (bit j set: variable j = -1)
 *   mode 1: D = npr + 1, vec = scale * (PAF, weight), M = 2 (plus mask, minus mask)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXV 64
#define MAXR 128
static int mode, nv, npr, target, scale;
static int vsize[MAXV], *Dg, *Cx, *Lv;
static int val[MAXV], suffix[MAXV + 1];
static int16_t *vbuf; static uint64_t *mbuf; static long cnt, cap;
static int D, M;

#define DG(j, k) Dg[(j) * npr + (k)]
#define CX(j, i, k) Cx[((size_t)(j) * nv + (i)) * npr + (k)]
#define LV(j, k) Lv[(j) * npr + (k)]

static void emit(const int *P, const int *X, int acc) {
    if (cnt == cap) {
        cap = cap ? 2 * cap : 1 << 16;
        vbuf = realloc(vbuf, sizeof(int16_t) * D * cap);
        mbuf = realloc(mbuf, sizeof(uint64_t) * M * cap);
        if (!vbuf || !mbuf) { fprintf(stderr, "out of memory at %ld entries\n", cnt); exit(1); }
    }
    int16_t *v = vbuf + (size_t)cnt * D;
    if (mode == 0) {
        for (int k = 0; k < npr; k++) { v[k] = (int16_t)P[k]; v[npr + k] = (int16_t)X[k]; }
        uint64_t m = 0;
        for (int j = 0; j < nv; j++) if (val[j] < 0) m |= 1ULL << j;
        mbuf[cnt] = m;
    } else {
        for (int k = 0; k < npr; k++) v[k] = (int16_t)(scale * P[k]);
        v[npr] = (int16_t)(scale * acc);
        uint64_t pl = 0, mi = 0;
        for (int j = 0; j < nv; j++) { if (val[j] > 0) pl |= 1ULL << j; else if (val[j] < 0) mi |= 1ULL << j; }
        mbuf[2 * cnt] = pl; mbuf[2 * cnt + 1] = mi;
    }
    cnt++;
}

static void dfs(int j, int acc, const int *P, const int *X) {
    if (mode == 0) {
        int rem = target - acc;
        if (rem > suffix[j] || rem < -suffix[j]) return;
        if (j == nv) { emit(P, X, acc); return; }
    } else if (j == nv) { emit(P, X, acc); return; }
    int P2[MAXR], X2[MAXR];
    if (mode == 1) {
        val[j] = 0;
        dfs(j + 1, acc, P, X);
        if (acc + vsize[j] > target) return;
    }
    for (int s = 0; s < 2; s++) {
        int v = s ? -1 : 1;
        for (int k = 0; k < npr; k++) {
            int d = DG(j, k);
            for (int i = 0; i < j; i++) if (val[i]) d += v * val[i] * CX(j, i, k);
            P2[k] = P[k] + d;
            X2[k] = X[k] + v * LV(j, k);
        }
        val[j] = v;
        dfs(j + 1, mode == 0 ? acc + v * vsize[j] : acc + vsize[j], P2, X2);
    }
    val[j] = 0;
}

static int rd(FILE *f) { int32_t x; if (fread(&x, 4, 1, f) != 1) { fprintf(stderr, "short table file\n"); exit(1); } return x; }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: gdfs tables.bin out.lst\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror(argv[1]); return 1; }
    mode = rd(f); nv = rd(f); npr = rd(f); target = rd(f); scale = rd(f);
    if (nv > MAXV || npr > MAXR) { fprintf(stderr, "nv %d or npr %d too large\n", nv, npr); return 1; }
    for (int j = 0; j < nv; j++) vsize[j] = rd(f);
    Dg = malloc(sizeof(int) * nv * npr); Cx = malloc(sizeof(int) * (size_t)nv * nv * npr); Lv = malloc(sizeof(int) * nv * npr);
    for (int i = 0; i < nv * npr; i++) Dg[i] = rd(f);
    for (size_t i = 0; i < (size_t)nv * nv * npr; i++) Cx[i] = rd(f);
    for (int i = 0; i < nv * npr; i++) Lv[i] = rd(f);
    fclose(f);
    suffix[nv] = 0;
    for (int j = nv - 1; j >= 0; j--) suffix[j] = suffix[j + 1] + abs(vsize[j]);
    D = mode == 0 ? 2 * npr : npr + 1;
    M = mode == 0 ? 1 : 2;
    int P0[MAXR] = {0}, X0[MAXR] = {0};
    dfs(0, 0, P0, X0);
    FILE *o = fopen(argv[2], "wb");
    int32_t hdr[2] = {(int32_t)cnt, D};
    fwrite(hdr, 4, 2, o);
    fwrite(vbuf, sizeof(int16_t), (size_t)cnt * D, o);
    fwrite(mbuf, sizeof(uint64_t), (size_t)cnt * M, o);
    fclose(o);
    fprintf(stderr, "gdfs: mode %d, %d variables, %ld entries\n", mode, nv, cnt);
    return 0;
}
