/*
 * Exhaustive low-energy enumeration of the members of one block of a cyclotomic-product structure
 * (attack/cp_lowe.py). A member is a choice of one pattern index per position r = 0 .. R-1 (position r has n_r
 * allowed patterns). Its key is a sum of pair terms: for r1 >= r2 the table Ti[r1][r2][i1][i2] (int32, D values)
 * and the whitened table Tz[r1][r2][i1][i2] (float, Z values) give the contribution of the patterns i1 at r1 and i2
 * at r2. The whitened key is z = z0 + sum of Tz; the member is kept when |z|^2 <= tau. Kept members are written as
 * records (uint64 mask, int16 key[D]), the mask holding the pattern indices in BITS bits per position.
 *
 * input (int32 unless noted): R D Z BITS; n_0 .. n_{R-1}; float tau; float z0[Z]; int32 key0[D];
 *        for r1 = 0..R-1, r2 = 0..r1: int32 Ti[n_r1][n_r2][D], float Tz[n_r1][n_r2][Z];
 *        optionally (when present): int32 TAG; for r = 0..R-1: int32 tag_r[n_r]. Then only members whose tags sum to
 *        TAG are considered (used for a fixed class sum when the positions are groups of classes).
 *   cc -O3 -march=native -fopenmp cp_enum.c -o cp_enum ;  ./cp_enum input.bin output.bin
 * Prints "KEPT k OF n" on stdout.
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXR 24
#define MAXZ 64
static int R, D, Z, BITS, nr[MAXR];
static float tau, z0[MAXZ];
static int32_t key0[MAXZ];
static int32_t *Ti[MAXR][MAXR];
static float *Tz[MAXR][MAXR];
static FILE *out;
static long kept = 0;
static int usetag = 0, TAG = 0, *tag[MAXR], tmin[MAXR + 1], tmax[MAXR + 1];   /* tmin/tmax[j]: range of the tags of positions j.. */
static int *lastby, *laststart, lastlo, lasthi;   /* patterns of the last position grouped by tag */

typedef struct { int idx[MAXR]; float zs[MAXR + 1][MAXZ]; char *buf; size_t nb, cap; int ts[MAXR + 1]; } State;

static void rd(void *p, size_t sz, size_t k, FILE *f) { if (fread(p, sz, k, f) != k) { fprintf(stderr, "cp_enum: short read\n"); exit(1); } }

static void flush(State *s) {
    #pragma omp critical
    { fwrite(s->buf, 8 + 2 * D, s->nb, out); kept += s->nb; }
    s->nb = 0;
}

static void leaf(State *s) {
    float en = 0; for (int k = 0; k < Z; k++) en += s->zs[R][k] * s->zs[R][k];
    if (en > tau) return;
    uint64_t mask = 0; int32_t key[MAXZ];
    for (int k = 0; k < D; k++) key[k] = key0[k];
    for (int r1 = 0; r1 < R; r1++) {
        mask |= (uint64_t)s->idx[r1] << (BITS * r1);
        for (int r2 = 0; r2 <= r1; r2++) {
            const int32_t *dd = Ti[r1][r2] + ((size_t)s->idx[r1] * nr[r2] + s->idx[r2]) * D;
            for (int k = 0; k < D; k++) key[k] += dd[k];
        }
    }
    if (s->nb == s->cap) flush(s);
    char *rec = s->buf + s->nb * (8 + 2 * D); memcpy(rec, &mask, 8);
    for (int k = 0; k < D; k++) { int16_t v = (int16_t)key[k]; memcpy(rec + 8 + 2 * k, &v, 2); }
    s->nb++;
}

static void dfs(State *s, int j) {
    if (j == R) { if (!usetag || s->ts[R] == TAG) leaf(s); return; }
    if (usetag && (s->ts[j] + tmin[j] > TAG || s->ts[j] + tmax[j] < TAG)) return;
    int lo = 0, hi = nr[j]; const int *list = NULL;
    if (usetag && j == R - 1) {                                  /* only the patterns that complete the class sum */
        int need = TAG - s->ts[j];
        if (need < lastlo || need > lasthi) return;
        list = lastby; lo = laststart[need - lastlo]; hi = laststart[need - lastlo + 1];
    }
    for (int ii = lo; ii < hi; ii++) {
        int i = list ? list[ii] : ii;
        s->idx[j] = i;
        if (usetag) s->ts[j + 1] = s->ts[j] + tag[j][i];
        float *zn = s->zs[j + 1]; const float *zo = s->zs[j];
        for (int k = 0; k < Z; k++) zn[k] = zo[k];
        for (int r = 0; r <= j; r++) {
            const float *dd = Tz[j][r] + ((size_t)i * nr[r] + s->idx[r]) * Z;
            for (int k = 0; k < Z; k++) zn[k] += dd[k];
        }
        dfs(s, j + 1);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: cp_enum input.bin output.bin\n"); return 2; }
    FILE *f = fopen(argv[1], "rb"); if (!f) { perror("input"); return 1; }
    int32_t hdr[4]; rd(hdr, 4, 4, f); R = hdr[0]; D = hdr[1]; Z = hdr[2]; BITS = hdr[3];
    if (R < 2 || R > MAXR || D > MAXZ || Z > MAXZ || R * BITS > 64) { fprintf(stderr, "cp_enum: sizes\n"); return 1; }
    rd(nr, 4, R, f); rd(&tau, 4, 1, f); rd(z0, 4, Z, f); rd(key0, 4, D, f);
    for (int r1 = 0; r1 < R; r1++)
        for (int r2 = 0; r2 <= r1; r2++) {
            size_t cnt = (size_t)nr[r1] * nr[r2];
            Ti[r1][r2] = malloc(cnt * D * 4); Tz[r1][r2] = malloc(cnt * Z * 4);
            rd(Ti[r1][r2], 4, cnt * D, f); rd(Tz[r1][r2], 4, cnt * Z, f);
        }
    if (fread(&TAG, 4, 1, f) == 1) {
        usetag = 1;
        for (int r = 0; r < R; r++) { tag[r] = malloc(4 * nr[r]); rd(tag[r], 4, nr[r], f); }
        tmin[R] = tmax[R] = 0;
        {   /* bucket the last position's patterns by tag */
            int r = R - 1; lastlo = 1 << 30; lasthi = -(1 << 30);
            for (int i = 0; i < nr[r]; i++) { if (tag[r][i] < lastlo) lastlo = tag[r][i]; if (tag[r][i] > lasthi) lasthi = tag[r][i]; }
            int nb_ = lasthi - lastlo + 1; laststart = calloc(nb_ + 1, sizeof(int)); lastby = malloc(sizeof(int) * nr[r]);
            for (int i = 0; i < nr[r]; i++) laststart[tag[r][i] - lastlo + 1]++;
            for (int t = 0; t < nb_; t++) laststart[t + 1] += laststart[t];
            int *fill = calloc(nb_, sizeof(int));
            for (int i = 0; i < nr[r]; i++) { int t = tag[r][i] - lastlo; lastby[laststart[t] + fill[t]++] = i; }
            free(fill);
        }
        for (int r = R - 1; r >= 0; r--) {
            int lo = 1 << 30, hi = -(1 << 30);
            for (int i = 0; i < nr[r]; i++) { if (tag[r][i] < lo) lo = tag[r][i]; if (tag[r][i] > hi) hi = tag[r][i]; }
            tmin[r] = tmin[r + 1] + lo; tmax[r] = tmax[r + 1] + hi;
        }
    }
    fclose(f);
    out = fopen(argv[2], "wb"); if (!out) { perror("output"); return 1; }
    long total = 1; for (int r = 0; r < R; r++) total *= nr[r];
    int n01 = nr[0] * nr[1];
    #pragma omp parallel
    {
        State *s = malloc(sizeof(State)); s->cap = 1 << 16; s->nb = 0; s->buf = malloc(s->cap * (8 + 2 * D));
        #pragma omp for schedule(dynamic, 1)
        for (int t = 0; t < n01; t++) {
            s->idx[0] = t / nr[1]; s->idx[1] = t % nr[1];
            s->ts[0] = 0;
            if (usetag) {
                s->ts[1] = tag[0][s->idx[0]]; s->ts[2] = s->ts[1] + tag[1][s->idx[1]];
                if (s->ts[2] + tmin[2] > TAG || s->ts[2] + tmax[2] < TAG) continue;
            }
            for (int k = 0; k < Z; k++) s->zs[0][k] = z0[k];
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < Z; k++) s->zs[j + 1][k] = s->zs[j][k];
                for (int r = 0; r <= j; r++) {
                    const float *dd = Tz[j][r] + ((size_t)s->idx[j] * nr[r] + s->idx[r]) * Z;
                    for (int k = 0; k < Z; k++) s->zs[j + 1][k] += dd[k];
                }
            }
            dfs(s, 2);
        }
        if (s->nb) flush(s);
        free(s->buf); free(s);
    }
    fclose(out);
    printf("KEPT %ld OF %ld\n", kept, total);
    return 0;
}
