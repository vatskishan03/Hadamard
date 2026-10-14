/*
 * Low-energy enumeration of one block (attack/lowe_sz.py): every sign vector y on O orbits with row sum r is visited
 * by the depth-first search of kstream.c / tsz.c (incremental class values v_c(y) = sum_{i,j} S_c[i][j] y_i y_j,
 * exact reachability pruning), and kept when its energy
 *
 *     |W (v - m)|^2 <= tau          (W a whitening matrix, m the conditional mean of the block's key)
 *
 * is at most tau. Kept members are written as uint64 masks (bit o = 1 when y_o = +1) and int16 class values.
 *
 *     lowe_enum INPUT OUTPUT [THREADS [FRAC]]
 * INPUT: int32 O, C, r, maxkeep; int32 sizes[O]; int32 S[C][O][O]; float m[C]; float W[C][C]; float tau.
 * OUTPUT: int64 n; uint64 masks[n]; int16 values[n][C]; float energy[n]. If more than maxkeep members pass, the
 * maxkeep of lowest energy are written. Prints "STATS leaves L kept K tau T seconds S".
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int O, C, CP;
static int *SZ, *S;
static int32_t *TM, *D4;
static unsigned char *REACH; static int64_t RW;
static float *MC, *WM, TAU;

#define SI(c, i, j) S[((int64_t)(c) * O + (i)) * O + (j)]

typedef struct { int y[64]; int32_t *v; int32_t *R; uint64_t code; int64_t sum; } state;

static void st_init(state *s) {
    s->R = aligned_alloc(64, sizeof(int32_t) * (((int64_t)O * CP + 15) / 16) * 16);
    s->v = aligned_alloc(64, sizeof(int32_t) * CP);
    for (int o = 0; o < O; o++) s->y[o] = -1;
    s->code = 0; s->sum = 0;
    for (int o = 0; o < O; o++) s->sum -= SZ[o];
    for (int c = 0; c < CP; c++) {
        int64_t v = 0;
        if (c < C) for (int i = 0; i < O; i++) for (int j = 0; j < O; j++) v += SI(c, i, j);
        s->v[c] = (int32_t)v;
    }
    for (int o = 0; o < O; o++) for (int c = 0; c < CP; c++) {
        int64_t r = 0;
        if (c < C) for (int j = 0; j < O; j++) r -= SI(c, o, j) + SI(c, j, o);
        s->R[o * CP + c] = (int32_t)r;
    }
}

static inline void flip(state *s, int o) {
    int32_t yo = s->y[o];
    int32_t *restrict v = s->v; const int32_t *restrict Ro = s->R + (int64_t)o * CP; const int32_t *restrict d4 = D4 + (int64_t)o * CP;
    for (int c = 0; c < CP; c++) v[c] += -2 * yo * Ro[c] + d4[c];
    int32_t *restrict R = s->R + (int64_t)o * CP; const int32_t *restrict T = TM + ((int64_t)o * O + o) * CP;
    int n = (O - o) * CP;
    if (yo < 0) { for (int i = 0; i < n; i++) R[i] += T[i]; }
    else { for (int i = 0; i < n; i++) R[i] -= T[i]; }
    s->y[o] = -yo; s->code ^= (1ULL << o); s->sum += 2 * (int64_t)(-yo) * SZ[o];
}

typedef struct { uint64_t *mask; int16_t *val; float *en; int64_t n, cap; int64_t leaves; } buf_t;

static inline float energy(const int32_t *v) {
    float d[256], e = 0;
    for (int c = 0; c < C; c++) d[c] = (float)v[c] - MC[c];
    for (int i = 0; i < C; i++) {
        float u = 0;
        const float *w = WM + (int64_t)i * C;
        for (int c = 0; c < C; c++) u += w[c] * d[c];
        e += u * u;
        if (e > TAU) return e;
    }
    return e;
}

static void leaf(state *s, buf_t *b) {
    b->leaves++;
    float e = energy(s->v);
    if (e > TAU) return;
    if (b->n == b->cap) {
        b->cap = b->cap ? 2 * b->cap : 1 << 16;
        b->mask = realloc(b->mask, sizeof(uint64_t) * b->cap);
        b->val = realloc(b->val, sizeof(int16_t) * b->cap * C);
        b->en = realloc(b->en, sizeof(float) * b->cap);
    }
    b->mask[b->n] = s->code; b->en[b->n] = e;
    for (int c = 0; c < C; c++) b->val[b->n * C + c] = (int16_t)s->v[c];
    b->n++;
}

static void dfs(state *s, int pos, int64_t target, buf_t *b) {
    int64_t need = target - s->sum;
    if (need == 0) { leaf(s, b); return; }
    if (pos >= O || need < 0 || need >= RW || !REACH[(int64_t)pos * RW + need]) return;
    flip(s, pos);
    dfs(s, pos + 1, target, b);
    flip(s, pos);
    dfs(s, pos + 1, target, b);
}

static double FRAC = 1.0;
static inline int keep_prefix(int64_t pfx) {
    if (FRAC >= 1.0) return 1;
    uint64_t h = (uint64_t)pfx * 0x9E3779B97F4A7C15ULL + 0x5851F42D4C957F2DULL;
    h ^= h >> 31; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 29;
    return (double)(h >> 11) / 9007199254740992.0 < FRAC;
}

static int cmpf_base; static float *cmpf_e;
static int cmpf(const void *a, const void *b) {
    float x = cmpf_e[*(const int64_t *)a], y = cmpf_e[*(const int64_t *)b];
    return x < y ? -1 : x > y;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: lowe_enum INPUT OUTPUT [THREADS [FRAC]]\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror("input"); return 2; }
    int hdr[4];
    if (fread(hdr, sizeof(int), 4, f) != 4) return 2;
    O = hdr[0]; C = hdr[1]; int r = hdr[2]; int64_t maxkeep = hdr[3];
    if (O > 63 || C > 255) { fprintf(stderr, "limits\n"); return 2; }
    SZ = malloc(sizeof(int) * O); S = malloc(sizeof(int) * (size_t)C * O * O);
    MC = malloc(sizeof(float) * C); WM = malloc(sizeof(float) * C * C);
    if (fread(SZ, sizeof(int), O, f) != (size_t)O || fread(S, sizeof(int), (size_t)C * O * O, f) != (size_t)C * O * O ||
        fread(MC, sizeof(float), C, f) != (size_t)C || fread(WM, sizeof(float), (size_t)C * C, f) != (size_t)C * C ||
        fread(&TAU, sizeof(float), 1, f) != 1) { fprintf(stderr, "short input\n"); return 2; }
    fclose(f);
    int threads = argc > 3 ? atoi(argv[3]) : omp_get_max_threads();
    if (argc > 4) FRAC = atof(argv[4]);
    CP = (C + 15) / 16 * 16;
    TM = aligned_alloc(64, sizeof(int32_t) * (((int64_t)O * O * CP + 15) / 16) * 16);
    D4 = aligned_alloc(64, sizeof(int32_t) * (((int64_t)O * CP + 15) / 16) * 16);
    for (int o = 0; o < O; o++) for (int p = 0; p < O; p++) for (int c = 0; c < CP; c++)
        TM[((int64_t)o * O + p) * CP + c] = c < C ? 2 * (SI(c, p, o) + SI(c, o, p)) : 0;
    for (int o = 0; o < O; o++) for (int c = 0; c < CP; c++) D4[o * CP + c] = c < C ? 4 * SI(c, o, o) : 0;
    int64_t tot = 0;
    for (int o = 0; o < O; o++) tot += 2 * SZ[o];
    RW = tot + 1;
    REACH = calloc((size_t)(O + 1) * RW, 1);
    REACH[(int64_t)O * RW] = 1;
    for (int pos = O - 1; pos >= 0; pos--)
        for (int64_t n = 0; n < RW; n++)
            REACH[(int64_t)pos * RW + n] = REACH[(int64_t)(pos + 1) * RW + n] | (n >= 2 * SZ[pos] ? REACH[(int64_t)(pos + 1) * RW + n - 2 * SZ[pos]] : 0);
    double t0 = omp_get_wtime();
    buf_t *bufs = calloc(threads, sizeof(buf_t));
    int d = O < 14 ? O / 2 : 12;
    int64_t np = 1LL << d;
    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        state s; st_init(&s);
        #pragma omp for schedule(dynamic, 1)
        for (int64_t pfx = 0; pfx < np; pfx++) {
            if (!keep_prefix(pfx)) continue;
            for (int o = 0; o < d; o++) if ((pfx >> o) & 1) flip(&s, o);
            int64_t need = (int64_t)r - s.sum;
            if (need >= 0 && need < RW && REACH[(int64_t)d * RW + need]) dfs(&s, d, r, &bufs[tid]);
            for (int o = d - 1; o >= 0; o--) if ((pfx >> o) & 1) flip(&s, o);
        }
        free(s.R); free(s.v);
    }
    int64_t n = 0, leaves = 0;
    for (int t = 0; t < threads; t++) { n += bufs[t].n; leaves += bufs[t].leaves; }
    uint64_t *mask = malloc(sizeof(uint64_t) * (n + 1)); int16_t *val = malloc(sizeof(int16_t) * (n + 1) * C);
    float *en = malloc(sizeof(float) * (n + 1));
    int64_t off = 0;
    for (int t = 0; t < threads; t++) {
        memcpy(mask + off, bufs[t].mask, sizeof(uint64_t) * bufs[t].n);
        memcpy(val + off * C, bufs[t].val, sizeof(int16_t) * bufs[t].n * C);
        memcpy(en + off, bufs[t].en, sizeof(float) * bufs[t].n);
        off += bufs[t].n; free(bufs[t].mask); free(bufs[t].val); free(bufs[t].en);
    }
    int64_t keep = n;
    int64_t *ord = NULL;
    if (n > maxkeep) {
        ord = malloc(sizeof(int64_t) * n);
        for (int64_t i = 0; i < n; i++) ord[i] = i;
        cmpf_e = en; qsort(ord, n, sizeof(int64_t), cmpf);
        keep = maxkeep;
    }
    FILE *g = fopen(argv[2], "wb");
    fwrite(&keep, sizeof(int64_t), 1, g);
    for (int64_t i = 0; i < keep; i++) { int64_t k = ord ? ord[i] : i; fwrite(mask + k, sizeof(uint64_t), 1, g); }
    for (int64_t i = 0; i < keep; i++) { int64_t k = ord ? ord[i] : i; fwrite(val + k * C, sizeof(int16_t), C, g); }
    for (int64_t i = 0; i < keep; i++) { int64_t k = ord ? ord[i] : i; fwrite(en + k, sizeof(float), 1, g); }
    fclose(g);
    printf("STATS leaves %lld passed %lld kept %lld tau %.3f seconds %.0f\n", (long long)leaves, (long long)n,
           (long long)keep, TAU, omp_get_wtime() - t0);
    (void)cmpf_base;
    return 0;
}
