/*
 * Two-list join for twisted Szöllősi quadruples (attack/tsz.py): sequences A (orbits of L, row sum rA) and C (orbits
 * of H, row sum rC) with class values a + c = 0 on all classes.
 *
 * Both lists are produced by the same depth-first search over the orbits (row sum fixed, exact reachability pruning),
 * with incremental class values v_m(y) = sum_{i,j} T_m[i][j] y_i y_j: flipping y_o changes v_m by
 * -2 y_o R[o][m] + 4 T_m[o][o], and R[p][m] = sum_j (T_m[p][j] + T_m[j][p]) y_j by -2 y_o (T_m[p][o] + T_m[o][p]).
 * A is stored in a hash table keyed by a random linear form of its value vector; every C probes the key of -c (the
 * key is linear, so key(-c) = -key(c)) and candidate matches are checked on the full vector.
 *
 *     tsz INPUT [THREADS [FRAC]]       prints  SOL codeA codeC  and  STATS listA probes matches solutions seconds
 * FRAC: fraction of the depth-first prefixes of C enumerated (a pseudo-random subset).
 * INPUT: int32 OL, OC, CM, rA, rC, T, sizesL[OL], sizesC[OC], TL[CM][OL][OL], TC[CM][OC][OC]; a + c = T on every class.
 * With OL = 0 the A side is an explicit list instead: int64 n, uint64 codes[n], int32 values[n][CM] (after TC).
 * Codes: bit o is 1 when y_o = +1.
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int O, C, CP;
    int *SZ;                 /* orbit sizes */
    int *S;                  /* S[(c * O + i) * O + j] */
    int32_t *TM;             /* TM[(o * O + p) * CP + c] = 2 (S_c[p][o] + S_c[o][p]) */
    int32_t *D4;             /* D4[o * CP + c] = 4 S_c[o][o] */
    int64_t suffix[65];
    unsigned char *REACH; int64_t RW;
} prob;

static int CM, TGT;
static uint64_t RK[256], KTGT;

static void prep(prob *P) {
    int O = P->O, C = P->C;
    P->CP = (C + 15) / 16 * 16;
    int CP = P->CP;
    P->TM = aligned_alloc(64, sizeof(int32_t) * (((int64_t)O * O * CP + 15) / 16) * 16);
    P->D4 = aligned_alloc(64, sizeof(int32_t) * (((int64_t)O * CP + 15) / 16) * 16);
    for (int o = 0; o < O; o++) for (int p = 0; p < O; p++) for (int c = 0; c < CP; c++)
        P->TM[((int64_t)o * O + p) * CP + c] = c < C ? 2 * (P->S[((int64_t)c * O + p) * O + o] + P->S[((int64_t)c * O + o) * O + p]) : 0;
    for (int o = 0; o < O; o++) for (int c = 0; c < CP; c++) P->D4[o * CP + c] = c < C ? 4 * P->S[((int64_t)c * O + o) * O + o] : 0;
    int64_t tot = 0;
    for (int o = 0; o < O; o++) tot += 2 * P->SZ[o];
    P->RW = tot + 1;
    P->REACH = calloc((size_t)(O + 1) * P->RW, 1);
    P->REACH[(int64_t)O * P->RW] = 1;
    for (int pos = O - 1; pos >= 0; pos--)
        for (int64_t n = 0; n < P->RW; n++)
            P->REACH[(int64_t)pos * P->RW + n] = P->REACH[(int64_t)(pos + 1) * P->RW + n] |
                (n >= 2 * P->SZ[pos] ? P->REACH[(int64_t)(pos + 1) * P->RW + n - 2 * P->SZ[pos]] : 0);
}

typedef struct { int y[64]; int32_t *v; int32_t *R; uint64_t code; int64_t sum; } state;

static void st_init(const prob *P, state *s) {
    int O = P->O, C = P->C, CP = P->CP;
    s->R = aligned_alloc(64, sizeof(int32_t) * (((int64_t)O * CP + 15) / 16) * 16);
    s->v = aligned_alloc(64, sizeof(int32_t) * CP);
    for (int o = 0; o < O; o++) s->y[o] = -1;
    s->code = 0; s->sum = 0;
    for (int o = 0; o < O; o++) s->sum -= P->SZ[o];
    for (int c = 0; c < CP; c++) {
        int64_t v = 0;
        if (c < C) for (int i = 0; i < O; i++) for (int j = 0; j < O; j++) v += P->S[((int64_t)c * O + i) * O + j];
        s->v[c] = (int32_t)v;
    }
    for (int o = 0; o < O; o++) for (int c = 0; c < CP; c++) {
        int64_t r = 0;
        if (c < C) for (int j = 0; j < O; j++) r -= P->S[((int64_t)c * O + o) * O + j] + P->S[((int64_t)c * O + j) * O + o];
        s->R[o * CP + c] = (int32_t)r;
    }
}

static inline void flip(const prob *P, state *s, int o) {
    int O = P->O, CP = P->CP;
    int32_t yo = s->y[o];
    int32_t *restrict v = s->v; const int32_t *restrict Ro = s->R + (int64_t)o * CP; const int32_t *restrict d4 = P->D4 + (int64_t)o * CP;
    for (int c = 0; c < CP; c++) v[c] += -2 * yo * Ro[c] + d4[c];
    int32_t *restrict R = s->R + (int64_t)o * CP; const int32_t *restrict T = P->TM + ((int64_t)o * O + o) * CP;
    int n = (O - o) * CP;
    if (yo < 0) { for (int i = 0; i < n; i++) R[i] += T[i]; }
    else { for (int i = 0; i < n; i++) R[i] -= T[i]; }
    s->y[o] = -yo; s->code ^= (1ULL << o); s->sum += 2 * (int64_t)(-yo) * P->SZ[o];
}

typedef void (*leaf_fn)(const prob *, state *, void *);

static void dfs(const prob *P, state *s, int pos, int64_t target, leaf_fn fn, void *ctx) {
    int64_t need = target - s->sum;
    if (need == 0) { fn(P, s, ctx); return; }
    if (pos >= P->O || need < 0 || need >= P->RW || !P->REACH[(int64_t)pos * P->RW + need]) return;
    flip(P, s, pos);
    dfs(P, s, pos + 1, target, fn, ctx);
    flip(P, s, pos);
    dfs(P, s, pos + 1, target, fn, ctx);
}

static inline uint64_t vkey(const int32_t *v) {
    uint64_t h = 0;
    for (int m = 0; m < CM; m++) h += RK[m] * (uint64_t)(int64_t)v[m];
    return h;
}

/* ------------------------------------------------------------------ list A */
typedef struct { uint64_t key, code; } ent;
typedef struct { ent *buf; int64_t n, cap; } buf_t;

static void leafA(const prob *P, state *s, void *vctx) {
    buf_t *b = (buf_t *)vctx;
    if (b->n == b->cap) { b->cap = b->cap ? 2 * b->cap : 1 << 16; b->buf = realloc(b->buf, sizeof(ent) * b->cap); }
    b->buf[b->n].key = vkey(s->v); b->buf[b->n].code = s->code; b->n++;
    (void)P;
}

static int cmp_ent(const void *p, const void *q) {
    uint64_t a = ((const ent *)p)->key, b = ((const ent *)q)->key;
    return a < b ? -1 : a > b;
}

static ent *LA; static int64_t NLA;
static int32_t *XV;          /* explicit A side: values, indexed by the code's position (code field = position) */
static uint64_t *XC;
static uint64_t *HK; static int64_t *HS; static uint64_t HMASK;    /* open addressing: key + 1 -> start in LA */

static void build_hash(void) {
    uint64_t size = 1;
    while (size < 2ULL * (uint64_t)NLA + 16) size <<= 1;
    HK = calloc(size, sizeof(uint64_t)); HS = malloc(sizeof(int64_t) * size); HMASK = size - 1;
    for (int64_t i = 0; i < NLA;) {
        int64_t j = i;
        while (j < NLA && LA[j].key == LA[i].key) j++;
        uint64_t k = LA[i].key + 1, pos = (k * 0x9E3779B97F4A7C15ULL) >> 17 & HMASK;
        while (HK[pos]) pos = (pos + 1) & HMASK;
        HK[pos] = k; HS[pos] = i;
        i = j;
    }
}

static inline int64_t lookup(uint64_t key) {
    uint64_t k = key + 1, pos = (k * 0x9E3779B97F4A7C15ULL) >> 17 & HMASK;
    while (HK[pos]) { if (HK[pos] == k) return HS[pos]; pos = (pos + 1) & HMASK; }
    return -1;
}

static prob PA, PC;

static void values_of(const prob *P, uint64_t code, int64_t *v) {
    int y[64];
    for (int o = 0; o < P->O; o++) y[o] = (code >> o) & 1 ? 1 : -1;
    for (int c = 0; c < P->C; c++) {
        int64_t s = 0;
        for (int i = 0; i < P->O; i++) for (int j = 0; j < P->O; j++) s += (int64_t)P->S[((int64_t)c * P->O + i) * P->O + j] * y[i] * y[j];
        v[c] = s;
    }
}

typedef struct { int64_t probes, matches, sols; } cctx;

static void leafC(const prob *P, state *s, void *vctx) {
    cctx *x = (cctx *)vctx;
    x->probes++;
    uint64_t key = KTGT - vkey(s->v);
    int64_t i = lookup(key);
    if (i < 0) return;
    for (; i < NLA && LA[i].key == key; i++) {
        x->matches++;
        int64_t va[256];
        if (XV) { for (int m = 0; m < CM; m++) va[m] = XV[(int64_t)LA[i].code * CM + m]; }
        else values_of(&PA, LA[i].code, va);
        int ok = 1;
        for (int m = 0; m < CM && ok; m++) ok = (va[m] + s->v[m] == TGT);
        if (ok) {
            x->sols++;
#pragma omp critical
            { printf("SOL %llu %llu\n", (unsigned long long)(XV ? XC[LA[i].code] : LA[i].code), (unsigned long long)s->code); fflush(stdout); }
        }
    }
    (void)P;
}

static double FRAC = 1.0;
static inline int keep_prefix(int64_t pfx) {
    if (FRAC >= 1.0) return 1;
    uint64_t h = (uint64_t)pfx * 0x9E3779B97F4A7C15ULL + 0xD1B54A32D192ED03ULL;
    h ^= h >> 31; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 29;
    return (double)(h >> 11) / 9007199254740992.0 < FRAC;
}

/* enumerate all sign vectors of P with row sum target, split over 2^d prefixes */
static void run(const prob *P, int64_t target, int threads, leaf_fn fn, void **ctxs, int sample) {
    int d = P->O < 14 ? P->O / 2 : 12;
    int64_t np = 1LL << d;
    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        state s; st_init(P, &s);
        #pragma omp for schedule(dynamic, 1)
        for (int64_t pfx = 0; pfx < np; pfx++) {
            if (sample && !keep_prefix(pfx)) continue;
            for (int o = 0; o < d; o++) if ((pfx >> o) & 1) flip(P, &s, o);
            int64_t need = target - s.sum;
            if (need >= 0 && need < P->RW && P->REACH[(int64_t)d * P->RW + need]) dfs(P, &s, d, target, fn, ctxs[tid]);
            for (int o = d - 1; o >= 0; o--) if ((pfx >> o) & 1) flip(P, &s, o);
        }
        free(s.R); free(s.v);
    }
}

static int rd(FILE *f, void *p, size_t sz, size_t n) { return fread(p, sz, n, f) == n; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: tsz INPUT [THREADS [FRAC]]\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror("input"); return 2; }
    int hdr[6];
    if (!rd(f, hdr, sizeof(int), 6)) return 2;
    PA.O = hdr[0]; PC.O = hdr[1]; CM = hdr[2]; int rA = hdr[3], rC = hdr[4]; TGT = hdr[5];
    PA.C = PC.C = CM;
    if (PA.O > 63 || PC.O > 63 || CM > 255) { fprintf(stderr, "limits\n"); return 2; }
    PA.SZ = malloc(sizeof(int) * (PA.O + 1)); PC.SZ = malloc(sizeof(int) * PC.O);
    PA.S = malloc(sizeof(int) * ((size_t)CM * PA.O * PA.O + 1)); PC.S = malloc(sizeof(int) * (size_t)CM * PC.O * PC.O);
    if (!rd(f, PA.SZ, sizeof(int), PA.O) || !rd(f, PC.SZ, sizeof(int), PC.O) ||
        !rd(f, PA.S, sizeof(int), (size_t)CM * PA.O * PA.O) || !rd(f, PC.S, sizeof(int), (size_t)CM * PC.O * PC.O)) {
        fprintf(stderr, "short input\n"); return 2;
    }
    int64_t nx = 0;
    if (PA.O == 0) {
        if (!rd(f, &nx, sizeof(int64_t), 1)) return 2;
        XC = malloc(sizeof(uint64_t) * (nx + 1)); XV = malloc(sizeof(int32_t) * ((size_t)nx * CM + 1));
        if (!rd(f, XC, sizeof(uint64_t), nx) || !rd(f, XV, sizeof(int32_t), (size_t)nx * CM)) { fprintf(stderr, "short list\n"); return 2; }
    }
    fclose(f);
    int threads = argc > 2 ? atoi(argv[2]) : omp_get_max_threads();
    if (argc > 3) FRAC = atof(argv[3]);
    if (PA.O) prep(&PA);
    prep(&PC);
    uint64_t x = 0x243F6A8885A308D3ULL;
    for (int m = 0; m < CM; m++) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; RK[m] = x | 1; }
    KTGT = 0;
    for (int m = 0; m < CM; m++) KTGT += RK[m] * (uint64_t)(int64_t)TGT;
    double t0 = omp_get_wtime();
    /* list A (single pass, merged per thread) */
    buf_t *bufs = calloc(threads, sizeof(buf_t));
    void **ctx = malloc(sizeof(void *) * threads);
    for (int t = 0; t < threads; t++) ctx[t] = &bufs[t];
    if (PA.O) {
        int64_t tA = rA;
        run(&PA, tA, threads, leafA, ctx, 0);
        NLA = 0;
        for (int t = 0; t < threads; t++) NLA += bufs[t].n;
        LA = malloc(sizeof(ent) * (NLA + 1));
        int64_t off = 0;
        for (int t = 0; t < threads; t++) { memcpy(LA + off, bufs[t].buf, sizeof(ent) * bufs[t].n); off += bufs[t].n; free(bufs[t].buf); }
    } else {
        NLA = nx; LA = malloc(sizeof(ent) * (NLA + 1));
        for (int64_t i = 0; i < nx; i++) { LA[i].key = vkey(XV + i * CM); LA[i].code = (uint64_t)i; }
    }
    qsort(LA, NLA, sizeof(ent), cmp_ent);
    build_hash();
    fprintf(stderr, "list A: %lld (%.0fs)\n", (long long)NLA, omp_get_wtime() - t0);
    /* stream C */
    cctx *cc = calloc(threads, sizeof(cctx));
    for (int t = 0; t < threads; t++) ctx[t] = &cc[t];
    run(&PC, (int64_t)rC, threads, leafC, ctx, 1);
    int64_t probes = 0, matches = 0, sols = 0;
    for (int t = 0; t < threads; t++) { probes += cc[t].probes; matches += cc[t].matches; sols += cc[t].sols; }
    printf("STATS listA %lld probes %lld matches %lld solutions %lld seconds %.0f\n", (long long)NLA, (long long)probes,
           (long long)matches, (long long)sols, omp_get_wtime() - t0);
    return 0;
}
