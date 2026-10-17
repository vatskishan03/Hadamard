/*
 * Streaming k-tree join for Szöllősi structures whose two arbitrary blocks are too large to list (search/kstream.py).
 *
 * Four blocks, class values a (list A), b (list B), c (block C, streamed), d (block D, streamed), with
 *
 *     a + b + c + d = T   on all classes.
 *
 * C and D are enumerated completely by a depth-first search over their orbits (row sums fixed), with incremental
 * class values: v_c(y) = sum_{i,j} S_c[i][j] y_i y_j and R[o][c] = sum_j (S_c[o][j] + S_c[j][o]) y_j, so flipping y_o
 * changes v_c by -2 y_o R[o][c] + 4 S_c[o][o] and every R[o'][c] by -2 y_o (S_c[o'][o] + S_c[o][o']).
 *
 * Pass 1: for every C and every target t (values of a + c on the K1 classes) the A entries with a|K1 = t - c|K1 are
 * looked up (direct-indexed table) and (fingerprint of (a + c) on the other classes, A index, C code) is stored.
 * Pass 2: for every D and target t the B entries with b|K1 = T|K1 - t - d|K1 are looked up and the fingerprint of
 * T - (b + d) on the other classes is searched in the sorted pass-1 table; every match is verified on all classes.
 *
 *     kstream INPUT MAXM1 [THREADS [FRAC1 [FRAC2]]]   prints  SOL a b ccode dcode  and  STATS ...
 * FRAC1, FRAC2: fractions of the DFS prefixes enumerated in the two passes (a pseudo-random subset).
 *
 * INPUT: int32 O, C, k, nT, rC, rD, nA, nB, then int32 sizes[O], int32 S[C][O][O], int32 K1[k], int32 targets[nT][k],
 * int32 T[C], int16 A[nA][C], int16 B[nB][C]. Codes: bit o of the code is 1 when y_o = +1.
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int O, C, K, NT, RC, RD, NA, NB;
static int *SZ, *S, *K1, *TG, *TT;
static int16_t *AV, *BV;
static int REST[256], NREST;
static int inK1[256];

#define SI(c, i, j) S[((int64_t)(c) * O + (i)) * O + (j)]

typedef struct { uint64_t fp; uint64_t code; } ent;   /* code: A index << 40 | C code (O <= 40) */

static inline uint64_t mix(uint64_t h, uint64_t x) {
    h ^= x + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
    h *= 0xff51afd7ed558ccdULL; h ^= h >> 33;
    return h;
}

typedef struct { uint64_t key; int32_t start, cnt; } hent;
typedef struct { hent *h; uint64_t mask; int32_t *idx; uint64_t *bloom; } htab;
#define BLOOM_BITS 26
static htab TA, TB;

/* linear key: sum_q RK[q] * w_q (mod 2^64), so key(t - v) = key(t) - key(v); the table stores key | 1 is avoided by
   storing key + 1 and never 0 */
static uint64_t RK[16];
static inline uint64_t vkey(const int *w) {
    uint64_t h = 0;
    for (int q = 0; q < K; q++) h += RK[q] * (uint64_t)(int64_t)w[q];
    return h;
}

typedef struct { uint64_t k; int32_t i; } kp;
static int cmpkp(const void *x, const void *y) { uint64_t p = ((const kp *)x)->k, q = ((const kp *)y)->k; return p < q ? -1 : p > q; }

static void build_table(int16_t *V, int n, htab *T) {
    uint64_t *key = malloc(sizeof(uint64_t) * n);
    int32_t *ord = malloc(sizeof(int32_t) * n);
    int w[16];
    for (int i = 0; i < n; i++) {
        for (int q = 0; q < K; q++) w[q] = V[(int64_t)i * C + K1[q]];
        key[i] = vkey(w); ord[i] = i;
    }
    /* sort indices by key */
    kp *a = malloc(sizeof(kp) * n);
    for (int i = 0; i < n; i++) { a[i].k = key[i]; a[i].i = i; }
    qsort(a, n, sizeof(kp), cmpkp);
    T->idx = malloc(sizeof(int32_t) * (n + 1));
    for (int i = 0; i < n; i++) T->idx[i] = a[i].i;
    uint64_t size = 1;
    while (size < 4ULL * (uint64_t)n + 16) size <<= 1;
    T->h = calloc(size, sizeof(hent)); T->mask = size - 1;
    T->bloom = calloc((1ULL << BLOOM_BITS) / 64, sizeof(uint64_t));
    for (int i = 0; i < n; i++) { uint64_t b = ((a[i].k + 1) * 0xD6E8FEB86659FD93ULL) >> (64 - BLOOM_BITS); T->bloom[b >> 6] |= 1ULL << (b & 63); }
    for (int i = 0; i < n;) {
        int j = i;
        while (j < n && a[j].k == a[i].k) j++;
        uint64_t kk = a[i].k + 1, pos = (kk * 0x9E3779B97F4A7C15ULL) >> 20 & T->mask;
        while (T->h[pos].key) pos = (pos + 1) & T->mask;
        T->h[pos].key = kk; T->h[pos].start = i; T->h[pos].cnt = j - i;
        i = j;
    }
    free(a); free(key); free(ord);
}

static inline const hent *lookup_key(const htab *T, uint64_t key) {
    uint64_t k = key + 1, b = (k * 0xD6E8FEB86659FD93ULL) >> (64 - BLOOM_BITS);
    if (!((T->bloom[b >> 6] >> (b & 63)) & 1)) return NULL;
    uint64_t pos = (k * 0x9E3779B97F4A7C15ULL) >> 20 & T->mask;
    while (T->h[pos].key) {
        if (T->h[pos].key == k) return &T->h[pos];
        pos = (pos + 1) & T->mask;
    }
    return NULL;
}
static uint64_t *KT;      /* key(targets) */

/* per-thread DFS state; TM[o][p*CP + c] = 2 (S_c[p][o] + S_c[o][p]) (int16), D4[o][c] = 4 S_c[o][o] (int32);
   R[p][c] = sum_j (S_c[p][j] + S_c[j][p]) y_j is kept up to date only for the orbits p the search may still flip */
static int CP;                       /* C rounded up to a multiple of 16 */
static int16_t *TM; static int32_t *D4;

typedef struct {
    int y[64]; int32_t *v; int16_t *R;    /* v[CP], R[p*CP + c] */
    uint64_t code; int64_t sum;
} state;

static void st_init(state *s) {
    s->R = aligned_alloc(64, sizeof(int16_t) * ((O * CP + 31) / 32) * 32);
    s->v = aligned_alloc(64, sizeof(int32_t) * ((CP + 15) / 16) * 16);
    for (int o = 0; o < O; o++) s->y[o] = -1;
    s->code = 0; s->sum = 0;
    for (int o = 0; o < O; o++) s->sum -= SZ[o];
    for (int c = 0; c < CP; c++) {
        int64_t v = 0;
        if (c < C) for (int i = 0; i < O; i++) for (int j = 0; j < O; j++) v += SI(c, i, j);   /* all y = -1 */
        s->v[c] = (int32_t)v;
    }
    for (int o = 0; o < O; o++) for (int c = 0; c < CP; c++) {
        int64_t r = 0;
        if (c < C) for (int j = 0; j < O; j++) r -= SI(c, o, j) + SI(c, j, o);
        s->R[o * CP + c] = (int16_t)r;
    }
}

static inline void flip(state *s, int o) {
    int32_t yo = s->y[o];
    int32_t *restrict v = s->v; const int16_t *Ro = s->R + (int64_t)o * CP; const int32_t *restrict d4 = D4 + (int64_t)o * CP;
    for (int c = 0; c < CP; c++) v[c] += -2 * yo * (int32_t)Ro[c] + d4[c];
    int16_t *restrict R = s->R + (int64_t)o * CP; const int16_t *restrict T = TM + ((int64_t)o * O + o) * CP;
    int n = (O - o) * CP;                 /* rows p >= o: the row of o itself changes through its diagonal term */
    if (yo < 0) { for (int i = 0; i < n; i++) R[i] += T[i]; }
    else { for (int i = 0; i < n; i++) R[i] -= T[i]; }
    s->y[o] = -yo; s->code ^= (1ULL << o); s->sum += 2 * (int64_t)(-yo) * SZ[o];
}

static int64_t suffix[65];
static unsigned char *REACH;          /* REACH[pos * RW + n]: n = sum of 2*SZ over some subset of orbits pos.. */
static int64_t RW;

static void build_reach(void) {
    int64_t tot = 0;
    for (int o = 0; o < O; o++) tot += 2 * SZ[o];
    RW = tot + 1;
    free(REACH);
    REACH = calloc((size_t)(O + 1) * RW, 1);
    REACH[(int64_t)O * RW + 0] = 1;
    for (int pos = O - 1; pos >= 0; pos--)
        for (int64_t n = 0; n < RW; n++)
            REACH[(int64_t)pos * RW + n] = REACH[(int64_t)(pos + 1) * RW + n] |
                (n >= 2 * SZ[pos] ? REACH[(int64_t)(pos + 1) * RW + n - 2 * SZ[pos]] : 0);
}

typedef void (*leaf_fn)(state *, void *);

static void dfs(state *s, int pos, int64_t target, leaf_fn fn, void *ctx) {
    int64_t need = target - s->sum;                      /* sum of 2*SZ over orbits still to flip */
    if (need == 0) { fn(s, ctx); return; }
    if (pos >= O || need < 0 || need >= RW || !REACH[(int64_t)pos * RW + need]) return;
    /* flip pos */
    flip(s, pos);
    dfs(s, pos + 1, target, fn, ctx);
    flip(s, pos);
    /* keep pos at -1 */
    dfs(s, pos + 1, target, fn, ctx);
}

/* ---------------------------------------------------------------- pass 1 */
typedef struct { ent *buf; int64_t n, cap; int64_t maxm1; int overflow; } p1ctx;
static int64_t M1MAX;

static void leaf1(state *s, void *vctx) {
    p1ctx *x = (p1ctx *)vctx;
    if (x->overflow) return;
    int want[16];
    uint64_t kv = 0;
    for (int q = 0; q < K; q++) kv += RK[q] * (uint64_t)(int64_t)s->v[K1[q]];
    for (int t = 0; t < NT; t++) {
        const hent *e = lookup_key(&TA, KT[t] - kv);
        if (!e) continue;
        for (int q = 0; q < K; q++) want[q] = TG[t * K + q] - (int)s->v[K1[q]];
        for (int u = 0; u < e->cnt; u++) {
            int ai = TA.idx[e->start + u];
            const int16_t *a0 = AV + (int64_t)ai * C;
            int okk = 1;
            for (int q = 0; q < K; q++) if (a0[K1[q]] != want[q]) { okk = 0; break; }
            if (!okk) continue;
            const int16_t *a = AV + (int64_t)ai * C;
            uint64_t h = 1469598103934665603ULL;
            for (int r = 0; r < NREST; r++) h = mix(h, (uint64_t)(int64_t)(a[REST[r]] + s->v[REST[r]]));
            if (x->n == x->cap) {
                x->cap = x->cap ? 2 * x->cap : 1 << 16;
                x->buf = realloc(x->buf, sizeof(ent) * x->cap);
            }
            x->buf[x->n].fp = h; x->buf[x->n].code = ((uint64_t)ai << 40) | s->code; x->n++;
            if (x->n > x->maxm1) { x->overflow = 1; return; }
        }
    }
}

static int cmp_ent(const void *p, const void *q) {
    uint64_t a = ((const ent *)p)->fp, b = ((const ent *)q)->fp;
    return a < b ? -1 : a > b;
}

/* ---------------------------------------------------------------- pass 2 */
static ent *M1; static int64_t NM1;
typedef struct { int64_t probes, matches, sols; } p2ctx;

static void class_values(uint64_t code, int64_t *v) {
    int y[64];
    for (int o = 0; o < O; o++) y[o] = (code >> o) & 1 ? 1 : -1;
    for (int c = 0; c < C; c++) {
        int64_t s = 0;
        for (int i = 0; i < O; i++) for (int j = 0; j < O; j++) s += SI(c, i, j) * y[i] * y[j];
        v[c] = s;
    }
}

static void leaf2(state *s, void *vctx) {
    p2ctx *x = (p2ctx *)vctx;
    int want[16];
    uint64_t kv = 0, kT = 0;
    for (int q = 0; q < K; q++) { kv += RK[q] * (uint64_t)(int64_t)s->v[K1[q]]; kT += RK[q] * (uint64_t)(int64_t)TT[K1[q]]; }
    for (int t = 0; t < NT; t++) {
        const hent *e = lookup_key(&TB, kT - KT[t] - kv);
        if (!e) continue;
        for (int q = 0; q < K; q++) want[q] = TT[K1[q]] - TG[t * K + q] - (int)s->v[K1[q]];
        for (int u = 0; u < e->cnt; u++) {
            int bi = TB.idx[e->start + u];
            const int16_t *b0 = BV + (int64_t)bi * C;
            int okk = 1;
            for (int q = 0; q < K; q++) if (b0[K1[q]] != want[q]) { okk = 0; break; }
            if (!okk) continue;
            const int16_t *b = BV + (int64_t)bi * C;
            uint64_t h = 1469598103934665603ULL;
            for (int r = 0; r < NREST; r++) h = mix(h, (uint64_t)(int64_t)(TT[REST[r]] - b[REST[r]] - s->v[REST[r]]));
            x->probes++;
            int64_t lo = 0, hi = NM1;
            while (lo < hi) { int64_t mid = (lo + hi) / 2; if (M1[mid].fp < h) lo = mid + 1; else hi = mid; }
            for (int64_t m = lo; m < NM1 && M1[m].fp == h; m++) {
                x->matches++;
                int ai = (int)(M1[m].code >> 40); uint64_t cc = M1[m].code & ((1ULL << 40) - 1);
                int64_t vc[256]; class_values(cc, vc);
                const int16_t *a = AV + (int64_t)ai * C;
                int ok = 1;
                for (int c = 0; c < C && ok; c++) ok = (a[c] + b[c] + vc[c] + s->v[c] == TT[c]);
                if (ok) {
                    x->sols++;
#pragma omp critical
                    { printf("SOL %d %d %llu %llu\n", ai, bi, (unsigned long long)cc, (unsigned long long)s->code); fflush(stdout); }
                }
            }
        }
    }
}

/* run the DFS for one row sum over 2^d prefixes in parallel */
static double FRAC[3] = {1.0, 1.0, 1.0};    /* fraction of DFS prefixes enumerated in pass 1 and pass 2 */

static inline int keep_prefix(int64_t pfx, int pass) {
    if (FRAC[pass] >= 1.0) return 1;
    uint64_t h = (uint64_t)pfx * 0x9E3779B97F4A7C15ULL + (uint64_t)pass * 0xD1B54A32D192ED03ULL;
    h ^= h >> 31; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 29;
    return (double)(h >> 11) / 9007199254740992.0 < FRAC[pass];
}

static void run_pass(int rowsum, int pass, int threads, void **ctxs) {
    int64_t total = 0;
    for (int o = 0; o < O; o++) total += SZ[o];
    int64_t target = rowsum;                                   /* final sum */
    for (int o = O - 1; o >= 0; o--) suffix[o] = suffix[o + 1] + SZ[o];
    int d = O < 12 ? O / 2 : 12;
    int64_t np = 1LL << d;
    build_reach();
    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        state s; st_init(&s);
        #pragma omp for schedule(dynamic, 1)
        for (int64_t pfx = 0; pfx < np; pfx++) {
            if (!keep_prefix(pfx, pass)) continue;
            /* apply prefix choices on orbits 0..d-1 */
            for (int o = 0; o < d; o++) if ((pfx >> o) & 1) flip(&s, o);
            int64_t need = target - s.sum;
            if (need >= 0 && need < RW && REACH[(int64_t)d * RW + need]) {
                if (pass == 1) dfs(&s, d, target, leaf1, ctxs[tid]);
                else dfs(&s, d, target, leaf2, ctxs[tid]);
            }
            for (int o = d - 1; o >= 0; o--) if ((pfx >> o) & 1) flip(&s, o);     /* undo in reverse order */
            if (pass == 1 && ((p1ctx *)ctxs[tid])->n > M1MAX / threads + 1) ((p1ctx *)ctxs[tid])->overflow = 1;
        }
        free(s.R); free(s.v);
    }
    (void)total;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: kstream INPUT MAXM1 [THREADS]\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    int hdr[8];
    if (fread(hdr, sizeof(int), 8, f) != 8) return 2;
    O = hdr[0]; C = hdr[1]; K = hdr[2]; NT = hdr[3]; RC = hdr[4]; RD = hdr[5]; NA = hdr[6]; NB = hdr[7];
    if (O > 40 || C > 255 || K > 15) { fprintf(stderr, "limits\n"); return 2; }
    SZ = malloc(sizeof(int) * O); S = malloc(sizeof(int) * C * O * O); K1 = malloc(sizeof(int) * K);
    TG = malloc(sizeof(int) * NT * K); TT = malloc(sizeof(int) * C);
    AV = malloc(sizeof(int16_t) * (int64_t)NA * C); BV = malloc(sizeof(int16_t) * (int64_t)NB * C);
    if (fread(SZ, sizeof(int), O, f) != (size_t)O || fread(S, sizeof(int), (size_t)C * O * O, f) != (size_t)C * O * O ||
        fread(K1, sizeof(int), K, f) != (size_t)K || fread(TG, sizeof(int), (size_t)NT * K, f) != (size_t)NT * K ||
        fread(TT, sizeof(int), C, f) != (size_t)C || fread(AV, sizeof(int16_t), (size_t)NA * C, f) != (size_t)NA * C ||
        fread(BV, sizeof(int16_t), (size_t)NB * C, f) != (size_t)NB * C) { fprintf(stderr, "short input\n"); return 2; }
    fclose(f);
    M1MAX = atoll(argv[2]);
    int threads = argc > 3 ? atoi(argv[3]) : omp_get_max_threads();
    if (argc > 4) FRAC[1] = atof(argv[4]);
    if (argc > 5) FRAC[2] = atof(argv[5]);
    memset(inK1, 0, sizeof(inK1));
    for (int q = 0; q < K; q++) inK1[K1[q]] = 1;
    NREST = 0;
    for (int c = 0; c < C; c++) if (!inK1[c]) REST[NREST++] = c;
    CP = (C + 15) / 16 * 16;
    TM = aligned_alloc(64, sizeof(int16_t) * ((O * O * CP + 31) / 32) * 32);
    D4 = aligned_alloc(64, sizeof(int32_t) * ((O * CP + 15) / 16) * 16);
    for (int o = 0; o < O; o++) for (int p = 0; p < O; p++) for (int c = 0; c < CP; c++)
        TM[((int64_t)o * O + p) * CP + c] = c < C ? 2 * (SI(c, p, o) + SI(c, o, p)) : 0;
    for (int o = 0; o < O; o++) for (int c = 0; c < CP; c++) D4[o * CP + c] = c < C ? 4 * SI(c, o, o) : 0;
    {
        uint64_t x = 0x243F6A8885A308D3ULL;
        for (int q = 0; q < K; q++) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; RK[q] = x | 1; }
        KT = malloc(sizeof(uint64_t) * NT);
        for (int t = 0; t < NT; t++) KT[t] = vkey(TG + t * K);
    }
    build_table(AV, NA, &TA);
    build_table(BV, NB, &TB);
    /* pass 1 */
    p1ctx *c1 = calloc(threads, sizeof(p1ctx));
    void **ctx = malloc(sizeof(void *) * threads);
    for (int t = 0; t < threads; t++) { ctx[t] = &c1[t]; c1[t].maxm1 = M1MAX / threads + 1; }
    double t0 = omp_get_wtime();
    run_pass(RC, 1, threads, ctx);
    NM1 = 0; int ovf = 0;
    for (int t = 0; t < threads; t++) { NM1 += c1[t].n; ovf |= c1[t].overflow; }
    fprintf(stderr, "pass 1: %lld entries (%.0fs)%s\n", (long long)NM1, omp_get_wtime() - t0, ovf ? " OVERFLOW" : "");
    printf("STATS pass1 %lld\n", (long long)NM1); fflush(stdout);
    if (ovf) { printf("OVERFLOW\n"); return 3; }
    M1 = malloc(sizeof(ent) * (NM1 + 1));
    int64_t off = 0;
    for (int t = 0; t < threads; t++) { memcpy(M1 + off, c1[t].buf, sizeof(ent) * c1[t].n); off += c1[t].n; free(c1[t].buf); }
    qsort(M1, NM1, sizeof(ent), cmp_ent);
    fprintf(stderr, "sorted (%.0fs)\n", omp_get_wtime() - t0);
    /* pass 2 */
    p2ctx *c2 = calloc(threads, sizeof(p2ctx));
    for (int t = 0; t < threads; t++) ctx[t] = &c2[t];
    run_pass(RD, 2, threads, ctx);
    int64_t probes = 0, matches = 0, sols = 0;
    for (int t = 0; t < threads; t++) { probes += c2[t].probes; matches += c2[t].matches; sols += c2[t].sols; }
    printf("STATS pass2 probes %lld matches %lld solutions %lld seconds %.0f\n", (long long)probes, (long long)matches,
           (long long)sols, omp_get_wtime() - t0);
    return 0;
}
