/* pgp.c -- exhaustive search for periodic Golay pairs over Z_2 x Z_w that are invariant under a multiplier group H.
 *
 *   w prime, w = 1 (mod 4); H = the subgroup of order h of Z_w^* (h odd, h | w-1); m = (w-1)/h orbits on Z_w \ {0}.
 *   A pair (a, b) of +-1 sequences on Gamma = Z_2 x Z_w, constant on the orbits {(e, hx)} of H, is a periodic Golay
 *   pair when PAF_a + PAF_b = 2 |Gamma| delta.  Write a = (a0, a1) (the two layers) and
 *       A = (a0 + a1)/2,  C = (a0 - a1)/2,  B = (b0 + b1)/2,  D = (b0 - b1)/2.
 *   A, C are ternary on Z_w with complementary supports (so are B, D) and the pair is Golay iff
 *       A A* + B B* = w delta   and   C C* + D D* = w delta   (in Z[Z_w]).
 *   A "row" is a pair (X, Y) of H-invariant ternary sequences with X X* + Y Y* = w delta.  A Golay pair is two rows
 *   (A, B), (C, D) whose supports are complementary (supp C = Z_w \ supp A, supp D = Z_w \ supp B).
 *
 *   Stage 1 (`rows`)  : enumerates all rows (X, Y) up to the multiplier rotations, by a join on the 16-vector of
 *                       autocorrelations within buckets given by the parity key  K(S)(k) = #{x : x, x+s_k in S} mod 2
 *                       (a function of the support only; K(supp X) = K(supp Y) for every row).
 *   Stage 2 (`join`)  : the support join: pairs of rows with complementary supports; prints Golay pairs.
 *
 *   cc -O3 -march=native -pthread -o pgp pgp.c -lm
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <stdatomic.h>

typedef int8_t v16 __attribute__((vector_size(16)));
typedef float v4f __attribute__((vector_size(16)));
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;

/* ------------------------------------------------------------------ parameters and tables */
static int P, HH, M, H2, G0;
static int *dlog, *cls;
static u8 Ncnt[16][32][32];
static v16 Wv[32][32], Dv[32], Uv[32];
static float ER[32][16], EI[32][16];
static u16 Cm[32][32], Dm[32];
static u32 FULL;                     /* (1 << M) - 1 */
static u8 *canon_tab;                /* canon_tab[key]: key minimal among its H2 rotations */
static int PSD_ON = 1;
static int REVCANON = 1;            /* keep only X <= X* (reversal = rotation by m/2) in each component */
static float PSD_EPS = 0.05f;

static int is_prime(int n) { if (n < 2) return 0; for (int d = 2; d * d <= n; d++) if (n % d == 0) return 0; return 1; }
static int primroot(int p) {
    for (int a = 2; a < p; a++) {
        int ok = 1;
        for (int q = 2; q < p; q++) if ((p - 1) % q == 0 && is_prime(q)) {
            long x = 1; for (int i = 0; i < (p - 1) / q; i++) x = x * a % p;
            if (x == 1) { ok = 0; break; }
        }
        if (ok) return a;
    }
    return -1;
}
static inline u32 rotr_m(u32 x, int s) {       /* bit j of the result = bit (j+s) mod M of x */
    s %= M; if (s < 0) s += M; if (!s) return x;
    return ((x >> s) | (x << (M - s))) & FULL;
}
static inline u32 kshift(u32 key, int s) {      /* bit k of the result = bit (k+s) mod H2 of key */
    u32 f = (1u << H2) - 1; s %= H2; if (s < 0) s += H2; if (!s) return key;
    return ((key >> s) | (key << (H2 - s))) & f;
}
static inline u32 fold(u32 S) { return (S & ((1u << H2) - 1)) ^ (S >> H2); }

static void die(const char *msg) { fprintf(stderr, "pgp: %s\n", msg); exit(2); }

static void init(int p, int h) {
    P = p; HH = h;
    if (!is_prime(p) || p % 4 != 1) die("w must be a prime = 1 (mod 4)");
    if (h % 2 == 0 || (p - 1) % h) die("h must be odd and divide w - 1");
    M = (p - 1) / h; H2 = M / 2; FULL = M == 32 ? 0xFFFFFFFFu : ((1u << M) - 1);
    if (M > 32 || H2 > 16) die("need m = (w-1)/h <= 32");
    G0 = primroot(p);
    dlog = calloc(p, sizeof(int)); cls = calloc(p, sizeof(int));
    long x = 1;
    for (int i = 0; i < p - 1; i++) { dlog[x] = i; cls[x] = i % M; x = x * G0 % p; }
    cls[0] = -1;
    memset(Ncnt, 0, sizeof Ncnt);
    for (int k = 0; k < H2; k++) {
        long s = 1; for (int i = 0; i < k; i++) s = s * G0 % p;
        for (int x2 = 1; x2 < p; x2++) { int y = (x2 + s) % p; if (y == 0) continue; Ncnt[k][cls[x2]][cls[y]]++; }
    }
    memset(Wv, 0, sizeof Wv); memset(Dv, 0, sizeof Dv); memset(Uv, 0, sizeof Uv);
    memset(Cm, 0, sizeof Cm); memset(Dm, 0, sizeof Dm);
    for (int i = 0; i < M; i++) {
        for (int k = 0; k < H2; k++) {
            Dv[i][k] = Ncnt[k][i][i];
            if (Ncnt[k][i][i] & 1) Dm[i] |= 1u << k;
            for (int j = 0; j < M; j++) if (j != i) {
                Wv[i][j][k] = Ncnt[k][i][j] + Ncnt[k][j][i];
                if ((Ncnt[k][i][j] + Ncnt[k][j][i]) & 1) Cm[i][j] |= 1u << k;
            }
        }
        Uv[i][i % H2] = 1;
    }
    /* Gaussian periods eta_j = sum_{t<h} zeta^(g^(j + m t)); F_k = X_inf + sum_j X_j eta_{j+k} */
    double *er = malloc(sizeof(double) * M), *ei = malloc(sizeof(double) * M);
    for (int j = 0; j < M; j++) {
        double re = 0, im = 0;
        for (int t = 0; t < HH; t++) {
            int e = j + M * t; long y = 1; for (int i = 0; i < e; i++) y = y * G0 % p;
            double ang = 2 * M_PI * (double)y / p; re += cos(ang); im += sin(ang);
        }
        er[j] = re; ei[j] = im;
    }
    memset(ER, 0, sizeof ER); memset(EI, 0, sizeof EI);
    for (int j = 0; j < M; j++) for (int k = 0; k < H2; k++) { ER[j][k] = (float)er[(j + k) % M]; EI[j][k] = (float)ei[(j + k) % M]; }
    free(er); free(ei);
    canon_tab = calloc(1u << H2, 1);
    for (u32 key = 0; key < (1u << H2); key++) {
        u32 mn = key; for (int s = 1; s < H2; s++) { u32 t = kshift(key, s); if (t < mn) mn = t; }
        canon_tab[key] = (mn == key);
    }
}

/* ------------------------------------------------------------------ reference (direct) computations, for tests */
/* expand an orbit-form sequence (xinf, x_0..x_{m-1}) to Z_w */
static void expand(int xinf, const int *xo, int *out) {
    out[0] = xinf;
    for (int x = 1; x < P; x++) out[x] = xo[cls[x]];
}
static void paf_direct(const int *seq, int *paf /* indexed by shift s = 0..w-1 */) {
    for (int s = 0; s < P; s++) { int t = 0; for (int x = 0; x < P; x++) t += seq[x] * seq[(x + s) % P]; paf[s] = t; }
}
/* PAF vector at the class representatives s_k = g^k, k < H2, from orbit form */
static void paf_classes_direct(int xinf, const int *xo, int *out) {
    int seq[512], paf[512]; expand(xinf, xo, seq); paf_direct(seq, paf);
    for (int k = 0; k < H2; k++) { long s = 1; for (int i = 0; i < k; i++) s = s * G0 % P; out[k] = paf[s]; }
}
static u32 key_direct(int e, u32 S) {     /* parity key from the definition */
    int xo[32], xinf = e; for (int j = 0; j < M; j++) xo[j] = (S >> j) & 1;
    int v[16]; paf_classes_direct(xinf, xo, v);
    u32 key = 0; for (int k = 0; k < H2; k++) key |= (u32)(v[k] & 1) << k;
    return key;
}
/* key through the tables Dm, Cm, fold */
static u32 key_tab(int e, u32 S) {
    u32 key = 0;
    for (int i = 0; i < M; i++) if ((S >> i) & 1) {
        key ^= Dm[i];
        for (int j = 0; j < i; j++) if ((S >> j) & 1) key ^= Cm[i][j];
    }
    if (e) key ^= fold(S);
    return key;
}

/* ------------------------------------------------------------------ types of rows */
typedef struct { int id, x, y, exi, eyi, sx, sy, ax, ay, swapcanon; } type_t;
static type_t *TYPES; static int NTYPES;
static i64 binom[40][40];

static void build_types(int swap_reduce) {
    for (int n = 0; n < 40; n++) { binom[n][0] = 1; for (int k = 1; k <= n; k++) binom[n][k] = binom[n - 1][k - 1] + (k <= n - 1 ? binom[n - 1][k] : 0); }
    TYPES = calloc(4096, sizeof(type_t)); NTYPES = 0;
    for (int x = 1; x * x < P; x++) {
        int y2 = P - x * x, y = (int)lround(sqrt((double)y2)); if (y * y != y2 || y < 1) continue;
        {
            int xs = x, ys = y;
            for (int exi = -1; exi <= 1; exi++) {
                if ((xs - exi) % HH) continue;
                int sx = (xs - exi) / HH;
                for (int eyi = -1; eyi <= 1; eyi++) {
                    if ((ys - eyi) % HH) continue;
                    int sy = (ys - eyi) / HH;
                    for (int ax = abs(sx); ax <= M; ax++) {
                        if (((ax + sx) & 1)) continue;
                        int rest = P - (exi != 0) - (eyi != 0) - HH * ax;
                        if (rest < 0 || rest % HH) continue;
                        int ay = rest / HH;
                        if (ay < abs(sy) || ay > M || ((ay + sy) & 1)) continue;
                        type_t t = {NTYPES, xs, ys, exi, eyi, sx, sy, ax, ay, 0};
                        /* swap-canonical: tuple (x,exi,sx,ax | y,eyi,sy,ay) <= its swap */
                        int a1[4] = {xs, exi, sx, ax}, a2[4] = {ys, eyi, sy, ay}, c = 0;
                        for (int i = 0; i < 4 && !c; i++) c = (a1[i] > a2[i]) - (a1[i] < a2[i]);
                        t.swapcanon = c <= 0;
                        if (swap_reduce && !t.swapcanon) continue;
                        TYPES[NTYPES++] = t; TYPES[NTYPES - 1].id = NTYPES - 1;
                    }
                }
            }
        }
    }
}
static void print_types(void) {
    for (int i = 0; i < NTYPES; i++) {
        type_t *t = &TYPES[i];
        printf("type %d: sums (%d,%d)  X_inf=%d sigma_X=%d alpha_X=%d   Y_inf=%d sigma_Y=%d alpha_Y=%d  %s\n", t->id, t->x, t->y, t->exi, t->sx, t->ax, t->eyi, t->sy, t->ay,
               t->swapcanon ? "swap-canonical" : "");
    }
}

/* ------------------------------------------------------------------ support index: supports of size alpha by key */
typedef struct { u32 *mask; u32 *start; i64 n; } supidx_t;   /* start has 2^H2 + 1 entries */
static supidx_t IDX[2][33]; static int idx_built[2][33];
static int ALLKEYS = 0;

typedef struct { u64 *buf; i64 n, cap; } vec64;
static void v64push(vec64 *v, u64 x) { if (v->n == v->cap) { v->cap = v->cap ? 2 * v->cap : 1024; v->buf = realloc(v->buf, v->cap * 8); if (!v->buf) die("oom"); } v->buf[v->n++] = x; }

/* one pass over all 2^M supports for the sizes in alphamask[e]; keys of supports kept: canonical (or all) */
static void build_index(int e, u64 alphamask /* bit alpha */) {
    u16 *Bq_lo = 0, *Bq_hi = 0;
    int hb = H2;
    u32 nlo = 1u << hb;
    Bq_lo = malloc(nlo * 2); Bq_hi = malloc(nlo * 2);
    for (u32 s = 0; s < nlo; s++) {
        u32 klo = 0, khi = 0;
        for (int i = 0; i < hb; i++) if ((s >> i) & 1) {
            klo ^= Dm[i]; khi ^= Dm[i + hb];
            for (int j = 0; j < i; j++) if ((s >> j) & 1) { klo ^= Cm[i][j]; khi ^= Cm[i + hb][j + hb]; }
        }
        Bq_lo[s] = klo; Bq_hi[s] = khi;
    }
    vec64 *out = calloc(33, sizeof(vec64));
    u8 *pc = malloc(nlo); for (u32 s = 0; s < nlo; s++) pc[s] = __builtin_popcount(s);
    for (u32 shi = 0; shi < nlo; shi++) {
        u32 R[16]; memset(R, 0, sizeof R);
        for (int i = 0; i < hb; i++) { u32 r = 0; for (int j = 0; j < hb; j++) if ((shi >> j) & 1) r ^= Cm[i][j + hb]; R[i] = r; }
        u32 base = Bq_hi[shi]; int pch = pc[shi];
        u32 L = 0, glo = 0;                       /* Gray code over s_lo */
        for (u32 it = 0; it < nlo; it++) {
            if (it) { u32 bit = __builtin_ctz(it); glo ^= 1u << bit; L ^= R[bit]; }
            int a = pch + pc[glo];
            if (!((alphamask >> a) & 1ull)) continue;
            u32 key = base ^ Bq_lo[glo] ^ L;
            u32 S = glo | (shi << hb);
            if (e) key ^= fold(S);
            if (REVCANON && S > rotr_m(S, H2)) continue;
            if (ALLKEYS || canon_tab[key]) v64push(&out[a], ((u64)key << 32) | S);
        }
    }
    for (int a = 0; a <= 32; a++) if ((alphamask >> a) & 1ull) {
        supidx_t *ix = &IDX[e][a]; u32 nk = 1u << hb;
        ix->start = calloc(nk + 1, 4); ix->n = out[a].n; ix->mask = malloc(out[a].n * 4 + 4);
        for (i64 i = 0; i < out[a].n; i++) ix->start[(out[a].buf[i] >> 32) + 1]++;
        for (u32 k = 0; k < nk; k++) ix->start[k + 1] += ix->start[k];
        u32 *fill = malloc(nk * 4); memcpy(fill, ix->start, nk * 4);
        for (i64 i = 0; i < out[a].n; i++) ix->mask[fill[out[a].buf[i] >> 32]++] = (u32)out[a].buf[i];
        free(fill); free(out[a].buf); idx_built[e][a] = 1;
    }
    free(out); free(pc); free(Bq_lo); free(Bq_hi);
}

/* ------------------------------------------------------------------ enumeration of sign patterns on a support */
typedef struct { v16 key; u32 pos, neg; } leaf_t;
typedef void (*leaf_fn)(void *ctx, u32 pos, u32 neg, v16 paf);

typedef struct { u64 leaves, psdpass; } enum_stats;

static inline int psd_ok(const v4f *fr, const v4f *fi, float lim) {
    /* all |F_k|^2 <= lim for k < H2 (unused lanes are zero) */
    v4f l = {lim, lim, lim, lim};
    for (int q = 0; q < 4; q++) {
        v4f s = fr[q] * fr[q] + fi[q] * fi[q];
        if (s[0] > l[0] || s[1] > l[1] || s[2] > l[2] || s[3] > l[3]) return 0;
    }
    return 1;
}

static void enum_support(const int *pos, int a, int np, int xinf, int symsup, leaf_fn emit, void *ctx, enum_stats *st) {
    v16 Acc[33][32]; v16 paf[33]; v4f fr[33][4], fi[33][4];
    int cnt[33]; u32 pm[33], nm[33]; int state[33];
    float lim = (float)P + PSD_EPS;
    for (int j = 0; j < a; j++) { v16 u = Uv[pos[j]]; Acc[0][j] = xinf > 0 ? u : (xinf < 0 ? -u : (v16){0}); }
    paf[0] = (v16){0};
    for (int q = 0; q < 4; q++) { fr[0][q] = (v4f){0, 0, 0, 0}; fi[0][q] = (v4f){0, 0, 0, 0}; }
    if (xinf) for (int k = 0; k < H2; k++) ((float *)fr[0])[k] = (float)xinf;
    cnt[0] = 0; pm[0] = nm[0] = 0; state[0] = 0;
    int d = 0;
    while (d >= 0) {
        if (d == a) {
            st->leaves++;
            if (symsup) {                      /* reversal-symmetric support: keep X <= X* (pos mask first, then neg mask) */
                u32 rp = rotr_m(pm[a], H2), rn = rotr_m(nm[a], H2);
                if (pm[a] > rp || (pm[a] == rp && nm[a] > rn)) { d--; continue; }
            }
            if (!PSD_ON || psd_ok(fr[a], fi[a], lim)) { st->psdpass++; emit(ctx, pm[a], nm[a], paf[a]); }
            d--; continue;
        }
        if (state[d] == 2) { d--; continue; }
        int s = state[d] == 0 ? 1 : -1; state[d]++;
        int c = cnt[d] + (s > 0), neg = d + 1 - c;
        if (c > np || neg > a - np) continue;
        int t = pos[d];
        paf[d + 1] = paf[d] + Dv[t] + (s > 0 ? Acc[d][d] : -Acc[d][d]);
        if (s > 0) { for (int j = d + 1; j < a; j++) Acc[d + 1][j] = Acc[d][j] + Wv[t][pos[j]]; }
        else       { for (int j = d + 1; j < a; j++) Acc[d + 1][j] = Acc[d][j] - Wv[t][pos[j]]; }
        const v4f *er = (const v4f *)ER[t], *ei = (const v4f *)EI[t];
        if (s > 0) for (int q = 0; q < 4; q++) { fr[d + 1][q] = fr[d][q] + er[q]; fi[d + 1][q] = fi[d][q] + ei[q]; }
        else       for (int q = 0; q < 4; q++) { fr[d + 1][q] = fr[d][q] - er[q]; fi[d + 1][q] = fi[d][q] - ei[q]; }
        cnt[d + 1] = c; pm[d + 1] = pm[d] | (s > 0 ? 1u << t : 0); nm[d + 1] = nm[d] | (s < 0 ? 1u << t : 0);
        d++; state[d] = 0;
    }
}

/* ------------------------------------------------------------------ hash of stored leaves */
typedef struct { u64 k0, k1; u32 pos, neg; } ent_t;
typedef struct { ent_t *e; i64 n, cap; u32 *tab; u64 tmask; int bits; } store_t;

static inline u64 hash16(u64 k0, u64 k1, int bits) {
    u64 x = k0 * 0x9E3779B97F4A7C15ull ^ ((k1 << 29 | k1 >> 35) * 0xC2B2AE3D27D4EB4Full);
    x ^= x >> 31; x *= 0xD6E8FEB86659FD93ull; x ^= x >> 32;
    return x >> (64 - bits);
}
static void store_push(store_t *s, u32 pos, u32 neg, v16 paf) {
    if (s->n == s->cap) { s->cap = s->cap ? 2 * s->cap : 4096; s->e = realloc(s->e, s->cap * sizeof(ent_t)); if (!s->e) die("oom"); }
    ent_t *x = &s->e[s->n++]; memcpy(&x->k0, &paf, 8); memcpy(&x->k1, ((char *)&paf) + 8, 8); x->pos = pos; x->neg = neg;
}
static void store_build(store_t *s) {
    int bits = 4; while ((1ll << bits) < 2 * s->n + 2) bits++;
    s->bits = bits; s->tmask = (1ull << bits) - 1; s->tab = malloc(((size_t)1 << bits) * 4); memset(s->tab, 0xFF, ((size_t)1 << bits) * 4);
    for (i64 i = 0; i < s->n; i++) {
        u64 h = hash16(s->e[i].k0, s->e[i].k1, bits);
        while (s->tab[h] != 0xFFFFFFFFu) h = (h + 1) & s->tmask;
        s->tab[h] = (u32)i;
    }
}
static void store_free(store_t *s) { free(s->e); free(s->tab); memset(s, 0, sizeof *s); }

/* ------------------------------------------------------------------ rows */
typedef struct __attribute__((packed)) { u32 xp, xn, yp, yn; u16 type, pad; } row_t;      /* 20 bytes */

typedef struct {
    int side;                 /* 0: stored side is X, 1: stored side is Y */
    store_t *st;
    type_t *t;
    row_t *rows; i64 nrows, caprows;
    u64 nmatch_chk, nprobe;
} bucket_ctx;

static void leaf_store(void *c, u32 pos, u32 neg, v16 paf) { store_push((store_t *)c, pos, neg, paf); }
static void leaf_probe(void *c, u32 pos, u32 neg, v16 paf) {
    bucket_ctx *b = c; store_t *s = b->st;
    v16 q = -paf; u64 k0, k1; memcpy(&k0, &q, 8); memcpy(&k1, ((char *)&q) + 8, 8);
    b->nprobe++;
    u64 h = hash16(k0, k1, s->bits);
    while (s->tab[h] != 0xFFFFFFFFu) {
        ent_t *e = &s->e[s->tab[h]];
        if (e->k0 == k0 && e->k1 == k1) {
            if (b->nrows == b->caprows) { b->caprows = b->caprows ? 2 * b->caprows : 256; b->rows = realloc(b->rows, b->caprows * sizeof(row_t)); }
            row_t *r = &b->rows[b->nrows++]; r->type = b->t->id; r->pad = 0;
            if (b->side == 0) { r->xp = e->pos; r->xn = e->neg; r->yp = pos; r->yn = neg; }
            else              { r->yp = e->pos; r->yn = e->neg; r->xp = pos; r->xn = neg; }
        }
        h = (h + 1) & s->tmask;
    }
}

/* ------------------------------------------------------------------ stage 1 driver */
typedef struct { int tid; u32 key; double work; } bucket_t;
static bucket_t *BUCKETS; static i64 NBUCKETS;
static atomic_llong NEXTB;
static pthread_mutex_t OUTMTX = PTHREAD_MUTEX_INITIALIZER;
static char OUTDIR[1024];
static FILE *DONELOG;
static atomic_ullong G_leaves, G_pass, G_stored, G_probe, G_rows, G_buckets;
static int NTHREADS = 1;
static double T0;
static int QUIET = 0;

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

static int decode_positions(u32 mask, int *pos) { int n = 0; for (int j = 0; j < M; j++) if ((mask >> j) & 1) pos[n++] = j; return n; }

static double thread_cpu(void) { struct timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static void process_bucket(const bucket_t *bk, i64 *cap_exceeded) {
    double tc0 = thread_cpu();
    type_t *t = &TYPES[bk->tid];
    int ex = t->exi != 0, ey = t->eyi != 0;
    supidx_t *sx = &IDX[ex][t->ax], *sy = &IDX[ey][t->ay];
    u32 kx0 = sx->start[bk->key], kx1 = sx->start[bk->key + 1], ky0 = sy->start[bk->key], ky1 = sy->start[bk->key + 1];
    if (kx0 == kx1 || ky0 == ky1) return;
    int npx = (t->ax + t->sx) / 2, npy = (t->ay + t->sy) / 2;
    double nx = (double)(kx1 - kx0) * binom[t->ax][npx], ny = (double)(ky1 - ky0) * binom[t->ay][npy];
    int side = nx <= ny ? 0 : 1;               /* store the side with fewer candidates */
    store_t st; memset(&st, 0, sizeof st);
    enum_stats es = {0, 0}, ep = {0, 0};
    int pos[32];
    /* store side */
    u32 s0 = side == 0 ? kx0 : ky0, s1 = side == 0 ? kx1 : ky1;
    supidx_t *si = side == 0 ? sx : sy; int sa = side == 0 ? t->ax : t->ay, snp = side == 0 ? npx : npy, sinf = side == 0 ? t->exi : t->eyi;
    for (u32 i = s0; i < s1; i++) { int n = decode_positions(si->mask[i], pos); enum_support(pos, n, snp, sinf, REVCANON && si->mask[i] == rotr_m(si->mask[i], H2), leaf_store, &st, &es); }
    store_build(&st);
    bucket_ctx bc = {0}; bc.side = side; bc.st = &st; bc.t = t;
    u32 p0 = side == 0 ? ky0 : kx0, p1 = side == 0 ? ky1 : kx1;
    supidx_t *pi = side == 0 ? sy : sx; int pa = side == 0 ? t->ay : t->ax, pnp = side == 0 ? npy : npx, pinf = side == 0 ? t->eyi : t->exi;
    for (u32 i = p0; i < p1; i++) { int n = decode_positions(pi->mask[i], pos); enum_support(pos, n, pnp, pinf, REVCANON && pi->mask[i] == rotr_m(pi->mask[i], H2), leaf_probe, &bc, &ep); }
    (void)sa; (void)pa;
    G_leaves += es.leaves + ep.leaves; G_pass += es.psdpass + ep.psdpass; G_stored += st.n; G_probe += bc.nprobe; G_rows += bc.nrows; G_buckets++;
    pthread_mutex_lock(&OUTMTX);
    if (bc.nrows > 0) {
        char fn[1200]; snprintf(fn, sizeof fn, "%s/rows/t%d_k%u.bin", OUTDIR, bk->tid, bk->key);
        FILE *f = fopen(fn, "wb"); if (!f) die("cannot write rows file");
        fwrite(bc.rows, sizeof(row_t), bc.nrows, f); fclose(f);
    }
    fprintf(DONELOG, "%d %u %lld %llu %llu %lld %.4f\n", bk->tid, bk->key, (long long)bc.nrows, (unsigned long long)(es.leaves + ep.leaves),
            (unsigned long long)(es.psdpass + ep.psdpass), (long long)st.n, thread_cpu() - tc0);
    fflush(DONELOG);
    pthread_mutex_unlock(&OUTMTX);
    free(bc.rows); store_free(&st);
}

static void *worker(void *arg) {
    (void)arg; i64 ce = 0;
    for (;;) {
        i64 i = atomic_fetch_add(&NEXTB, 1);
        if (i >= NBUCKETS) break;
        process_bucket(&BUCKETS[i], &ce);
        if (!QUIET && (i % 200 == 0)) fprintf(stderr, "[stage1] bucket %lld/%lld  leaves %.3g  rows %llu  t=%.1fs\n", (long long)i, (long long)NBUCKETS, (double)G_leaves, (unsigned long long)G_rows, now() - T0);
    }
    return 0;
}
static int cmp_work(const void *a, const void *b) { double x = ((bucket_t *)a)->work, y = ((bucket_t *)b)->work; return (x < y) - (x > y); }

static void load_done(const char *fn, i64 **doneset, i64 *nd) {
    FILE *f = fopen(fn, "r"); *nd = 0; *doneset = 0; if (!f) return;
    int t; unsigned key; long long a, b, c, d; double tt; i64 cap = 0;
    while (fscanf(f, "%d %u %lld %lld %lld %lld %lf", &t, &key, &a, &b, &c, &d, &tt) == 7) {
        if (*nd == cap) { cap = cap ? 2 * cap : 1024; *doneset = realloc(*doneset, cap * 8); }
        (*doneset)[(*nd)++] = ((i64)t << 32) | key;
    }
    fclose(f);
}
static int cmp_i64(const void *a, const void *b) { i64 x = *(i64 *)a, y = *(i64 *)b; return (x > y) - (x < y); }

/* argument parsing helpers */
static const char *argval(int argc, char **argv, const char *name, const char *def) {
    for (int i = 0; i < argc - 1; i++) if (!strcmp(argv[i], name)) return argv[i + 1];
    return def;
}
static int argflag(int argc, char **argv, const char *name) { for (int i = 0; i < argc; i++) if (!strcmp(argv[i], name)) return 1; return 0; }

static void select_types(int argc, char **argv, int *sel) {
    const char *ts = argval(argc, argv, "--types", 0);
    for (int i = 0; i < NTYPES; i++) sel[i] = 1;
    if (ts) {
        for (int i = 0; i < NTYPES; i++) sel[i] = 0;
        char *c = strdup(ts), *tok = strtok(c, ",");
        while (tok) { int a, b; if (sscanf(tok, "%d-%d", &a, &b) == 2) { for (int i = a; i <= b && i < NTYPES; i++) sel[i] = 1; } else { a = atoi(tok); if (a < NTYPES) sel[a] = 1; } tok = strtok(0, ","); }
        free(c);
    }
}

static int cmd_rows(int argc, char **argv) {
    /* pgp rows w h OUTDIR [--threads T] [--all-keys] [--swap-reduce] [--types ...] [--no-psd] [--limit N] [--stride S] */
    if (argc < 5) die("usage: pgp rows w h OUTDIR [options]");
    int w = atoi(argv[2]), h = atoi(argv[3]); snprintf(OUTDIR, sizeof OUTDIR, "%s", argv[4]);
    init(w, h);
    ALLKEYS = argflag(argc, argv, "--all-keys"); PSD_ON = !argflag(argc, argv, "--no-psd"); REVCANON = !argflag(argc, argv, "--no-rev");
    NTHREADS = atoi(argval(argc, argv, "--threads", "1")); QUIET = argflag(argc, argv, "--quiet");
    double maxwork = atof(argval(argc, argv, "--maxwork", "0")), minwork = atof(argval(argc, argv, "--minwork", "0"));
    i64 limit = atoll(argval(argc, argv, "--limit", "0")), stride = atoll(argval(argc, argv, "--stride", "1")), offset = atoll(argval(argc, argv, "--offset", "0"));
    build_types(argflag(argc, argv, "--swap-reduce"));
    int *sel = calloc(NTYPES + 1, sizeof(int)); select_types(argc, argv, sel);
    char d[1200]; mkdir(OUTDIR, 0755); snprintf(d, sizeof d, "%s/rows", OUTDIR); mkdir(d, 0755);
    T0 = now();
    u64 am[2] = {0, 0};
    for (int i = 0; i < NTYPES; i++) if (sel[i]) { am[TYPES[i].exi != 0] |= 1ull << TYPES[i].ax; am[TYPES[i].eyi != 0] |= 1ull << TYPES[i].ay; }
    for (int e = 0; e < 2; e++) if (am[e]) { build_index(e, am[e]); if (!QUIET) fprintf(stderr, "[stage1] support index e=%d built (t=%.1fs)\n", e, now() - T0); }
    /* buckets */
    i64 *done; i64 nd; snprintf(d, sizeof d, "%s/done.log", OUTDIR); load_done(d, &done, &nd); qsort(done, nd, 8, cmp_i64);
    BUCKETS = malloc(sizeof(bucket_t) * (size_t)NTYPES * (1u << H2)); NBUCKETS = 0; i64 seen = 0;
    (void)seen;
    for (int i = 0; i < NTYPES; i++) if (sel[i]) {
        type_t *t = &TYPES[i]; supidx_t *sx = &IDX[t->exi != 0][t->ax], *sy = &IDX[t->eyi != 0][t->ay];
        int npx = (t->ax + t->sx) / 2, npy = (t->ay + t->sy) / 2;
        for (u32 key = 0; key < (1u << H2); key++) {
            u32 nx = sx->start[key + 1] - sx->start[key], ny = sy->start[key + 1] - sy->start[key];
            if (!nx || !ny) continue;
            i64 id = ((i64)i << 32) | key;
            if (nd && bsearch(&id, done, nd, 8, cmp_i64)) continue;
            { double wk = (double)nx * binom[t->ax][npx] + (double)ny * binom[t->ay][npy]; if ((maxwork > 0 && wk > maxwork) || wk < minwork) continue; }
            { i64 ib = seen++; if (ib < offset || (ib - offset) % stride != 0) continue; }
            BUCKETS[NBUCKETS].tid = i; BUCKETS[NBUCKETS].key = key;
            BUCKETS[NBUCKETS].work = (double)nx * binom[t->ax][npx] + (double)ny * binom[t->ay][npy]; NBUCKETS++;
            if (limit && NBUCKETS >= limit) goto bdone;
        }
    }
bdone:
    qsort(BUCKETS, NBUCKETS, sizeof(bucket_t), cmp_work);
    if (!QUIET) fprintf(stderr, "[stage1] %lld buckets to do (%lld already done), types %d, m=%d h2=%d\n", (long long)NBUCKETS, (long long)nd, NTYPES, M, H2);
    DONELOG = fopen(d, "a"); if (!DONELOG) die("cannot open done.log");
    pthread_t th[256]; atomic_store(&NEXTB, 0);
    for (int i = 0; i < NTHREADS; i++) pthread_create(&th[i], 0, worker, 0);
    for (int i = 0; i < NTHREADS; i++) pthread_join(th[i], 0);
    fclose(DONELOG);
    double el = now() - T0;
    printf("stage1 done: buckets %llu leaves %llu psd_pass %llu stored %llu probes %llu rows %llu wall %.2fs\n", (unsigned long long)G_buckets, (unsigned long long)G_leaves,
           (unsigned long long)G_pass, (unsigned long long)G_stored, (unsigned long long)G_probe, (unsigned long long)G_rows, el);
    return 0;
}

/* ------------------------------------------------------------------ stage 2: the support join */
typedef struct { u32 xs, ys; u32 row; } hent_t;     /* support pair -> row index (xs, ys = support masks) */
static row_t *ROWS; static i64 NROWS;
static void load_rows(const char *dir, int *typesel) {
    char d[1200]; snprintf(d, sizeof d, "%s/rows", dir);
    DIR *dp = opendir(d); if (!dp) die("no rows directory");
    struct dirent *de; i64 cap = 0; NROWS = 0; ROWS = 0;
    while ((de = readdir(dp))) {
        int t; unsigned k; if (sscanf(de->d_name, "t%d_k%u.bin", &t, &k) != 2) continue;
        if (typesel && !typesel[t]) continue;
        char fn[1400]; snprintf(fn, sizeof fn, "%s/%s", d, de->d_name);
        FILE *f = fopen(fn, "rb"); fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
        i64 n = sz / sizeof(row_t);
        if (NROWS + n > cap) { cap = (NROWS + n) * 2; ROWS = realloc(ROWS, cap * sizeof(row_t)); }
        if (fread(ROWS + NROWS, sizeof(row_t), n, f) != (size_t)n) die("short read");
        NROWS += n; fclose(f);
    }
    closedir(dp);
}
static int cmp_row(const void *a, const void *b) { return memcmp(a, b, sizeof(row_t)); }

static u64 hsup(u32 xs, u32 ys, int tag) { u64 x = ((u64)xs << 32 | ys) ^ ((u64)tag * 0x9E3779B97F4A7C15ull); x ^= x >> 33; x *= 0xFF51AFD7ED558CCDull; x ^= x >> 33; x *= 0xC4CEB9FE1A85EC53ull; x ^= x >> 33; return x; }

typedef struct { u64 hk; u32 xs, ys; u32 row; } he_t;
static int cmp_he(const void *a, const void *b) { u64 x = ((const he_t *)a)->hk, y = ((const he_t *)b)->hk; return (x > y) - (x < y); }

static int cmd_join(int argc, char **argv) {
    /* pgp join w h OUTDIR [--swap-reduce] [--types ...]  -> prints solutions: 4 hex pairs per line */
    if (argc < 5) die("usage: pgp join w h OUTDIR [options]");
    int w = atoi(argv[2]), h = atoi(argv[3]); snprintf(OUTDIR, sizeof OUTDIR, "%s", argv[4]);
    init(w, h);
    int swap_reduce = argflag(argc, argv, "--swap-reduce");
    build_types(swap_reduce);
    int *sel = calloc(NTYPES + 1, sizeof(int)); select_types(argc, argv, sel);
    load_rows(OUTDIR, sel);
    qsort(ROWS, NROWS, sizeof(row_t), cmp_row);
    /* dedupe identical rows */
    i64 u = 0; for (i64 i = 0; i < NROWS; i++) if (!u || memcmp(&ROWS[i], &ROWS[u - 1], sizeof(row_t))) ROWS[u++] = ROWS[i];
    NROWS = u;
    fprintf(stderr, "[join] %lld rows loaded\n", (long long)NROWS);
    /* hash: (type-tag, xsupp, ysupp) -> list of rows, sorted array of (hashkey,row) */
    he_t *H = malloc(sizeof(he_t) * (NROWS + 1));
    for (i64 i = 0; i < NROWS; i++) { row_t *r = &ROWS[i]; u32 xs = r->xp | r->xn, ys = r->yp | r->yn; H[i].hk = hsup(xs, ys, 0); H[i].xs = xs; H[i].ys = ys; H[i].row = (u32)i; }
    qsort(H, NROWS, sizeof(he_t), cmp_he);
    i64 nsol = 0;
    /* r1 over all rows of selected types (all rows if not swap-reduced; canonical types if swap-reduced);
       r2 over stored rows, rotated (all m rotations) / reversed (4 variants) and, if swap-reduced, also swapped */
    for (i64 i = 0; i < NROWS; i++) {
        row_t *r1 = &ROWS[i]; type_t *t1 = &TYPES[r1->type];
        u32 xs = r1->xp | r1->xn, ys = r1->yp | r1->yn;
        u32 tc = (~xs) & FULL, td = (~ys) & FULL;           /* target supports of C and D */
        int c_inf = t1->exi == 0, d_inf = t1->eyi == 0;
        for (int order = 0; order < (swap_reduce ? 2 : 1); order++) {
            /* order 0: stored row (C, D); order 1: stored row (D, C) = swap */
            u32 q1 = order ? td : tc, q2 = order ? tc : td;
            for (int rot = 0; rot < M; rot++) for (int dd = 0; dd < 2; dd++) {
                const int ra = 0, rb = dd;
                /* a stored row g(r2') = target  <=>  r2' = g^{-1}(target); g = rotate-by-rot then reversals.
                   We look for stored supports S' with  rotr(rev(S'), rot)... simply test all transforms of the target. */
                u32 s1 = rotr_m(q1, rot + (ra ? H2 : 0)), s2 = rotr_m(q2, rot + (rb ? H2 : 0));
                u64 hk = hsup(s1, s2, 0);
                i64 lo = 0, hi = NROWS;
                while (lo < hi) { i64 mid = (lo + hi) / 2; if (H[mid].hk < hk) lo = mid + 1; else hi = mid; }
                for (i64 j = lo; j < NROWS && H[j].hk == hk; j++) {
                    if (H[j].xs != s1 || H[j].ys != s2) continue;
                    row_t *r2 = &ROWS[H[j].row]; type_t *t2 = &TYPES[r2->type];
                    /* infinity compatibility: (C_inf, D_inf) of the matched row */
                    int ci = order ? (t2->eyi != 0) : (t2->exi != 0), di = order ? (t2->exi != 0) : (t2->eyi != 0);
                    if (ci != c_inf || di != d_inf) continue;
                    /* rebuild C and D (sign masks) by undoing the transform: stored (X', Y'); C = rot^{-1} of ... */
                    u32 cx_p, cx_n, dx_p, dx_n;                     /* stored first, second comps */
                    u32 f1p = r2->xp, f1n = r2->xn, f2p = r2->yp, f2n = r2->yn;
                    /* g inverse of the lookup transform: the stored comps were found as rotr(q, rot + H2*r) so the true comp
                       is rotr(stored, -(rot + H2*r)) */
                    u32 a1p = rotr_m(f1p, -(rot + (ra ? H2 : 0))), a1n = rotr_m(f1n, -(rot + (ra ? H2 : 0)));
                    u32 a2p = rotr_m(f2p, -(rot + (rb ? H2 : 0))), a2n = rotr_m(f2n, -(rot + (rb ? H2 : 0)));
                    (void)cx_p; (void)cx_n; (void)dx_p; (void)dx_n;
                    u32 Cp, Cn, Dp, Dn;
                    if (order == 0) { Cp = a1p; Cn = a1n; Dp = a2p; Dn = a2n; } else { Dp = a1p; Dn = a1n; Cp = a2p; Cn = a2n; }
                    /* but the signs of the infinity point of C and D are fixed by the type of the matched row */
                    int cinfv = order ? t2->eyi : t2->exi, dinfv = order ? t2->exi : t2->eyi;
                    printf("%u %u %u %u %u %u %u %u %d %d %d %d %d %d\n", r1->xp, r1->xn, r1->yp, r1->yn, Cp, Cn, Dp, Dn, t1->exi, t1->eyi, cinfv, dinfv, (int)r1->type, (int)r2->type);
                    nsol++;
                }
            }
        }
    }
    fprintf(stderr, "[join] %lld (row,row) matches printed\n", (long long)nsol);
    return 0;
}


/* ------------------------------------------------------------------ plan: exact leaf counts per type without enumerating */
static u32 (*CNT[2])[1 << 16];          /* CNT[e][alpha][key] = number of supports of size alpha with that key */
static void count_keys(int e, u64 alphamask) {
    int hb = H2; u32 nlo = 1u << hb;
    u16 *Bq_lo = malloc(nlo * 2), *Bq_hi = malloc(nlo * 2);
    for (u32 s = 0; s < nlo; s++) {
        u32 klo = 0, khi = 0;
        for (int i = 0; i < hb; i++) if ((s >> i) & 1) {
            klo ^= Dm[i]; khi ^= Dm[i + hb];
            for (int j = 0; j < i; j++) if ((s >> j) & 1) { klo ^= Cm[i][j]; khi ^= Cm[i + hb][j + hb]; }
        }
        Bq_lo[s] = klo; Bq_hi[s] = khi;
    }
    u8 *pc = malloc(nlo); for (u32 s = 0; s < nlo; s++) pc[s] = __builtin_popcount(s);
    if (!CNT[e]) CNT[e] = calloc(33, sizeof(u32) * (1 << 16));
    u64 lo_by_pc_mask = alphamask;
    for (u32 shi = 0; shi < nlo; shi++) {
        u32 R[16]; memset(R, 0, sizeof R);
        for (int i = 0; i < hb; i++) { u32 r = 0; for (int j = 0; j < hb; j++) if ((shi >> j) & 1) r ^= Cm[i][j + hb]; R[i] = r; }
        u32 base = Bq_hi[shi]; int pch = pc[shi]; u32 L = 0, glo = 0;
        for (u32 it = 0; it < nlo; it++) {
            if (it) { u32 bit = __builtin_ctz(it); glo ^= 1u << bit; L ^= R[bit]; }
            int a = pch + pc[glo];
            if (!((lo_by_pc_mask >> a) & 1ull)) continue;
            u32 S = glo | (shi << hb);
            u32 key = base ^ Bq_lo[glo] ^ L; if (e) key ^= fold(S);
            if (REVCANON && S > rotr_m(S, H2)) continue;
            CNT[e][a][key]++;
        }
    }
    free(Bq_lo); free(Bq_hi); free(pc);
}

static int cmd_plan(int argc, char **argv) {
    /* pgp plan w h [--swap-reduce] [--all-keys] : exact numbers of candidate sequences ("leaves") that stage 1 enumerates */
    init(atoi(argv[2]), atoi(argv[3]));
    int allk = argflag(argc, argv, "--all-keys"); REVCANON = !argflag(argc, argv, "--no-rev");
    build_types(argflag(argc, argv, "--swap-reduce"));
    u64 am[2] = {0, 0};
    for (int i = 0; i < NTYPES; i++) { am[TYPES[i].exi != 0] |= 1ull << TYPES[i].ax; am[TYPES[i].eyi != 0] |= 1ull << TYPES[i].ay; }
    double t0 = now();
    for (int e = 0; e < 2; e++) if (am[e]) { count_keys(e, am[e]); fprintf(stderr, "[plan] key counts e=%d done (%.1fs)\n", e, now() - t0); }
    double tot = 0, totstore = 0, totprobe = 0, maxb = 0, maxstore = 0, nb = 0, supx_all = 0, supx_canon = 0;
    printf("type  sums  Xinf sX aX | Yinf sY aY | buckets   leaves(X+Y, both nonempty)   store-side leaves   canonical-support fraction\n");
    for (int i = 0; i < NTYPES; i++) {
        type_t *t = &TYPES[i]; int ex = t->exi != 0, ey = t->eyi != 0;
        int npx = (t->ax + t->sx) / 2, npy = (t->ay + t->sy) / 2;
        double lt = 0, ls = 0, lp = 0, nbk = 0, sall = 0, scan = 0;
        for (u32 key = 0; key < (1u << H2); key++) {
            double nx = CNT[ex][t->ax][key], ny = CNT[ey][t->ay][key];
            sall += nx + ny; if (allk || canon_tab[key]) { scan += nx + ny; } else continue;
            if (nx == 0 || ny == 0) continue;
            double lx = nx * binom[t->ax][npx], ly = ny * binom[t->ay][npy];
            lt += lx + ly; ls += lx < ly ? lx : ly; lp += lx < ly ? ly : lx; nbk++;
            if (lx + ly > maxb) maxb = lx + ly;
            if ((lx < ly ? lx : ly) > maxstore) maxstore = lx < ly ? lx : ly;
        }
        printf("%3d  (%d,%d)  %2d %2d %2d | %2d %2d %2d | %8.0f  %.4e  %.4e  %.4f\n", i, t->x, t->y, t->exi, t->sx, t->ax, t->eyi, t->sy, t->ay, nbk, lt, ls, scan / (sall > 0 ? sall : 1));
        tot += lt; totstore += ls; totprobe += lp; nb += nbk; supx_all += sall; supx_canon += scan;
    }
    printf("TOTAL leaves %.6e   (store side %.6e, probe side %.6e)   buckets %.0f   max bucket leaves %.4e   max store-side leaves per bucket %.4e\n", tot, totstore, totprobe, nb, maxb, maxstore);
    printf("(leaves are the sign patterns enumerated before the PSD filter; 'canonical' = key minimal among its %d rotations)\n", H2);
    return 0;
}

/* ------------------------------------------------------------------ small commands */
static int cmd_info(int argc, char **argv) {
    init(atoi(argv[2]), atoi(argv[3]));
    build_types(argflag(argc, argv, "--swap-reduce"));
    printf("w=%d h=%d m=%d h2=%d g0=%d  (|Gamma| = %d)\n", P, HH, M, H2, G0, 2 * P);
    print_types();
    return 0;
}

static int cmd_selftest(int argc, char **argv) {
    init(atoi(argv[2]), atoi(argv[3]));
    srand(12345); REVCANON = 0;
    int bad = 0;
    /* (1) class counts: sum_{i,j} N_k(i,j) = w - 2 */
    for (int k = 0; k < H2; k++) { int s = 0; for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) s += Ncnt[k][i][j]; if (s != P - 2) { printf("N sum bad k=%d: %d\n", k, s); bad++; } }
    /* (2) PAF incremental vs direct, parity key tab vs direct, key rotation covariance, Fourier vs direct DFT */
    for (int it = 0; it < 300; it++) {
        int xo[32], xinf = (rand() % 3) - 1; for (int j = 0; j < M; j++) xo[j] = (rand() % 3) - 1;
        int v[16]; paf_classes_direct(xinf, xo, v);
        /* via tables */
        int vt[16];
        for (int k = 0; k < H2; k++) {
            int s = 0;
            for (int i = 0; i < M; i++) { s += Dv[i][k] * xo[i] * xo[i]; for (int j = i + 1; j < M; j++) s += Wv[i][j][k] * xo[i] * xo[j]; }
            s += xinf * (xo[k] + xo[k + H2]); vt[k] = s;
            if (s != v[k]) { bad++; if (bad < 10) printf("PAF table mismatch it=%d k=%d %d vs %d\n", it, k, s, v[k]); }
        }
        /* Fourier */
        for (int k = 0; k < H2; k++) {
            double re = xinf, im = 0;
            long gk = 1; for (int i = 0; i < k; i++) gk = gk * G0 % P;
            for (int x = 1; x < P; x++) { double ang = 2 * M_PI * (double)(gk * x % P) / P; re += xo[cls[x]] * cos(ang); im += xo[cls[x]] * sin(ang); }
            double re2 = xinf, im2 = 0; for (int j = 0; j < M; j++) { re2 += xo[j] * ER[j][k]; im2 += xo[j] * EI[j][k]; }
            if (fabs(re - re2) > 1e-3 || fabs(im - im2) > 1e-3) { bad++; if (bad < 10) printf("Fourier mismatch k=%d\n", k); }
        }
        /* keys */
        u32 S = 0; for (int j = 0; j < M; j++) if (xo[j]) S |= 1u << j;
        int e = xinf != 0;
        u32 kd = key_direct(e, S), kt = key_tab(e, S);
        if (kd != kt) { bad++; if (bad < 10) printf("key mismatch %u %u\n", kd, kt); }
        for (int s = 0; s < M; s++) { u32 S2 = rotr_m(S, s); if (key_direct(e, S2) != kshift(kd, s)) { bad++; if (bad < 10) printf("key rotation covariance fails s=%d\n", s); } }
        /* rotation covariance of PAF: rotate xo */
        for (int s = 0; s < 3; s++) {
            int xr[32]; for (int j = 0; j < M; j++) xr[j] = xo[(j + s) % M]; int vr[16]; paf_classes_direct(xinf, xr, vr);
            for (int k = 0; k < H2; k++) if (vr[k] != v[(k + s) % H2]) { bad++; if (bad < 10) printf("PAF rotation covariance fails\n"); }
        }
    }
    /* (3) index vs direct keys */
    ALLKEYS = 1;
    if (M <= 16) {
        u64 am = 0; for (int a = 0; a <= M; a++) am |= 1ull << a;
        for (int e = 0; e < 2; e++) {
            build_index(e, am); i64 tot = 0;
            for (int a = 0; a <= M; a++) { supidx_t *ix = &IDX[e][a]; tot += ix->n;
                for (u32 key = 0; key < (1u << H2); key++) for (u32 i = ix->start[key]; i < ix->start[key + 1]; i++) {
                    u32 S = ix->mask[i]; if (__builtin_popcount(S) != a || key_direct(e, S) != key) { bad++; if (bad < 10) printf("index key bad\n"); }
                } }
            if (tot != (1ll << M)) { bad++; printf("index size bad %lld\n", (long long)tot); }
        }
    }
    printf("selftest w=%d h=%d: %s (%d problems)\n", P, HH, bad ? "FAIL" : "ok", bad);
    return bad != 0;
}

int main(int argc, char **argv) {
    if (argc < 2) die("commands: info selftest rows join");
    if (!strcmp(argv[1], "info")) return cmd_info(argc, argv);
    if (!strcmp(argv[1], "selftest")) return cmd_selftest(argc, argv);
    if (!strcmp(argv[1], "rows")) return cmd_rows(argc, argv);
    if (!strcmp(argv[1], "join")) return cmd_join(argc, argv);
    if (!strcmp(argv[1], "plan")) return cmd_plan(argc, argv);
    die("unknown command");
    return 1;
}
