/*
 * p3.c -- three-list join for propus quadruples (A, B, B, D) on Z_v with a multiplier group H.
 *
 *   PAF_A + 2 PAF_B + PAF_D = 4 v delta,  A constant on <H,-1>-orbits (symmetric), B and D constant on H-orbits.
 *   Class form: for each of the C classes c of shifts,  F_A(c) + 2 F_B(c) + F_D(c) = TGT(c)   (TGT = 0 for a search).
 *
 * Lists
 *   A  all sign vectors on the O_A orbits with row sum rA (optionally PSD-filtered),
 *   B  all sign vectors on the O_B orbits with row sum rB, filtered by (a) canonical form under the multiplier
 *      group Z_v^* modulo H (lossless: a solution maps to one with canonical B), (b) the power-spectrum bound
 *      |B^(w)|^2 <= 2v (necessary: P_A + 2 P_B + P_D = 4v), (c) optionally a low-energy cut,
 *   D  all sign vectors with row sum rD, filtered by canonical form under the reversal x(g) -> x(-g) of D alone (lossless:
 *      PAF is reversal invariant), |D^(w)|^2 <= 4v and optionally a low-energy cut.
 * Join
 *   Every list member carries its class-value vector F (C int16) and the 64-bit linear hash h = sum_c R_c F_c.  A solution
 *   needs h_D = H(TGT) - h_A - 2 h_B and, on the NS "split" classes S, F_D(S) = TGT(S) - F_A(S) - 2 F_B(S).  D is grouped
 *   by its F on S (dense key kappa; F = v mod 4 always, so (F - v mod 4)/4 is an exact lattice coordinate).  For each
 *   D group (a "sigma") a small open-addressing table of the D hashes is built; every a in A looks up the B group that
 *   matches (group of 2F_B(S) = TGT - sigma - F_A(S)) and every member of that B group probes the table.  Every probe
 *   hit is re-verified exactly on all C classes from the three orbit masks.  The join is exact: no solution is lost
 *   by grouping, hashing (false positives are re-checked) or the lossless list filters.
 *
 * Work is split over sigma groups (threads, --shard i/n for several machines, checkpoint file, resume).
 *
 *   cc -O3 -march=native -ffast-math -pthread p3.c -o p3 -lm
 *   ./p3 PROBLEM [options]        (see usage())
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <time.h>
#include <unistd.h>

#define MAXO 64
#define MAXC 32
typedef unsigned long long u64;

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void die(const char *m) { fprintf(stderr, "p3: %s\n", m); exit(2); }

/* ---------------------------------------------------------------------------------------------------------------- */
static int V, C, CP, R0, OFF, W;           /* R0 = v mod 4; key window: q = (F - R0)/4 + OFF in [0, W) */
static int NS, SC[4];                      /* split classes */
static u64 STRIDE[4], NKAPPA;
static u64 RH[MAXC], HTGT;
static int TGT[MAXC];
static float THP[MAXC][MAXC] __attribute__((aligned(64)));   /* PSD matrix, padded rows */

typedef struct {
    int O, r;
    int n[MAXO];
    int tail[MAXO + 1];                    /* sum of sizes of orbits k.. */
    int16_t c0[MAXC];
    int16_t *Tv;                           /* [O][O][CP]: T[c][i][j] as a vector over c */
    uint8_t *reach;                        /* [(O+1)][2V+1]: can orbits k.. produce sum s */
    int psd; float bound;
    int energy; double tau; double m[MAXC]; double Wm[MAXC * MAXC];
} Blk;
static Blk BA, BB, BD;
static int NPERM, NBYTES;                  /* canonical form of B under the multiplier group */
static u64 (*PTB)[8][256];                 /* PTB[t][byte][value] */
static u64 PTD[8][256];                    /* reversal on D */
static int CANONB = 0, CANOND = 0;

static int16_t *alloc16(size_t n) { size_t b = ((n * 2 + 63) / 64) * 64; int16_t *p = aligned_alloc(64, b ? b : 64); if (!p) die("oom"); memset(p, 0, b ? b : 64); return p; }

/* ---------------------------------------------------------------------------------------------------------------- */
static FILE *rf;
static int rdi(void) { int x; if (fscanf(rf, "%d", &x) != 1) die("problem file: int expected"); return x; }
static double rdd(void) { double x; if (fscanf(rf, "%lf", &x) != 1) die("problem file: double expected"); return x; }
static u64 rdu(void) { u64 x; if (fscanf(rf, "%llu", &x) != 1) die("problem file: u64 expected"); return x; }

static void make_reach(Blk *B) {
    int O = B->O; B->tail[O] = 0;
    for (int k = O - 1; k >= 0; k--) B->tail[k] = B->tail[k + 1] + B->n[k];
    size_t w = 2 * (size_t)V + 1;
    B->reach = calloc((O + 1) * w, 1);
    B->reach[(size_t)O * w + V] = 1;
    for (int k = O - 1; k >= 0; k--)
        for (int s = -V; s <= V; s++) {
            int a = s - B->n[k], b = s + B->n[k];
            int ok = 0;
            if (a >= -V && a <= V && B->reach[(size_t)(k + 1) * w + a + V]) ok = 1;
            if (b >= -V && b <= V && B->reach[(size_t)(k + 1) * w + b + V]) ok = 1;
            B->reach[(size_t)k * w + s + V] = ok;
        }
}

static void read_block(Blk *B, int O, int r) {
    B->O = O; B->r = r;
    for (int i = 0; i < O; i++) B->n[i] = rdi();
}
static void read_const_T(Blk *B, int shareT, const Blk *src) {
    int O = B->O;
    for (int c = 0; c < C; c++) B->c0[c] = (int16_t)rdi();
    if (shareT) { B->Tv = src->Tv; return; }
    B->Tv = alloc16((size_t)O * O * CP);
    for (int c = 0; c < C; c++) for (int i = 0; i < O; i++) for (int j = 0; j < O; j++) B->Tv[((size_t)i * O + j) * CP + c] = (int16_t)rdi();
}

static void read_problem(const char *path) {
    rf = fopen(path, "r"); if (!rf) die("cannot open problem file");
    char tag[32]; if (fscanf(rf, "%31s", tag) != 1 || strcmp(tag, "P3PROB")) die("bad problem file");
    if (rdi() != 1) die("problem version");
    V = rdi(); C = rdi(); int OA = rdi(), OB = rdi();
    if (C > MAXC) die("C > 32"); CP = ((C + 15) / 16) * 16; R0 = V & 3;
    int rA = rdi(), rB = rdi(), rD = rdi();
    for (int c = 0; c < C; c++) TGT[c] = rdi();
    read_block(&BA, OA, rA); read_block(&BB, OB, rB);
    BD.O = OB; BD.r = rD; memcpy(BD.n, BB.n, sizeof BB.n);
    read_const_T(&BA, 0, NULL);
    read_const_T(&BB, 0, NULL);
    BD.Tv = BB.Tv; memcpy(BD.c0, BB.c0, sizeof BB.c0);
    for (int cp = 0; cp < C; cp++) for (int c = 0; c < C; c++) THP[cp][c] = (float)rdd();
    BA.bound = (float)rdd(); BB.bound = (float)rdd(); BD.bound = (float)rdd();
    BA.psd = BA.bound > 0; BB.psd = BB.bound > 0; BD.psd = BD.bound > 0;
    CANONB = rdi(); CANOND = rdi();
    NPERM = rdi();
    NBYTES = (OB + 7) / 8;
    PTB = calloc(NPERM ? NPERM : 1, sizeof *PTB);
    for (int t = 0; t < NPERM; t++) {
        int pi[MAXO]; for (int i = 0; i < OB; i++) pi[i] = rdi();
        for (int b = 0; b < NBYTES; b++) for (int x = 0; x < 256; x++) {
            u64 im = 0;
            for (int bit = 0; bit < 8; bit++) {          /* source bit (8b+bit) set => image bits i with pi[i] = 8b+bit */
                if (!((x >> bit) & 1)) continue;
                int src = 8 * b + bit; if (src >= OB) continue;
                for (int i = 0; i < OB; i++) if (pi[i] == src) im |= 1ULL << i;
            }
            PTB[t][b][x] = im;
        }
    }
    { int pi[MAXO]; for (int i = 0; i < OB; i++) pi[i] = rdi();
      for (int b = 0; b < NBYTES; b++) for (int x = 0; x < 256; x++) {
          u64 im = 0;
          for (int bit = 0; bit < 8; bit++) {
              if (!((x >> bit) & 1)) continue;
              int src = 8 * b + bit; if (src >= OB) continue;
              for (int i = 0; i < OB; i++) if (pi[i] == src) im |= 1ULL << i;
          }
          PTD[b][x] = im;
      } }
    for (int c = 0; c < C; c++) RH[c] = rdu();
    for (int k = 0; k < 2; k++) {
        Blk *B = k == 0 ? &BB : &BD;
        B->energy = rdi(); B->tau = rdd();
        for (int c = 0; c < C; c++) B->m[c] = rdd();
        for (int i = 0; i < C * C; i++) B->Wm[i] = rdd();
    }
    fclose(rf);
    make_reach(&BA); make_reach(&BB); BD.reach = BB.reach;
    for (int i = 0; i <= BB.O; i++) BD.tail[i] = BB.tail[i];
    HTGT = 0; for (int c = 0; c < C; c++) HTGT += RH[c] * (u64)(int64_t)TGT[c];
    OFF = (V + 3) / 4 + 1; W = OFF + V / 4 + 2;
}

/* ---------------------------------------------------------------------------------------------------------------- */
static inline u64 hash_of(const int16_t *val) { u64 h = 0; for (int c = 0; c < C; c++) h += RH[c] * (u64)(int64_t)val[c]; return h; }
static inline u64 kappa_of(const int16_t *val) {
    u64 k = 0;
    for (int i = 0; i < NS; i++) k += (u64)(((val[SC[i]] - R0) >> 2) + OFF) * STRIDE[i];
    return k;
}
static inline int psd_ok(const int16_t *val, float bound) {
    float f[MAXC] __attribute__((aligned(64)));
    for (int c = 0; c < CP; c++) f[c] = (float)val[c];
    const float lim = bound + 0.5f - (float)V;
    for (int cp = 0; cp < C; cp++) {
        float s = 0; const float *t = THP[cp];
        for (int c = 0; c < C; c++) s += t[c] * f[c];
        if (s > lim) return 0;
    }
    return 1;
}
static inline int energy_ok(const Blk *B, const int16_t *val) {
    double d[MAXC], e = 0;
    for (int c = 0; c < C; c++) d[c] = (double)val[c] - B->m[c];
    for (int i = 0; i < C; i++) {
        double u = 0; const double *w = B->Wm + (size_t)i * C;
        for (int c = 0; c < C; c++) u += w[c] * d[c];
        e += u * u; if (e > B->tau) return 0;
    }
    return 1;
}
static inline int canonB_ok(u64 mask) {
    for (int t = 0; t < NPERM; t++) {
        u64 im = 0; for (int b = 0; b < NBYTES; b++) im |= PTB[t][b][(mask >> (8 * b)) & 255];
        if (im < mask) return 0;
    }
    return 1;
}
static inline int canonD_ok(u64 mask) {
    u64 im = 0; for (int b = 0; b < NBYTES; b++) im |= PTD[b][(mask >> (8 * b)) & 255];
    return im >= mask;
}

/* class values from an orbit mask (bit i = 1 <=> y_i = +1) */
static void vals_from_mask(const Blk *B, u64 mask, int16_t *val) {
    int O = B->O;
    int16_t u[MAXO][MAXC] __attribute__((aligned(64)));
    memset(u, 0, sizeof u);
    for (int c = 0; c < CP; c++) val[c] = c < C ? B->c0[c] : 0;
    for (int k = 0; k < O; k++) {
        int s = ((mask >> k) & 1) ? 1 : -1;
        for (int c = 0; c < CP; c++) val[c] += s * u[k][c];
        for (int j = k + 1; j < O; j++) { const int16_t *t = B->Tv + ((size_t)k * O + j) * CP; for (int c = 0; c < CP; c++) u[j][c] += s * t[c]; }
    }
}

/* ---------------------------------------------------------------------------------------------------------------- */
typedef struct { u64 n; u64 *hash, *mask; u64 *start; int16_t *fs; } List;   /* start[kappa]; fs = F on S (A only) */

enum { KA = 0, KB = 1, KD = 2 };
static int SHARD_I = 0, SHARD_N = 1, NTHREADS = 1;
static int SLICE_S = 0, SLICE_N = 1, PASS1ONLY = 0;   /* benchmarking: every SLICE_N-th DFS task of B and D only; counts only */
static u64 *OVR_MASK[3]; static u64 OVR_N[3];     /* override lists */

typedef struct {
    const Blk *B; int kind, pass, tid;
    uint32_t *cnt; u64 *pos; List *L;
    int16_t *U, *val;
    u64 leaves, kept;
} Ctx;

static inline void leaf(Ctx *cx, u64 mask, const int16_t *val) {
    cx->leaves++;
    const Blk *B = cx->B;
    if (cx->kind == KB && CANONB && !canonB_ok(mask)) return;
    if (cx->kind == KD && CANOND && !canonD_ok(mask)) return;
    if (B->psd && !psd_ok(val, B->bound)) return;
    if (B->energy && !energy_ok(B, val)) return;
    u64 k = kappa_of(val);
    if (cx->kind == KD && SHARD_N > 1 && (k % SHARD_N) != (u64)SHARD_I) return;
    cx->kept++;
    if (cx->pass == 1) cx->cnt[k]++;
    else { u64 p = cx->pos[k]++; cx->L->hash[p] = hash_of(val); cx->L->mask[p] = mask; }
}

static void dfs(Ctx *cx, int k, int rem, u64 mask) {
    const Blk *B = cx->B; int O = B->O;
    const int16_t *val = cx->val + (size_t)k * CP;
    const int16_t *U = cx->U + (size_t)k * O * CP;
    int16_t *val2 = cx->val + (size_t)(k + 1) * CP;
    size_t w = 2 * (size_t)V + 1;
    for (int sidx = 0; sidx < 2; sidx++) {
        int s = sidx ? -1 : 1;
        int rem2 = rem - s * B->n[k];
        if (rem2 > B->tail[k + 1] || rem2 < -B->tail[k + 1]) continue;
        if (!B->reach[(size_t)(k + 1) * w + rem2 + V]) continue;
        const int16_t *uk = U + (size_t)k * CP;
        if (s > 0) for (int c = 0; c < CP; c++) val2[c] = val[c] + uk[c];
        else       for (int c = 0; c < CP; c++) val2[c] = val[c] - uk[c];
        u64 m2 = mask | ((u64)(s > 0) << k);
        if (k == O - 1) { leaf(cx, m2, val2); continue; }
        int16_t *U2 = cx->U + (size_t)(k + 1) * O * CP;
        for (int j = k + 1; j < O; j++) {
            const int16_t *t = B->Tv + ((size_t)k * O + j) * CP; const int16_t *a = U + (size_t)j * CP; int16_t *o = U2 + (size_t)j * CP;
            if (s > 0) for (int c = 0; c < CP; c++) o[c] = a[c] + t[c];
            else       for (int c = 0; c < CP; c++) o[c] = a[c] - t[c];
        }
        dfs(cx, k + 1, rem2, m2);
    }
}

typedef struct { u64 mask; int rem, depth; } Task;
static void collect(const Blk *B, int k, int rem, u64 mask, int depth, Task *T, int *nt, int cap) {
    if (k == depth) { if (*nt < cap) { T[*nt].mask = mask; T[*nt].rem = rem; T[*nt].depth = depth; } (*nt)++; return; }
    size_t w = 2 * (size_t)V + 1;
    for (int sidx = 0; sidx < 2; sidx++) {
        int s = sidx ? -1 : 1; int rem2 = rem - s * B->n[k];
        if (rem2 > B->tail[k + 1] || rem2 < -B->tail[k + 1]) continue;
        if (!B->reach[(size_t)(k + 1) * w + rem2 + V]) continue;
        collect(B, k + 1, rem2, mask | ((u64)(s > 0) << k), depth, T, nt, cap);
    }
}

static Task *TASKS; static int NTASKS;
static void make_tasks(const Blk *B, int want) {
    int d = 0, nt = 0;
    for (d = 1; d < B->O - 1; d++) { nt = 0; collect(B, 0, B->r, 0, d, NULL, &nt, 0); if (nt >= want) break; }
    if (d >= B->O) d = B->O - 1;
    if (d < 1) d = 1;
    nt = 0; collect(B, 0, B->r, 0, d, NULL, &nt, 0);
    TASKS = malloc(sizeof(Task) * (nt ? nt : 1)); int n2 = 0;
    collect(B, 0, B->r, 0, d, TASKS, &n2, nt);
    NTASKS = n2;
}

static void run_task(Ctx *cx, const Task *t) {
    const Blk *B = cx->B; int O = B->O;
    memset(cx->U, 0, sizeof(int16_t) * (size_t)O * CP);
    for (int c = 0; c < CP; c++) cx->val[c] = c < C ? B->c0[c] : 0;
    for (int k = 0; k < t->depth; k++) {      /* replay the prefix */
        int s = ((t->mask >> k) & 1) ? 1 : -1;
        const int16_t *val = cx->val + (size_t)k * CP; int16_t *val2 = cx->val + (size_t)(k + 1) * CP;
        const int16_t *U = cx->U + (size_t)k * O * CP; int16_t *U2 = cx->U + (size_t)(k + 1) * O * CP;
        const int16_t *uk = U + (size_t)k * CP;
        for (int c = 0; c < CP; c++) val2[c] = val[c] + s * uk[c];
        for (int j = k + 1; j < O; j++) { const int16_t *tt = B->Tv + ((size_t)k * O + j) * CP; for (int c = 0; c < CP; c++) U2[(size_t)j * CP + c] = U[(size_t)j * CP + c] + s * tt[c]; }
    }
    dfs(cx, t->depth, t->rem, t->mask & ((1ULL << t->depth) - 1));
}

typedef struct { Ctx cx; int nth; } Arg;
static void *worker_list(void *p) {
    Arg *a = p; Ctx *cx = &a->cx; const Blk *B = cx->B;
    cx->U = alloc16((size_t)(B->O + 1) * B->O * CP); cx->val = alloc16((size_t)(B->O + 2) * CP);
    if (OVR_MASK[cx->kind]) {                               /* override list: tasks are chunks of 4096 masks */
        u64 n = OVR_N[cx->kind]; int nch = (int)((n + 4095) / 4096);
        for (int ch = cx->tid; ch < nch; ch += a->nth) {
            u64 lo = (u64)ch * 4096, hi = lo + 4096 < n ? lo + 4096 : n;
            for (u64 i = lo; i < hi; i++) { vals_from_mask(B, OVR_MASK[cx->kind][i], cx->val); leaf(cx, OVR_MASK[cx->kind][i], cx->val); }
        }
    } else {
        for (int i = cx->tid; i < NTASKS; i += a->nth) { if (SLICE_N > 1 && cx->kind != KA && (i % SLICE_N) != SLICE_S) continue; run_task(cx, &TASKS[i]); }
    }
    free(cx->U); free(cx->val);
    return NULL;
}

static void build_list(const Blk *B, int kind, List *L, const char *name) {
    double t0 = now();
    if (!OVR_MASK[kind]) make_tasks(B, 64 * NTHREADS);
    else NTASKS = (int)((OVR_N[kind] + 4095) / 4096);
    int T = NTHREADS; if (T > NTASKS) T = NTASKS > 0 ? NTASKS : 1;
    Arg *args = calloc(T, sizeof(Arg)); pthread_t *th = calloc(T, sizeof(pthread_t));
    for (int pass = 1; pass <= 2; pass++) {
        if (pass == 2) {
            u64 tot = 0; L->start = calloc(NKAPPA + 1, sizeof(u64));
            for (u64 k = 0; k < NKAPPA; k++) { u64 s = 0; for (int t = 0; t < T; t++) s += args[t].cx.cnt[k]; L->start[k] = tot; tot += s; }
            L->start[NKAPPA] = tot; L->n = tot;
            L->hash = malloc(8 * (tot ? tot : 1)); L->mask = malloc(8 * (tot ? tot : 1));
            if (!L->hash || !L->mask) die("oom (list)");
            for (int t = 0; t < T; t++) {
                args[t].cx.pos = malloc(8 * NKAPPA); if (!args[t].cx.pos) die("oom (pos)");
            }
            for (u64 k = 0; k < NKAPPA; k++) { u64 off = L->start[k]; for (int t = 0; t < T; t++) { args[t].cx.pos[k] = off; off += args[t].cx.cnt[k]; } }
            for (int t = 0; t < T; t++) { free(args[t].cx.cnt); args[t].cx.cnt = NULL; }
        }
        for (int t = 0; t < T; t++) {
            args[t].nth = T; args[t].cx.B = B; args[t].cx.kind = kind; args[t].cx.pass = pass; args[t].cx.tid = t; args[t].cx.L = L;
            args[t].cx.leaves = args[t].cx.kept = 0;
            if (pass == 1) { args[t].cx.cnt = calloc(NKAPPA, 4); if (!args[t].cx.cnt) die("oom (cnt)"); }
            pthread_create(&th[t], NULL, worker_list, &args[t]);
        }
        for (int t = 0; t < T; t++) pthread_join(th[t], NULL);
        if (pass == 1) {
            u64 lv = 0, kp = 0; for (int t = 0; t < T; t++) { lv += args[t].cx.leaves; kp += args[t].cx.kept; }
            fprintf(stderr, "  %s: pass 1: %llu leaves, %llu kept [%.1fs]\n", name, lv, kp, now() - t0);
            printf("COUNT %s %llu %llu\n", name, lv, kp); fflush(stdout);
            if (PASS1ONLY) { for (int t = 0; t < T; t++) free(args[t].cx.cnt); free(args); free(th); L->n = 0; L->start = NULL; return; }
        }
    }
    for (int t = 0; t < T; t++) free(args[t].cx.pos);
    free(args); free(th);
    fprintf(stderr, "  %s: list of %llu members [%.1fs]\n", name, L->n, now() - t0);
}

/* ---------------------------------------------------------------------------------------------------------------- */
static List LA, LB, LD;
static int16_t *FAS;                       /* A: class values on S, [nA][NS] */
static u64 HITS = 0, MATCHES = 0, BADM = 0;
static pthread_mutex_t MX = PTHREAD_MUTEX_INITIALIZER;
static FILE *HITF = NULL;

static void report(u64 ma, u64 mb, u64 md) {
    int16_t fa[MAXC], fb[MAXC], fd[MAXC];
    vals_from_mask(&BA, ma, fa); vals_from_mask(&BB, mb, fb); vals_from_mask(&BD, md, fd);
    int ok = 1; for (int c = 0; c < C; c++) if (fa[c] + 2 * fb[c] + fd[c] != TGT[c]) ok = 0;
    pthread_mutex_lock(&MX);
    MATCHES++;
    if (ok) {
        HITS++;
        printf("HIT %llx %llx %llx\n", ma, mb, md); fflush(stdout);
        if (HITF) { fprintf(HITF, "%llx %llx %llx\n", ma, mb, md); fflush(HITF); }
    } else BADM++;
    pthread_mutex_unlock(&MX);
}

typedef struct {
    u64 *glist; u64 ng; int bs; u64 nbatch;
    atomic_ullong next; uint8_t *done; FILE *ck; int countonly;
    atomic_ullong pairs, probes;
    double t0; atomic_ullong ndone;
} JoinG;
static JoinG JG;

static void *worker_join(void *p) {
    (void)p;
    u64 tcap = 0; u64 *tkey = NULL; uint32_t *tidx = NULL;
    u64 pairs = 0;
    for (;;) {
        u64 b = atomic_fetch_add(&JG.next, 1);
        if (b >= JG.nbatch) break;
        if (JG.done[b]) continue;
        u64 g0 = b * JG.bs, g1 = g0 + JG.bs < JG.ng ? g0 + JG.bs : JG.ng;
        for (u64 gi = g0; gi < g1; gi++) {
            u64 kap = JG.glist[gi];
            u64 d0 = LD.start[kap], d1 = LD.start[kap + 1], nd = d1 - d0;
            int sig[4]; { u64 q = kap; for (int i = 0; i < NS; i++) { int qi = (int)(q % W); q /= W; sig[i] = 4 * (qi - OFF) + R0; } }
            u64 size = 64; while (size < 4 * nd) size <<= 1;
            u64 mask = size - 1; int sh = 64 - __builtin_ctzll(size);
            if (!JG.countonly) {
                if (size > tcap) { free(tkey); free(tidx); tkey = malloc(8 * size); tidx = malloc(4 * size); if (!tkey || !tidx) die("oom (table)"); tcap = size; }
                memset(tkey, 0, 8 * size);
                for (u64 i = 0; i < nd; i++) {
                    u64 key = LD.hash[d0 + i]; u64 s = (key * 0x9E3779B97F4A7C15ULL) >> sh;
                    if (key == 0) die("D hash 0: choose another seed");
                    while (tkey[s]) s = (s + 1) & mask;
                    tkey[s] = key; tidx[s] = (uint32_t)i;
                }
            }
            for (u64 a = 0; a < LA.n; a++) {
                u64 kb = 0; int ok = 1;
                for (int i = 0; i < NS; i++) {
                    int t = TGT[SC[i]] - sig[i] - FAS[a * NS + i];       /* = 2 F_B(S_i) */
                    if (t & 1) { ok = 0; break; }
                    int fb = t >> 1; if ((fb - R0) & 3) { ok = 0; break; }
                    int q = ((fb - R0) >> 2) + OFF; if (q < 0 || q >= W) { ok = 0; break; }
                    kb += (u64)q * STRIDE[i];
                }
                if (!ok) continue;
                u64 b0 = LB.start[kb], b1 = LB.start[kb + 1];
                if (b0 == b1) continue;
                pairs += b1 - b0;
                if (JG.countonly) continue;
                u64 base = HTGT - LA.hash[a];
                for (u64 j = b0; j < b1; j++) {
                    u64 key = base - 2 * LB.hash[j];
                    u64 s = (key * 0x9E3779B97F4A7C15ULL) >> sh;
                    u64 tk;
                    while ((tk = tkey[s]) != 0) {
                        if (tk == key) report(LA.mask[a], LB.mask[j], LD.mask[d0 + tidx[s]]);
                        s = (s + 1) & mask;
                    }
                }
            }
        }
        pthread_mutex_lock(&MX);
        JG.done[b] = 1;
        if (JG.ck) { fprintf(JG.ck, "%llu\n", b); fflush(JG.ck); }
        pthread_mutex_unlock(&MX);
        u64 nd = atomic_fetch_add(&JG.ndone, 1) + 1;
        if ((nd & 63) == 0) fprintf(stderr, "  join: %llu/%llu batches [%.0fs]\n", nd, JG.nbatch, now() - JG.t0);
    }
    atomic_fetch_add(&JG.pairs, pairs);
    free(tkey); free(tidx);
    return NULL;
}

static u64 fnv(const char *buf, size_t n, u64 h) { for (size_t i = 0; i < n; i++) { h ^= (unsigned char)buf[i]; h *= 1099511628211ULL; } return h; }

static u64 *load_masks(const char *path, u64 *n) {
    FILE *f = fopen(path, "rb"); if (!f) die("cannot open list file");
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    *n = sz / 8; u64 *m = malloc(8 * (*n ? *n : 1));
    if (fread(m, 8, *n, f) != *n) die("short read"); fclose(f); return m;
}

static void usage(void) {
    fprintf(stderr, "usage: p3 PROBLEM [-t threads] [-ns k] [-shard i/n] [-ckpt file] [-resume] [-hits file] [-count]\n"
                    "          [-alist f] [-blist f] [-dlist f]   (binary uint64 orbit masks; bit i = 1 <=> y_i = +1)\n"
                    "          [-listonly]  (build and report the lists, no join)  [-bs batch] [-S c0,c1,c2]\n"
                    "          [-countlists] (pass 1 only: leaf and kept counts per list)  [-slice s/n] (BENCHMARK ONLY: every n-th DFS task of B, D)\n");
    exit(1);
}

int main(int argc, char **argv) {
    if (argc < 2) usage();
    const char *prob = argv[1], *ckpath = NULL, *hitpath = NULL, *al = NULL, *bl = NULL, *dl = NULL;
    int resume = 0, countonly = 0, listonly = 0; NS = 3; JG.bs = 16;
    int sset = 0;
    NTHREADS = (int)sysconf(_SC_NPROCESSORS_ONLN);
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) NTHREADS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-ns") && i + 1 < argc) NS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-shard") && i + 1 < argc) { if (sscanf(argv[++i], "%d/%d", &SHARD_I, &SHARD_N) != 2) usage(); }
        else if (!strcmp(argv[i], "-ckpt") && i + 1 < argc) ckpath = argv[++i];
        else if (!strcmp(argv[i], "-resume")) resume = 1;
        else if (!strcmp(argv[i], "-hits") && i + 1 < argc) hitpath = argv[++i];
        else if (!strcmp(argv[i], "-count")) countonly = 1;
        else if (!strcmp(argv[i], "-listonly")) listonly = 1;
        else if (!strcmp(argv[i], "-countlists")) { PASS1ONLY = 1; listonly = 1; }
        else if (!strcmp(argv[i], "-slice") && i + 1 < argc) { if (sscanf(argv[++i], "%d/%d", &SLICE_S, &SLICE_N) != 2) usage(); }
        else if (!strcmp(argv[i], "-bs") && i + 1 < argc) JG.bs = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-alist") && i + 1 < argc) al = argv[++i];
        else if (!strcmp(argv[i], "-blist") && i + 1 < argc) bl = argv[++i];
        else if (!strcmp(argv[i], "-dlist") && i + 1 < argc) dl = argv[++i];
        else if (!strcmp(argv[i], "-S") && i + 1 < argc) { char *s = argv[++i]; NS = 0; for (char *tok = strtok(s, ","); tok && NS < 4; tok = strtok(NULL, ",")) SC[NS++] = atoi(tok); sset = 1; }
        else usage();
    }
    if (NTHREADS < 1) NTHREADS = 1;
    read_problem(prob);
    if (NS < 1 || NS > 3) die("-ns must be 1, 2 or 3");
    if (!sset) for (int i = 0; i < NS; i++) SC[i] = i;
    for (int i = 0; i < NS; i++) if (SC[i] < 0 || SC[i] >= C) die("bad split class");
    STRIDE[0] = 1; for (int i = 1; i < 4; i++) STRIDE[i] = STRIDE[i - 1] * W;
    NKAPPA = 1; for (int i = 0; i < NS; i++) NKAPPA *= W;
    if (SHARD_I < 0 || SHARD_I >= SHARD_N) die("bad shard");
    if (al) OVR_MASK[KA] = load_masks(al, &OVR_N[KA]);
    if (bl) OVR_MASK[KB] = load_masks(bl, &OVR_N[KB]);
    if (dl) OVR_MASK[KD] = load_masks(dl, &OVR_N[KD]);
    fprintf(stderr, "p3: v=%d C=%d OA=%d OB=%d rows=(%d,%d,%d) NS=%d W=%d threads=%d shard %d/%d canonB=%d(perms %d) canonD=%d psd=(%d,%d,%d) energy=(%d,%d)\n",
            V, C, BA.O, BB.O, BA.r, BB.r, BD.r, NS, W, NTHREADS, SHARD_I, SHARD_N, CANONB, NPERM, CANOND, BA.psd, BB.psd, BD.psd, BB.energy, BD.energy);
    double t0 = now();
    int sA = SHARD_N; SHARD_N = 1;      /* A and B are never sharded */
    build_list(&BA, KA, &LA, "A");
    build_list(&BB, KB, &LB, "B");
    SHARD_N = sA;
    FAS = malloc(2 * (LA.n ? LA.n : 1) * NS);
    for (u64 a = 0; a < LA.n; a++) { int16_t v[MAXC]; vals_from_mask(&BA, LA.mask[a], v); for (int i = 0; i < NS; i++) FAS[a * NS + i] = v[SC[i]]; }
    /* LA must be a flat list: re-order is irrelevant (only the order of entries matters), kappa grouping of A is unused */
    build_list(&BD, KD, &LD, "D");
    fprintf(stderr, "p3: lists A=%llu B=%llu D=%llu [%.1fs]  (D keys with 0 hash: checked in join)\n", LA.n, LB.n, LD.n, now() - t0);
    printf("LISTS %llu %llu %llu\n", LA.n, LB.n, LD.n); fflush(stdout);
    if (listonly) return 0;

    /* sigma groups: nonempty D groups of this shard */
    u64 ng = 0; for (u64 k = 0; k < NKAPPA; k++) if (LD.start[k + 1] > LD.start[k]) ng++;
    JG.glist = malloc(8 * (ng ? ng : 1)); u64 i2 = 0;
    for (u64 k = 0; k < NKAPPA; k++) if (LD.start[k + 1] > LD.start[k]) JG.glist[i2++] = k;
    JG.ng = ng; JG.nbatch = (ng + JG.bs - 1) / JG.bs; JG.done = calloc(JG.nbatch ? JG.nbatch : 1, 1); JG.countonly = countonly;
    u64 fp = 1469598103934665603ULL; { FILE *f = fopen(prob, "rb"); char buf[1 << 16]; size_t n; while ((n = fread(buf, 1, sizeof buf, f)) > 0) fp = fnv(buf, n, fp); fclose(f); }
    { char o[256]; snprintf(o, sizeof o, "ns%d S%d,%d,%d shard%d/%d bs%d A%llu B%llu D%llu", NS, SC[0], SC[1], SC[2], SHARD_I, SHARD_N, JG.bs, LA.n, LB.n, LD.n); fp = fnv(o, strlen(o), fp); }
    JG.ck = NULL;
    if (ckpath && !countonly) {
        if (resume) {
            FILE *f = fopen(ckpath, "r");
            if (f) {
                unsigned long long hdr = 0; char tag[16];
                if (fscanf(f, "%15s %llx", tag, &hdr) == 2 && !strcmp(tag, "P3CKPT")) {
                    if (hdr != fp) die("checkpoint does not match this problem/options");
                    unsigned long long b; u64 nd = 0; while (fscanf(f, "%llu", &b) == 1) if (b < JG.nbatch) { JG.done[b] = 1; nd++; }
                    fprintf(stderr, "p3: resume: %llu/%llu batches already done\n", nd, JG.nbatch);
                } else die("bad checkpoint file");
                fclose(f);
                JG.ck = fopen(ckpath, "a");
            }
        }
        if (!JG.ck) { JG.ck = fopen(ckpath, "w"); if (!JG.ck) die("cannot write checkpoint"); fprintf(JG.ck, "P3CKPT %llx\n", fp); fflush(JG.ck); }
    }
    if (hitpath) HITF = fopen(hitpath, "a");
    atomic_store(&JG.next, 0); atomic_store(&JG.pairs, 0); atomic_store(&JG.ndone, 0);
    JG.t0 = now();
    fprintf(stderr, "p3: join: %llu sigma groups in %llu batches\n", ng, JG.nbatch);
    double tj = now();
    int T = NTHREADS; pthread_t *th = calloc(T, sizeof(pthread_t));
    for (int t = 0; t < T; t++) pthread_create(&th[t], NULL, worker_join, NULL);
    for (int t = 0; t < T; t++) pthread_join(th[t], NULL);
    double dt = now() - tj;
    printf("DONE hits=%llu matches=%llu badmatches=%llu pairs=%llu join_seconds=%.2f total_seconds=%.2f\n", HITS, MATCHES, BADM, (u64)atomic_load(&JG.pairs), dt, now() - t0);
    fprintf(stderr, "p3: done: %llu pairs in %.1fs (%.3g pairs/s over %d threads)\n", (u64)atomic_load(&JG.pairs), dt, (double)atomic_load(&JG.pairs) / (dt > 0 ? dt : 1), T);
    if (JG.ck) fclose(JG.ck);
    if (HITF) fclose(HITF);
    return 0;
}
