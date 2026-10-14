/* Propus on Z_3 x GF(p) (v = 3p, p prime, p != 3) with H = {h : h = 1 mod 3, h mod p an e-th power}: symmetric
 * Hadamard matrices of order 4v from (A, B, B, D), A symmetric, PAF_A + 2 PAF_B + PAF_D = 4v delta.
 *
 * Coordinates g <-> (x, y) = (g mod 3, g mod p). An H-invariant X is three rows x_0, x_1, x_2 over GF(p), each
 * constant on {0} and on the e cosets C_c = g^c <g^e> of the e-th powers (1 + e bits). A symmetric means
 * a_0(-y) = a_0(y) and a_2(y) = a_1(-y) (y -> -y maps C_c to C_{c + e/2}; e is even, -1 is not an e-th power).
 * With P(u) the PAF of a row and CC(u, w)(t) = sum_y u(y) w(y + t), the classes of <H, -1> on Z_v \ 0 are
 *     (0, k)  k = 0 .. e/2 - 1 :  P(x_0) + P(x_1) + P(x_2)               at t in C_k
 *     (1, j)  j = 0 .. e       :  CC(x_0, x_1) + CC(x_1, x_2) + CC(x_2, x_0)  at t = 0 (j = 0) or t in C_{j-1}
 * and the equation is F_A + 2 F_B + F_D = 0 in every class.
 *
 * The run enumerates the Z_3 compressions (row sums) compatible with the compressed equation, with the row sums of
 * B and of D each taken minimal under S_3 (the translations x -> x + 1 of B, resp. D, and its reversal permute its
 * rows), and A canonical under the multipliers y -> g^j y (which fix all row sums). PSD filters
 * |B^(chi)|^2 <= 2v, |D^(chi)|^2 <= 4v, A^(chi)^2 <= 4v; the B survivors are hashed on their class values and the D
 * survivors streamed against A (prow), or (prow2 join) the D survivors hashed in parts that fit in memory and probed
 * with every pair (A, B). Every hit is rebuilt on Z_v and checked before it is printed.
 *
 *   cc -O3 -march=native -o prow prow.c -lm
 *   ./prow p e all [count | list | run START STEP PARTCAP]      (run: triples START, START+STEP, ...; D part <= PARTCAP)
 *   ./prow p e test SEED
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>

static int p, e, f, v, g, NB, NM, C, E2;
static int *ind;                 /* discrete log base g mod e: coset index of y != 0 */
static int (*Nk)[17][17];        /* Nk[j][a][b] = #{y : cls(y) = a, cls(y + t_j) = b}, t_0 = 0, t_{1+k} = g^k */
static int *rsum;                /* row sums */
static signed char *Pt;          /* Pt[m * (e/2) + k]: PAF of row m at t in C_k */
static double complex *hat;      /* hat[m * e + k]: sum_y row(y) zeta^(g^k y) */
static double bA, bB, bD;
static int canonA = 1, tgtC[32];      /* class targets (0 in a search, planted values in test mode) */
static long *tgtV;                    /* PAF-sum target on Z_v for verify() */

static inline int bit(int m, int i) { return 1 - 2 * ((m >> i) & 1); }
static inline int cls(int y) { return y == 0 ? 0 : 1 + ind[y] % e; }
static inline int rowval(int m, int y) { return bit(m, cls(y)); }
static int refl(int m) {            /* y -> -y: coset c -> c + e/2 */
    int o = m & 1;
    for (int c = 0; c < e; c++) if (m >> (1 + c) & 1) o |= 1 << (1 + (c + e / 2) % e);
    return o;
}
static int rotc(int m, int j) {     /* row'(y) = row(g^j y): coset c of row' = coset c + j of row */
    int o = m & 1;
    for (int c = 0; c < e; c++) if (m >> (1 + (c + j) % e) & 1) o |= 1 << (1 + c);
    return o;
}

static int primroot(int q) {
    for (int a = 2; a < q; a++) {
        int ok = 1, n = q - 1;
        for (int d = 2; d <= n; d++) if (n % d == 0) {
            int is_p = 1; for (int s = 2; s * s <= d; s++) if (d % s == 0) { is_p = 0; break; }
            if (!is_p) continue;
            long x = 1; for (int i = 0; i < n / d; i++) x = x * a % q;
            if (x == 1) { ok = 0; break; }
        }
        if (ok) return a;
    }
    return -1;
}

static void init(void) {
    f = (p - 1) / e; v = 3 * p; NB = 1 + e; NM = 1 << NB; C = e / 2 + 1 + e; E2 = e / 2;
    g = primroot(p);
    ind = calloc(p, sizeof(int));
    long x = 1; for (int i = 0; i < p - 1; i++) { ind[x] = i; x = x * g % p; }
    Nk = calloc(1 + e, sizeof *Nk);
    for (int j = 0; j <= e; j++) {
        long t = 0; if (j > 0) { t = 1; for (int i = 0; i < j - 1; i++) t = t * g % p; }
        for (int y = 0; y < p; y++) Nk[j][cls(y)][cls((y + t) % p)]++;
    }
    rsum = malloc(sizeof(int) * NM); Pt = malloc((size_t)NM * E2); hat = malloc(sizeof(double complex) * (size_t)NM * e);
    double complex *eta = malloc(sizeof(double complex) * e * e);   /* eta[k * e + c] = sum_{y in C_c} zeta^(g^k y) */
    for (int k = 0; k < e; k++) {
        long gk = 1; for (int i = 0; i < k; i++) gk = gk * g % p;
        for (int c = 0; c < e; c++) eta[k * e + c] = 0;
        for (int y = 1; y < p; y++) eta[k * e + ind[y] % e] += cexp(2 * M_PI * I * (double)(gk * y % p) / p);
    }
    for (int m = 0; m < NM; m++) {
        int s = bit(m, 0); for (int c = 0; c < e; c++) s += f * bit(m, 1 + c); rsum[m] = s;
        for (int k = 0; k < E2; k++) {          /* PAF at t = g^k: sum over class pairs */
            int a = 0; for (int i = 0; i < NB; i++) for (int j = 0; j < NB; j++) a += Nk[1 + k][i][j] * bit(m, i) * bit(m, j);
            Pt[(size_t)m * E2 + k] = a;
        }
        for (int k = 0; k < e; k++) {
            double complex z = bit(m, 0);
            for (int c = 0; c < e; c++) z += bit(m, 1 + c) * eta[k * e + c];
            hat[(size_t)m * e + k] = z;
        }
    }
}

static inline void ccv(int u, int w, int *out) {     /* CC(u, w) at t_j, j = 0..e */
    for (int j = 0; j <= e; j++) {
        int a = 0;
        for (int i = 0; i < NB; i++) { int ui = bit(u, i); for (int k = 0; k < NB; k++) a += Nk[j][i][k] * ui * bit(w, k); }
        out[j] = a;
    }
}

/* ---- sequences: three row lists and pair tables ---- */
typedef struct { int S[3], n[3]; int *L[3]; signed char *T01, *T12, *T20; double bound; } Seq;

static int *rows_with_sum(int S, int *n) {
    int *L = malloc(sizeof(int) * NM), k = 0;
    for (int m = 0; m < NM; m++) if (rsum[m] == S) L[k++] = m;
    *n = k; return L;
}
static signed char *pairtab(int *La, int na, int *Lb, int nb) {
    signed char *T = malloc((size_t)na * nb * (1 + e)); int out[17];
    for (int i = 0; i < na; i++) for (int j = 0; j < nb; j++) {
        ccv(La[i], Lb[j], out);
        for (int t = 0; t <= e; t++) T[((size_t)i * nb + j) * (1 + e) + t] = out[t];
    }
    return T;
}
static void seq_init(Seq *q, int S0, int S1, int S2, double bound) {
    q->S[0] = S0; q->S[1] = S1; q->S[2] = S2; q->bound = bound;
    for (int i = 0; i < 3; i++) q->L[i] = rows_with_sum(q->S[i], &q->n[i]);
    q->T01 = pairtab(q->L[0], q->n[0], q->L[1], q->n[1]);
    q->T12 = pairtab(q->L[1], q->n[1], q->L[2], q->n[2]);
    q->T20 = pairtab(q->L[2], q->n[2], q->L[0], q->n[0]);
}
static void seq_free(Seq *q) { for (int i = 0; i < 3; i++) free(q->L[i]); free(q->T01); free(q->T12); free(q->T20); }

static inline void classvals(const Seq *q, int i0, int i1, int i2, signed char *F) {
    int a = q->L[0][i0], b = q->L[1][i1], c = q->L[2][i2], W = 1 + e;
    const signed char *x01 = q->T01 + ((size_t)i0 * q->n[1] + i1) * W;
    const signed char *x12 = q->T12 + ((size_t)i1 * q->n[2] + i2) * W;
    const signed char *x20 = q->T20 + ((size_t)i2 * q->n[0] + i0) * W;
    for (int k = 0; k < E2; k++) F[k] = Pt[(size_t)a * E2 + k] + Pt[(size_t)b * E2 + k] + Pt[(size_t)c * E2 + k];
    for (int j = 0; j <= e; j++) F[E2 + j] = x01[j] + x12[j] + x20[j];
}

typedef void (*cb_t)(const Seq *, int, int, int, void *);
static const double complex W3[3] = {1.0, -0.5 + 0.86602540378443864676 * I, -0.5 - 0.86602540378443864676 * I};

static long long seq_iter(const Seq *q, cb_t cb, void *ctx) {
    long long cnt = 0;
    double complex u[3][16];
    for (int i0 = 0; i0 < q->n[0]; i0++) for (int i1 = 0; i1 < q->n[1]; i1++) {
        int a = q->L[0][i0], b = q->L[1][i1];
        for (int k = 0; k < 3; k++) for (int j = 0; j < e; j++)
            u[k][j] = hat[(size_t)a * e + j] + W3[k] * hat[(size_t)b * e + j];
        for (int i2 = 0; i2 < q->n[2]; i2++) {
            int c = q->L[2][i2], ok = 1;
            for (int k = 0; k < 3 && ok; k++) {
                const double complex w2 = W3[(2 * k) % 3];
                for (int j = 0; j < e; j++) {
                    double complex z = u[k][j] + w2 * hat[(size_t)c * e + j];
                    if (creal(z) * creal(z) + cimag(z) * cimag(z) > q->bound + 1e-6) { ok = 0; break; }
                }
            }
            if (!ok) continue;
            cnt++;
            if (cb) cb(q, i0, i1, i2, ctx);
        }
    }
    return cnt;
}

/* ---- hashing ---- */
#define MAXC 32
static inline uint64_t mix(uint64_t x) { x ^= x >> 31; x *= 0x7fb5d329728ea185ULL; x ^= x >> 27; x *= 0x81dadef4bc2dd44dULL; x ^= x >> 33; return x; }
static inline uint64_t hkey(const signed char *F) {
    uint64_t w[4] = {0, 0, 0, 0}; memcpy(w, F, C);
    return mix(w[0] ^ mix(w[1] ^ mix(w[2] ^ mix(w[3] + 0x9e3779b97f4a7c15ULL))));
}
typedef struct { int i0, i1, i2; signed char F[MAXC]; } Ent;
typedef struct { Ent *e; long long n, cap; uint32_t *tab; uint64_t mask; } Store;
static void store_cb(const Seq *q, int i0, int i1, int i2, void *ctx) {
    Store *s = ctx;
    if (s->n == s->cap) { s->cap = s->cap ? 2 * s->cap : 1 << 20; s->e = realloc(s->e, sizeof(Ent) * s->cap); }
    Ent *x = &s->e[s->n++]; x->i0 = i0; x->i1 = i1; x->i2 = i2; classvals(q, i0, i1, i2, x->F);
}
static void store_hash(Store *s) {
    uint64_t sz = 1; while (sz < 2 * (uint64_t)s->n + 16) sz <<= 1;
    s->tab = malloc(sizeof(uint32_t) * sz); memset(s->tab, 0xff, sizeof(uint32_t) * sz); s->mask = sz - 1;
    for (long long i = 0; i < s->n; i++) {
        uint64_t h = hkey(s->e[i].F) & s->mask;
        while (s->tab[h] != 0xffffffffu) h = (h + 1) & s->mask;
        s->tab[h] = (uint32_t)i;
    }
}

/* ---- A ---- */
typedef struct { int a0, a1; signed char F[MAXC]; } AEnt;
static AEnt *Al; static int nA;
static void build_A(int S0, int S1) {
    int n0, n1, *L0 = rows_with_sum(S0, &n0), *L1 = rows_with_sum(S1, &n1);
    Al = malloc(sizeof(AEnt) * ((size_t)n0 * n1 + 1)); nA = 0;
    for (int i = 0; i < n0; i++) {
        int a0 = L0[i]; if (refl(a0) != a0) continue;
        for (int j = 0; j < n1; j++) {
            int a1 = L1[j], a2 = refl(a1), keep = 1;
            for (int r = 1; r < e && keep && canonA; r++) {           /* multiplier-canonical (a0, a1) */
                int b0 = rotc(a0, r), b1 = rotc(a1, r);
                if (b0 < a0 || (b0 == a0 && b1 < a1)) keep = 0;
            }
            if (!keep) continue;
            for (int k = 0; k < 3 && keep; k++) for (int jj = 0; jj < e; jj++) {
                double complex z = hat[(size_t)a0 * e + jj] + W3[k] * hat[(size_t)a1 * e + jj] + W3[(2 * k) % 3] * hat[(size_t)a2 * e + jj];
                if (creal(z) * creal(z) + cimag(z) * cimag(z) > bA + 1e-6) { keep = 0; break; }
            }
            if (!keep) continue;
            AEnt *x = &Al[nA++]; x->a0 = a0; x->a1 = a1;
            int c01[17], c12[17], c20[17];
            ccv(a0, a1, c01); ccv(a1, a2, c12); ccv(a2, a0, c20);
            for (int k = 0; k < E2; k++) x->F[k] = Pt[(size_t)a0 * E2 + k] + Pt[(size_t)a1 * E2 + k] + Pt[(size_t)a2 * E2 + k];
            for (int jj = 0; jj <= e; jj++) x->F[E2 + jj] = c01[jj] + c12[jj] + c20[jj];
        }
    }
    free(L0); free(L1);
}

/* ---- verification on Z_v ---- */
static void tov(int x0, int x1, int x2, int *X) {
    int r[3] = {x0, x1, x2};
    for (int gg = 0; gg < v; gg++) X[gg] = rowval(r[gg % 3], gg % p);
}
static int verify(const int *A, const int *B, const int *D) {
    for (int s = 0; s < v; s++) {
        long t = 0;
        for (int gg = 0; gg < v; gg++) t += A[gg] * A[(gg + s) % v] + 2 * B[gg] * B[(gg + s) % v] + D[gg] * D[(gg + s) % v];
        if (t != tgtV[s]) return 0;
    }
    return 1;
}
static long long nhits;
static void report(const AEnt *a, int b0, int b1, int b2, int d0, int d1, int d2) {
    int *A = malloc(sizeof(int) * v), *B = malloc(sizeof(int) * v), *D = malloc(sizeof(int) * v);
    tov(a->a0, a->a1, refl(a->a1), A); tov(b0, b1, b2, B); tov(d0, d1, d2, D);
    int ok = verify(A, B, D);
    printf("HIT %s A ", ok ? "verified" : "WRONG");
    for (int gg = 0; gg < v; gg++) putchar(A[gg] > 0 ? '+' : '-');
    printf(" B "); for (int gg = 0; gg < v; gg++) putchar(B[gg] > 0 ? '+' : '-');
    printf(" D "); for (int gg = 0; gg < v; gg++) putchar(D[gg] > 0 ? '+' : '-');
    printf("\n"); fflush(stdout);
    if (ok) nhits++;
    free(A); free(B); free(D);
}

typedef struct { Store *hs; const Seq *hq; long long probes; } JCtx;
static void join_cb(const Seq *q, int i0, int i1, int i2, void *vctx) {     /* streamed D, hashed B */
    JCtx *c = vctx;
    signed char F[MAXC], need[MAXC];
    classvals(q, i0, i1, i2, F);
    for (int ia = 0; ia < nA; ia++) {
        const signed char *FA = Al[ia].F; int ok = 1;
        memset(need, 0, MAXC);
        for (int k = 0; k < C; k++) {
            int s = tgtC[k] - FA[k] - F[k];
            if (s & 1) { ok = 0; break; }
            s /= 2; if (s < -127 || s > 127) { ok = 0; break; }
            need[k] = (signed char)s;
        }
        if (!ok) continue;
        c->probes++;
        uint64_t h = hkey(need) & c->hs->mask;
        while (c->hs->tab[h] != 0xffffffffu) {
            const Ent *x = &c->hs->e[c->hs->tab[h]];
            if (!memcmp(x->F, need, C))
                report(&Al[ia], c->hq->L[0][x->i0], c->hq->L[1][x->i1], c->hq->L[2][x->i2], q->L[0][i0], q->L[1][i1], q->L[2][i2]);
            h = (h + 1) & c->hs->mask;
        }
    }
}

/* ---- partitioned join: D hashed in parts, probed by all pairs (A, B) ---- */
typedef struct { uint64_t h; uint16_t i0, i1, i2, used; } DEnt;
typedef struct { DEnt *t; uint64_t mask; long long n; int part, nparts; } DTab;
static inline int partof(uint64_t h, int nparts) { return (int)((h >> 40) % (uint64_t)nparts); }
static void dtab_cb(const Seq *q, int i0, int i1, int i2, void *ctx) {
    DTab *d = ctx; signed char F[MAXC];
    memset(F, 0, MAXC); classvals(q, i0, i1, i2, F);
    uint64_t h = hkey(F);
    if (partof(h, d->nparts) != d->part) return;
    uint64_t s = h & d->mask;
    while (d->t[s].used) s = (s + 1) & d->mask;
    d->t[s].h = h; d->t[s].i0 = i0; d->t[s].i1 = i1; d->t[s].i2 = i2; d->t[s].used = 1; d->n++;
}
static long long join2(Seq *qb, Seq *qd, long long nd, long long partcap) {
    Store sb = {0};
    seq_iter(qb, store_cb, &sb);
    int nparts = (int)((nd + partcap - 1) / partcap); if (nparts < 1) nparts = 1;
    uint64_t sz = 1; while (sz < 2 * (uint64_t)((nd + nparts - 1) / nparts) + 16) sz <<= 1;
    DTab d; d.t = malloc(sizeof(DEnt) * sz); d.mask = sz - 1; d.nparts = nparts;
    long long probes = 0;
    for (int part = 0; part < nparts; part++) {
        memset(d.t, 0, sizeof(DEnt) * sz); d.n = 0; d.part = part;
        seq_iter(qd, dtab_cb, &d);
        for (int ia = 0; ia < nA; ia++) {
            const signed char *FA = Al[ia].F;
            for (long long ib = 0; ib < sb.n; ib++) {
                const signed char *FB = sb.e[ib].F; signed char need[MAXC]; int ok = 1;
                memset(need, 0, MAXC);
                for (int k = 0; k < C; k++) { int s = tgtC[k] - FA[k] - 2 * FB[k]; if (s < -127 || s > 127) { ok = 0; break; } need[k] = (signed char)s; }
                if (!ok) continue;
                uint64_t h = hkey(need);
                if (partof(h, nparts) != part) continue;
                probes++;
                uint64_t s = h & d.mask;
                while (d.t[s].used) {
                    if (d.t[s].h == h) {
                        signed char FD[MAXC]; memset(FD, 0, MAXC); classvals(qd, d.t[s].i0, d.t[s].i1, d.t[s].i2, FD);
                        if (!memcmp(FD, need, C))
                            report(&Al[ia], qb->L[0][sb.e[ib].i0], qb->L[1][sb.e[ib].i1], qb->L[2][sb.e[ib].i2],
                                   qd->L[0][d.t[s].i0], qd->L[1][d.t[s].i1], qd->L[2][d.t[s].i2]);
                    }
                    s = (s + 1) & d.mask;
                }
            }
        }
        printf("    part %d/%d: D entries %lld, probes so far %lld\n", part + 1, nparts, d.n, probes); fflush(stdout);
    }
    free(d.t); free(sb.e);
    return probes;
}

/* ---- join3: D partitioned by an additive key of its s = 0 classes, so each pass enumerates only its own part.
 * Every row PAF is 1 mod 4, so c(r) = ((P_0(r) - 1)/4 + 7 (P_1(r) - 1)/4) mod nparts is additive over the three rows and
 * the part of a D is c(d_0) + c(d_1) + c(d_2) mod nparts; the part of a needed class vector is computed from its
 * entries 0 and 1 in the same way. */
static inline int pmod(int a, int n) { a %= n; return a < 0 ? a + n : a; }
static inline int rowkey(int m, int np) { return pmod((Pt[(size_t)m * E2] - 1) / 4 + 7 * ((Pt[(size_t)m * E2 + 1] - 1) / 4), np); }
static long long join3(Seq *qb, Seq *qd, long long nd, long long partcap) {
    Store sb = {0};
    seq_iter(qb, store_cb, &sb);
    int np = (int)((nd + partcap - 1) / partcap); if (np < 1) np = 1;
    int *c0 = malloc(sizeof(int) * qd->n[0]), *c1 = malloc(sizeof(int) * qd->n[1]);
    for (int i = 0; i < qd->n[0]; i++) c0[i] = rowkey(qd->L[0][i], np);
    for (int i = 0; i < qd->n[1]; i++) c1[i] = rowkey(qd->L[1][i], np);
    int **bk = malloc(sizeof(int *) * np), *bn = calloc(np, sizeof(int));     /* third rows bucketed by key */
    for (int b = 0; b < np; b++) bk[b] = malloc(sizeof(int) * qd->n[2]);
    for (int i = 0; i < qd->n[2]; i++) { int k = rowkey(qd->L[2][i], np); bk[k][bn[k]++] = i; }
    uint64_t sz = 1; while (sz < 3 * (uint64_t)((nd + np - 1) / np) + 16) sz <<= 1;
    DTab d; d.t = malloc(sizeof(DEnt) * sz); d.mask = sz - 1; d.nparts = 1; d.part = 0;
    long long probes = 0;
    double complex u[3][16];
    for (int part = 0; part < np; part++) {
        memset(d.t, 0, sizeof(DEnt) * sz); d.n = 0;
        for (int i0 = 0; i0 < qd->n[0]; i0++) for (int i1 = 0; i1 < qd->n[1]; i1++) {
            int want = pmod(part - c0[i0] - c1[i1], np), a = qd->L[0][i0], b = qd->L[1][i1];
            if (!bn[want]) continue;
            for (int k = 0; k < 3; k++) for (int j = 0; j < e; j++) u[k][j] = hat[(size_t)a * e + j] + W3[k] * hat[(size_t)b * e + j];
            for (int x = 0; x < bn[want]; x++) {
                int i2 = bk[want][x], c = qd->L[2][i2], ok = 1;
                for (int k = 0; k < 3 && ok; k++) {
                    const double complex w2 = W3[(2 * k) % 3];
                    for (int j = 0; j < e; j++) {
                        double complex z = u[k][j] + w2 * hat[(size_t)c * e + j];
                        if (creal(z) * creal(z) + cimag(z) * cimag(z) > qd->bound + 1e-6) { ok = 0; break; }
                    }
                }
                if (!ok) continue;
                if ((uint64_t)d.n * 2 >= sz) { fprintf(stderr, "part %d overflows its table\n", part); exit(2); }
                dtab_cb(qd, i0, i1, i2, &d);
            }
        }
        for (int ia = 0; ia < nA; ia++) {
            const signed char *FA = Al[ia].F;
            for (long long ib = 0; ib < sb.n; ib++) {
                const signed char *FB = sb.e[ib].F;
                int n0 = tgtC[0] - FA[0] - 2 * FB[0], n1 = tgtC[1] - FA[1] - 2 * FB[1];
                if (pmod(n0, 4) != 3 || pmod(n1, 4) != 3) continue;
                if (pmod((n0 - 3) / 4 + 7 * ((n1 - 3) / 4), np) != part) continue;
                signed char need[MAXC]; int ok = 1;
                memset(need, 0, MAXC);
                for (int k = 0; k < C; k++) { int s = tgtC[k] - FA[k] - 2 * FB[k]; if (s < -127 || s > 127) { ok = 0; break; } need[k] = (signed char)s; }
                if (!ok) continue;
                probes++;
                uint64_t h = hkey(need), s = h & d.mask;
                while (d.t[s].used) {
                    if (d.t[s].h == h) {
                        signed char FD[MAXC]; memset(FD, 0, MAXC); classvals(qd, d.t[s].i0, d.t[s].i1, d.t[s].i2, FD);
                        if (!memcmp(FD, need, C))
                            report(&Al[ia], qb->L[0][sb.e[ib].i0], qb->L[1][sb.e[ib].i1], qb->L[2][sb.e[ib].i2],
                                   qd->L[0][d.t[s].i0], qd->L[1][d.t[s].i1], qd->L[2][d.t[s].i2]);
                    }
                    s = (s + 1) & d.mask;
                }
            }
        }
        printf("    part %d/%d: D entries %lld, probes so far %lld\n", part + 1, np, d.n, probes); fflush(stdout);
    }
    free(d.t); free(sb.e); free(c0); free(c1); for (int b = 0; b < np; b++) free(bk[b]); free(bk); free(bn);
    return probes;
}

/* ---- compatible Z_3 compressions ---- */
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: prow p e [ra,rb,rd|all] [count]\n"); return 1; }
    p = atoi(argv[1]); e = atoi(argv[2]);
    int want[3] = {0, 0, 0}, all = argc < 4 || !strcmp(argv[3], "all");
    if (!all) sscanf(argv[3], "%d,%d,%d", &want[0], &want[1], &want[2]);
    if ((p - 1) % e || e % 2 || ((p - 1) / e) % 2 == 0 || 1 + e > 16) { fprintf(stderr, "need e | p-1, e even, (p-1)/e odd, e <= 15\n"); return 1; }
    init();
    bA = 4.0 * v; bB = 2.0 * v; bD = 4.0 * v;
    tgtV = calloc(v, sizeof(long)); tgtV[0] = 4 * v;
    if (argc > 3 && !strcmp(argv[3], "test")) {     /* plant a triple with rows of extreme sums; PSD and A-canonical off */
        uint64_t st = argc > 4 ? strtoull(argv[4], 0, 10) : 1;
        int rows[9];
        for (int i = 0; i < 9; i++) for (;;) {
            st = mix(st + 0x9e3779b97f4a7c15ULL); int m = (int)(st % NM);
            int S = 0; for (int c = 0; c < e; c++) S += bit(m, 1 + c);
            if (S < e - 4 && S > -(e - 4)) continue;
            if (i == 0 && refl(m) != m) continue;            /* a0 symmetric */
            rows[i] = m; break;
        }
        int a0 = rows[0], a1 = rows[1], a2 = refl(a1);
        int *A = malloc(sizeof(int) * v), *B = malloc(sizeof(int) * v), *D = malloc(sizeof(int) * v);
        tov(a0, a1, a2, A); tov(rows[3], rows[4], rows[5], B); tov(rows[6], rows[7], rows[8], D);
        for (int s = 0; s < v; s++) { long t = 0; for (int gg = 0; gg < v; gg++) t += A[gg] * A[(gg + s) % v] + 2 * B[gg] * B[(gg + s) % v] + D[gg] * D[(gg + s) % v]; tgtV[s] = t; }
        /* class targets from the Z_v PAF sum at class representatives */
        for (int k = 0; k < E2; k++) { long t = 1; for (int i = 0; i < k; i++) t = t * g % p; int gg; for (gg = 0; gg < v; gg++) if (gg % 3 == 0 && gg % p == t) break; tgtC[k] = tgtV[gg]; }
        for (int j = 0; j <= e; j++) { long t = 0; if (j) { t = 1; for (int i = 0; i < j - 1; i++) t = t * g % p; } int gg; for (gg = 0; gg < v; gg++) if (gg % 3 == 1 && gg % p == t) break; tgtC[E2 + j] = tgtV[gg]; }
        bA = bB = bD = 1e18; canonA = 0;
        int R[3] = {rsum[rows[3]], rsum[rows[4]], rsum[rows[5]]}, Q[3] = {rsum[rows[6]], rsum[rows[7]], rsum[rows[8]]};
        for (int i = 0; i < 3; i++) for (int j = i + 1; j < 3; j++) { if (R[j] < R[i]) { int x = R[i]; R[i] = R[j]; R[j] = x; } if (Q[j] < Q[i]) { int x = Q[i]; Q[i] = Q[j]; Q[j] = x; } }
        printf("planted rows A %d %d B %d %d %d D %d %d %d; sums A (%d,%d) B (%d,%d,%d) D (%d,%d,%d)\n", a0, a1, rows[3], rows[4], rows[5], rows[6], rows[7], rows[8],
               rsum[a0], rsum[a1], R[0], R[1], R[2], Q[0], Q[1], Q[2]);
        Seq qb, qd; seq_init(&qb, R[0], R[1], R[2], bB); seq_init(&qd, Q[0], Q[1], Q[2], bD);
        build_A(rsum[a0], rsum[a1]);
        long long nd = seq_iter(&qd, NULL, NULL);
        long long pr = join3(&qb, &qd, nd, nd / 3 + 1);
        printf("test: |A| %d, D %lld, probes %lld, verified hits %lld\n", nA, nd, pr, nhits);
        return 0;
    }
    printf("v = %d = 3 * %d, e = %d, f = %d, g = %d, classes %d\n", v, p, e, f, g, C); fflush(stdout);
    /* attainable row sums */
    int *ok = calloc(2 * p + 1, sizeof(int));
    for (int m = 0; m < NM; m++) ok[rsum[m] + p] = 1;
    int nv = 0, vals[256]; for (int s = -p; s <= p; s++) if (ok[s + p]) vals[nv++] = s;
    /* A0 must be attainable by a symmetric row */
    int *oksym = calloc(2 * p + 1, sizeof(int));
    for (int m = 0; m < NM; m++) if (refl(m) == m) oksym[rsum[m] + p] = 1;
    const char *mode = argc > 4 ? argv[4] : "list";
    int start = argc > 5 ? atoi(argv[5]) : 0, step = argc > 6 ? atoi(argv[6]) : 1;
    long long partcap = argc > 7 ? atoll(argv[7]) : 60000000LL;
    int nt = 0; static int T[100000][8];
    for (int ia0 = 0; ia0 < nv; ia0++) for (int ia1 = 0; ia1 < nv; ia1++) {
        int A0 = vals[ia0], A1 = vals[ia1], ra = A0 + 2 * A1;
        if (ra <= 0 || !oksym[A0 + p]) continue;
        int sA0 = A0 * A0 + 2 * A1 * A1, sA1 = 2 * A0 * A1 + A1 * A1;
        for (int ib0 = 0; ib0 < nv; ib0++) for (int ib1 = ib0; ib1 < nv; ib1++) for (int ib2 = ib1; ib2 < nv; ib2++) {
            int R0 = vals[ib0], R1 = vals[ib1], R2 = vals[ib2], rb = R0 + R1 + R2;
            int sB0 = R0 * R0 + R1 * R1 + R2 * R2, sB1 = R0 * R1 + R1 * R2 + R2 * R0;
            if (rb < 0 || sA0 + 2 * sB0 > 4 * v) continue;
            for (int id0 = 0; id0 < nv; id0++) for (int id1 = id0; id1 < nv; id1++) for (int id2 = id1; id2 < nv; id2++) {
                int Q0 = vals[id0], Q1 = vals[id1], Q2 = vals[id2], rd = Q0 + Q1 + Q2;
                if (rd < 0) continue;
                if (sA0 + 2 * sB0 + Q0 * Q0 + Q1 * Q1 + Q2 * Q2 != 4 * v) continue;
                if (sA1 + 2 * sB1 + Q0 * Q1 + Q1 * Q2 + Q2 * Q0 != 0) continue;
                if (!all && (ra != want[0] || rb != want[1] || rd != want[2])) continue;
                int *t = T[nt++]; t[0] = A0; t[1] = A1; t[2] = R0; t[3] = R1; t[4] = R2; t[5] = Q0; t[6] = Q1; t[7] = Q2;
            }
        }
    }
    printf("%d compatible Z_3 triples (B, D row sums sorted, ra, rb, rd > 0)\n", nt); fflush(stdout);
    double work = 0;
    for (int it = start; it < nt; it += step) {
        int *t = T[it];
        printf("triple %d: rows (%d,%d,%d) A (%d,%d,%d) B (%d,%d,%d) D (%d,%d,%d)", it, t[0] + 2 * t[1], t[2] + t[3] + t[4], t[5] + t[6] + t[7],
               t[0], t[1], t[1], t[2], t[3], t[4], t[5], t[6], t[7]);
        if (!strcmp(mode, "list")) { printf("\n"); continue; }
        Seq qb, qd;
        seq_init(&qb, t[2], t[3], t[4], bB); seq_init(&qd, t[5], t[6], t[7], bD);
        build_A(t[0], t[1]);
        long long nb = nA ? seq_iter(&qb, NULL, NULL) : 0, nd = (nA && nb) ? seq_iter(&qd, NULL, NULL) : 0;
        printf(": |A| %d, B %lld, D %lld, pairs %.3g\n", nA, nb, nd, (double)nA * nb); fflush(stdout);
        work += (double)nA * nb;
        if (!strcmp(mode, "run") && nA && nb && nd) {
            long long pr = join3(&qb, &qd, nd, partcap);
            printf("  done triple %d: probes %lld, verified hits so far %lld\n", it, pr, nhits); fflush(stdout);
        }
        free(Al); seq_free(&qb); seq_free(&qd);
    }
    printf("DONE %d triples, pairs %.3g, verified hits %lld\n", nt, work, nhits);
    return 0;
}
