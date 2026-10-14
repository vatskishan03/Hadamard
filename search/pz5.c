/* Propus on Z_65 = Z_5 x Z_13 with H = <14> (H acts as -1 on Z_5 and trivially on Z_13): symmetric Hadamard
 * matrices of order 260 from (A, B, B, D), A symmetric, PAF_A + 2 PAF_B + PAF_D = 260 delta.
 *
 * An H-invariant X is three +-1 rows x0, x1, x2 of length 13: X(a, y) = x_{|a|}(y) (|0| = 0, |+-1| = 1,
 * |+-2| = 2); A's rows are symmetric in y. With C_ii the PAF of row i and S_ij(t) = CC_ij(t) + CC_ij(-t)
 * (S_ij(0) = 2 CC_ij(0)), the 20 classes of <H, -1> on Z_65 \ 0 are
 *     F(0, t) = C_00 + 2 C_11 + 2 C_22              t = 1..6
 *     F(1, t) = S_01 + S_12 + C_22                  t = 0..6
 *     F(2, t) = S_02 + S_12 + C_11                  t = 0..6
 * and the equation is F_A + 2 F_B + F_D = 0 in every class.
 *
 * One run handles one compatible Z_5 compression triple (row sums P_i of A, R_i of B, Q_i of D).
 * Symmetries used: B and D may each be translated by Z_13 and reflected (y -> -y) independently (translations by
 * the fixed points of H and the reversal keep H-invariance and the PAF), so their first rows are taken canonical
 * under the dihedral group D_13; A is taken canonical under the multipliers y -> mu y (mu = 1..6), which fix the
 * Z_5 compression. PSD filters: |B^(chi)|^2 <= 130, |D^(chi)|^2 <= 260, A^(chi)^2 <= 260 at every character.
 * The smaller of the B and D survivor lists is hashed on its class values; the other is streamed against A.
 * Every hit is rebuilt on Z_65 and its PAF sum checked before it is printed.
 *
 *   cc -O3 -march=native -o p65 p65.c -lm
 *   ./p65 P0 P1 P2 R0 R1 R2 Q0 Q1 Q2 [hashcap]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>

static int N = 13, M = 1 << 13, QB = 0;   /* row length (Z_N, 5 does not divide N); QB = 5: Z_5-bordered targets */
#define MAXM (1 << 16)
#define NC 20

static int rsum[MAXM], paf[MAXM][7];
static int target[NC];
static int V;                          /* 5 N */                 /* class targets (0 for the search; planted values in test mode) */
static long tgt65[5 * 16];                 /* PAF-sum target on Z_65 used by verify() */
static double boundA = 260, boundB = 130, boundD = 260;
static int canonA = 1;                 /* multiplier-canonical A (valid for multiplier-invariant targets only) */
static double complex dft[MAXM][7];
static unsigned char canon[MAXM], issym[MAXM];

static int rot(int m, int t) { int o = 0; for (int y = 0; y < N; y++) if (m >> ((y + t) % N) & 1) o |= 1 << y; return o; }
static int refl(int m) { int o = 0; for (int y = 0; y < N; y++) if (m >> ((N - y) % N) & 1) o |= 1 << y; return o; }
static int mulm(int m, int mu) { int o = 0; for (int y = 0; y < N; y++) if (m >> ((mu * y) % N) & 1) o |= 1 << y; return o; }
static inline int val(int m, int y) { return 1 - 2 * ((m >> y) & 1); }
static inline int cc(int a, int b, int t) { return N - 2 * __builtin_popcount((unsigned)(a ^ rot(b, t))); }

static void init(void) {
    for (int m = 0; m < M; m++) {
        int s = 0; for (int y = 0; y < N; y++) s += val(m, y); rsum[m] = s;
        for (int t = 0; t < 7; t++) { int p = 0; for (int y = 0; y < N; y++) p += val(m, y) * val(m, (y + t) % N); paf[m][t] = p; }
        for (int j = 0; j < 7; j++) { double complex z = 0; for (int y = 0; y < N; y++) z += val(m, y) * cexp(2 * M_PI * I * j * y / N); dft[m][j] = z; }
        int best = m, r = refl(m);
        for (int t = 0; t < N; t++) { int a = rot(m, t), b = rot(r, t); if (a < best) best = a; if (b < best) best = b; }
        canon[m] = best == m; issym[m] = r == m;
    }
}

static int *bysum(int S, int canonical, int *n) {
    int *L = malloc(sizeof(int) * M); int k = 0;
    for (int m = 0; m < M; m++) if (rsum[m] == S && (!canonical || canon[m])) L[k++] = m;
    *n = k; return L;
}

/* S table for pairs (a in La, b in Lb): S[(ia * nb + ib) * 7 + t] */
static signed char *stab(int *La, int na, int *Lb, int nb) {
    signed char *T = malloc((size_t)na * nb * 7);
    for (int i = 0; i < na; i++) for (int j = 0; j < nb; j++) {
        signed char *o = T + ((size_t)i * nb + j) * 7;
        o[0] = 2 * cc(La[i], Lb[j], 0);
        for (int t = 1; t < 7; t++) o[t] = cc(La[i], Lb[j], t) + cc(La[i], Lb[j], N - t);
    }
    return T;
}

static const double CA[3] = {2.0, 0.6180339887498949, -1.6180339887498949};
static const double CB[3] = {2.0, -1.6180339887498949, 0.6180339887498949};

typedef struct { int S[3]; int *L[3]; int n[3]; signed char *s01, *s12, *s02; double bound; } Seq;

static void seq_init(Seq *q, int S0, int S1, int S2, double bound) {
    q->S[0] = S0; q->S[1] = S1; q->S[2] = S2; q->bound = bound;
    q->L[0] = bysum(S0, 1, &q->n[0]); q->L[1] = bysum(S1, 0, &q->n[1]); q->L[2] = bysum(S2, 0, &q->n[2]);
    q->s01 = stab(q->L[0], q->n[0], q->L[1], q->n[1]);
    q->s12 = stab(q->L[1], q->n[1], q->L[2], q->n[2]);
    q->s02 = stab(q->L[0], q->n[0], q->L[2], q->n[2]);
}

static inline void classvals(const Seq *q, int i0, int i1, int i2, signed char *F) {
    int a = q->L[0][i0], b = q->L[1][i1], c = q->L[2][i2];
    const signed char *x01 = q->s01 + ((size_t)i0 * q->n[1] + i1) * 7;
    const signed char *x12 = q->s12 + ((size_t)i1 * q->n[2] + i2) * 7;
    const signed char *x02 = q->s02 + ((size_t)i0 * q->n[2] + i2) * 7;
    for (int t = 1; t < 7; t++) F[t - 1] = paf[a][t] + 2 * paf[b][t] + 2 * paf[c][t];
    for (int t = 0; t < 7; t++) F[6 + t] = x01[t] + x12[t] + paf[c][t];
    for (int t = 0; t < 7; t++) F[13 + t] = x02[t] + x12[t] + paf[b][t];
}

/* iterate the PSD survivors; callback gets indices */
typedef void (*cb_t)(const Seq *, int, int, int, void *);
static long long seq_iter(const Seq *q, cb_t cb, void *ctx) {
    long long cnt = 0;
    double complex u[3][7];
    for (int i0 = 0; i0 < q->n[0]; i0++) for (int i1 = 0; i1 < q->n[1]; i1++) {
        int a = q->L[0][i0], b = q->L[1][i1];
        for (int k = 0; k < 3; k++) for (int j = 0; j < 7; j++) u[k][j] = dft[a][j] + CA[k] * dft[b][j];
        for (int i2 = 0; i2 < q->n[2]; i2++) {
            int c = q->L[2][i2], ok = 1;
            for (int k = 0; k < 3 && ok; k++) for (int j = (k == 0); j < 7; j++) {
                double complex z = u[k][j] + CB[k] * dft[c][j];
                if (creal(z) * creal(z) + cimag(z) * cimag(z) > q->bound + 1e-6) { ok = 0; break; }
            }
            if (!ok) continue;
            cnt++;
            if (cb) cb(q, i0, i1, i2, ctx);
        }
    }
    return cnt;
}

/* ---- hash of class-value vectors ---- */
static inline uint64_t mix(uint64_t x) { x ^= x >> 31; x *= 0x7fb5d329728ea185ULL; x ^= x >> 27; x *= 0x81dadef4bc2dd44dULL; x ^= x >> 33; return x; }
static inline uint64_t hkey(const signed char *F) {
    uint64_t w[3] = {0, 0, 0}; memcpy(w, F, NC);
    return mix(w[0] ^ mix(w[1] ^ mix(w[2] + 0x9e3779b97f4a7c15ULL)));
}

typedef struct { int i0, i1, i2; signed char F[NC]; } Ent;
typedef struct { Ent *e; long long n, cap; uint32_t *tab; uint64_t mask; } Store;

static void store_cb(const Seq *q, int i0, int i1, int i2, void *ctx) {
    Store *s = ctx;
    if (s->n == s->cap) { s->cap = s->cap ? 2 * s->cap : 1 << 20; s->e = realloc(s->e, sizeof(Ent) * s->cap); }
    Ent *e = &s->e[s->n++]; e->i0 = i0; e->i1 = i1; e->i2 = i2; classvals(q, i0, i1, i2, e->F);
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

/* ---- A list ---- */
typedef struct { int p[3]; signed char F[NC]; } AEnt;
static AEnt *Alist; static int nA;

static void build_A(int P0, int P1, int P2) {
    int *L[3], n[3];
    int Ps[3] = {P0, P1, P2};
    for (int i = 0; i < 3; i++) { L[i] = malloc(sizeof(int) * 1024); n[i] = 0; for (int m = 0; m < M; m++) if (issym[m] && rsum[m] == Ps[i]) L[i][n[i]++] = m; }
    Alist = malloc(sizeof(AEnt) * (size_t)(n[0] * n[1] * n[2] + 1)); nA = 0;
    for (int a = 0; a < n[0]; a++) for (int b = 0; b < n[1]; b++) for (int c = 0; c < n[2]; c++) {
        int p0 = L[0][a], p1 = L[1][b], p2 = L[2][c], keep = 1;
        for (int mu = 2; mu <= N / 2 && keep && canonA; mu++) {
            int gg = mu, nn = N; while (nn) { int t = gg % nn; gg = nn; nn = t; }
            if (gg != 1) continue;                                  /* multipliers must be units of Z_N */
            int q0 = mulm(p0, mu), q1 = mulm(p1, mu), q2 = mulm(p2, mu);
            if (q0 < p0 || (q0 == p0 && (q1 < p1 || (q1 == p1 && q2 < p2)))) keep = 0;
        }
        if (!keep) continue;
        for (int k = 0; k < 3 && keep; k++) for (int j = (k == 0); j < 7; j++) {
            double z = creal(dft[p0][j] + CA[k] * dft[p1][j] + CB[k] * dft[p2][j]);
            if (z * z > boundA + 1e-6) { keep = 0; break; }
        }
        if (!keep) continue;
        AEnt *e = &Alist[nA++]; e->p[0] = p0; e->p[1] = p1; e->p[2] = p2;
        /* class values of A: S_ij from cc directly */
        signed char *F = e->F;
        int S01[7], S12[7], S02[7];
        S01[0] = 2 * cc(p0, p1, 0); S12[0] = 2 * cc(p1, p2, 0); S02[0] = 2 * cc(p0, p2, 0);
        for (int t = 1; t < 7; t++) { S01[t] = cc(p0, p1, t) + cc(p0, p1, N - t); S12[t] = cc(p1, p2, t) + cc(p1, p2, N - t); S02[t] = cc(p0, p2, t) + cc(p0, p2, N - t); }
        for (int t = 1; t < 7; t++) F[t - 1] = paf[p0][t] + 2 * paf[p1][t] + 2 * paf[p2][t];
        for (int t = 0; t < 7; t++) F[6 + t] = S01[t] + S12[t] + paf[p2][t];
        for (int t = 0; t < 7; t++) F[13 + t] = S02[t] + S12[t] + paf[p1][t];
    }
}

/* ---- verification on Z_65 ---- */
static void to65(int r0, int r1, int r2, int *X) {
    int rr[3] = {r0, r1, r2};
    for (int g = 0; g < V; g++) { int a = g % 5, y = g % N; int ia = a <= 2 ? a : 5 - a; X[g] = val(rr[ia], y); }
}
static int verify(const int *A, const int *B, const int *D) {
    for (int s = 0; s < V; s++) {
        long t = 0;
        for (int g = 0; g < V; g++) t += A[g] * A[(g + s) % V] + 2 * B[g] * B[(g + s) % V] + D[g] * D[(g + s) % V];
        if (t != tgt65[s]) return 0;
    }
    return 1;
}

static long long nhits = 0;
static void report(const AEnt *a, const Seq *qb, const Ent *eb, const Seq *qd, int d0, int d1, int d2) {
    int A[80], B[80], D[80];
    to65(a->p[0], a->p[1], a->p[2], A);
    to65(qb->L[0][eb->i0], qb->L[1][eb->i1], qb->L[2][eb->i2], B);
    to65(qd->L[0][d0], qd->L[1][d1], qd->L[2][d2], D);
    int ok = verify(A, B, D);
    printf("HIT %s A ", ok ? "verified" : "WRONG");
    for (int g = 0; g < V; g++) putchar(A[g] > 0 ? '+' : '-');
    printf(" B "); for (int g = 0; g < V; g++) putchar(B[g] > 0 ? '+' : '-');
    printf(" D "); for (int g = 0; g < V; g++) putchar(D[g] > 0 ? '+' : '-');
    printf("\n"); fflush(stdout);
    if (ok) nhits++;
}

/* streaming contexts */
typedef struct { Store *hs; const Seq *hq; long long probes; int hashedB; const Seq *sq; } JCtx;

static void join_cb(const Seq *q, int i0, int i1, int i2, void *vctx) {
    JCtx *c = vctx;
    signed char F[NC], need[NC];
    classvals(q, i0, i1, i2, F);
    for (int ia = 0; ia < nA; ia++) {
        const signed char *FA = Alist[ia].F; int ok = 1;
        for (int k = 0; k < NC; k++) {
            int v;
            if (c->hashedB) { int s = target[k] - FA[k] - F[k]; if (s & 1) { ok = 0; break; } v = s / 2; }   /* streamed D */
            else v = target[k] - FA[k] - 2 * F[k];                                                  /* streamed B */
            if (v < -127 || v > 127) { ok = 0; break; }
            need[k] = (signed char)v;
        }
        if (!ok) continue;
        c->probes++;
        uint64_t h = hkey(need) & c->hs->mask;
        while (c->hs->tab[h] != 0xffffffffu) {
            const Ent *e = &c->hs->e[c->hs->tab[h]];
            if (!memcmp(e->F, need, NC)) {
                if (c->hashedB) report(&Alist[ia], c->hq, e, q, i0, i1, i2);
                else { Ent eb = {i0, i1, i2}; memcpy(eb.F, F, NC); report(&Alist[ia], q, &eb, c->hq, e->i0, e->i1, e->i2); }
            }
            h = (h + 1) & c->hs->mask;
        }
    }
}

static void classvals_rows(int x0, int x1, int x2, int *F) {
    int r[3] = {x0, x1, x2}, S[3][7];
    int pr[3][2] = {{0, 1}, {1, 2}, {0, 2}};
    for (int k = 0; k < 3; k++) {
        int a = r[pr[k][0]], b = r[pr[k][1]];
        S[k][0] = 2 * cc(a, b, 0);
        for (int t = 1; t < 7; t++) S[k][t] = cc(a, b, t) + cc(a, b, N - t);
    }
    for (int t = 1; t < 7; t++) F[t - 1] = paf[x0][t] + 2 * paf[x1][t] + 2 * paf[x2][t];
    for (int t = 0; t < 7; t++) F[6 + t] = S[0][t] + S[1][t] + paf[x2][t];
    for (int t = 0; t < 7; t++) F[13 + t] = S[2][t] + S[1][t] + paf[x1][t];
}

static int randrow(int S, int sym, uint64_t *st) {
    for (;;) {
        *st = mix(*st + 0x9e3779b97f4a7c15ULL);
        int m = (int)(*st & (M - 1));
        if (sym) { int o = 0; for (int y = 0; y < N; y++) if (m >> (y <= 6 ? y : N - y) & 1) o |= 1 << y; m = o; }
        if (rsum[m] == S) return m;
    }
}

int main(int argc, char **argv) {
    if (argc > 2 && !strcmp(argv[1], "-N")) { N = atoi(argv[2]); M = 1 << N; argv += 2; argc -= 2; }
    if (argc > 2 && !strcmp(argv[1], "-q")) { QB = atoi(argv[2]); argv += 2; argc -= 2; }
    V = 5 * N;
    if (N != 12 && N != 13) { fprintf(stderr, "N must be 12 or 13 (7 shift classes)\n"); return 1; }
    int test = argc > 1 && !strcmp(argv[1], "test");
    if (argc < 10 + test) { fprintf(stderr, "usage: p65 [-N n] [-q 5] [test] P0 P1 P2 R0 R1 R2 Q0 Q1 Q2 [seed]\n"); return 1; }
    int v[9]; for (int i = 0; i < 9; i++) v[i] = atoi(argv[i + 1 + test]);
    init();
    tgt65[0] = 4 * V;
    if (QB) {          /* Z_5-bordered (kernel of Z_V -> Z_5 is {x5 = 0}): -4q on it off 0, order 4(V + q) */
        tgt65[0] = 4 * V;
        for (int s = 1; s < V; s++) tgt65[s] = (s % 5 == 0) ? -4 * QB : 0;
        for (int k = 0; k < 6; k++) target[k] = -4 * QB;      /* classes (0, t), t = 1..6 */
    }
    if (test) {                         /* plant a random triple, target = its class values, PSD filters off */
        uint64_t st = argc > 11 ? strtoull(argv[11], 0, 10) : 1;
        int p[3], r[3], q[3], FA[NC], FB[NC], FD[NC];
        for (int i = 0; i < 3; i++) { p[i] = randrow(v[i], 1, &st); r[i] = randrow(v[3 + i], 0, &st); q[i] = randrow(v[6 + i], 0, &st); }
        classvals_rows(p[0], p[1], p[2], FA); classvals_rows(r[0], r[1], r[2], FB); classvals_rows(q[0], q[1], q[2], FD);
        for (int k = 0; k < NC; k++) target[k] = FA[k] + 2 * FB[k] + FD[k];
        int A[80], B[80], D[80];
        to65(p[0], p[1], p[2], A); to65(r[0], r[1], r[2], B); to65(q[0], q[1], q[2], D);
        for (int s = 0; s < V; s++) { long t = 0; for (int g = 0; g < V; g++) t += A[g] * A[(g + s) % V] + 2 * B[g] * B[(g + s) % V] + D[g] * D[(g + s) % V]; tgt65[s] = t; }
        boundA = boundB = boundD = 1e9; canonA = 0;   /* a planted target is not multiplier-invariant */
        printf("planted A %d %d %d B %d %d %d D %d %d %d\n", p[0], p[1], p[2], r[0], r[1], r[2], q[0], q[1], q[2]);
    }
    build_A(v[0], v[1], v[2]);
    Seq qb, qd;
    seq_init(&qb, v[3], v[4], v[5], boundB);
    seq_init(&qd, v[6], v[7], v[8], boundD);
    long long nb = seq_iter(&qb, NULL, NULL), nd = seq_iter(&qd, NULL, NULL);
    printf("A (%d,%d,%d) B (%d,%d,%d) D (%d,%d,%d): |A| %d, rows B %d/%d/%d D %d/%d/%d, PSD survivors B %lld D %lld\n",
           v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], nA, qb.n[0], qb.n[1], qb.n[2], qd.n[0], qd.n[1], qd.n[2], nb, nd);
    fflush(stdout);
    int hashB = nb <= nd;
    Store s = {0};
    const Seq *hq = hashB ? &qb : &qd, *sq = hashB ? &qd : &qb;
    seq_iter(hq, store_cb, &s);
    store_hash(&s);
    JCtx c = {&s, hq, 0, hashB, sq};
    seq_iter(sq, join_cb, &c);
    printf("DONE hashed %s (%lld), probes %lld, verified hits %lld\n", hashB ? "B" : "D", s.n, c.probes, nhits);
    return 0;
}
