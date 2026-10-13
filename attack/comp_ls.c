/* Lift search for cyclic sequence families on Z_v, v = p q with gcd(p, q) = 1, with prescribed compressions.
 *
 * Every block X_k (+-1 sequence on Z_v = Z_p x Z_q, type f free, s symmetric, k skew) has its compressions to Z_p
 * (row sums of the p x q array) and to Z_q (column sums) fixed at values that satisfy the compressed equations; the
 * search moves only through arrays with these margins, by 2 x 2 switches (rows a1, a2, columns b1, b2 with an
 * alternating pattern; for symmetric and skew blocks together with the mirror switch), so the compressed equations
 * stay satisfied and the tabu search works on the remaining autocorrelation classes:
 *     E = sum_{s=1}^{(v-1)/2} (sum_k w_k PAF_{X_k}(s) - T)^2.
 *
 *   comp_ls FILE seed maxsec [R] [TT] [nsol]
 * FILE: "v p q nb T" then, per block, "type weight" followed by p row sums (index a = x mod p) and q column sums
 * (index b = x mod q). The initial arrays are random with these margins (a randomized greedy fill, then random
 * switches). Prints "SOL <X_1> ... <X_nb>" and "STAT ...".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

#define MAXB 6
static int v, p, q, h, nb, T;
static char typ[MAXB];
static int w[MAXB];
static int *X[MAXB];
static int *rowsum[MAXB], *colsum[MAXB];
static int *P, *dP;
static int *cell;                          /* cell[a * q + b] = x with x = a mod p, x = b mod q */

typedef struct { int k, n, u[8]; } Move;
static Move *mv; static int nmv, cap;

static void full_paf(void) {
    for (int s = 1; s <= h; s++) {
        long t = 0;
        for (int k = 0; k < nb; k++) {
            long u = 0;
            for (int y = 0; y < v; y++) u += X[k][y] * X[k][(y + s) % v];
            t += w[k] * u;
        }
        P[s] = (int)t;
    }
}
static long energy(void) { long e = 0; for (int s = 1; s <= h; s++) e += (long)(P[s] - T) * (P[s] - T); return e; }

static int fold(int s) { s %= v; if (s < 0) s += v; return s > h ? v - s : s; }

/* change of P when the cells u[0..n-1] of block k flip (exact, with the old values) */
static void delta(const Move *m, int *out) {
    int k = m->k, *a = X[k], wk = w[k];
    for (int s = 1; s <= h; s++) {
        int f = 0;
        for (int i = 0; i < m->n; i++) {
            int x = m->u[i];
            int xp = x + s; if (xp >= v) xp -= v;
            int xm = x - s; if (xm < 0) xm += v;
            f += -2 * a[x] * (a[xm] + a[xp]);
        }
        out[s] = wk * f;
    }
    for (int i = 0; i < m->n; i++)              /* second order: pairs of flipped cells at distance s */
        for (int j = 0; j < m->n; j++) {
            if (i == j) continue;
            int s = m->u[j] - m->u[i];
            int fs = fold(s);
            if (!fs) continue;
            /* sum_y d(y) d(y+s) counts the ordered pair (u_i, u_j) at shift s = u_j - u_i; PAF(s) = PAF(-s), and P
               stores one value per +-s, so each ordered pair contributes once to the folded shift */
            out[fs] += wk * 4 * a[m->u[i]] * a[m->u[j]] / 2;
        }
}

static void apply(const Move *m) { for (int i = 0; i < m->n; i++) X[m->k][m->u[i]] = -X[m->k][m->u[i]]; }

static int neg(int x) { return x ? v - x : 0; }

/* all currently valid switches */
static void moves(void) {
    nmv = 0;
    for (int k = 0; k < nb; k++) {
        int *a = X[k];
        for (int a1 = 0; a1 < p; a1++) for (int a2 = a1 + 1; a2 < p; a2++)
        for (int b1 = 0; b1 < q; b1++) for (int b2 = b1 + 1; b2 < q; b2++) {
            int x11 = cell[a1 * q + b1], x12 = cell[a1 * q + b2], x21 = cell[a2 * q + b1], x22 = cell[a2 * q + b2];
            if (!(a[x11] == a[x22] && a[x12] == a[x21] && a[x11] == -a[x12])) continue;
            Move m; m.k = k; m.n = 4; m.u[0] = x11; m.u[1] = x12; m.u[2] = x21; m.u[3] = x22;
            if (typ[k] != 'f') {
                int mir[4] = {neg(x11), neg(x12), neg(x21), neg(x22)};
                int same = 0, any = 0;
                for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) if (mir[i] == m.u[j]) any++;
                int self = (mir[0] == x22 && mir[3] == x11 && mir[1] == x21 && mir[2] == x12);
                if (self) {                              /* S = -S: flip S alone (symmetric blocks only) */
                    if (typ[k] == 'k') continue;
                    same = 1;
                } else if (any) continue;                /* partial overlap: skip */
                if (typ[k] == 'k') { int bad = 0; for (int i = 0; i < 4; i++) if (!m.u[i]) bad = 1; if (bad) continue; }
                if (!same) { for (int i = 0; i < 4; i++) m.u[4 + i] = mir[i]; m.n = 8; }
            }
            if (nmv == cap) { cap = cap ? 2 * cap : 4096; mv = realloc(mv, cap * sizeof(Move)); }
            mv[nmv++] = m;
        }
    }
}

/* margin error of block k */
static long merr(int k) {
    long e = 0;
    for (int a1 = 0; a1 < p; a1++) { int r = 0; for (int b = 0; b < q; b++) r += X[k][cell[a1 * q + b]]; e += (long)(r - rowsum[k][a1]) * (r - rowsum[k][a1]); }
    for (int b = 0; b < q; b++) { int c = 0; for (int a1 = 0; a1 < p; a1++) c += X[k][cell[a1 * q + b]]; e += (long)(c - colsum[k][b]) * (c - colsum[k][b]); }
    return e;
}

/* random array of the block's type with the given margins: local search on the margin error (flips of one cell, or
   of a mirror pair for symmetric and skew blocks), greedy with random ties and random kicks */
static void fix_margins(int k) {
    int *a = X[k];
    for (int tries = 0; tries < 200; tries++) {
        if (typ[k] == 'f') { for (int x = 0; x < v; x++) a[x] = (rnd() & 1) ? 1 : -1; }
        else if (typ[k] == 's') { for (int x = 0; x <= h; x++) { a[x] = (rnd() & 1) ? 1 : -1; if (x) a[v - x] = a[x]; } }
        else { a[0] = 1; for (int x = 1; x <= h; x++) { a[x] = (rnd() & 1) ? 1 : -1; a[v - x] = -a[x]; } }
        int R[64], C[256];
        for (int a1 = 0; a1 < p; a1++) { R[a1] = 0; for (int b = 0; b < q; b++) R[a1] += a[cell[a1 * q + b]]; }
        for (int b = 0; b < q; b++) { C[b] = 0; for (int a1 = 0; a1 < p; a1++) C[b] += a[cell[a1 * q + b]]; }
        long e = merr(k);
        int lo = typ[k] == 'f' ? 0 : (typ[k] == 's' ? 0 : 1), hi = typ[k] == 'f' ? v - 1 : h;
        for (long it = 0; it < 200000 && e > 0; it++) {
            long bd = 1L << 60; int bx[512], nbx = 0;
            for (int x = lo; x <= hi; x++) {
                int cells[2] = {x, (typ[k] == 'f' || x == 0) ? -1 : v - x};
                long d = 0;
                int dR[64] = {0}, dC[256] = {0};
                for (int i = 0; i < 2; i++) { if (cells[i] < 0) continue; int y = cells[i]; dR[y % p] -= 2 * a[y]; dC[y % q] -= 2 * a[y]; }
                for (int i = 0; i < 2; i++) {
                    if (cells[i] < 0) continue; int y = cells[i];
                    int ra = y % p, cb = y % q;
                    if (dR[ra]) { d += (long)(R[ra] + dR[ra] - rowsum[k][ra]) * (R[ra] + dR[ra] - rowsum[k][ra]) - (long)(R[ra] - rowsum[k][ra]) * (R[ra] - rowsum[k][ra]); dR[ra] = 0; }
                    if (dC[cb]) { d += (long)(C[cb] + dC[cb] - colsum[k][cb]) * (C[cb] + dC[cb] - colsum[k][cb]) - (long)(C[cb] - colsum[k][cb]) * (C[cb] - colsum[k][cb]); dC[cb] = 0; }
                }
                if (d < bd) { bd = d; nbx = 0; bx[nbx++] = x; } else if (d == bd && nbx < 512) bx[nbx++] = x;
            }
            int x = (bd >= 0 && (rnd() % 8) == 0) ? lo + (int)(rnd() % (hi - lo + 1)) : bx[rnd() % nbx];
            int cells[2] = {x, (typ[k] == 'f' || x == 0) ? -1 : v - x};
            for (int i = 0; i < 2; i++) { if (cells[i] < 0) continue; int y = cells[i]; R[y % p] -= 2 * a[y]; C[y % q] -= 2 * a[y]; a[y] = -a[y]; }
            e = merr(k);
        }
        if (e == 0) return;
    }
    fprintf(stderr, "no array of type %c with the margins of block %d\n", typ[k], k); exit(3);
}

static void randomize(void) {
    for (int k = 0; k < nb; k++) fix_margins(k);
    for (int r = 0; r < 200; r++) { moves(); if (!nmv) break; apply(&mv[rnd() % nmv]); }
    full_paf();
}

static void put(const int *x) { for (int i = 0; i < v; i++) putchar(x[i] > 0 ? '+' : '-'); }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: comp_ls FILE seed maxsec [R] [TT] [nsol]\n"); return 1; }
    FILE *fh = fopen(argv[1], "r"); if (!fh) { perror("open"); return 1; }
    if (fscanf(fh, "%d %d %d %d %d", &v, &p, &q, &nb, &T) != 5) return 1;
    h = (v - 1) / 2;
    for (int k = 0; k < nb; k++) {
        char t[4]; if (fscanf(fh, "%3s %d", t, &w[k]) != 2) return 1; typ[k] = t[0];
        rowsum[k] = calloc(p, sizeof(int)); colsum[k] = calloc(q, sizeof(int)); X[k] = calloc(v, sizeof(int));
        for (int a = 0; a < p; a++) if (fscanf(fh, "%d", &rowsum[k][a]) != 1) return 1;
        for (int b = 0; b < q; b++) if (fscanf(fh, "%d", &colsum[k][b]) != 1) return 1;
    }
    fclose(fh);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[2]) * 0x2545F4914F6CDD1DULL; if (!rs) rs = 1;
    double maxsec = atof(argv[3]);
    long R = argc > 4 ? atol(argv[4]) : 20000;
    int TT = argc > 5 ? atoi(argv[5]) : 0;
    int nsol = argc > 6 ? atoi(argv[6]) : 1;
    cell = calloc(p * q, sizeof(int));
    for (int x = 0; x < v; x++) cell[(x % p) * q + (x % q)] = x;
    P = calloc(h + 1, sizeof(int)); dP = calloc(h + 1, sizeof(int));
    /* tabu on cells: a move is tabu if any of its cells was flipped within the tenure */
    long *tabu = calloc(nb * v, sizeof(long));
    clock_t t0 = clock();
    long step = 0, restarts = 0, found = 0, sincebest = 0, E, best;
    randomize(); E = energy(); best = E;
    if (TT <= 0) TT = 8;
    int *cand = NULL; int ccap = 0;
    while (found < nsol) {
        if ((step & 63) == 0 && (double)(clock() - t0) / CLOCKS_PER_SEC > maxsec) break;
        moves();
        if (ccap < nmv) { ccap = nmv; cand = realloc(cand, ccap * sizeof(int)); }
        long bdE = 1L << 62; int nc = 0;
        for (int i = 0; i < nmv; i++) {
            delta(&mv[i], dP);
            long dE = 0;
            for (int s = 1; s <= h; s++) dE += (long)dP[s] * (2L * (P[s] - T) + dP[s]);
            int tb = 0;
            for (int j = 0; j < mv[i].n; j++) if (tabu[mv[i].k * v + mv[i].u[j]] > step) { tb = 1; break; }
            if (tb && E + dE >= best) continue;
            if (dE < bdE) { bdE = dE; nc = 0; cand[nc++] = i; }
            else if (dE == bdE) cand[nc++] = i;
        }
        if (!nc) { step++; sincebest++; if (sincebest > R) goto restart; continue; }
        {
            Move m = mv[cand[rnd() % nc]];
            delta(&m, dP);
            for (int s = 1; s <= h; s++) P[s] += dP[s];
            apply(&m);
            E += bdE;
            long until = step + TT + (long)(rnd() % (TT + 1));
            for (int j = 0; j < m.n; j++) tabu[m.k * v + m.u[j]] = until;
        }
        step++;
        if (E < best) { best = E; sincebest = 0; } else sincebest++;
        if (E == 0) {
            full_paf();
            if (energy() != 0) { fprintf(stderr, "internal error: energy mismatch\n"); return 2; }
            printf("SOL");
            for (int k = 0; k < nb; k++) { putchar(' '); put(X[k]); }
            printf(" step %ld restart %ld\n", step, restarts);
            fflush(stdout);
            found++;
            goto restart;
        }
        if (sincebest > R) goto restart;
        continue;
restart:
        randomize(); E = energy(); best = E; sincebest = 0; restarts++;
        memset(tabu, 0, nb * v * sizeof(long));
    }
    printf("STAT v %d steps %ld restarts %ld found %ld best %ld sec %.1f\n", v, step, restarts, found, best,
           (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
