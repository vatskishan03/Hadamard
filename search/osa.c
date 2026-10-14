/* osa.c: simulated annealing over orbit variables for the c (or y) parts of the orthogonal designs of
 * search/odm_g.py (Goethals-Seidel, Theorem 7) and search/symod_g.py (Balonin), given the b (x) parts.
 *
 * usage: osa problem.bin seed restarts steps T0 T1 out.txt
 *
 * problem.bin (int32): nb, npr, then for each block: mult, nv, target, vsize[nv], Dg[nv][npr], Cx[nv][nv][npr],
 *   Lv[nv][npr]; then wk[npr] (the number of shifts represented by each representative).
 * Block b with values val in {+-1}^nv has row sum sum_j val_j vsize_j (held at target: moves flip one variable of
 * size 0 or swap two variables of equal size and opposite values), PAF_b[k] = sum_j Dg[j][k]
 * + sum_{i<j} val_i val_j Cx[j][i][k] and X_b[k] = sum_j val_j Lv[j][k].  The energy is
 *   E = sum_k wk[k] ((sum_b mult_b PAF_b[k])^2 + (sum_b mult_b X_b[k])^2),
 * zero exactly for a design.  out.txt receives "SOL" and the values of every block when E = 0.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXB 4
#define MAXV 160
#define MAXR 128
static int nb, npr, mult[MAXB], nv[MAXB], target[MAXB], vsize[MAXB][MAXV];
static int *Dg[MAXB], *Cx[MAXB], *Lv[MAXB];
static int wk[MAXR];
static int val[MAXB][MAXV];
static long P[MAXR], X[MAXR];          /* totals over blocks (with multiplicities) */
static uint64_t rng;

static inline uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static inline double ur(void) { return (xr() >> 11) * (1.0 / 9007199254740992.0); }
#define CX(b, j, i, k) Cx[b][((size_t)(j) * nv[b] + (i)) * npr + (k)]
#define DG(b, j, k) Dg[b][(j) * npr + (k)]
#define LV(b, j, k) Lv[b][(j) * npr + (k)]

static int rd(FILE *f) { int32_t x; if (fread(&x, 4, 1, f) != 1) { fprintf(stderr, "short problem file\n"); exit(1); } return x; }

static void totals(void) {
    memset(P, 0, sizeof P); memset(X, 0, sizeof X);
    for (int b = 0; b < nb; b++) for (int k = 0; k < npr; k++) {
        long p = 0, x = 0;
        for (int j = 0; j < nv[b]; j++) {
            p += DG(b, j, k);
            for (int i = 0; i < j; i++) p += (long)val[b][i] * val[b][j] * CX(b, j, i, k);
            x += (long)val[b][j] * LV(b, j, k);
        }
        P[k] += mult[b] * p; X[k] += mult[b] * x;
    }
}

static long energy_of(const long *p, const long *x) {
    long e = 0;
    for (int k = 0; k < npr; k++) e += (long)wk[k] * (p[k] * p[k] + x[k] * x[k]);
    return e;
}

/* change of the block-b PAF and X vectors when the variables in J (n of them) are negated */
static void delta(int b, const int *J, int n, long *dp, long *dx) {
    for (int k = 0; k < npr; k++) { dp[k] = 0; dx[k] = 0; }
    int inJ[MAXV]; memset(inJ, 0, sizeof(int) * nv[b]);
    for (int a = 0; a < n; a++) inJ[J[a]] = 1;
    for (int a = 0; a < n; a++) {
        int j = J[a];
        for (int k = 0; k < npr; k++) dx[k] += -2L * val[b][j] * LV(b, j, k);
        for (int i = 0; i < nv[b]; i++) {
            if (i == j || inJ[i]) continue;                 /* pairs inside J keep their product */
            int hi = i > j ? i : j, lo = i > j ? j : i;
            for (int k = 0; k < npr; k++) dp[k] += -2L * val[b][i] * val[b][j] * CX(b, hi, lo, k);
        }
    }
    for (int k = 0; k < npr; k++) { dp[k] *= mult[b]; dx[k] *= mult[b]; }
}

static void init_block(int b) {
    /* random values with the target row sum: greedy on a random order, then repair */
    for (int tries = 0; tries < 100000; tries++) {
        int rs = 0;
        for (int j = 0; j < nv[b]; j++) { val[b][j] = (xr() & 1) ? 1 : -1; rs += val[b][j] * vsize[b][j]; }
        for (int it = 0; it < 10 * nv[b] && rs != target[b]; it++) {
            int j = xr() % nv[b];
            int nrs = rs - 2 * val[b][j] * vsize[b][j];
            if (abs(nrs - target[b]) < abs(rs - target[b])) { val[b][j] = -val[b][j]; rs = nrs; }
        }
        if (rs == target[b]) return;
    }
    fprintf(stderr, "block %d: row sum %d not reached\n", b, target[b]);
    exit(2);
}

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage: osa problem.bin seed restarts steps T0 T1 out.txt\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    nb = rd(f); npr = rd(f);
    for (int b = 0; b < nb; b++) {
        mult[b] = rd(f); nv[b] = rd(f); target[b] = rd(f);
        if (nv[b] > MAXV) { fprintf(stderr, "too many variables\n"); return 1; }
        for (int j = 0; j < nv[b]; j++) vsize[b][j] = rd(f);
        Dg[b] = malloc(sizeof(int) * nv[b] * npr); Cx[b] = malloc(sizeof(int) * (size_t)nv[b] * nv[b] * npr);
        Lv[b] = malloc(sizeof(int) * nv[b] * npr);
        for (int i = 0; i < nv[b] * npr; i++) Dg[b][i] = rd(f);
        for (size_t i = 0; i < (size_t)nv[b] * nv[b] * npr; i++) Cx[b][i] = rd(f);
        for (int i = 0; i < nv[b] * npr; i++) Lv[b][i] = rd(f);
    }
    for (int k = 0; k < npr; k++) wk[k] = rd(f);
    fclose(f);
    rng = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[2]) * 0x2545F4914F6CDD1DULL;
    int restarts = atoi(argv[3]); long steps = atol(argv[4]);
    double T0 = atof(argv[5]), T1 = atof(argv[6]);
    long best_all = -1;
    for (int r = 0; r < restarts; r++) {
        for (int b = 0; b < nb; b++) init_block(b);
        totals();
        long E = energy_of(P, X), best = E;
        long dp[MAXR], dx[MAXR], np_[MAXR], nx[MAXR];
        for (long s = 0; s < steps && E > 0; s++) {
            double T = T0 * pow(T1 / T0, (double)s / steps);
            int b = xr() % nb;
            int J[2], n;
            int j1 = xr() % nv[b];
            if (vsize[b][j1] == 0) { J[0] = j1; n = 1; }
            else {
                int j2 = -1;
                for (int tries = 0; tries < 8; tries++) {
                    int c = xr() % nv[b];
                    if (c != j1 && vsize[b][c] == vsize[b][j1] && val[b][c] != val[b][j1]) { j2 = c; break; }
                }
                if (j2 < 0) continue;
                J[0] = j1; J[1] = j2; n = 2;
            }
            delta(b, J, n, dp, dx);
            for (int k = 0; k < npr; k++) { np_[k] = P[k] + dp[k]; nx[k] = X[k] + dx[k]; }
            long E2 = energy_of(np_, nx);
            if (E2 <= E || ur() < exp((E - E2) / T)) {
                for (int a = 0; a < n; a++) val[b][J[a]] = -val[b][J[a]];
                memcpy(P, np_, sizeof(long) * npr); memcpy(X, nx, sizeof(long) * npr);
                E = E2;
                if (E < best) best = E;
            }
        }
        if (best_all < 0 || best < best_all) best_all = best;
        fprintf(stderr, "restart %d: final E %ld, best %ld\n", r, E, best);
        if (E == 0) {
            totals();
            if (energy_of(P, X) != 0) { fprintf(stderr, "internal error: recomputed energy nonzero\n"); return 3; }
            FILE *o = fopen(argv[7], "w");
            fprintf(o, "SOL restart %d\n", r);
            for (int b = 0; b < nb; b++) {
                for (int j = 0; j < nv[b]; j++) fprintf(o, "%d%c", val[b][j], j == nv[b] - 1 ? '\n' : ' ');
            }
            fclose(o);
            printf("SOL restart %d\n", r);
            return 0;
        }
    }
    printf("NONE best %ld\n", best_all);
    return 0;
}
