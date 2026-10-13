/* Tabu search for cyclic sequence families with multiplier symmetry: blocks constant on the orbits of a multiplier
 * group H on Z_v (structure file from orbit_ls.py), sum_k w_k PAF_{X_k}(s) = T for s != 0.
 * Block types: f (free on the H-orbits), s (symmetric: orbit and mirror orbit equal), k (skew: X(0) = 1, mirror orbit
 * opposite). Energy: sum over the classes of nonzero shifts (orbits of <H, -1>) of |class| (P(s_c) - T)^2.
 * Moves: flip one orbit (free) or an orbit with its mirror (symmetric, skew). Tabu with aspiration; restart after R
 * steps without a new best.
 *
 *   orbit_ls FILE seed maxsec [R] [TT] [nsol]
 * FILE: "v nb T" / per block "type weight" / "O" / per orbit "size e1 e2 ..." / "mirror m_0 .. m_{O-1}" / "C" / per
 * class "size rep". Prints "SOL <X_1> ... <X_nb>" (+/- strings over Z_v) and "STAT ...".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

#define MAXB 6
static int v, nb, T, O, C;
static char typ[MAXB];
static int w[MAXB];
static int *X[MAXB];
static int *osz, **orb, *mir, *csz, *crep;
static long *P, *dP;
static int *stamp, stampv;

typedef struct { int k, o1, o2; } Move;      /* flip orbit o1 (and o2 >= 0) of block k */
static Move *mv; static int M;

static void full_paf(void) {
    for (int c = 0; c < C; c++) {
        int s = crep[c]; long t = 0;
        for (int k = 0; k < nb; k++) { long u = 0; for (int y = 0; y < v; y++) u += X[k][y] * X[k][(y + s) % v]; t += w[k] * u; }
        P[c] = t;
    }
}
static long energy(void) { long e = 0; for (int c = 0; c < C; c++) e += (long)csz[c] * (P[c] - T) * (P[c] - T); return e; }

static void delta(const Move *m, long *out) {
    int k = m->k, *a = X[k], wk = w[k];
    int n = 0, cells[4096];
    stampv++;
    for (int i = 0; i < osz[m->o1]; i++) { cells[n++] = orb[m->o1][i]; stamp[orb[m->o1][i]] = stampv; }
    if (m->o2 >= 0) for (int i = 0; i < osz[m->o2]; i++) { cells[n++] = orb[m->o2][i]; stamp[orb[m->o2][i]] = stampv; }
    for (int c = 0; c < C; c++) {
        int s = crep[c]; long f = 0;
        for (int i = 0; i < n; i++) {
            int x = cells[i];
            int xp = x + s; if (xp >= v) xp -= v;
            int xm = x - s; if (xm < 0) xm += v;
            f += -2 * a[x] * (a[xm] + a[xp]);
            if (stamp[xp] == stampv) f += 4 * a[x] * a[xp];
        }
        out[c] = wk * f;
    }
}

static void apply(const Move *m) {
    int *a = X[m->k];
    for (int i = 0; i < osz[m->o1]; i++) a[orb[m->o1][i]] = -a[orb[m->o1][i]];
    if (m->o2 >= 0) for (int i = 0; i < osz[m->o2]; i++) a[orb[m->o2][i]] = -a[orb[m->o2][i]];
}

static int zero_orbit;
static void randomize(void) {
    for (int k = 0; k < nb; k++) {
        int *a = X[k];
        for (int o = 0; o < O; o++) {
            if (typ[k] != 'f' && mir[o] < o) continue;
            int val = (rnd() & 1) ? 1 : -1;
            if (typ[k] == 'k' && o == zero_orbit) val = 1;
            for (int i = 0; i < osz[o]; i++) a[orb[o][i]] = val;
            if (typ[k] != 'f' && mir[o] != o) { int vm = typ[k] == 's' ? val : -val; for (int i = 0; i < osz[mir[o]]; i++) a[orb[mir[o]][i]] = vm; }
        }
    }
    full_paf();
}

static void put(const int *x) { for (int i = 0; i < v; i++) putchar(x[i] > 0 ? '+' : '-'); }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: orbit_ls FILE seed maxsec [R] [TT] [nsol]\n"); return 1; }
    FILE *fh = fopen(argv[1], "r"); if (!fh) { perror("open"); return 1; }
    if (fscanf(fh, "%d %d %d", &v, &nb, &T) != 3) return 1;
    for (int k = 0; k < nb; k++) { char t[4]; if (fscanf(fh, "%3s %d", t, &w[k]) != 2) return 1; typ[k] = t[0]; X[k] = calloc(v, sizeof(int)); }
    if (fscanf(fh, "%d", &O) != 1) return 1;
    osz = calloc(O, sizeof(int)); orb = calloc(O, sizeof(int *)); mir = calloc(O, sizeof(int));
    for (int o = 0; o < O; o++) {
        if (fscanf(fh, "%d", &osz[o]) != 1) return 1;
        orb[o] = calloc(osz[o], sizeof(int));
        for (int i = 0; i < osz[o]; i++) if (fscanf(fh, "%d", &orb[o][i]) != 1) return 1;
        if (osz[o] == 1 && orb[o][0] == 0) zero_orbit = o;
    }
    char word[16]; if (fscanf(fh, "%15s", word) != 1) return 1;
    for (int o = 0; o < O; o++) if (fscanf(fh, "%d", &mir[o]) != 1) return 1;
    if (fscanf(fh, "%d", &C) != 1) return 1;
    csz = calloc(C, sizeof(int)); crep = calloc(C, sizeof(int));
    for (int c = 0; c < C; c++) if (fscanf(fh, "%d %d", &csz[c], &crep[c]) != 2) return 1;
    fclose(fh);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[2]) * 0x2545F4914F6CDD1DULL; if (!rs) rs = 1;
    double maxsec = atof(argv[3]);
    long R = argc > 4 ? atol(argv[4]) : 20000;
    int TT = argc > 5 ? atoi(argv[5]) : 0;
    int nsol = argc > 6 ? atoi(argv[6]) : 1;
    mv = calloc(nb * O, sizeof(Move)); M = 0;
    for (int k = 0; k < nb; k++)
        for (int o = 0; o < O; o++) {
            if (typ[k] == 'f') { mv[M].k = k; mv[M].o1 = o; mv[M].o2 = -1; M++; continue; }
            if (mir[o] < o) continue;
            if (typ[k] == 'k' && o == zero_orbit) continue;
            if (typ[k] == 'k' && mir[o] == o) { fprintf(stderr, "skew block with a self-mirror orbit\n"); return 1; }
            mv[M].k = k; mv[M].o1 = o; mv[M].o2 = mir[o] == o ? -1 : mir[o]; M++;
        }
    if (TT <= 0) TT = M / 8 + 1;
    P = calloc(C, sizeof(long)); dP = calloc(C, sizeof(long));
    stamp = calloc(v, sizeof(int));
    long *tabu = calloc(M, sizeof(long)); int *cand = calloc(M, sizeof(int));
    clock_t t0 = clock();
    long step = 0, restarts = 0, found = 0, sincebest = 0, E, best, gbest;
    randomize(); E = energy(); best = E; gbest = E;
    while (found < nsol) {
        if ((step & 255) == 0 && (double)(clock() - t0) / CLOCKS_PER_SEC > maxsec) break;
        long bdE = 1L << 62; int nc = 0;
        for (int m = 0; m < M; m++) {
            delta(&mv[m], dP);
            long dE = 0;
            for (int c = 0; c < C; c++) dE += (long)csz[c] * dP[c] * (2L * (P[c] - T) + dP[c]);
            if (tabu[m] > step && E + dE >= best) continue;
            if (dE < bdE) { bdE = dE; nc = 0; cand[nc++] = m; }
            else if (dE == bdE) cand[nc++] = m;
        }
        if (!nc) { step++; continue; }
        int m = cand[rnd() % nc];
        delta(&mv[m], dP);
        for (int c = 0; c < C; c++) P[c] += dP[c];
        apply(&mv[m]);
        E += bdE;
        tabu[m] = step + TT + (long)(rnd() % (TT + 1));
        step++;
        if (E < best) { best = E; sincebest = 0; if (E < gbest) gbest = E; } else sincebest++;
        if (E == 0) {
            full_paf();
            if (energy() != 0) { fprintf(stderr, "internal error: energy mismatch\n"); return 2; }
            printf("SOL");
            for (int k = 0; k < nb; k++) { putchar(' '); put(X[k]); }
            printf(" step %ld restart %ld\n", step, restarts);
            fflush(stdout);
            found++;
            randomize(); E = energy(); best = E; sincebest = 0; restarts++;
            memset(tabu, 0, M * sizeof(long));
            continue;
        }
        if (sincebest > R) {
            randomize(); E = energy(); best = E; sincebest = 0; restarts++;
            memset(tabu, 0, M * sizeof(long));
        }
    }
    printf("STAT v %d O %d C %d moves %d steps %ld restarts %ld found %ld best %ld sec %.1f\n", v, O, C, M, step, restarts,
           found, gbest, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
