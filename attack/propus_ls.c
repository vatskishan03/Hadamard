/* Tabu search for cyclic propus triples without multiplier symmetry: +-1 sequences A, B, D of odd length v with
 * A symmetric (A(-x) = A(x)) and PAF_A + 2 PAF_B + PAF_D = 4v delta, which give a symmetric Hadamard matrix of
 * order 4v in the array of verify.py (construction "propus": blocks A, B, B, D).
 *
 * Energy: E = sum_{s=1}^{(v-1)/2} P(s)^2 with P = PAF_A + 2 PAF_B + PAF_D. Moves: flip one entry of B or D, or one
 * symmetric pair {x, -x} of A (or A(0)). Every step evaluates all moves exactly (O(v^2)), takes the best one that is
 * not tabu (a move is tabu for TT..2TT steps after it is made, unless it gives a new best energy), and restarts from
 * a random point after R steps without a new best.
 *
 *   propus_ls v seed maxsec [R] [TT] [nsol]
 * prints "SOL <A> <B> <D>" (strings of + and -) for each solution found, "STAT ..." at the end.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

static int v, h;
static int *a, *b, *d, *P, *dP;

static void full_paf(void) {
    for (int s = 1; s <= h; s++) {
        long t = 0;
        for (int y = 0; y < v; y++) {
            int z = (y + s) % v;
            t += a[y] * a[z] + 2 * b[y] * b[z] + d[y] * d[z];
        }
        P[s] = (int)t;
    }
}

static long energy(void) { long e = 0; for (int s = 1; s <= h; s++) e += (long)P[s] * P[s]; return e; }

/* change of P for move m (0..h: A pair x = m; h+1 .. h+v: B entry; h+v+1 .. h+2v: D entry) */
static void delta(int m, int *out) {
    if (m <= h) {
        int x = m, al = a[x];
        if (x == 0) {
            for (int s = 1; s <= h; s++) out[s] = -4 * al * a[s];
        } else {
            int t2 = (2 * x) % v; if (t2 > h) t2 = v - t2;
            for (int s = 1; s <= h; s++) {
                int p = x + s; if (p >= v) p -= v;
                int q = x - s; if (q < 0) q += v;
                out[s] = -4 * al * (a[p] + a[q]) + (s == t2 ? 4 : 0);
            }
        }
    } else if (m <= h + v) {
        int x = m - h - 1, bx = b[x];
        for (int s = 1; s <= h; s++) {
            int p = x + s; if (p >= v) p -= v;
            int q = x - s; if (q < 0) q += v;
            out[s] = -4 * bx * (b[p] + b[q]);
        }
    } else {
        int x = m - h - v - 1, dx = d[x];
        for (int s = 1; s <= h; s++) {
            int p = x + s; if (p >= v) p -= v;
            int q = x - s; if (q < 0) q += v;
            out[s] = -2 * dx * (d[p] + d[q]);
        }
    }
}

static void apply(int m) {
    if (m <= h) { a[m] = -a[m]; if (m) a[v - m] = -a[v - m]; }
    else if (m <= h + v) b[m - h - 1] = -b[m - h - 1];
    else d[m - h - v - 1] = -d[m - h - v - 1];
}

static void randomize(void) {
    for (int x = 0; x <= h; x++) { a[x] = (rnd() & 1) ? 1 : -1; if (x) a[v - x] = a[x]; }
    for (int x = 0; x < v; x++) { b[x] = (rnd() & 1) ? 1 : -1; d[x] = (rnd() & 1) ? 1 : -1; }
    full_paf();
}

static void put(const int *x) { for (int i = 0; i < v; i++) putchar(x[i] > 0 ? '+' : '-'); }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: propus_ls v seed maxsec [R] [TT] [nsol]\n"); return 1; }
    v = atoi(argv[1]); rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[2]) * 0x2545F4914F6CDD1DULL; if (!rs) rs = 1;
    double maxsec = atof(argv[3]);
    long R = argc > 4 ? atol(argv[4]) : 20000;
    int TT = argc > 5 ? atoi(argv[5]) : 0;
    int nsol = argc > 6 ? atoi(argv[6]) : 1;
    if (v < 3 || v % 2 == 0) { fprintf(stderr, "v must be odd\n"); return 1; }
    h = (v - 1) / 2;
    int M = h + 1 + 2 * v;
    if (TT <= 0) TT = M / 8 + 1;
    a = calloc(v, sizeof(int)); b = calloc(v, sizeof(int)); d = calloc(v, sizeof(int));
    P = calloc(h + 1, sizeof(int)); dP = calloc(h + 1, sizeof(int));
    long *tabu = calloc(M, sizeof(long));
    int *cand = calloc(M, sizeof(int));
    clock_t t0 = clock();
    long step = 0, restarts = 0, found = 0, sincebest = 0;
    long E, best;
    randomize(); E = energy(); best = E;
    while (found < nsol) {
        if ((step & 1023) == 0 && (double)(clock() - t0) / CLOCKS_PER_SEC > maxsec) break;
        long bdE = 1L << 62; int nc = 0;
        for (int m = 0; m < M; m++) {
            delta(m, dP);
            long dE = 0;
            for (int s = 1; s <= h; s++) dE += (long)dP[s] * (2L * P[s] + dP[s]);
            if (tabu[m] > step && E + dE >= best) continue;
            if (dE < bdE) { bdE = dE; nc = 0; cand[nc++] = m; }
            else if (dE == bdE) cand[nc++] = m;
        }
        if (!nc) { step++; continue; }
        int m = cand[rnd() % nc];
        delta(m, dP);
        for (int s = 1; s <= h; s++) P[s] += dP[s];
        apply(m);
        E += bdE;
        tabu[m] = step + TT + (long)(rnd() % (TT + 1));
        step++;
        if (E < best) { best = E; sincebest = 0; } else sincebest++;
        if (E == 0) {
            full_paf();
            if (energy() != 0) { fprintf(stderr, "internal error: energy mismatch\n"); return 2; }
            printf("SOL "); put(a); putchar(' '); put(b); putchar(' '); put(d); printf(" step %ld restart %ld\n", step, restarts);
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
    printf("STAT v %d steps %ld restarts %ld found %ld sec %.1f\n", v, step, restarts, found,
           (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
