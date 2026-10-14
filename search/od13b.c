/* od13b.c: simulated annealing for a skew-type OD(4t; 1, 3, 4t - 4) in Goethals-Seidel form (see od13.c), with
 * the row sums of the x3 parts held fixed.
 *
 * Blocks: A1 = x1 I + x3 W0 (W0 skew), Ai = x2 I + x3 Wi (i = 1, 2, 3; Wi[0] = 0, +-1 elsewhere).  Write
 * Wi = Si + Ki with Si symmetric and Ki skew; on each pair {s, -s} exactly one of Si(s), Ki(s) is non-zero.  The
 * x2 x3 condition is S1 + S2 + S3 = 0, so the vector (S1(s), S2(s), S3(s)) of a pair is (0, 0, 0) or a permutation of
 * (1, -1, 0).  Row sums: ri = 2 sum_pairs Si; the PAF condition forces r1^2 + r2^2 + r3^2 = 4t - 4 together with
 * r1 + r2 + r3 = 0, i.e. a^2 + ab + b^2 = (t - 1)/2 for (r1, r2) = (2a, 2b).  The search keeps (sum_pairs Si) equal to
 * (a, b, -a - b), given on the command line (default: the first solution found).
 *
 * Moves: flip a pair of W0; flip the skew sign of a zero-S entry; swap the configurations of two pairs; move a unit
 * of S between two pairs (keeping both valid).
 *
 * Usage: od13b t seed [restarts] [steps] [T0] [T1] [a b]
 * Output: "SOL" and the four blocks as integer codes (1 = x1, 2 = x2, +-3 = +-x3).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

static int t, h;
static int W[4][512];
static long paf[4][512];
static int S[3][256];      /* S[i][p] for pair p (1..h) */
static int K[3][256];      /* skew sign of the entry at p when S = 0 (W[i][p] = K, W[i][t-p] = -K) */

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static void full_paf(int i) {
    for (int s = 0; s < t; s++) {
        long a = 0;
        for (int g = 0; g < t; g++) a += W[i][g] * W[i][(g + s) % t];
        paf[i][s] = a;
    }
}

static void setv(int i, int a, int v) {
    int d = v - W[i][a];
    if (!d) return;
    for (int s = 1; s < t; s++) paf[i][s] += (long)d * (W[i][(a + s) % t] + W[i][(a - s + t) % t]);
    W[i][a] = v;
}

static long energy(void) {
    long e = 0;
    for (int s = 1; s <= h; s++) {
        long z = paf[0][s] + paf[1][s] + paf[2][s] + paf[3][s];
        e += z * z;
    }
    return e;
}

/* write the entries of sequence i (1..3) at pair p from S/K */
static void sync(int i, int p) {
    int s = S[i - 1][p], k = K[i - 1][p];
    if (s) { setv(i, p, s); setv(i, t - p, s); }
    else { setv(i, p, k); setv(i, t - p, -k); }
}

static int valid(int p) {
    int z = 0, sum = 0;
    for (int i = 0; i < 3; i++) { if (S[i][p] == 0) z++; sum += S[i][p]; if (S[i][p] > 1 || S[i][p] < -1) return 0; }
    return sum == 0 && (z == 3 || z == 1);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: od13b t seed [restarts] [steps] [T0] [T1] [a b]\n"); return 1; }
    t = atoi(argv[1]);
    rs = strtoull(argv[2], 0, 10) * 0x9E3779B97F4A7C15ULL + 7777;
    long maxr = argc > 3 ? atol(argv[3]) : 1000000;
    long steps = argc > 4 ? atol(argv[4]) : 5000000;
    double T0 = argc > 5 ? atof(argv[5]) : 60.0, T1 = argc > 6 ? atof(argv[6]) : 2.0;
    if (t % 2 == 0 || t < 3 || t > 511) { fprintf(stderr, "t must be odd, 3..511\n"); return 1; }
    h = (t - 1) / 2;
    int A = 0, B = 0, found = 0;
    if (argc > 8) { A = atoi(argv[7]); B = atoi(argv[8]); found = (A * A + A * B + B * B == h); }
    else {
        for (int a = 0; a <= 64 && !found; a++)
            for (int b = a; b <= 64 && !found; b++)
                if (a * a + a * b + b * b == h) { A = a; B = b; found = 1; }
    }
    if (!found) { fprintf(stderr, "no a, b with a^2 + ab + b^2 = %d: no such design in this form\n", h); return 3; }
    int target[3] = {A, B, -A - B};
    fprintf(stderr, "t=%d targets (sum S) = (%d, %d, %d), row sums (%d, %d, %d)\n", t, A, B, -A - B, 2 * A, 2 * B,
            -2 * A - 2 * B);
    for (long r = 0; r < maxr; r++) {
        memset(W, 0, sizeof W);
        memset(S, 0, sizeof S);
        /* initial S: place +-1 units greedily on random pairs */
        int need[3] = {target[0], target[1], target[2]};
        for (int p = 1; p <= h; p++) for (int i = 0; i < 3; i++) K[i][p] = (rnd() & 1) ? 1 : -1;
        int guard = 0;
        while ((need[0] || need[1] || need[2]) && guard++ < 100000) {
            int ip = -1, im = -1;
            for (int i = 0; i < 3; i++) { if (need[i] > 0 && (ip < 0 || rnd() & 1)) ip = i; }
            for (int i = 0; i < 3; i++) { if (need[i] < 0 && (im < 0 || rnd() & 1)) im = i; }
            if (ip < 0 || im < 0) break;
            int p = 1 + rnd() % h;
            if (S[0][p] || S[1][p] || S[2][p]) continue;
            S[ip][p] = 1; S[im][p] = -1; need[ip]--; need[im]++;
        }
        if (need[0] || need[1] || need[2]) { fprintf(stderr, "init failed\n"); return 4; }
        for (int p = 1; p <= h; p++) {
            int sg = (rnd() & 1) ? 1 : -1;
            W[0][p] = sg; W[0][t - p] = -sg;
            for (int i = 1; i <= 3; i++) {
                int s = S[i - 1][p], k = K[i - 1][p];
                if (s) { W[i][p] = s; W[i][t - p] = s; } else { W[i][p] = k; W[i][t - p] = -k; }
            }
        }
        for (int i = 0; i < 4; i++) full_paf(i);
        long E = energy(), best = E;
        double T = T0;
        for (long st = 0; st < steps && E; st++) {
            if (!(st & 1023)) T = T0 * pow(T1 / T0, (double)st / steps);
            int mv = rnd() % 8;
            long E2;
            if (mv < 2) {                    /* flip W0 pair */
                int p = 1 + rnd() % h, a = W[0][p];
                setv(0, p, -a); setv(0, t - p, a);
                E2 = energy();
                if (E2 <= E || urand() < exp((E - E2) / T)) E = E2;
                else { setv(0, p, a); setv(0, t - p, -a); }
            } else if (mv < 5) {             /* flip a skew sign */
                int p = 1 + rnd() % h, i = rnd() % 3;
                if (S[i][p]) continue;
                K[i][p] = -K[i][p]; sync(i + 1, p);
                E2 = energy();
                if (E2 <= E || urand() < exp((E - E2) / T)) E = E2;
                else { K[i][p] = -K[i][p]; sync(i + 1, p); }
            } else if (mv < 6) {             /* swap two pairs */
                int p = 1 + rnd() % h, q = 1 + rnd() % h;
                if (p == q) continue;
                int sp[3], kp[3];
                for (int i = 0; i < 3; i++) { sp[i] = S[i][p]; kp[i] = K[i][p]; S[i][p] = S[i][q]; K[i][p] = K[i][q]; S[i][q] = sp[i]; K[i][q] = kp[i]; }
                for (int i = 1; i <= 3; i++) { sync(i, p); sync(i, q); }
                E2 = energy();
                if (E2 <= E || urand() < exp((E - E2) / T)) E = E2;
                else {
                    for (int i = 0; i < 3; i++) { int a = S[i][p], b = K[i][p]; S[i][p] = S[i][q]; K[i][p] = K[i][q]; S[i][q] = a; K[i][q] = b; }
                    for (int i = 1; i <= 3; i++) { sync(i, p); sync(i, q); }
                }
            } else {                         /* move one unit of S from pair p to pair q in sequences i, j */
                int p = 1 + rnd() % h, q = 1 + rnd() % h, i = rnd() % 3, j = rnd() % 3;
                if (p == q || i == j) continue;
                /* at p: S_i -= 1, S_j += 1; at q: S_i += 1, S_j -= 1 */
                int os[2][3];
                for (int c = 0; c < 3; c++) { os[0][c] = S[c][p]; os[1][c] = S[c][q]; }
                S[i][p]--; S[j][p]++; S[i][q]++; S[j][q]--;
                if (!valid(p) || !valid(q)) {
                    for (int c = 0; c < 3; c++) { S[c][p] = os[0][c]; S[c][q] = os[1][c]; }
                    continue;
                }
                for (int c = 1; c <= 3; c++) { sync(c, p); sync(c, q); }
                E2 = energy();
                if (E2 <= E || urand() < exp((E - E2) / T)) E = E2;
                else {
                    for (int c = 0; c < 3; c++) { S[c][p] = os[0][c]; S[c][q] = os[1][c]; }
                    for (int c = 1; c <= 3; c++) { sync(c, p); sync(c, q); }
                }
            }
            if (E < best) best = E;
        }
        if (!E) {
            /* final check from scratch */
            for (int i = 0; i < 4; i++) full_paf(i);
            if (energy()) { fprintf(stderr, "internal error: energy mismatch\n"); return 5; }
            printf("SOL t=%d restart=%ld\n", t, r);
            for (int i = 0; i < 4; i++)
                for (int g = 0; g < t; g++) {
                    int c = g == 0 ? (i == 0 ? 1 : 2) : 3 * W[i][g];
                    printf("%d%c", c, g == t - 1 ? '\n' : ' ');
                }
            fflush(stdout);
            return 0;
        }
        if (r % 5 == 0) fprintf(stderr, "restart %ld final E=%ld best=%ld\n", r, E, best);
    }
    return 2;
}
