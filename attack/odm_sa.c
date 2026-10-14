/* odm_sa.c: simulated annealing for a skew-type OD(4t; 1, m, 4t - m - 1) in Goethals-Seidel form with a given b part.
 *
 * Blocks (circulant, length t): A1 = a e0 + b u1 + c w1, Ai = b ui + c wi (i = 2, 3, 4).  The b part (u1..u4, total
 * weight m) is given on the command line and must satisfy sum_i PAF(ui) = m delta, u1 skew (checked).  The c parts wi
 * are +-1 off the supports of the a and b parts, w1 skew.  The design is an OD iff
 *   (cc) sum_i PAF(wi)(s) = 0 for s != 0,
 *   (bc) X(s) = sum_i sum_a ui(a) [wi(a + s) + wi(a - s)] = 0 for all s.
 * The row sums of w2, w3, w4 are held at the given targets (moves: flip a +1 and a -1 of one block, or a pair of the
 * skew block w1).  Energy = sum_{s=1..h} cc(s)^2 + lambda sum_{s=1..h} X(s)^2.
 *
 * Usage: odm_sa t seed restarts steps T0 T1 lambda "i:pos:sign,..." r2 r3 r4
 *   e.g. odm_sa 25 1 1000 2000000 40 1 2 "2:0:1,2:5:1,3:0:1,3:5:-1,4:0:1" 3 7 -6
 * Output: "SOL" and the four blocks as codes (1 = a, +-2 = +-b, +-3 = +-c).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define MAXT 1024
static int t, h;
static int U[4][MAXT], W[4][MAXT], fixedpos[4][MAXT];
static long paf[4][MAXT], X[MAXT];
static int usup[4][64], nus[4];
static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static void setv(int i, int g, int v) {
    int d = v - W[i][g];
    if (!d) return;
    for (int s = 1; s < t; s++) paf[i][s] += (long)d * (W[i][(g + s) % t] + W[i][(g - s + t) % t]);
    W[i][g] = v;
    for (int k = 0; k < nus[i]; k++) {
        int a = usup[i][k], ua = U[i][a];
        X[((g - a) % t + t) % t] += (long)ua * d;
        X[((a - g) % t + t) % t] += (long)ua * d;
    }
}

static double lam;
static double energy(void) {
    double e = 0;
    for (int s = 1; s <= h; s++) {
        long z = paf[0][s] + paf[1][s] + paf[2][s] + paf[3][s];
        e += (double)z * z + lam * (double)X[s] * X[s];
    }
    return e;
}

int main(int argc, char **argv) {
    if (argc < 12) { fprintf(stderr, "usage: odm_sa t seed restarts steps T0 T1 lambda template r2 r3 r4\n"); return 1; }
    t = atoi(argv[1]); h = (t - 1) / 2;
    rs = strtoull(argv[2], 0, 10) * 0x9E3779B97F4A7C15ULL + 99991;
    long maxr = atol(argv[3]), steps = atol(argv[4]);
    double T0 = atof(argv[5]), T1 = atof(argv[6]);
    lam = atof(argv[7]);
    char *tp = strdup(argv[8]);
    int tgt[4] = {0, atoi(argv[9]), atoi(argv[10]), atoi(argv[11])};
    if (t % 2 == 0 || t > MAXT) { fprintf(stderr, "t must be odd and <= %d\n", MAXT); return 1; }
    memset(U, 0, sizeof U);
    for (char *tok = strtok(tp, ","); tok; tok = strtok(0, ",")) {
        int i, p, sg;
        if (sscanf(tok, "%d:%d:%d", &i, &p, &sg) != 3 || i < 1 || i > 4) { fprintf(stderr, "bad template %s\n", tok); return 1; }
        U[i - 1][((p % t) + t) % t] = sg;
    }
    int m = 0;
    for (int i = 0; i < 4; i++) {
        nus[i] = 0;
        for (int g = 0; g < t; g++) if (U[i][g]) { usup[i][nus[i]++] = g; m++; }
    }
    /* check b part: sum PAF = m delta, u1 skew, u1(0) = 0 */
    for (int s = 1; s < t; s++) {
        long z = 0;
        for (int i = 0; i < 4; i++) for (int g = 0; g < t; g++) z += U[i][g] * U[i][(g + s) % t];
        if (z) { fprintf(stderr, "b part: PAF sum at %d is %ld\n", s, z); return 1; }
    }
    if (U[0][0]) { fprintf(stderr, "u1(0) must be 0 (a sits there)\n"); return 1; }
    for (int g = 1; g < t; g++) if (U[0][g] != -U[0][t - g]) { fprintf(stderr, "u1 not skew\n"); return 1; }
    for (int i = 0; i < 4; i++) for (int g = 0; g < t; g++) fixedpos[i][g] = (U[i][g] != 0) || (i == 0 && g == 0);
    int free_[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++) for (int g = 0; g < t; g++) free_[i] += !fixedpos[i][g];
    for (int i = 1; i < 4; i++)
        if ((free_[i] - tgt[i]) % 2 || abs(tgt[i]) > free_[i]) { fprintf(stderr, "row sum %d impossible for block %d (%d free)\n", tgt[i], i + 1, free_[i]); return 1; }
    long ssq = (long)tgt[1] * tgt[1] + (long)tgt[2] * tgt[2] + (long)tgt[3] * tgt[3];
    if (ssq != 4L * t - m - 1) fprintf(stderr, "warning: row sums squared %ld != %ld\n", ssq, 4L * t - m - 1);
    fprintf(stderr, "t=%d m=%d c-weight %d, targets (0, %d, %d, %d)\n", t, m, 4 * t - m - 1, tgt[1], tgt[2], tgt[3]);
    for (long r = 0; r < maxr; r++) {
        memset(W, 0, sizeof W);
        memset(X, 0, sizeof X);
        for (int p = 1; p <= h; p++) if (!fixedpos[0][p]) { int sg = rnd() & 1 ? 1 : -1; W[0][p] = sg; W[0][t - p] = -sg; }
        for (int i = 1; i < 4; i++) {
            int pos[MAXT], n = 0;
            for (int g = 0; g < t; g++) if (!fixedpos[i][g]) pos[n++] = g;
            int nplus = (n + tgt[i]) / 2;
            for (int k = n - 1; k > 0; k--) { int j = rnd() % (k + 1), x = pos[k]; pos[k] = pos[j]; pos[j] = x; }
            for (int k = 0; k < n; k++) W[i][pos[k]] = k < nplus ? 1 : -1;
        }
        for (int i = 0; i < 4; i++)
            for (int s = 0; s < t; s++) { long a = 0; for (int g = 0; g < t; g++) a += W[i][g] * W[i][(g + s) % t]; paf[i][s] = a; }
        for (int s = 0; s < t; s++) {
            long a = 0;
            for (int i = 0; i < 4; i++) for (int k = 0; k < nus[i]; k++) { int q = usup[i][k]; a += U[i][q] * (W[i][(q + s) % t] + W[i][(q - s + t) % t]); }
            X[s] = a;
        }
        double E = energy(), best = E, T = T0;
        for (long st = 0; st < steps && E > 0; st++) {
            if (!(st & 1023)) T = T0 * pow(T1 / T0, (double)st / steps);
            int i = rnd() % 4;
            if (i == 0) {
                int p = 1 + rnd() % h;
                if (fixedpos[0][p]) continue;
                int a = W[0][p];
                setv(0, p, -a); setv(0, t - p, a);
                double E2 = energy();
                if (E2 <= E || urand() < exp((E - E2) / T)) E = E2;
                else { setv(0, p, a); setv(0, t - p, -a); }
            } else {
                int g1 = rnd() % t, g2 = rnd() % t;
                if (fixedpos[i][g1] || fixedpos[i][g2] || W[i][g1] == W[i][g2]) continue;
                int a = W[i][g1], b = W[i][g2];
                setv(i, g1, b); setv(i, g2, a);
                double E2 = energy();
                if (E2 <= E || urand() < exp((E - E2) / T)) E = E2;
                else { setv(i, g1, a); setv(i, g2, b); }
            }
            if (E < best) best = E;
        }
        if (E == 0) {
            printf("SOL t=%d m=%d restart=%ld\n", t, m, r);
            for (int i = 0; i < 4; i++)
                for (int g = 0; g < t; g++) {
                    int c = (i == 0 && g == 0) ? 1 : U[i][g] ? 2 * U[i][g] : 3 * W[i][g];
                    printf("%d%c", c, g == t - 1 ? '\n' : ' ');
                }
            fflush(stdout);
            return 0;
        }
        if (r % 10 == 0) fprintf(stderr, "restart %ld final %.0f best %.0f\n", r, E, best);
    }
    return 2;
}
