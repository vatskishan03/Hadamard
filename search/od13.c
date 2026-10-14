/* od13.c: simulated annealing for a skew-type orthogonal design of type (1, 3, 4t - 4) in order 4t (t odd)
 * in Goethals-Seidel form.
 *
 * Four circulant blocks of length t.  A1 = x1 I + x3 W0 with W0 skew (W0[0] = 0, W0[-g] = -W0[g]);
 * A2, A3, A4 = x2 I + x3 Wi with Wi[0] = 0 and Wi = +-1 elsewhere.  The Goethals-Seidel array of A1..A4 is an
 * orthogonal design of type (1, 3, 4t - 4) iff
 *   (a) sum_i PAF(Wi)(s) = 0 for s != 0                       (x3 x3 terms),
 *   (b) sum_{i=1..3} (Wi[s] + Wi[-s]) = 0 for every s != 0     (x2 x3 terms; the x1 x3 terms vanish because W0 is
 *                                                                skew, the x2 x2 terms because each x2 part is a
 *                                                                single point).
 * (b) is kept invariant: for each pair {s, -s} the six values (W1[s], W1[-s], W2[s], W2[-s], W3[s], W3[-s]) are one
 * of the 20 patterns whose symmetric parts cancel.  The energy is sum_{s=1..(t-1)/2} (sum_i PAF(Wi)(s))^2.
 *
 * Usage: od13 t seed [max_restarts] [steps_per_restart]
 * Output: one line "SOL" followed by the four blocks as integer codes (1 = x1, 2 = x2, +-3 = +-x3).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

static int t, h;
static int W[4][256];
static long paf[4][256];
static int pat[20][6], npat = 0;
static int cfg[128];

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

/* change W[i][a] to v, updating paf[i] */
static void setv(int i, int a, int v) {
    int d = v - W[i][a];
    if (!d) return;
    for (int s = 1; s < t; s++) paf[i][s] += (long)d * (W[i][(a + s) % t] + W[i][(a - s + t) % t]);
    W[i][a] = v;
    long sq = 0; for (int g = 0; g < t; g++) sq += W[i][g] * W[i][g];
    paf[i][0] = sq;
}

static long energy(void) {
    long e = 0;
    for (int s = 1; s <= h; s++) {
        long z = paf[0][s] + paf[1][s] + paf[2][s] + paf[3][s];
        e += z * z;
    }
    return e;
}

static void apply_cfg(int p, int c) {
    for (int i = 0; i < 3; i++) {
        setv(i + 1, p, pat[c][2 * i]);
        setv(i + 1, t - p, pat[c][2 * i + 1]);
    }
    cfg[p] = c;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: od13 t seed [max_restarts] [steps]\n"); return 1; }
    t = atoi(argv[1]);
    rs = strtoull(argv[2], 0, 10) * 0x9E3779B97F4A7C15ULL + 12345;
    long maxr = argc > 3 ? atol(argv[3]) : 1000000;
    long steps = argc > 4 ? atol(argv[4]) : 2000000;
    if (t % 2 == 0 || t < 3 || t > 255) { fprintf(stderr, "t must be odd, 3..255\n"); return 1; }
    h = (t - 1) / 2;
    for (int m = 0; m < 64; m++) {
        int v[6], sym = 0;
        for (int k = 0; k < 6; k++) v[k] = (m >> k) & 1 ? 1 : -1;
        for (int i = 0; i < 3; i++) sym += v[2 * i] + v[2 * i + 1];
        if (sym == 0) { memcpy(pat[npat], v, sizeof v); npat++; }
    }
    if (npat != 20) { fprintf(stderr, "pattern count %d\n", npat); return 1; }
    for (long r = 0; r < maxr; r++) {
        memset(W, 0, sizeof W);
        for (int p = 1; p <= h; p++) {
            int sg = (rnd() & 1) ? 1 : -1;
            W[0][p] = sg; W[0][t - p] = -sg;
            int c = rnd() % 20;
            for (int i = 0; i < 3; i++) { W[i + 1][p] = pat[c][2 * i]; W[i + 1][t - p] = pat[c][2 * i + 1]; }
            cfg[p] = c;
        }
        for (int i = 0; i < 4; i++) full_paf(i);
        long E = energy();
        double T0 = 40.0, T1 = 0.5;
        for (long st = 0; st < steps && E; st++) {
            double T = T0 * pow(T1 / T0, (double)st / steps);
            int p = 1 + rnd() % h;
            long E2;
            if (rnd() % 4 == 0) {
                int a = W[0][p];
                setv(0, p, -a); setv(0, t - p, a);
                E2 = energy();
                if (E2 <= E || urand() < exp((E - E2) / T)) E = E2;
                else { setv(0, p, a); setv(0, t - p, -a); }
            } else {
                int old = cfg[p], c = rnd() % 19;
                if (c >= old) c++;
                apply_cfg(p, c);
                E2 = energy();
                if (E2 <= E || urand() < exp((E - E2) / T)) E = E2;
                else apply_cfg(p, old);
            }
        }
        if (!E) {
            printf("SOL t=%d restart=%ld\n", t, r);
            for (int i = 0; i < 4; i++) {
                for (int g = 0; g < t; g++) {
                    int c = g == 0 ? (i == 0 ? 1 : 2) : 3 * W[i][g];
                    printf("%d%c", c, g == t - 1 ? '\n' : ' ');
                }
            }
            fflush(stdout);
            return 0;
        }
        if (r % 10 == 0) { fprintf(stderr, "restart %ld best-final E=%ld\n", r, E); }
    }
    return 2;
}
