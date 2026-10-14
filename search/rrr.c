/*
 * RRR (Elser's relaxed reflect-reflect / difference map) for four +-1 sequences on Z_m with a flat summed power
 * spectrum (search/rrr.py describes the problems):
 *
 *     sum_k |X_k(w)|^2 = TGT  for every w != 0,   row sums r_k fixed,   blocks with SYM[k] = 1 symmetric.
 *
 * A: +-1 vectors with the row sums (symmetric blocks constant on {t, -t}); projection = keep the largest entries.
 * B: per frequency w != 0, (X_1(w), .., X_4(w)) on the sphere of radius sqrt(TGT) in C^4; projection = rescaling.
 * x <- x + beta (P_A(2 P_B(x) - x) - P_B(x)). Every CHECK iterations P_A(x) is tested exactly (integer PAFs).
 *
 *     rrr m TGT TOFF r1,r2,r3,r4 s1,s2,s3,s4 beta iters seed nruns [threads]
 * TOFF: required value of sum_k PAF_k(g) for g != 0 (0 plain, -4 bordered). Runs nruns independent starts (seeds
 * seed, seed+1, ...) in parallel; prints "SOL run it" and the four sequences of the first solution, and
 * "STATS solved/runs iterations".
 */
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int m, H;                  /* H = (m - 1) / 2 for odd m, m / 2 for even m */
static double TGT, TGT0; static int TOFF, FREE;
static int R[4], SYM[4];
static double *CS, *SN;           /* CS[w * m + t] = cos(2 pi w t / m), w = 0..H */

static uint64_t rs_next(uint64_t *s) { *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17; return *s; }
static double gauss(uint64_t *s) {
    double u = ((rs_next(s) >> 11) + 0.5) / 9007199254740992.0, v = ((rs_next(s) >> 11) + 0.5) / 9007199254740992.0;
    return sqrt(-2 * log(u)) * cos(2 * M_PI * v);
}

/* forward: re[w] = sum x cos, im[w] = -sum x sin, w = 0..H */
static void dft(const double *x, double *re, double *im) {
    for (int w = 0; w <= H; w++) {
        const double *c = CS + (int64_t)w * m, *s = SN + (int64_t)w * m;
        double a = 0, b = 0;
        for (int t = 0; t < m; t++) { a += x[t] * c[t]; b -= x[t] * s[t]; }
        re[w] = a; im[w] = b;
    }
}
static void idft(const double *re, const double *im, double *x) {
    int even = (m % 2 == 0);
    for (int t = 0; t < m; t++) {
        double v = re[0];
        for (int w = 1; w <= H; w++) {
            double f = (even && w == H) ? 1.0 : 2.0;
            v += f * (re[w] * CS[(int64_t)w * m + t] - im[w] * SN[(int64_t)w * m + t]);
        }
        x[t] = v / m;
    }
}

static void projB(double x[4][1024], double out[4][1024]) {
    double re[4][1024], im[4][1024];
    for (int k = 0; k < 4; k++) dft(x[k], re[k], im[k]);
    for (int w = FREE ? 0 : 1; w <= H; w++) {
        double n2 = 0;
        double tg = w == 0 ? TGT0 : TGT;
        for (int k = 0; k < 4; k++) n2 += re[k][w] * re[k][w] + im[k][w] * im[k][w];
        double sc = n2 > 1e-24 ? sqrt(tg / n2) : 1.0;
        for (int k = 0; k < 4; k++) { re[k][w] *= sc; im[k][w] *= sc; }
    }
    for (int k = 0; k < 4; k++) idft(re[k], im[k], out[k]);
}

static int cmpd_desc(const void *a, const void *b) {
    double x = ((const double *)a)[0], y = ((const double *)b)[0];
    return x < y ? 1 : x > y ? -1 : 0;
}

/* nearest +-1 vector with row sum r (symmetric: constant on {t, -t}) */
static void projA1(const double *v, int r, int sym, double *out) {
    static __thread double buf[2 * 1024];
    if (FREE) {
        if (!sym) { for (int t = 0; t < m; t++) out[t] = v[t] >= 0 ? 1.0 : -1.0; return; }
        out[0] = v[0] >= 0 ? 1.0 : -1.0;
        if (m % 2 == 0) out[m / 2] = v[m / 2] >= 0 ? 1.0 : -1.0;
        for (int t = 1; t < m - t; t++) { double a = v[t] + v[m - t]; out[t] = out[m - t] = a >= 0 ? 1.0 : -1.0; }
        return;
    }
    if (!sym) {
        int npos = (m + r) / 2;
        for (int t = 0; t < m; t++) { buf[2 * t] = v[t]; buf[2 * t + 1] = t; }
        qsort(buf, m, 2 * sizeof(double), cmpd_desc);
        for (int i = 0; i < m; i++) out[(int)buf[2 * i + 1]] = i < npos ? 1.0 : -1.0;
        return;
    }
    /* orbits: {0} (size 1), {t, m - t} (size 2), {m/2} (size 1, m even) */
    int n1 = 0, n2 = 0; double o1[2][2], o2v[1024]; int o2t[1024];
    o1[n1][0] = v[0]; o1[n1][1] = 0; n1++;
    if (m % 2 == 0) { o1[n1][0] = v[m / 2]; o1[n1][1] = m / 2; n1++; }
    for (int t = 1; t < m - t; t++) { o2v[n2] = (v[t] + v[m - t]) / 2; o2t[n2] = t; n2++; }
    int total = n1 + 2 * n2, need = (r + total) / 2;     /* total size of the + orbits */
    /* sort size-2 orbits by value */
    for (int i = 0; i < n2; i++) { buf[2 * i] = o2v[i]; buf[2 * i + 1] = o2t[i]; }
    qsort(buf, n2, 2 * sizeof(double), cmpd_desc);
    if (n1 == 2 && o1[1][0] > o1[0][0]) { double a = o1[0][0], b = o1[0][1]; o1[0][0] = o1[1][0]; o1[0][1] = o1[1][1]; o1[1][0] = a; o1[1][1] = b; }
    double best = -1e300; int ba = -1, bb = -1;
    for (int a = 0; a <= n1; a++) {
        int b2 = need - a;
        if (b2 < 0 || (b2 & 1) || b2 / 2 > n2) continue;
        int b = b2 / 2; double sc = 0;
        for (int i = 0; i < n1; i++) sc += (i < a ? 1 : -1) * o1[i][0];
        for (int i = 0; i < n2; i++) sc += 2 * (i < b ? 1 : -1) * buf[2 * i];
        if (sc > best) { best = sc; ba = a; bb = b; }
    }
    if (ba < 0) { ba = 0; bb = 0; }
    for (int i = 0; i < n1; i++) out[(int)o1[i][1]] = i < ba ? 1.0 : -1.0;
    for (int i = 0; i < n2; i++) { int t = (int)buf[2 * i + 1]; double s = i < bb ? 1.0 : -1.0; out[t] = s; out[m - t] = s; }
}

static int exact_ok(double y[4][1024]) {
    for (int g = 1; g <= H; g++) {
        long s = 0;
        for (int k = 0; k < 4; k++) for (int t = 0; t < m; t++) s += (long)y[k][t] * (long)y[k][(t + g) % m];
        if (s != TOFF) return 0;
    }
    return 1;
}

static long error(double y[4][1024]) {
    long e = 0;
    for (int g = 1; g <= H; g++) {
        long s = 0;
        for (int k = 0; k < 4; k++) for (int t = 0; t < m; t++) s += (long)y[k][t] * (long)y[k][(t + g) % m];
        e += (s - TOFF) * (s - TOFF);
    }
    return e;
}

int main(int argc, char **argv) {
    if (argc < 10) { fprintf(stderr, "usage: rrr m TGT TOFF r1,r2,r3,r4 s1,s2,s3,s4 beta iters seed nruns [threads]\n"); return 1; }
    m = atoi(argv[1]); TGT = atof(argv[2]); TOFF = atoi(argv[3]);
    if (!strcmp(argv[4], "free")) { FREE = 1; } else sscanf(argv[4], "%d,%d,%d,%d", &R[0], &R[1], &R[2], &R[3]);
    sscanf(argv[5], "%d,%d,%d,%d", &SYM[0], &SYM[1], &SYM[2], &SYM[3]);
    double beta = atof(argv[6]); long iters = atol(argv[7]); uint64_t seed = strtoull(argv[8], 0, 10);
    int nruns = atoi(argv[9]); int threads = argc > 10 ? atoi(argv[10]) : omp_get_max_threads();
    if (m > 1000) return 2;
    TGT0 = TOFF == 0 ? TGT : 4.0;                 /* sum of squared row sums: 4m plain, 4 bordered */
    H = m / 2;
    CS = malloc(sizeof(double) * (H + 1) * m); SN = malloc(sizeof(double) * (H + 1) * m);
    for (int w = 0; w <= H; w++) for (int t = 0; t < m; t++) {
        double a = 2 * M_PI * (double)((long)w * t % m) / m;
        CS[(int64_t)w * m + t] = cos(a); SN[(int64_t)w * m + t] = sin(a);
    }
    int solved = 0; long totit = 0; volatile int stop = 0;
    double t0 = omp_get_wtime();
    #pragma omp parallel for num_threads(threads) schedule(dynamic, 1) reduction(+:solved, totit)
    for (int run = 0; run < nruns; run++) {
        if (stop) continue;
        uint64_t s = seed * 0x9E3779B97F4A7C15ULL + run * 0xD1B54A32D192ED03ULL + 1;
        static __thread double x[4][1024], pb[4][1024], z[4][1024], pa[4][1024], y[4][1024];
        for (int k = 0; k < 4; k++) for (int t = 0; t < m; t++) x[k][t] = gauss(&s);
        for (int k = 0; k < 4; k++) if (SYM[k]) for (int t = 1; t < m; t++) if (t < m - t) { double a = (x[k][t] + x[k][m - t]) / 2; x[k][t] = x[k][m - t] = a; }
        long it;
        for (it = 0; it < iters && !stop; it++) {
            projB(x, pb);
            for (int k = 0; k < 4; k++) for (int t = 0; t < m; t++) z[k][t] = 2 * pb[k][t] - x[k][t];
            for (int k = 0; k < 4; k++) projA1(z[k], R[k], SYM[k], pa[k]);
            for (int k = 0; k < 4; k++) for (int t = 0; t < m; t++) x[k][t] += beta * (pa[k][t] - pb[k][t]);
            if (it % 4 == 0) {
                for (int k = 0; k < 4; k++) memcpy(y[k], pa[k], sizeof(double) * m);   /* P_A(2 P_B(x) - x) */
                if (getenv("RRR_TRACE") && it % 20000 == 0) { fprintf(stderr, "run %d it %ld error %ld\n", run, it, error(y)); }
                if (exact_ok(y)) {
                    #pragma omp critical
                    {
                        if (!stop) {
                            printf("SOL run %d it %ld (%.1fs)\n", run, it, omp_get_wtime() - t0);
                            for (int k = 0; k < 4; k++) { for (int t = 0; t < m; t++) putchar(y[k][t] > 0 ? '+' : '-'); putchar('\n'); }
                            fflush(stdout);
                        }
                        stop = 1;
                    }
                    solved++;
                    break;
                }
            }
        }
        totit += it;
    }
    printf("STATS solved %d of %d runs, %ld iterations, %.1fs\n", solved, nruns, totit, omp_get_wtime() - t0);
    (void)error;
    return 0;
}
