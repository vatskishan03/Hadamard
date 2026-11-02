/*
 * lg_pt.c -- list generator for the skew k-tree: parallel tempering over the +-1 class signs of one block of a
 * cyclotomic product over a prime field, collecting DISTINCT members whose whitened energy |W (key - c)|^2 is below tau.
 * The output plugs into cp_lowe's kt4 pipeline: records (uint64 class mask, int16 key[D]); the class mask has bit c set
 * when the sign of class c is -1 (B, C, D blocks: f at 0 is not in the mask; skew block A: x_c = -x_(c+h)).
 * Only one of a pair of reversal twins (class mask and its rotation by h = e/2, same PAF key) is kept.
 *
 * Variables y_i in {+1,-1}, i < m. key_c(y) = const_c + sum_i lin[i][c] y_i + sum_{i<j} S[i][j][c] y_i y_j.
 * mode 0: moves are swaps of a +1 with a -1 (class sum fixed); mode 1: single flips (skew block A).
 *   cc -O3 -march=native -ffast-math -fopenmp lg_pt.c -o lg_pt -lm
 *   lg_pt in.bin out.bin tau tau_out seconds target threads seed R Tmin Tmax nprop maxdistinct
 * input: int32 hdr[12] = m D Z mode npos e skewA 0..; double W[Z][D]; double ck[D]; int32 const[D]; int32 lin[m][D];
 *        int16 SI[m][m][D] (symmetric in the first two indices, zero diagonal)
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXM 64
#define MAXD 64
#define MAXZ 64
#define MAXR 32

static int m, D, Z, mode, npos, e, skewA, hh;
static double Wd[MAXZ][MAXD], ckd[MAXD];
static int32_t cst[MAXD], lin[MAXM][MAXD];
static int16_t *SI;                 /* [m][m][D] */
static float *WS;                   /* [m][m][Z] */
static float WL[MAXM][MAXZ];
static float tau;
static float lam = 0.0f;
static uint64_t fullmask;

typedef struct {
    int8_t y[MAXM];
    int32_t a[MAXM][MAXD];
    float za[MAXM][MAXZ];
    float z[MAXZ];
    float E;
    int nplus, nminus, plus[MAXM], minus[MAXM], pos[MAXM];
    uint64_t mask;                  /* bit i set when y_i = -1 */
} Rep;

typedef struct { uint64_t s[4]; } Rng;
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnext(Rng *r) {
    uint64_t *s = r->s; uint64_t res = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return res;
}
static void rseed(Rng *r, uint64_t seed) {
    uint64_t x = seed * 0x9E3779B97F4A7C15ULL + 12345;
    for (int i = 0; i < 4; i++) { x += 0x9E3779B97F4A7C15ULL; uint64_t z = x; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; r->s[i] = z ^ (z >> 31); }
}
static inline double runif(Rng *r) { return ((rnext(r) >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
static inline uint32_t rbelow(Rng *r, uint32_t n) { return (uint32_t)(((rnext(r) >> 32) * (uint64_t)n) >> 32); }

/* hash set of canonical class masks with payload */
static uint64_t *hk; static float *hE; static int16_t *hkey; static uint64_t hcap;
static volatile long ndistinct = 0, nvisits_tot = 0;
static long maxdistinct;

static inline uint64_t canon_of(uint64_t cm) {
    uint64_t r = ((cm >> hh) | (cm << hh)) & fullmask;
    return r < cm ? r : cm;
}
static inline uint64_t class_mask(uint64_t ymask) {
    if (!skewA) return ymask;
    uint64_t lo = ymask & ((1ULL << hh) - 1), hi = (~ymask) & ((1ULL << hh) - 1);
    return lo | (hi << hh);
}
static inline uint64_t hash64(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }

/* exact key and whitened energy of a state from its integer fields */
static double exact_energy(const Rep *R, int32_t *key, double *zd) {
    for (int c = 0; c < D; c++) {
        long s = 0;
        for (int i = 0; i < m; i++) s += (long)R->y[i] * (R->a[i][c] + lin[i][c]);
        key[c] = cst[c] + (int32_t)(s / 2);
    }
    double E = 0;
    for (int k = 0; k < Z; k++) { double s = 0; for (int c = 0; c < D; c++) s += Wd[k][c] * (key[c] - ckd[c]); zd[k] = s; E += s * s; }
    return E;
}

static void init_fields(Rep *R) {
    for (int i = 0; i < m; i++) {
        for (int c = 0; c < D; c++) {
            int32_t s = lin[i][c];
            for (int j = 0; j < m; j++) s += (int32_t)SI[((size_t)i * m + j) * D + c] * R->y[j];
            R->a[i][c] = s;
        }
        for (int k = 0; k < Z; k++) { double s = 0; for (int c = 0; c < D; c++) s += Wd[k][c] * R->a[i][c]; R->za[i][k] = (float)s; }
    }
    int32_t key[MAXD]; double zd[MAXZ];
    double E = exact_energy(R, key, zd);
    for (int k = 0; k < Z; k++) R->z[k] = (float)zd[k];
    R->E = (float)E;
}

static void set_lists(Rep *R) {
    R->nplus = R->nminus = 0; R->mask = 0;
    for (int i = 0; i < m; i++) {
        if (R->y[i] > 0) { R->pos[i] = R->nplus; R->plus[R->nplus++] = i; }
        else { R->pos[i] = R->nminus; R->minus[R->nminus++] = i; R->mask |= 1ULL << i; }
    }
}

static void random_state(Rep *R, Rng *g) {
    if (mode == 0 || mode == 2) {
        for (int i = 0; i < m; i++) R->y[i] = -1;
        int cnt = 0;
        while (cnt < npos) { int i = rbelow(g, m); if (R->y[i] < 0) { R->y[i] = 1; cnt++; } }
    } else for (int i = 0; i < m; i++) R->y[i] = (rnext(g) >> 40) & 1 ? 1 : -1;
    set_lists(R); init_fields(R);
}

static inline void flip_update(Rep *R, int k) {
    const float coef = -2.0f * R->y[k]; const int ci = -2 * R->y[k];
    const int16_t *sk = SI + (size_t)k * m * D; const float *wk = WS + (size_t)k * m * Z;
    for (int j = 0; j < m; j++) {
        int32_t *aj = R->a[j]; const int16_t *sj = sk + (size_t)j * D;
        for (int c = 0; c < D; c++) aj[c] += ci * sj[c];
        float *zj = R->za[j]; const float *wj = wk + (size_t)j * Z;
        for (int q = 0; q < Z; q++) zj[q] += coef * wj[q];
    }
    for (int q = 0; q < Z; q++) R->z[q] += coef * R->za[k][q];
    if (mode != 1) {                                               /* maintain the plus / minus lists */
        if (R->y[k] > 0) { int p = R->pos[k]; int last = R->minus[R->nminus - 1]; /* move k from plus to minus */
            R->plus[p] = R->plus[R->nplus - 1]; R->pos[R->plus[p]] = p; R->nplus--;
            R->pos[k] = R->nminus; R->minus[R->nminus++] = k; (void)last; }
        else { int p = R->pos[k];
            R->minus[p] = R->minus[R->nminus - 1]; R->pos[R->minus[p]] = p; R->nminus--;
            R->pos[k] = R->nplus; R->plus[R->nplus++] = k; }
    }
    R->y[k] = -R->y[k]; R->mask ^= 1ULL << k;
}

static inline float energy_of(const float *z) { float s = 0; for (int q = 0; q < Z; q++) s += z[q] * z[q]; return s; }

static int record(Rep *R, long *novel) {
    int32_t key[MAXD]; double zd[MAXZ];
    double E = exact_energy(R, key, zd);
    if (fabs(E - R->E) > 1e-2) {                                  /* resynchronise the float state */
        for (int k = 0; k < Z; k++) R->z[k] = (float)zd[k];
        R->E = (float)E;
        for (int i = 0; i < m; i++) for (int k = 0; k < Z; k++) { double s = 0; for (int c = 0; c < D; c++) s += Wd[k][c] * R->a[i][c]; R->za[i][k] = (float)s; }
    }
    if (E > tau) return 0;
    uint64_t cm = canon_of(class_mask(R->mask)) + 1;
    uint64_t slot = hash64(cm) & (hcap - 1);
    for (;;) {
        uint64_t cur = hk[slot];
        if (cur == cm) return 1;
        if (cur == 0) {
            if (__sync_bool_compare_and_swap(&hk[slot], 0ULL, cm)) {
                hE[slot] = (float)E;
                for (int c = 0; c < D; c++) hkey[slot * D + c] = (int16_t)key[c];
                __sync_fetch_and_add(&ndistinct, 1); (*novel)++;
                return 2;
            }
            continue;
        }
        slot = (slot + 1) & (hcap - 1);
    }
}

int main(int argc, char **argv) {
    if (argc < 15) { fprintf(stderr, "usage: lg_pt in out tau tau_out seconds target threads seed R Tmin Tmax nprop maxdistinct\n"); return 2; }
    FILE *f = fopen(argv[1], "rb"); if (!f) { perror("in"); return 1; }
    int32_t hdr[12]; if (fread(hdr, 4, 12, f) != 12) return 1;
    m = hdr[0]; D = hdr[1]; Z = hdr[2]; mode = hdr[3]; npos = hdr[4]; e = hdr[5]; skewA = hdr[6]; hh = e / 2;
    fullmask = e == 64 ? ~0ULL : ((1ULL << e) - 1);
    if (m > MAXM || D > MAXD || Z > MAXZ) { fprintf(stderr, "sizes\n"); return 1; }
    double Wt[MAXZ * MAXD]; if (fread(Wt, 8, (size_t)Z * D, f) != (size_t)Z * D) return 1;
    for (int k = 0; k < Z; k++) for (int c = 0; c < D; c++) Wd[k][c] = Wt[k * D + c];
    if (fread(ckd, 8, D, f) != (size_t)D) return 1;
    if (fread(cst, 4, D, f) != (size_t)D) return 1;
    if (fread(lin, 4, (size_t)m * D, f) != (size_t)m * D) return 1;
    int32_t (*lt)[MAXD] = lin;  (void)lt;
    /* lin was read densely as [m][D]; spread to the strided array */
    { int32_t tmp[MAXM * MAXD]; memcpy(tmp, lin, sizeof(int32_t) * (size_t)m * D); for (int i = 0; i < m; i++) for (int c = 0; c < D; c++) lin[i][c] = tmp[i * D + c]; }
    SI = malloc(sizeof(int16_t) * (size_t)m * m * D); if (fread(SI, 2, (size_t)m * m * D, f) != (size_t)m * m * D) return 1;
    fclose(f);
    tau = (float)atof(argv[3]); float tau_out = (float)atof(argv[4]); double seconds = atof(argv[5]); long target = atol(argv[6]);
    int nth = atoi(argv[7]); uint64_t seed = strtoull(argv[8], 0, 10); int R_ = atoi(argv[9]); double Tmin = atof(argv[10]), Tmax = atof(argv[11]);
    long nprop = atol(argv[12]); maxdistinct = atol(argv[13]);
    const char *statpath = argc > 14 ? argv[14] : NULL;
    if (argc > 15) mode = atoi(argv[15]);
    if (argc > 16) lam = (float)atof(argv[16]);
    WS = malloc(sizeof(float) * (size_t)m * m * Z);
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) for (int k = 0; k < Z; k++) {
        double s = 0; for (int c = 0; c < D; c++) s += Wd[k][c] * SI[((size_t)i * m + j) * D + c];
        WS[((size_t)i * m + j) * Z + k] = (float)s;
    }
    hcap = 1; while (hcap < (uint64_t)(maxdistinct * 1.5)) hcap <<= 1;
    hk = calloc(hcap, 8); hE = malloc(hcap * 4); hkey = malloc(hcap * 2 * (size_t)D);
    if (!hk || !hE || !hkey) { fprintf(stderr, "out of memory\n"); return 3; }
    double T[MAXR]; if (R_ > MAXR) R_ = MAXR;
    for (int r = 0; r < R_; r++) T[r] = R_ == 1 ? Tmin : Tmin * pow(Tmax / Tmin, (double)r / (R_ - 1));
    omp_set_num_threads(nth);
    long tot_prop = 0, tot_acc = 0, tot_vis = 0; double t0 = omp_get_wtime();
    volatile int stop = 0;
    FILE *sf = statpath ? fopen(statpath, "w") : NULL;
    if (argc > 15 && atoi(argv[15]) == 3) {                          /* exhaustive Gray-code enumeration of the skew block */
        mode = 1;
        Rep *R = aligned_alloc(64, (sizeof(Rep) + 63) / 64 * 64);
        for (int i = 0; i < m; i++) R->y[i] = 1;
        set_lists(R); init_fields(R);
        long novel = 0; R->E = energy_of(R->z);
        if (R->E <= tau * 1.0005f) record(R, &novel);
        uint64_t total = 1ULL << (m - 1);
        for (uint64_t gc = 1; gc < total; gc++) {
            int k = __builtin_ctzll(gc) + 1; flip_update(R, k); R->E = energy_of(R->z);
            if (R->E <= tau * 1.0005f) record(R, &novel);
        }
        tot_prop = (long)total;
        goto output;
    }
    #pragma omp parallel reduction(+:tot_prop,tot_acc,tot_vis)
    {
        int tid = omp_get_thread_num(); Rng g; rseed(&g, seed * 1000003ULL + tid);
        Rep **rep = malloc(sizeof(Rep *) * R_);
        for (int r = 0; r < R_; r++) { rep[r] = aligned_alloc(64, (sizeof(Rep) + 63) / 64 * 64); random_state(rep[r], &g); }
        long novel = 0; long acc_t[MAXR] = {0}, prop_t[MAXR] = {0};
        double lastlog = t0; long rounds = 0;
        while (!stop) {
            for (int r = 0; r < R_; r++) {
                Rep *S = rep[r]; const float Tr = (float)T[r]; const float cutoff = 25.0f * Tr;
                for (long it = 0; it < nprop; it++) {
                    float dz[MAXZ]; int i, j = -1; float dE;
                    if (mode == 0) {
                        i = S->plus[rbelow(&g, S->nplus)]; j = S->minus[rbelow(&g, S->nminus)];
                        const float *zi = S->za[i], *zj = S->za[j], *w = WS + ((size_t)i * m + j) * Z;
                        float d1 = 0, d2 = 0;
                        for (int q = 0; q < Z; q++) { float d = -2.0f * zi[q] + 2.0f * zj[q] - 4.0f * w[q]; dz[q] = d; d1 += S->z[q] * d; d2 += d * d; }
                        dE = 2.0f * d1 + d2;
                    } else {
                        i = rbelow(&g, m); const float c = -2.0f * S->y[i]; const float *zi = S->za[i];
                        float d1 = 0, d2 = 0;
                        for (int q = 0; q < Z; q++) { float d = c * zi[q]; dz[q] = d; d1 += S->z[q] * d; d2 += d * d; }
                        dE = 2.0f * d1 + d2;
                        if (mode == 2) { int np = S->nplus, np2 = np + (S->y[i] > 0 ? -1 : 1); float a0 = np - npos, a1 = np2 - npos; dE += lam * (a1 * a1 - a0 * a0); }
                    }
                    prop_t[r]++;
                    if (dE > 0) {
                        if (dE > cutoff) continue;
                        if (runif(&g) >= expf(-dE / Tr)) continue;
                    }
                    flip_update(S, i); if (mode == 0) flip_update(S, j);
                    S->E = energy_of(S->z);
                    acc_t[r]++;
                    if (S->E <= tau * 1.0005f && (mode != 2 || S->nplus == npos)) {
                        tot_vis++;
                        int rc = record(S, &novel);
                        (void)rc;
                    }
                }
            }
            for (int r = (rounds & 1); r + 1 < R_; r += 2) {          /* replica exchange */
                double ea = rep[r]->E, eb = rep[r + 1]->E;
                if (mode == 2) { ea += lam * (double)(rep[r]->nplus - npos) * (rep[r]->nplus - npos); eb += lam * (double)(rep[r + 1]->nplus - npos) * (rep[r + 1]->nplus - npos); }
                double d = (1.0 / T[r] - 1.0 / T[r + 1]) * (ea - eb);
                if (d >= 0 || runif(&g) < exp(d)) { Rep *t = rep[r]; rep[r] = rep[r + 1]; rep[r + 1] = t; }
            }
            rounds++;
            double now = omp_get_wtime();
            if (ndistinct >= target || now - t0 >= seconds || ndistinct >= maxdistinct - 1) stop = 1;
            if (tid == 0 && now - lastlog > 5.0) {
                lastlog = now;
                fprintf(stderr, "t=%.0f distinct=%ld E[r]:", now - t0, ndistinct);
                for (int r = 0; r < R_; r++) fprintf(stderr, " %.2f", rep[r]->E);
                fprintf(stderr, "\n");
                if (sf) { fprintf(sf, "%.1f %ld\n", now - t0, ndistinct); fflush(sf); }
            }
        }
        for (int r = 0; r < R_; r++) { tot_prop += prop_t[r]; tot_acc += acc_t[r]; }
        if (tid == 0) { for (int r = 0; r < R_; r++) fprintf(stderr, "rung %d T=%.3f accept %.3g\n", r, T[r], (double)acc_t[r] / (prop_t[r] + 1)); }
    }
output:;
    double el = omp_get_wtime() - t0;
    if (sf) { fprintf(sf, "%.1f %ld\n", el, ndistinct); fclose(sf); }
    FILE *out = fopen(argv[2], "wb"); long nout = 0;
    for (uint64_t s = 0; s < hcap; s++) {
        if (!hk[s] || hE[s] > tau_out) continue;
        uint64_t cm = hk[s] - 1; fwrite(&cm, 8, 1, out); fwrite(hkey + s * D, 2, D, out); nout++;
    }
    fclose(out);
    printf("DISTINCT %ld WRITTEN %ld VISITS %ld PROPOSALS %ld ACCEPTS %ld WALL %.1f\n", ndistinct, nout, tot_vis, tot_prop, tot_acc, el);
    return 0;
}
