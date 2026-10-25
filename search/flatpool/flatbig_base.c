// flat.c -- near-flat orbit-sign vectors for Goethals-Seidel quadruples over Z_v (v prime), multiplier group U.
//
// A sequence is (f0; s_0..s_{m-1}), s_j = sign on the coset rho^j U (m = (v-1)/|U| cosets, cyclic order), bit j of the
// mask set = sign -1.  Spectrum  fhat_c = f0 + sum_j s_j eta_{(j+c) mod m}, c in Z_C (C = m/2); PSD_c = |fhat_c|^2.
// Flatness  Phi = (1/C) sum_c (PSD_c / v - 1)^2   (Egan-Kotsireas-London-O Cathain's Phi_f).
// Multiplication by rho rotates the mask (cyclic shift of j); the key of the rotated vector is the cyclic shift of
// the key by k in Z_C; rotation by C (= x -> -x) leaves the key unchanged.  Necklace = canonical (min) rotation.
//
//   cc -O3 -march=native -o flat flat.c -lm
//   flat gen  tab.txt f0 npos thr tau nsteps seed out.bin [thin] [burn]    Metropolis on swaps (+ <-> -), logs states with Phi<=thr
//   flat keys tab.txt f0 masks.bin keys.bin                                  int16 keys kappa_e = (PAF_e - pm)/4, e in Z_C
//   flat phi  tab.txt f0 masks.bin                                           Phi histogram
// tab.txt: v m C ulen pm, then m reals eta_re, then m reals eta_im.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXM 256
#ifndef LW
#define LW 64
#endif
typedef unsigned __int128 u128;
static int V, M, C, ULEN, PM;
static double ETR[MAXM], ETI[MAXM];
static double ER[MAXM][LW] __attribute__((aligned(64))), EI[MAXM][LW] __attribute__((aligned(64)));
static double INVV;
static double TH[MAXM];

static double now(void) { struct timespec t; clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static void load_tab(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
    if (fscanf(f, "%d %d %d %d %d", &V, &M, &C, &ULEN, &PM) != 5) exit(1);
    if (M > MAXM || C > LW) { fprintf(stderr, "too large\n"); exit(1); }
    for (int k = 0; k < M; k++) if (fscanf(f, "%lf", &ETR[k]) != 1) exit(1);
    for (int k = 0; k < M; k++) if (fscanf(f, "%lf", &ETI[k]) != 1) exit(1);
    fclose(f);
    INVV = 1.0 / V;
    for (int j = 0; j < M; j++) for (int c = 0; c < LW; c++) {
        if (c < C) { ER[j][c] = ETR[(j + c) % M]; EI[j][c] = ETI[(j + c) % M]; } else { ER[j][c] = 0; EI[j][c] = 0; }
    }
    for (int k = 0; k < M; k++) TH[k] = 2 * ETR[k];
}

static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnd(void) {
    uint64_t r = rs[0] + rs[3], t = rs[1] << 17;
    rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45);
    return r;
}
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static void seed_rng(uint64_t s) {
    for (int i = 0; i < 4; i++) { s += 0x9E3779B97F4A7C15ULL; uint64_t z = s; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; rs[i] = z ^ (z >> 31); }
}

static inline u128 rotmask(u128 x, int k) {
    u128 full = (M == 128) ? ~(u128)0 : (((u128)1 << M) - 1);
    if (k == 0) return x;
    return ((x << k) | (x >> (M - k))) & full;
}
static u128 canon(u128 x) {
    u128 b = x;
    for (int k = 1; k < M; k++) { u128 y = rotmask(x, k); if (y < b) b = y; }
    return b;
}

// full spectrum from mask: fr, fi (lane arrays)
static void spec(u128 mask, double f0, double *fr, double *fi) {
    for (int c = 0; c < LW; c++) { fr[c] = (c < C) ? f0 : sqrt((double)V); fi[c] = 0; }
    for (int j = 0; j < M; j++) {
        double sg = (mask >> j & 1) ? -1.0 : 1.0;
        for (int c = 0; c < LW; c++) { fr[c] += sg * ER[j][c]; fi[c] += sg * EI[j][c]; }
    }
}
static double energy(const double *fr, const double *fi) {
    double E = 0;
    for (int c = 0; c < LW; c++) { double d = (fr[c] * fr[c] + fi[c] * fi[c]) * INVV - 1.0; E += d * d; }
    return E;
}

static void cmd_gen(int argc, char **argv) {
    // flat gen tab f0 npos thr tau nsteps seed out [thin] [burn]
    load_tab(argv[2]);
    double f0 = atof(argv[3]); int npos = atoi(argv[4]); double thr = atof(argv[5]); double tau = atof(argv[6]);
    long nsteps = (long)atof(argv[7]); uint64_t seed = strtoull(argv[8], 0, 10); const char *out = argv[9];
    long thin = argc > 10 ? atol(argv[10]) : 1; long burn = argc > 11 ? (long)atof(argv[11]) : 100000;
    seed_rng(seed);
    int sg[MAXM], P[MAXM], N[MAXM], np_ = 0, nn = 0, pos[MAXM];
    // random start with npos plus signs
    int perm[MAXM]; for (int i = 0; i < M; i++) perm[i] = i;
    for (int i = M - 1; i > 0; i--) { int j = rnd() % (i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
    for (int i = 0; i < M; i++) {
        int j = perm[i];
        if (i < npos) { sg[j] = 1; pos[j] = np_; P[np_++] = j; } else { sg[j] = -1; pos[j] = nn; N[nn++] = j; }
    }
    double fr[LW] __attribute__((aligned(64))), fi[LW] __attribute__((aligned(64)));
    double nr[LW] __attribute__((aligned(64))), ni[LW] __attribute__((aligned(64)));
    uint64_t mask = 0; for (int j = 0; j < M; j++) if (sg[j] < 0) mask |= 1ULL << j;
    spec(mask, f0, fr, fi);
    double E = energy(fr, fi), Ethr = thr * C;
    FILE *fo = fopen(out, "ab"); if (!fo) { perror(out); exit(1); }
    uint64_t *buf = malloc(8 * (1 << 16)); int nb = 0;
    long nrec = 0, nacc = 0, nin = 0, since = 0; double t0 = now(); double Emin = 1e9;
    long total = burn + nsteps;
    for (long it = 0; it < total; it++) {
        int a = rnd() % np_, b = rnd() % nn; int i = P[a], j = N[b];   // i: + -> -, j: - -> +
        const double *Ei_r = ER[i], *Ei_i = EI[i], *Ej_r = ER[j], *Ej_i = EI[j];
        double En = 0;
        for (int c = 0; c < LW; c++) {
            double xr = fr[c] + 2.0 * (Ej_r[c] - Ei_r[c]), xi = fi[c] + 2.0 * (Ej_i[c] - Ei_i[c]);
            nr[c] = xr; ni[c] = xi;
            double d = (xr * xr + xi * xi) * INVV - 1.0; En += d * d;
        }
        int acc;
        if (tau > 0) acc = (En <= E) || (urand() < exp((E - En) / tau));
        else acc = (En <= Ethr) || (En <= E);      /* tau <= 0: hard wall at Phi <= thr (uniform on the flat set) */
        if (acc) {
            memcpy(fr, nr, sizeof fr); memcpy(fi, ni, sizeof fi); E = En; nacc++;
            sg[i] = -1; sg[j] = 1; mask ^= (1ULL << i) | (1ULL << j);
            int pa = pos[i], pb = pos[j]; P[a] = j; N[b] = i; pos[j] = a; pos[i] = b; (void)pa; (void)pb;
            if (E < Emin) Emin = E;
            if (it >= burn && E <= Ethr) {
                nin++;
                if (++since >= thin) {
                    since = 0; buf[nb++] = canon(mask); nrec++;
                    if (nb == (1 << 16)) { fwrite(buf, 8, nb, fo); nb = 0; }
                }
            }
        }
        if ((it & 0xFFFFFFF) == 0xFFFFFFF) {   // recompute to kill drift
            spec(mask, f0, fr, fi); E = energy(fr, fi);
        }
    }
    if (nb) fwrite(buf, 8, nb, fo);
    fclose(fo);
    double dt = now() - t0;
    fprintf(stderr, "gen: steps %ld, accepted %ld (%.3f), states under thr %ld, recorded %ld, Emin/C %.4f, cpu %.2fs, %.2f Mprop/s\n",
            total, nacc, (double)nacc / total, nin, nrec, Emin / C, dt, total / dt / 1e6);
    printf("%ld %ld %.3f\n", total, nrec, dt);
}

// PAF classes (double) from the spectrum
static void paf_from_mask(uint64_t mask, double f0, double *paf) {
    double fr[LW], fi[LW]; spec(mask, f0, fr, fi);
    int S = 0; for (int j = 0; j < M; j++) S += (mask >> j & 1) ? -1 : 1;
    double r = f0 + ULEN * (double)S;
    double psd[LW]; for (int c = 0; c < C; c++) psd[c] = fr[c] * fr[c] + fi[c] * fi[c];
    for (int e = 0; e < C; e++) {
        double a = r * r;
        for (int c = 0; c < C; c++) a += psd[c] * TH[(c + e) % M];
        paf[e] = a * INVV;
    }
}

static void cmd_keys(int argc, char **argv) {
    load_tab(argv[2]); double f0 = atof(argv[3]);
    FILE *fi = fopen(argv[4], "rb"); FILE *fo = fopen(argv[5], "wb");
    if (!fi || !fo) { perror("open"); exit(1); }
    uint64_t *buf = malloc(8 * 4096); int16_t *ob = malloc(2 * C * 4096); size_t n;
    long bad = 0, tot = 0;
    while ((n = fread(buf, 8, 4096, fi)) > 0) {
        for (size_t t = 0; t < n; t++) {
            double paf[LW]; paf_from_mask(buf[t], f0, paf);
            for (int e = 0; e < C; e++) {
                long p = lround(paf[e]);
                if (fabs(paf[e] - p) > 1e-4 || ((p - PM) % 4 + 4) % 4) bad++;
                ob[t * C + e] = (int16_t)((p - PM) / 4);
            }
        }
        fwrite(ob, 2 * C, n, fo); tot += n;
    }
    fclose(fo); fclose(fi);
    fprintf(stderr, "keys: %ld vectors, %ld bad entries\n", tot, bad);
}

static void cmd_phi(int argc, char **argv) {
    load_tab(argv[2]); double f0 = atof(argv[3]);
    FILE *fi = fopen(argv[4], "rb"); if (!fi) { perror("open"); exit(1); }
    uint64_t x; long n = 0; double s1 = 0, mn = 1e9, mx = 0; long hist[40] = {0};
    while (fread(&x, 8, 1, fi) == 1) {
        double fr[LW], fi_[LW]; spec(x, f0, fr, fi_); double ph = energy(fr, fi_) / C;
        s1 += ph; n++; if (ph < mn) mn = ph; if (ph > mx) mx = ph; int b = (int)(ph * 100); if (b > 39) b = 39; hist[b]++;
    }
    printf("n=%ld mean Phi %.4f min %.4f max %.4f\n", n, s1 / n, mn, mx);
    for (int b = 0; b < 40; b++) if (hist[b]) printf("  Phi in [%.2f,%.2f): %ld\n", b / 100.0, (b + 1) / 100.0, hist[b]);
}


// flat rem tab f0 npos nbase seed : Phi histogram of swap neighbours of random bases (all npos*nneg swaps per base)
static void cmd_rem(int argc, char **argv) {
    load_tab(argv[2]); double f0 = atof(argv[3]); int npos = atoi(argv[4]); long nbase = (long)atof(argv[5]); seed_rng(strtoull(argv[6], 0, 10));
    int nneg = M - npos; long hist[400] = {0}; long tot = 0; double t0 = now();
    double fr[LW] __attribute__((aligned(64))), fi[LW] __attribute__((aligned(64)));
    for (long b = 0; b < nbase; b++) {
        int perm[MAXM]; for (int i = 0; i < M; i++) perm[i] = i;
        for (int i = M - 1; i > 0; i--) { int j = rnd() % (i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        uint64_t mask = 0; for (int i = npos; i < M; i++) mask |= 1ULL << perm[i];
        spec(mask, f0, fr, fi);
        for (int a = 0; a < npos; a++) for (int c2 = 0; c2 < nneg; c2++) {
            int i = perm[a], j = perm[npos + c2];   // i (+) -> -, j (-) -> +
            double En = 0;
            for (int c = 0; c < LW; c++) {
                double xr = fr[c] + 2.0 * (ER[j][c] - ER[i][c]), xi = fi[c] + 2.0 * (EI[j][c] - EI[i][c]);
                double d = (xr * xr + xi * xi) * INVV - 1.0; En += d * d;
            }
            int bin = (int)(En / C * 100); if (bin > 399) bin = 399; hist[bin]++; tot++;
        }
    }
    double dt = now() - t0;
    fprintf(stderr, "rem: %ld evaluations in %.1fs (%.2f M/s)\n", tot, dt, tot / dt / 1e6);
    long cum = 0;
    for (int b = 0; b < 400; b++) { cum += hist[b]; if (b < 40 || (b % 20 == 0)) printf("Phi<%.2f  %ld  p=%.3e\n", (b + 1) / 100.0, cum, (double)cum / tot); }
}


// flat gen2 tab f0 npos thr tau K nrestart seed out [tau_hi]
// restarts from random bases; each runs K Metropolis proposals (tau annealed geometrically from tau_hi to tau if tau_hi>tau);
// EVERY proposal with Phi<=thr is recorded (canonical necklace), deduplicated in a hash table; prints distinct counts.
static void cmd_gen2(int argc, char **argv) {
    load_tab(argv[2]);
    double f0 = atof(argv[3]); int npos = atoi(argv[4]); double thr = atof(argv[5]); double tau = atof(argv[6]);
    long K = (long)atof(argv[7]); long nrest = (long)atof(argv[8]); uint64_t seed = strtoull(argv[9], 0, 10); const char *out = argv[10];
    double tau_hi = argc > 11 ? atof(argv[11]) : tau;
    seed_rng(seed);
    int nneg = M - npos; double Ethr = thr * C;
    int HB = 24; uint64_t HS = 1ULL << HB; uint64_t *ht = calloc(HS, 8);
    FILE *fo = fopen(out, "ab"); if (!fo) { perror(out); exit(1); }
    long ndist = 0, nrec = 0, nprop = 0, nacc = 0; double t0 = now();
    double fr[LW] __attribute__((aligned(64))), fi[LW] __attribute__((aligned(64)));
    double nr[LW] __attribute__((aligned(64))), ni[LW] __attribute__((aligned(64)));
    double cool = (tau_hi > tau && K > 1) ? pow(tau / tau_hi, 1.0 / (K * 0.7)) : 1.0;
    for (long rr = 0; rr < nrest; rr++) {
        int P[MAXM], N[MAXM]; int perm[MAXM]; for (int i = 0; i < M; i++) perm[i] = i;
        for (int i = M - 1; i > 0; i--) { int j = rnd() % (i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        uint64_t mask = 0;
        for (int i = 0; i < npos; i++) P[i] = perm[i];
        for (int i = 0; i < nneg; i++) { N[i] = perm[npos + i]; mask |= 1ULL << N[i]; }
        spec(mask, f0, fr, fi);
        double E = energy(fr, fi), T = tau_hi;
        for (long it = 0; it < K; it++) {
            int a = rnd() % npos, b = rnd() % nneg; int i = P[a], j = N[b];
            const double *Ei_r = ER[i], *Ei_i = EI[i], *Ej_r = ER[j], *Ej_i = EI[j];
            double En = 0;
            for (int c = 0; c < LW; c++) {
                double xr = fr[c] + 2.0 * (Ej_r[c] - Ei_r[c]), xi = fi[c] + 2.0 * (Ej_i[c] - Ei_i[c]);
                nr[c] = xr; ni[c] = xi;
                double d = (xr * xr + xi * xi) * INVV - 1.0; En += d * d;
            }
            nprop++;
            if (En <= Ethr) {
                uint64_t nm = mask ^ ((1ULL << i) | (1ULL << j)); uint64_t cn = canon(nm); nrec++;
                uint64_t h = (cn * 0x9E3779B97F4A7C15ULL) >> (64 - HB);
                while (ht[h] && ht[h] != cn) h = (h + 1) & (HS - 1);
                if (!ht[h]) { ht[h] = cn; ndist++; fwrite(&cn, 8, 1, fo); if (ndist > (long)(HS * 0.7)) { fprintf(stderr, "hash full\n"); goto done; } }
            }
            int acc = (En <= E) || (T > 0 && urand() < exp((E - En) / T));
            if (acc) {
                memcpy(fr, nr, sizeof fr); memcpy(fi, ni, sizeof fi); E = En; nacc++;
                mask ^= (1ULL << i) | (1ULL << j); P[a] = j; N[b] = i;
            }
            if (T > tau) { T *= cool; if (T < tau) T = tau; }
        }
    }
done:
    fclose(fo);
    double dt = now() - t0;
    fprintf(stderr, "gen2: proposals %ld (%.1f M/s), acc %.5f, under-thr proposals %ld, DISTINCT necklaces %ld, cpu %.2fs => %.2f distinct/s\n",
            nprop, nprop / dt / 1e6, (double)nacc / nprop, nrec, ndist, dt, ndist / dt);
    printf("%ld %ld %.3f\n", nprop, ndist, dt);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: flat gen|keys|phi ...\n"); return 2; }
    if (!strcmp(argv[1], "gen")) cmd_gen(argc, argv);
    else if (!strcmp(argv[1], "keys")) cmd_keys(argc, argv);
    else if (!strcmp(argv[1], "phi")) cmd_phi(argc, argv);
    else if (!strcmp(argv[1], "rem")) cmd_rem(argc, argv);
    else if (!strcmp(argv[1], "gen2")) cmd_gen2(argc, argv);
    else return 2;
    return 0;
}
