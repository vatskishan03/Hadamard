// flatbig.c -- like flat.c (rem, gen2 only) for structures with up to 128 cosets and 64 classes (128-bit masks).
//   cc -O3 -march=native -DLW=64 -o flatbig flatbig.c -lm   (LW = lanes >= C, multiple of 8)
//   flatbig rem  tab f0 npos nbase seed
//   flatbig gen2 tab f0 npos thr tau K nrestart seed out [tau_hi]      (out: distinct canonical masks, 16 bytes each)
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
static double now(void) { struct timespec t; clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void load_tab(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
    if (fscanf(f, "%d %d %d %d %d", &V, &M, &C, &ULEN, &PM) != 5) exit(1);
    if (M > MAXM || C > LW) { fprintf(stderr, "too large (M=%d C=%d LW=%d)\n", M, C, LW); exit(1); }
    for (int k = 0; k < M; k++) if (fscanf(f, "%lf", &ETR[k]) != 1) exit(1);
    for (int k = 0; k < M; k++) if (fscanf(f, "%lf", &ETI[k]) != 1) exit(1);
    fclose(f); INVV = 1.0 / V;
    for (int j = 0; j < M; j++) for (int c = 0; c < LW; c++) {
        if (c < C) { ER[j][c] = ETR[(j + c) % M]; EI[j][c] = ETI[(j + c) % M]; } else { ER[j][c] = 0; EI[j][c] = 0; }
    }
}
static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnd(void) { uint64_t r = rs[0] + rs[3], t = rs[1] << 17; rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45); return r; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static void seed_rng(uint64_t s) { for (int i = 0; i < 4; i++) { s += 0x9E3779B97F4A7C15ULL; uint64_t z = s; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; rs[i] = z ^ (z >> 31); } }
static inline u128 rotmask(u128 x, int k) { u128 full = (M == 128) ? ~(u128)0 : (((u128)1 << M) - 1); if (k == 0) return x; return ((x << k) | (x >> (M - k))) & full; }
static u128 canon(u128 x) { u128 b = x; for (int k = 1; k < M; k++) { u128 y = rotmask(x, k); if (y < b) b = y; } return b; }
static void spec(u128 mask, double f0, double *fr, double *fi) {
    for (int c = 0; c < LW; c++) { fr[c] = (c < C) ? f0 : sqrt((double)V); fi[c] = 0; }
    for (int j = 0; j < M; j++) { double sg = ((mask >> j) & 1) ? -1.0 : 1.0; for (int c = 0; c < LW; c++) { fr[c] += sg * ER[j][c]; fi[c] += sg * EI[j][c]; } }
}
static double energy(const double *fr, const double *fi) { double E = 0; for (int c = 0; c < LW; c++) { double d = (fr[c] * fr[c] + fi[c] * fi[c]) * INVV - 1.0; E += d * d; } return E; }
static void cmd_rem(char **argv) {
    load_tab(argv[2]); double f0 = atof(argv[3]); int npos = atoi(argv[4]); long nbase = (long)atof(argv[5]); seed_rng(strtoull(argv[6], 0, 10));
    int nneg = M - npos; long hist[400] = {0}; long tot = 0; double t0 = now();
    double fr[LW] __attribute__((aligned(64))), fi[LW] __attribute__((aligned(64)));
    for (long b = 0; b < nbase; b++) {
        int perm[MAXM]; for (int i = 0; i < M; i++) perm[i] = i;
        for (int i = M - 1; i > 0; i--) { int j = rnd() % (i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        u128 mask = 0; for (int i = npos; i < M; i++) mask |= (u128)1 << perm[i];
        spec(mask, f0, fr, fi);
        for (int a = 0; a < npos; a++) for (int c2 = 0; c2 < nneg; c2++) {
            int i = perm[a], j = perm[npos + c2]; double En = 0;
            for (int c = 0; c < LW; c++) { double xr = fr[c] + 2.0 * (ER[j][c] - ER[i][c]), xi = fi[c] + 2.0 * (EI[j][c] - EI[i][c]); double d = (xr * xr + xi * xi) * INVV - 1.0; En += d * d; }
            int bin = (int)(En / C * 100); if (bin > 399) bin = 399; hist[bin]++; tot++;
        }
    }
    double dt = now() - t0; fprintf(stderr, "rem: %ld evaluations in %.1fs (%.2f M/s)\n", tot, dt, tot / dt / 1e6);
    long cum = 0; for (int b = 0; b < 400; b++) { cum += hist[b]; if (b < 60 || (b % 20 == 0)) printf("Phi<%.2f  %ld  p=%.3e\n", (b + 1) / 100.0, cum, (double)cum / tot); }
}
static void cmd_gen2(int argc, char **argv) {
    load_tab(argv[2]);
    double f0 = atof(argv[3]); int npos = atoi(argv[4]); double thr = atof(argv[5]); double tau = atof(argv[6]);
    long K = (long)atof(argv[7]); long nrest = (long)atof(argv[8]); uint64_t seed = strtoull(argv[9], 0, 10); const char *out = argv[10];
    double tau_hi = argc > 11 ? atof(argv[11]) : tau; seed_rng(seed);
    int nneg = M - npos; double Ethr = thr * C;
    int HB = 22; uint64_t HS = 1ULL << HB; u128 *ht = calloc(HS, 16);
    FILE *fo = fopen(out, "ab"); if (!fo) { perror(out); exit(1); }
    long ndist = 0, nrec = 0, nprop = 0, nacc = 0; double t0 = now();
    double fr[LW] __attribute__((aligned(64))), fi[LW] __attribute__((aligned(64)));
    double nr[LW] __attribute__((aligned(64))), ni[LW] __attribute__((aligned(64)));
    double cool = (tau_hi > tau && K > 1) ? pow(tau / tau_hi, 1.0 / (K * 0.7)) : 1.0;
    for (long rr = 0; rr < nrest; rr++) {
        int P[MAXM], N[MAXM]; int perm[MAXM]; for (int i = 0; i < M; i++) perm[i] = i;
        for (int i = M - 1; i > 0; i--) { int j = rnd() % (i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        u128 mask = 0; for (int i = 0; i < npos; i++) P[i] = perm[i];
        for (int i = 0; i < nneg; i++) { N[i] = perm[npos + i]; mask |= (u128)1 << N[i]; }
        spec(mask, f0, fr, fi);
        double E = energy(fr, fi), T = tau_hi;
        for (long it = 0; it < K; it++) {
            int a = rnd() % npos, b = rnd() % nneg; int i = P[a], j = N[b]; double En = 0;
            for (int c = 0; c < LW; c++) { double xr = fr[c] + 2.0 * (ER[j][c] - ER[i][c]), xi = fi[c] + 2.0 * (EI[j][c] - EI[i][c]); nr[c] = xr; ni[c] = xi; double d = (xr * xr + xi * xi) * INVV - 1.0; En += d * d; }
            nprop++;
            if (En <= Ethr) {
                u128 nm = mask ^ (((u128)1 << i) | ((u128)1 << j)); u128 cn = canon(nm); nrec++;
                uint64_t h = (((uint64_t)cn ^ (uint64_t)(cn >> 64) * 0x9E3779B97F4A7C15ULL) * 0xD6E8FEB86659FD93ULL) >> (64 - HB);
                while (ht[h] && ht[h] != cn) h = (h + 1) & (HS - 1);
                if (!ht[h]) { ht[h] = cn; ndist++; fwrite(&cn, 16, 1, fo); if (ndist > (long)(HS * 0.7)) { fprintf(stderr, "hash full\n"); goto done; } }
            }
            int acc = (En <= E) || (T > 0 && urand() < exp((E - En) / T));
            if (acc) { memcpy(fr, nr, sizeof fr); memcpy(fi, ni, sizeof fi); E = En; nacc++; mask ^= ((u128)1 << i) | ((u128)1 << j); P[a] = j; N[b] = i; }
            if (T > tau) { T *= cool; if (T < tau) T = tau; }
        }
    }
done:
    fclose(fo); double dt = now() - t0;
    fprintf(stderr, "gen2: proposals %ld (%.1f M/s), acc %.5f, under-thr proposals %ld, DISTINCT necklaces %ld, cpu %.2fs => %.2f distinct/s\n", nprop, nprop / dt / 1e6, (double)nacc / nprop, nrec, ndist, dt, ndist / dt);
    printf("%ld %ld %.3f\n", nprop, ndist, dt);
}
int main(int argc, char **argv) {
    if (argc < 3) return 2;
    if (!strcmp(argv[1], "rem")) cmd_rem(argv); else if (!strcmp(argv[1], "gen2")) cmd_gen2(argc, argv); else return 2;
    return 0;
}
