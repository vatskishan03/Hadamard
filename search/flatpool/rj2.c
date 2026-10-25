// rj.c -- rotation-virtual four-list k-tree for flat necklace pools (Goethals-Seidel over Z_v, multiplier group U).
//
// Four roles 0..3 (pairs (0,1) and (2,3)), each a list of necklaces with int8 keys K[x][e], e in Z_C.  The entries of a
// role are (x,k), k in Z_C (or k = 0 only when rot = 0), with key kappa[e] = K[x][(e+k) % C]  (rotation by k of the
// orbit-sign vector; see flat.c).  Find (e0,e1,e2,e3) with  sum_i kappa_i[e] = T  for every class e.
// k-tree: merge classes S = {0..ns-1}; for each split value sigma (ns integers, in order of decreasing solution density)
// the pairs (e0,e1) with kappa_0+kappa_1 = sigma on S and (e2,e3) with kappa_2+kappa_3 = T - sigma on S are formed by a
// merge pass over the S-sorted lists and joined on a 64-bit linear hash of the other classes (partitioned hash join);
// every hash match is verified exactly on all C classes.
//   cc -O3 -march=native -o rj rj.c -lm
//   rj INPUT NS NBUCKETS [MAXPAIRS=1.2e8] [FIRST=0] [LAST=NBUCKETS]
// INPUT: int32 C, T, rot[4], n[4]; then int8 K[n_i][C] for i = 0..3.
// Output: "SOL x0 k0 x1 k1 x2 k2 x3 k3" ; "BUCKET i sigma pairsAB pairsCD matches found cpu" on stderr ; "TOTAL ...".
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static int C, T0, NS, ROT[4], NC_;
static long N[4], M[4];                 // necklaces, entries
static int8_t *K[4];
static uint64_t R[64];
static uint64_t *SK[4];                 // sorted S-keys, as int64 W
static uint32_t *ID[4];                 // entry ids (x*C + k or x)
static uint64_t *HS[4];                 // hashes aligned with sorted order
static uint64_t hT;

static double now(void) { struct timespec t; clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static inline int keyval(int L, uint32_t id, int e) {
    if (ROT[L]) { uint32_t x = id / C, k = id % C; int ee = e + k; if (ee >= C) ee -= C; return K[L][(size_t)x * C + ee]; }
    return K[L][(size_t)id * C + e];
}

// W-key of an entry: sum_s kappa[s] * 64^s (signed)
static inline int64_t wkey(int L, uint32_t id) {
    int64_t w = 0;
    for (int s = NS - 1; s >= 0; s--) w = w * 64 + keyval(L, id, s);
    return w;
}
static inline uint64_t hkey(int L, uint32_t id) {
    uint64_t h = 0;
    for (int e = NS; e < NC_; e++) h += R[e] * (uint64_t)(int64_t)keyval(L, id, e);
    return h;
}

static void radix_sort(uint64_t *key, uint32_t *id, long n) {   // sort by key as signed int64 (flip sign bit)
    uint64_t *k2 = malloc(8 * n); uint32_t *i2 = malloc(4 * n);
    if (!k2 || !i2) { fprintf(stderr, "rj: out of memory (sort)\n"); exit(3); }
    for (long i = 0; i < n; i++) key[i] ^= 0x8000000000000000ULL;
    uint64_t *a = key, *b = k2; uint32_t *ia = id, *ib = i2;
    for (int pass = 0; pass < 4; pass++) {
        int sh = pass * 16; long *cnt = calloc(65537, sizeof(long));
        for (long i = 0; i < n; i++) cnt[((a[i] >> sh) & 0xFFFF) + 1]++;
        int skip = 0; for (int q = 0; q < 65536; q++) if (cnt[q + 1] == n) skip = 1;
        if (!skip) {
            for (int q = 0; q < 65536; q++) cnt[q + 1] += cnt[q];
            for (long i = 0; i < n; i++) { long p = cnt[(a[i] >> sh) & 0xFFFF]++; b[p] = a[i]; ib[p] = ia[i]; }
            uint64_t *t = a; a = b; b = t; uint32_t *ti = ia; ia = ib; ib = ti;
        }
        free(cnt);
    }
    if (a != key) { memcpy(key, a, 8 * n); memcpy(id, ia, 4 * n); }
    for (long i = 0; i < n; i++) key[i] ^= 0x8000000000000000ULL;
    free(k2); free(i2);
}

static inline uint64_t mixp(uint64_t x) { return x * 0xD6E8FEB86659FD93ULL; }

typedef struct { long g0, g1, h0, h1; } Run;     // run [g0,g1) in list a, [h0,h1) in list b

// merge pass: runs of list a with W = w and runs of list b with W = target - w
static long match_runs(int la, int lb, int64_t target, Run **runs, long *cap, long *pairs) {
    long nr = 0; *pairs = 0;
    long j = M[lb] - 1, i = 0;
    while (i < M[la] && j >= 0) {
        int64_t wa = (int64_t)SK[la][i], want = target - wa;
        long i1 = i + 1; while (i1 < M[la] && (int64_t)SK[la][i1] == wa) i1++;
        while (j >= 0 && (int64_t)SK[lb][j] > want) j--;
        if (j < 0) break;
        if ((int64_t)SK[lb][j] == want) {
            long j1 = j + 1; long j0 = j; while (j0 > 0 && (int64_t)SK[lb][j0 - 1] == want) j0--;
            if (nr == *cap) { *cap = *cap ? 2 * *cap : 1 << 16; *runs = realloc(*runs, sizeof(Run) * *cap); }
            (*runs)[nr++] = (Run){i, i1, j0, j1}; *pairs += (i1 - i) * (j1 - j0);
            j = j0 - 1;
        }
        i = i1;
    }
    return nr;
}

static long verify_and_report(uint32_t i0, uint32_t i1, uint32_t i2, uint32_t i3, long *found) {
    for (int e = 0; e < NC_; e++)
        if (keyval(0, i0, e) + keyval(1, i1, e) + keyval(2, i2, e) + keyval(3, i3, e) != T0) return 0;
    uint32_t id[4] = {i0, i1, i2, i3};
    printf("SOL");
    for (int L = 0; L < 4; L++) { if (ROT[L]) printf(" %u %u", id[L] / C, id[L] % C); else printf(" %u 0", id[L]); }
    printf("\n"); fflush(stdout); (*found)++;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: rj INPUT NS NBUCKETS [MAXPAIRS] [FIRST] [LAST]\n"); return 2; }
    FILE *f = fopen(argv[1], "rb"); if (!f) { perror(argv[1]); return 1; }
    int32_t hdr[10]; if (fread(hdr, 4, 10, f) != 10) return 1;
    C = hdr[0]; T0 = hdr[1]; for (int i = 0; i < 4; i++) { ROT[i] = hdr[2 + i]; N[i] = hdr[6 + i]; }
    NS = atoi(argv[2]); NC_ = getenv("RJ_NC") ? atoi(getenv("RJ_NC")) : C; long NB = (long)atof(argv[3]);
    long MAXPAIRS = argc > 4 ? (long)atof(argv[4]) : 120000000L;
    long FIRST = argc > 5 ? atol(argv[5]) : 0, LAST = argc > 6 ? atol(argv[6]) : NB;
    if (NS > 9) { fprintf(stderr, "NS <= 9\n"); return 1; }
    for (int L = 0; L < 4; L++) {
        K[L] = malloc((size_t)N[L] * C); if (!K[L]) { fprintf(stderr, "oom K\n"); return 3; }
        if (fread(K[L], 1, (size_t)N[L] * C, f) != (size_t)N[L] * C) { fprintf(stderr, "short file\n"); return 1; }
        M[L] = ROT[L] ? N[L] * C : N[L];
    }
    fclose(f);
    uint64_t sd = 0x243F6A8885A308D3ULL;
    for (int e = 0; e < C; e++) { sd ^= sd << 13; sd ^= sd >> 7; sd ^= sd << 17; R[e] = sd | 1; }
    hT = 0; for (int e = NS; e < NC_; e++) hT += R[e] * (uint64_t)(int64_t)T0;
    double t0 = now();
    // statistics of the S classes (over entries) for the order of the split values
    double mu[4][9], var[4][9];
    for (int L = 0; L < 4; L++) for (int s = 0; s < NS; s++) {
        double a = 0, b = 0; long cnt = 0;
        for (long x = 0; x < N[L]; x++) for (int e = 0; e < C; e++) { /* rotation-symmetrised: use every class */ double v = K[L][(size_t)x * C + e]; a += v; b += v * v; cnt++; }
        if (!ROT[L]) { a = 0; b = 0; cnt = 0; for (long x = 0; x < N[L]; x++) { double v = K[L][(size_t)x * C + s]; a += v; b += v * v; cnt++; } }
        mu[L][s] = a / cnt; var[L][s] = b / cnt - mu[L][s] * mu[L][s];
    }
    for (int L = 0; L < 4; L++) {
        SK[L] = malloc(8 * M[L]); ID[L] = malloc(4 * M[L]);
        if (!SK[L] || !ID[L]) { fprintf(stderr, "oom lists\n"); return 3; }
        for (long i = 0; i < M[L]; i++) { ID[L][i] = (uint32_t)i; SK[L][i] = (uint64_t)wkey(L, (uint32_t)i); }
        radix_sort(SK[L], ID[L], M[L]);
        HS[L] = malloc(8 * M[L]); if (!HS[L]) { fprintf(stderr, "oom hs\n"); return 3; }
        for (long i = 0; i < M[L]; i++) HS[L][i] = hkey(L, ID[L][i]);
        fprintf(stderr, "rj: list %d: %ld necklaces, %ld entries sorted [%.1fs]\n", L, N[L], M[L], now() - t0);
    }
    // split values: candidates in a box around the centre, ordered by solution density
    double muAB[9], vAB[9], muCD[9], vCD[9];
    for (int s = 0; s < NS; s++) { muAB[s] = mu[0][s] + mu[1][s]; vAB[s] = var[0][s] + var[1][s]; muCD[s] = mu[2][s] + mu[3][s]; vCD[s] = var[2][s] + var[3][s]; }
    int Rb = 1; long nc;
    while (1) { nc = 1; for (int s = 0; s < NS; s++) nc *= (2 * Rb + 1); if (nc >= 2 * LAST || Rb >= 6) break; Rb++; }
    int *cand = malloc(sizeof(int) * nc * NS); double *score = malloc(sizeof(double) * nc); long *ord = malloc(sizeof(long) * nc);
    // joint centre: sigma maximising N(sigma; muAB, vAB) N(T - sigma; muCD, vCD)  (diagonal)
    double cen[9]; for (int s = 0; s < NS; s++) cen[s] = (muAB[s] / vAB[s] + (T0 - muCD[s]) / vCD[s]) / (1 / vAB[s] + 1 / vCD[s]);
    for (long t = 0; t < nc; t++) {
        long q = t; double sc = 0;
        for (int s = 0; s < NS; s++) {
            int k = (int)(q % (2 * Rb + 1)) - Rb; q /= (2 * Rb + 1);
            int sig = (int)lround(cen[s]) + k; cand[t * NS + s] = sig;
            sc += (sig - muAB[s]) * (sig - muAB[s]) / vAB[s] + (T0 - sig - muCD[s]) * (T0 - sig - muCD[s]) / vCD[s];
        }
        score[t] = sc; ord[t] = t;
    }
    int cmpd(const void *a, const void *b) { double x = score[*(const long *)a], y = score[*(const long *)b]; return (x > y) - (x < y); }
    qsort(ord, nc, sizeof(long), cmpd);
    long found = 0; long long totAB = 0, totCD = 0, totM = 0; Run *runs = NULL; long rcap = 0; Run *runs2 = NULL; long rcap2 = 0;
    for (long bi = FIRST; bi < LAST && bi < nc; bi++) {
        const int *sig = cand + ord[bi] * NS;
        int64_t wsig = 0; for (int s = NS - 1; s >= 0; s--) wsig = wsig * 64 + sig[s];
        int64_t wtg = 0; for (int s = NS - 1; s >= 0; s--) wtg = wtg * 64 + (T0 - sig[s]);
        long pAB, pCD;
        long n1 = match_runs(0, 1, wsig, &runs, &rcap, &pAB);
        long n2 = match_runs(2, 3, wtg, &runs2, &rcap2, &pCD);
        double tb = now();
        if (pAB > MAXPAIRS || pCD > 4 * MAXPAIRS) { fprintf(stderr, "BUCKET %ld skipped (pairs %ld %ld)\n", bi, pAB, pCD); continue; }
        long fb = 0, matches = 0;
        if (pAB && pCD) {
            // partitioned hash join: AB pairs built, CD pairs probed
            int pbits = 0; while ((pAB >> pbits) > 16384 && pbits < 18) pbits++;
            long P = 1L << pbits; int sh = 64 - pbits;
            long *pc = calloc(P + 1, sizeof(long));
            for (long u = 0; u < n1; u++) for (long a = runs[u].g0; a < runs[u].g1; a++) for (long b = runs[u].h0; b < runs[u].h1; b++) {
                uint64_t x = HS[0][a] + HS[1][b]; pc[(pbits ? mixp(x) >> sh : 0) + 1]++;
            }
            for (long p = 0; p < P; p++) pc[p + 1] += pc[p];
            long *pos = malloc(sizeof(long) * P); memcpy(pos, pc, sizeof(long) * P);
            uint64_t *ek = malloc(8 * pAB); uint32_t *ea = malloc(4 * pAB), *eb = malloc(4 * pAB);
            if (!ek || !ea || !eb) { fprintf(stderr, "rj: out of memory (%ld pairs)\n", pAB); return 3; }
            for (long u = 0; u < n1; u++) for (long a = runs[u].g0; a < runs[u].g1; a++) for (long b = runs[u].h0; b < runs[u].h1; b++) {
                uint64_t x = HS[0][a] + HS[1][b]; long i = pos[pbits ? mixp(x) >> sh : 0]++; ek[i] = x; ea[i] = ID[0][a]; eb[i] = ID[1][b];
            }
            // tables per partition
            long *tst = malloc(sizeof(long) * (P + 1)); long tt = 0;
            for (long p = 0; p < P; p++) { long n_ = pc[p + 1] - pc[p], c_ = 1; while (c_ < 2 * n_ + 2) c_ <<= 1; tst[p] = tt; tt += c_; }
            tst[P] = tt;
            uint32_t *tab = calloc(tt, 4);
            if (!tab) { fprintf(stderr, "rj: out of memory (table)\n"); return 3; }
            for (long p = 0; p < P; p++) {
                uint32_t *tb_ = tab + tst[p]; uint64_t m_ = tst[p + 1] - tst[p] - 1;
                for (long i = pc[p]; i < pc[p + 1]; i++) {
                    uint64_t slot = (ek[i] * 0x9E3779B97F4A7C15ULL) >> 20 & m_;
                    while (tb_[slot]) slot = (slot + 1) & m_;
                    tb_[slot] = (uint32_t)(i - pc[p] + 1);
                }
            }
            // probe in chunks of CD runs
            const long CH = 1L << 22;
            uint64_t *ck = malloc(8 * (CH + 1)); uint32_t *cc_ = malloc(4 * (CH + 1)), *cd_ = malloc(4 * (CH + 1));
            long *cpc = malloc(sizeof(long) * (P + 2)); long *cpos = malloc(sizeof(long) * P);
            long u0 = 0;
            while (u0 < n2) {
                long u1 = u0, np_ = 0;
                while (u1 < n2) { long k_ = (runs2[u1].g1 - runs2[u1].g0) * (runs2[u1].h1 - runs2[u1].h0); if (np_ && np_ + k_ > CH) break; np_ += k_; u1++; }
                if (np_ > CH) { ck = realloc(ck, 8 * np_); cc_ = realloc(cc_, 4 * np_); cd_ = realloc(cd_, 4 * np_); }
                memset(cpc, 0, sizeof(long) * (P + 2));
                for (long u = u0; u < u1; u++) for (long c = runs2[u].g0; c < runs2[u].g1; c++) { uint64_t xc = hT - HS[2][c];
                    for (long d = runs2[u].h0; d < runs2[u].h1; d++) { uint64_t x = xc - HS[3][d]; cpc[(pbits ? mixp(x) >> sh : 0) + 1]++; } }
                for (long p = 0; p < P; p++) cpc[p + 1] += cpc[p];
                memcpy(cpos, cpc, sizeof(long) * P);
                for (long u = u0; u < u1; u++) for (long c = runs2[u].g0; c < runs2[u].g1; c++) { uint64_t xc = hT - HS[2][c];
                    for (long d = runs2[u].h0; d < runs2[u].h1; d++) { uint64_t x = xc - HS[3][d]; long i = cpos[pbits ? mixp(x) >> sh : 0]++; ck[i] = x; cc_[i] = ID[2][c]; cd_[i] = ID[3][d]; } }
                for (long p = 0; p < P; p++) {
                    const uint32_t *tb_ = tab + tst[p]; uint64_t m_ = tst[p + 1] - tst[p] - 1;
                    for (long i = cpc[p]; i < cpc[p + 1]; i++) {
                        uint64_t x = ck[i]; uint64_t slot = (x * 0x9E3779B97F4A7C15ULL) >> 20 & m_;
                        while (tb_[slot]) {
                            long j = pc[p] + tb_[slot] - 1;
                            if (ek[j] == x) { matches++; fb += verify_and_report(ea[j], eb[j], cc_[i], cd_[i], &found); }
                            slot = (slot + 1) & m_;
                        }
                    }
                }
                u0 = u1;
            }
            free(ck); free(cc_); free(cd_); free(cpc); free(cpos);
            free(tab); free(tst); free(ek); free(ea); free(eb); free(pc); free(pos);
        }
        totAB += pAB; totCD += pCD; totM += matches;
        fprintf(stderr, "BUCKET %ld sigma", bi); for (int s = 0; s < NS; s++) fprintf(stderr, " %d", sig[s]);
        fprintf(stderr, " pairsAB %ld pairsCD %ld matches %ld found %ld (total %ld) join %.1fs cpu %.0fs\n", pAB, pCD, matches, fb, found, now() - tb, now() - t0);
    }
    printf("TOTAL found %ld pairsAB %lld pairsCD %lld hashmatches %lld cpu %.1f\n", found, totAB, totCD, totM, now() - t0);
    return 0;
}
