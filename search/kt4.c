/*
 * Four-list k-tree (Wagner) for class-value vectors: find (a, b, c, d) with VA[a] + VB[b] + VC[c] + VD[d] = T on all
 * D classes, among sampled lists.
 *
 * For each split value sigma of the merge classes S (|S| <= 5, values within +-2047), the pairs (a, b) with (VA[a] + VB[b])_S = sigma and
 * the pairs (c, d) with (VC[c] + VD[d])_S = (T - sigma)_S are formed group by group: every list is sorted by its S-tuple, and
 * the matching groups are found by one merge pass (or by binary search when one list has few groups), so the pairs of two
 * groups read contiguous memory. The pairs are joined on a 64-bit linear hash of the other classes; hash matches are
 * re-checked on all classes. (The set of solutions found for a given sigma does not depend on this organisation.) The values of
 * sigma are taken in order of increasing distance from the centre mu_A + mu_B + (T - mu)/2 (on S), stepping by the
 * lattice of each merge class (the gcd of the differences of its values within each list, combined).
 *
 * input: int32 nA nB nC nD D; int16 T[D]; int16 VA[nA][D] VB[nB][D] VC[nC][D] VD[nD][D]
 *   cc -O3 -march=native -fopenmp kt4.c -o kt4 ;  ./kt4 input.bin S0,S1,.. NSIG [K]
 * NSIG may also be @file: the split values are then read from the file (one line of |S| integers per value, in the
 * order given) instead of being generated around the centre.
 * Prints "SOL a b c d" (at most K, default 1e9), progress on stderr, and "TOTAL found pairsAB pairsCD".
 */
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

static int D, NS, Sx[5];
static long MAXPAIRS = 300000000L;
static int16_t *T, *V[4];   /* point into the mapped input */
static long n[4];
static uint64_t Rh[256], hT;

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void rd(void *p, size_t sz, size_t k, FILE *f) { if (fread(p, sz, k, f) != k) { fprintf(stderr, "kt4: short read\n"); exit(1); } }

static inline uint64_t pk(int x) { return (uint64_t)((x + 2048) & 0xFFF); }
static inline uint64_t skey(const int16_t *v, const int *off) {          /* pack the S-values (12 bits each) */
    uint64_t k = 0;
    for (int s = 0; s < NS; s++) k = (k << 12) | pk(v[Sx[s]] + off[s]);
    return k;
}
static inline int isS(int c) { for (int s = 0; s < NS; s++) if (Sx[s] == c) return 1; return 0; }

typedef struct { long m; uint64_t *key; uint32_t *start; uint32_t *idx; } Index;      /* list sorted by S-key, CSR */

static int cmpk_base;
static uint64_t *cmp_keys;
static int cmpk(const void *a, const void *b) {
    uint64_t x = cmp_keys[*(const uint32_t *)a], y = cmp_keys[*(const uint32_t *)b];
    return x < y ? -1 : x > y;
}

static Index build_index(int L, const int *off) {
    Index I; long N = n[L];
    uint64_t *k = malloc(8 * N); uint32_t *o = malloc(sizeof(uint32_t) * N);
    for (long i = 0; i < N; i++) { k[i] = skey(V[L] + (size_t)D * i, off); o[i] = (uint32_t)i; }
    cmp_keys = k; qsort(o, N, sizeof(uint32_t), cmpk);
    I.idx = o; I.key = malloc(8 * (N + 1)); I.start = malloc(sizeof(uint32_t) * (N + 1)); I.m = 0;
    for (long i = 0; i < N; i++) {
        if (i == 0 || k[o[i]] != k[o[i - 1]]) { I.key[I.m] = k[o[i]]; I.start[I.m] = i; I.m++; }
    }
    I.start[I.m] = (uint32_t)N; free(k);
    I.key = realloc(I.key, 8 * (I.m + 1)); I.start = realloc(I.start, sizeof(uint32_t) * (I.m + 1));
    return I;
}

static long find(const Index *I, uint64_t key) {                         /* group index or -1 */
    long lo = 0, hi = I->m - 1;
    while (lo <= hi) { long mid = (lo + hi) >> 1; if (I->key[mid] < key) lo = mid + 1; else if (I->key[mid] > key) hi = mid - 1; else return mid; }
    return -1;
}


#ifndef RADIX_MIN
#define RADIX_MIN (1L << 18)
#endif
static inline uint64_t mixp(uint64_t x) { return x * 0xD6E8FEB86659FD93ULL; }
/* Partitioned hash join of the AB pairs (built) and the CD pairs (probed) of one split value. The pairs are scattered
 * into P = 2^pbits partitions by the top bits of a mixed hash; each partition's table holds about 2^15 entries, so
 * probes stay in cache. The CD pairs are processed in chunks of about 2^23. */
static void radix_join(const long *nm, const long *mp, long cnt, const Index *IA, const Index *IB, const Index *IC,
                       const Index *ID, uint64_t *const *hs, long long *ncd_out, long *fnd_out, long *found, long K) {
    int pbits = 0; while ((cnt >> pbits) > 32768 && pbits < 16) pbits++;
    const long P = 1L << pbits; const int sh = 64 - pbits;
    int nt = omp_get_max_threads();
    /* AB: count, offsets, fill */
    long *pc = calloc((size_t)nt * P, sizeof(long));
    #pragma omp parallel
    {
        long *my = pc + (size_t)omp_get_thread_num() * P;
        #pragma omp for schedule(static)
        for (long u = 0; u < nm[0]; u++) {
            long ga = mp[2 * u], gb = mp[2 * u + 1];
            for (long pa = IA->start[ga]; pa < IA->start[ga + 1]; pa++)
                for (long pb = IB->start[gb]; pb < IB->start[gb + 1]; pb++) {
                    uint64_t x = hs[0][pa] + hs[1][pb]; if (x == ~0ULL) x--;
                    my[pbits ? mixp(x) >> sh : 0]++;
                }
        }
    }
    long *off = malloc(sizeof(long) * (size_t)nt * P), *pstart = malloc(sizeof(long) * (P + 1)); long run = 0;
    for (long p = 0; p < P; p++) { pstart[p] = run; for (int t = 0; t < nt; t++) { off[(size_t)t * P + p] = run; run += pc[(size_t)t * P + p]; } }
    pstart[P] = run;
    uint64_t *ek = malloc(8 * run); uint32_t *ea = malloc(4 * run), *eb = malloc(4 * run);
    if (!ek || !ea || !eb) { fprintf(stderr, "kt4: out of memory (%ld pairs)\n", run); exit(3); }
    #pragma omp parallel
    {
        long *my = off + (size_t)omp_get_thread_num() * P;
        #pragma omp for schedule(static)
        for (long u = 0; u < nm[0]; u++) {
            long ga = mp[2 * u], gb = mp[2 * u + 1];
            for (long pa = IA->start[ga]; pa < IA->start[ga + 1]; pa++)
                for (long pb = IB->start[gb]; pb < IB->start[gb + 1]; pb++) {
                    uint64_t x = hs[0][pa] + hs[1][pb]; if (x == ~0ULL) x--;
                    long i = my[pbits ? mixp(x) >> sh : 0]++;
                    ek[i] = x; ea[i] = IA->idx[pa]; eb[i] = IB->idx[pb];
                }
        }
    }
    /* per-partition tables of entry indices (+1; 0 = empty) */
    long *tstart = malloc(sizeof(long) * (P + 1)); long tt = 0;
    for (long p = 0; p < P; p++) { long n_ = pstart[p + 1] - pstart[p], c_ = 1; while (c_ < 2 * n_ + 2) c_ <<= 1; tstart[p] = tt; tt += c_; }
    tstart[P] = tt;
    uint32_t *tab = calloc(tt, sizeof(uint32_t));
    if (!tab) { fprintf(stderr, "kt4: out of memory (table of %ld)\n", tt); exit(3); }
    #pragma omp parallel for schedule(dynamic, 1)
    for (long p = 0; p < P; p++) {
        uint32_t *tb = tab + tstart[p]; uint64_t m_ = tstart[p + 1] - tstart[p] - 1;
        for (long i = pstart[p]; i < pstart[p + 1]; i++) {
            uint64_t slot = (ek[i] * 0x9E3779B97F4A7C15ULL) & m_;
            while (tb[slot]) slot = (slot + 1) & m_;
            tb[slot] = (uint32_t)(i - pstart[p] + 1);
        }
    }
    /* CD in chunks */
    long long ncd = 0; long fnd = 0;
    const long CH = 1L << 23;
    long u0 = nm[0], uend = nm[0] + nm[1];
    uint64_t *ck = malloc(8 * (CH + (1L << 20))); uint32_t *cc_ = malloc(4 * (CH + (1L << 20))), *cd_ = malloc(4 * (CH + (1L << 20)));
    if (!ck || !cc_ || !cd_) { fprintf(stderr, "kt4: out of memory (CD chunk)\n"); exit(3); }
    while (u0 < uend) {
        /* chunk [u0, u1): at most about CH pairs (a single group pair larger than CH is taken alone) */
        long u1 = u0; long npair = 0;
        while (u1 < uend) {
            long gc = mp[2 * u1], gd = mp[2 * u1 + 1];
            long k_ = (long)(IC->start[gc + 1] - IC->start[gc]) * (ID->start[gd + 1] - ID->start[gd]);
            if (npair && npair + k_ > CH) break;
            npair += k_; u1++;
        }
        if (npair > CH + (1L << 20)) {                        /* a huge single group pair: enlarge the buffers */
            ck = realloc(ck, 8 * npair); cc_ = realloc(cc_, 4 * npair); cd_ = realloc(cd_, 4 * npair);
            if (!ck || !cc_ || !cd_) { fprintf(stderr, "kt4: out of memory (%ld CD pairs)\n", npair); exit(3); }
        }
        memset(pc, 0, sizeof(long) * (size_t)nt * P);
        #pragma omp parallel
        {
            long *my = pc + (size_t)omp_get_thread_num() * P;
            #pragma omp for schedule(static)
            for (long u = u0; u < u1; u++) {
                long gc = mp[2 * u], gd = mp[2 * u + 1];
                for (long pcx = IC->start[gc]; pcx < IC->start[gc + 1]; pcx++) {
                    uint64_t xc = hT - hs[2][pcx];
                    for (long pd = ID->start[gd]; pd < ID->start[gd + 1]; pd++) {
                        uint64_t x = xc - hs[3][pd]; if (x == ~0ULL) x--;
                        my[pbits ? mixp(x) >> sh : 0]++;
                    }
                }
            }
        }
        long *cstart = malloc(sizeof(long) * (P + 1)); run = 0;
        for (long p = 0; p < P; p++) { cstart[p] = run; for (int t = 0; t < nt; t++) { off[(size_t)t * P + p] = run; run += pc[(size_t)t * P + p]; } }
        cstart[P] = run;
        #pragma omp parallel
        {
            long *my = off + (size_t)omp_get_thread_num() * P;
            #pragma omp for schedule(static)
            for (long u = u0; u < u1; u++) {
                long gc = mp[2 * u], gd = mp[2 * u + 1];
                for (long pcx = IC->start[gc]; pcx < IC->start[gc + 1]; pcx++) {
                    uint64_t xc = hT - hs[2][pcx];
                    for (long pd = ID->start[gd]; pd < ID->start[gd + 1]; pd++) {
                        uint64_t x = xc - hs[3][pd]; if (x == ~0ULL) x--;
                        long i = my[pbits ? mixp(x) >> sh : 0]++;
                        ck[i] = x; cc_[i] = IC->idx[pcx]; cd_[i] = ID->idx[pd];
                    }
                }
            }
        }
        ncd += run;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:fnd)
        for (long p = 0; p < P; p++) {
            const uint32_t *tb = tab + tstart[p]; uint64_t m_ = tstart[p + 1] - tstart[p] - 1;
            const uint64_t *pk_ = ek + pstart[p];
            for (long i = cstart[p]; i < cstart[p + 1]; i++) {
                uint64_t x = ck[i];
                uint64_t slot = (x * 0x9E3779B97F4A7C15ULL) & m_;
                while (tb[slot]) {
                    long j = tb[slot] - 1;
                    if (pk_[j] == x) {
                        long a = ea[pstart[p] + j], b = eb[pstart[p] + j], c = cc_[i], d = cd_[i];
                        const int16_t *va = V[0] + (size_t)D * a, *vb = V[1] + (size_t)D * b, *vc = V[2] + (size_t)D * c, *vd = V[3] + (size_t)D * d;
                        int ok = 1; for (int q = 0; q < D; q++) if (va[q] + vb[q] + vc[q] + vd[q] != T[q]) { ok = 0; break; }
                        if (ok) {
                            fnd++;
                            #pragma omp critical
                            { if (*found < K) { printf("SOL %ld %ld %ld %ld\n", a, b, c, d); fflush(stdout); } (*found)++; }
                        }
                    }
                    slot = (slot + 1) & m_;
                }
            }
        }
        free(cstart);
        u0 = u1;
    }
    free(ck); free(cc_); free(cd_); free(tab); free(tstart); free(ek); free(ea); free(eb); free(off); free(pstart); free(pc);
    *ncd_out = ncd; *fnd_out = fnd;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: kt4 input.bin S0,S1,.. NSIG [K]\n"); return 2; }
    /* the lists are memory-mapped from the input file (read-only, paged in on demand) */
    int fdin = open(argv[1], O_RDONLY); if (fdin < 0) { perror(argv[1]); return 1; }
    struct stat stt; fstat(fdin, &stt);
    const char *mbase = mmap(NULL, stt.st_size, PROT_READ, MAP_PRIVATE, fdin, 0);
    if (mbase == MAP_FAILED) { perror("mmap"); return 1; }
    { const int32_t *hdr = (const int32_t *)mbase; for (int i = 0; i < 4; i++) n[i] = hdr[i]; D = hdr[4]; }
    T = (int16_t *)(mbase + 20);
    { int16_t *p = T + D; for (int i = 0; i < 4; i++) { V[i] = p; p += (size_t)D * n[i]; } }
    if ((size_t)(20 + 2 * D) + 2 * (size_t)D * (n[0] + n[1] + n[2] + n[3]) > (size_t)stt.st_size) { fprintf(stderr, "kt4: short file\n"); return 1; }
    for (int i = 0; i < 4; i++) if (n[i] >= 4294967295L) { fprintf(stderr, "kt4: list too long\n"); return 1; }
    NS = 0; for (char *p = strtok(argv[2], ","); p && NS < 5; p = strtok(NULL, ",")) Sx[NS++] = atoi(p);
    long NSIG = argv[3][0] == '@' ? 0 : atol(argv[3]); long K = argc > 4 ? atol(argv[4]) : 1000000000L;
    if (argc > 5) MAXPAIRS = atol(argv[5]);
    uint64_t seed = 0x243F6A8885A308D3ULL;
    for (int c = 0; c < D; c++) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; Rh[c] = isS(c) ? 0 : (seed | 1); }
    /* per-list hash over the non-S classes */
    uint64_t *h[4];
    for (int L = 0; L < 4; L++) {
        h[L] = malloc(8 * n[L]);
        #pragma omp parallel for
        for (long i = 0; i < n[L]; i++) { uint64_t x = 0; const int16_t *v = V[L] + (size_t)D * i; for (int c = 0; c < D; c++) x += Rh[c] * (uint64_t)(int64_t)v[c]; h[L][i] = x; }
    }
    hT = 0; for (int c = 0; c < D; c++) hT += Rh[c] * (uint64_t)(int64_t)T[c];
    /* means, lattices and offsets for the S classes */
    double mu[4][5]; int g[5], vmin[5], vmax[5];
    for (int s = 0; s < NS; s++) {
        g[s] = 0; vmin[s] = 1 << 30; vmax[s] = -(1 << 30);
        for (int L = 0; L < 4; L++) {
            double acc = 0; int first = V[L][Sx[s]];
            for (long i = 0; i < n[L]; i++) {
                int x = V[L][(size_t)D * i + Sx[s]]; acc += x;
                int d = abs(x - first); while (d) { int t = g[s] % d; g[s] = d; d = t; }
                if (x < vmin[s]) vmin[s] = x; if (x > vmax[s]) vmax[s] = x;
            }
            mu[L][s] = acc / n[L];
        }
        if (g[s] == 0) g[s] = 1;
    }
    int off[5]; for (int s = 0; s < NS; s++) {
        off[s] = 0;
        if (vmin[s] < -1000 || vmax[s] > 1000) { fprintf(stderr, "kt4: class %d out of range\n", Sx[s]); return 1; }
    }
    Index IB = build_index(1, off), ID = build_index(3, off);
    fprintf(stderr, "kt4: lists %ld %ld %ld %ld, D=%d, S=", n[0], n[1], n[2], n[3], D);
    for (int s = 0; s < NS; s++) fprintf(stderr, "%d(g=%d) ", Sx[s], g[s]);
    fprintf(stderr, ", groups B %ld D %ld\n", IB.m, ID.m);
    /* candidate sigmas: lattice points near the centre, by increasing normalised distance */
    double cen[5], sd[5];
    for (int s = 0; s < NS; s++) {
        double tot = 0; for (int L = 0; L < 4; L++) tot += mu[L][s];
        cen[s] = mu[0][s] + mu[1][s] + 0.5 * (T[Sx[s]] - tot);
        sd[s] = 1.0;
    }
    long nc; int *cand; long *ord;
    if (argv[3][0] == '@') {                                     /* explicit split values */
        FILE *sf = fopen(argv[3] + 1, "r"); if (!sf) { perror(argv[3] + 1); return 1; }
        long capc = 1024; cand = malloc(sizeof(int) * capc * NS); nc = 0;
        while (1) {
            if (nc == capc) { capc *= 2; cand = realloc(cand, sizeof(int) * capc * NS); }
            int ok = 1; for (int s = 0; s < NS; s++) if (fscanf(sf, "%d", &cand[nc * NS + s]) != 1) { ok = 0; break; }
            if (!ok) break;
            nc++;
        }
        fclose(sf); NSIG = nc; ord = malloc(sizeof(long) * (nc + 1)); for (long t = 0; t < nc; t++) ord[t] = t;
        goto search;
    }
    int R = 1; long ncand = 1; while (ncand < 4 * NSIG) { R++; ncand = 1; for (int s = 0; s < NS; s++) ncand *= (2 * R + 1); }
    nc = 1; for (int s = 0; s < NS; s++) nc *= (2 * R + 1);
    cand = malloc(sizeof(int) * nc * NS); double *dist = malloc(sizeof(double) * nc); ord = malloc(sizeof(long) * nc);
    /* residue of sigma_s: sum of a value of A and of B is congruent to (vA + vB) mod g */
    int base[5]; for (int s = 0; s < NS; s++) base[s] = V[0][Sx[s]] + V[1][Sx[s]];
    for (long t = 0; t < nc; t++) {
        long q = t; double dd = 0;
        for (int s = 0; s < NS; s++) {
            int k = (int)(q % (2 * R + 1)) - R; q /= (2 * R + 1);
            long c0 = lround((cen[s] - base[s]) / g[s]); int sig = base[s] + (int)(c0 + k) * g[s];
            cand[t * NS + s] = sig; dd += (sig - cen[s]) * (sig - cen[s]);
        }
        dist[t] = dd; ord[t] = t;
    }
    /* sort candidates by distance (simple insertion into an index array by qsort on dist) */
    { double *dk = dist; for (long i = 1; i < nc; i++) { long x = ord[i]; long j = i - 1; while (j >= 0 && dk[ord[j]] > dk[x]) { ord[j + 1] = ord[j]; j--; } ord[j + 1] = x; } }
  search:;
    /* All four lists indexed by their S-keys (groups in increasing key order). The S-key order is the lexicographic
     * order of the S-tuples, so for a split value sigma the partners of the A-groups in increasing order are the
     * B-groups with keys sigma - vA in decreasing order: one merge pass per side instead of a binary search per
     * member. The hashes are copied in index order, so the pairs of two groups read contiguous memory. */
    Index IA = build_index(0, off), IC = build_index(2, off);
    uint64_t *hs[4]; Index *IX[4] = {&IA, &IB, &IC, &ID};
    for (int L = 0; L < 4; L++) {
        hs[L] = malloc(8 * n[L]);
        for (long i = 0; i < n[L]; i++) hs[L][i] = h[L][IX[L]->idx[i]];
        free(h[L]); h[L] = NULL;                                  /* only the sorted copies are used below */
    }
    /* the S-tuple of each group (decoded from its key) */
    int16_t *tup[4];
    for (int L = 0; L < 4; L++) {
        tup[L] = malloc(sizeof(int16_t) * NS * (IX[L]->m + 1));
        for (long gi = 0; gi < IX[L]->m; gi++) {
            uint64_t k = IX[L]->key[gi];
            for (int s = NS - 1; s >= 0; s--) { tup[L][gi * NS + s] = (int16_t)((int)(k & 0xFFF) - 2048); k >>= 12; }
        }
    }
    long found = 0; long long totAB = 0, totCD = 0; double t0 = now();
    long *mp = NULL; long mcap = 0;                                   /* matched group pairs (g1, g2) */
    for (long ci = 0; ci < NSIG && ci < nc && found < K; ci++) {
        const int *sig = cand + ord[ci] * NS;
        int tgt[5]; for (int s = 0; s < NS; s++) tgt[s] = T[Sx[s]] - sig[s];
        /* match groups: (L1, L2) = (A, B) with vA + vB = sig, then (C, D) with vC + vD = T - sig */
        long nm[2] = {0, 0}; long start[2] = {0, 0}; long long npairs[2] = {0, 0};
        for (int side = 0; side < 2; side++) {
            const Index *I1 = IX[2 * side], *I2 = IX[2 * side + 1]; const int16_t *t1 = tup[2 * side], *t2 = tup[2 * side + 1];
            const int *want = side == 0 ? sig : tgt;
            long j = I2->m - 1; start[side] = nm[0];
            if ((double)I1->m * log2((double)I2->m + 2) < (double)(I1->m + I2->m)) {   /* few groups: binary search */
                for (long gi = 0; gi < I1->m; gi++) {
                    uint64_t k = 0; int ok = 1;
                    for (int s = 0; s < NS; s++) {
                        int part = want[s] - t1[gi * NS + s];
                        if (part < -2048 || part > 2047) { ok = 0; break; }
                        k = (k << 12) | pk(part);
                    }
                    if (!ok) continue;
                    long jj = find(I2, k); if (jj < 0) continue;
                    long idx_ = side == 0 ? nm[0] : nm[0] + nm[1];
                    if (2 * (idx_ + 1) > mcap) { mcap = mcap ? 2 * mcap : 1 << 20; mp = realloc(mp, sizeof(long) * mcap); }
                    mp[2 * idx_] = gi; mp[2 * idx_ + 1] = jj; nm[side]++;
                    npairs[side] += (long long)(I1->start[gi + 1] - I1->start[gi]) * (I2->start[jj + 1] - I2->start[jj]);
                }
                continue;
            }
            for (long gi = 0; gi < I1->m && j >= 0; gi++) {
                /* key of the partner: want - t1[gi]; advance j (decreasing) while key(t2[j]) > partner */
                int part[5]; int ok = 1;
                for (int s = 0; s < NS; s++) { part[s] = want[s] - t1[gi * NS + s]; if (part[s] < -2048 || part[s] > 2047) ok = 0; }
                if (!ok) continue;
                while (j >= 0) {
                    int cmp = 0;
                    for (int s = 0; s < NS; s++) { int x = t2[j * NS + s]; if (x != part[s]) { cmp = x > part[s] ? 1 : -1; break; } }
                    if (cmp > 0) { j--; continue; }
                    if (cmp == 0) {
                        long idx_ = side == 0 ? nm[0] : nm[0] + nm[1];
                        if (2 * (idx_ + 1) > mcap) { mcap = mcap ? 2 * mcap : 1 << 20; mp = realloc(mp, sizeof(long) * mcap); }
                        mp[2 * idx_] = gi; mp[2 * idx_ + 1] = j; nm[side]++;
                        npairs[side] += (long long)(I1->start[gi + 1] - I1->start[gi]) * (I2->start[j + 1] - I2->start[j]);
                    }
                    break;
                }
            }
        }
        long cnt = npairs[0];
        if (cnt > MAXPAIRS) { fprintf(stderr, "kt4: sigma %ld: %ld pairs exceed the cap, skipped\n", ci + 1, cnt); continue; }
        long long ncd = 0; long fnd = 0;
        if (cnt >= RADIX_MIN) {                                       /* partitioned join: every table fits in cache */
            radix_join(nm, mp, cnt, &IA, &IB, &IC, &ID, hs, &ncd, &fnd, &found, K);
        } else {
            long cap = 1; while (cap < 2 * cnt + 2) cap <<= 1;
            uint64_t *hk = malloc(8 * cap); uint32_t *ha = malloc(4 * cap), *hb = malloc(4 * cap);
            if (!hk || !ha || !hb) { fprintf(stderr, "kt4: out of memory at sigma %ld\n", ci + 1); return 3; }
            memset(hk, 0xff, 8 * cap);
            uint64_t mask = cap - 1;
            /* threads only when there is work: with many split values of a few pairs each, the fork/join of the
             * parallel regions (and spinning threads on a loaded machine) cost far more than the pairs */
            #pragma omp parallel for schedule(dynamic, 64) if(cnt > 65536)
            for (long u = 0; u < nm[0]; u++) {
                long ga = mp[2 * u], gb = mp[2 * u + 1];
                for (long pa = IA.start[ga]; pa < IA.start[ga + 1]; pa++) {
                    uint64_t xa = hs[0][pa];
                    for (long pb = IB.start[gb]; pb < IB.start[gb + 1]; pb++) {
                        uint64_t x = xa + hs[1][pb]; if (x == ~0ULL) x--;
                        uint64_t slot = (x * 0x9E3779B97F4A7C15ULL) & mask;
                        while (1) {
                            uint64_t cur = __atomic_load_n(&hk[slot], __ATOMIC_RELAXED);
                            if (cur == ~0ULL && __atomic_compare_exchange_n(&hk[slot], &cur, x, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                                ha[slot] = (uint32_t)IA.idx[pa]; hb[slot] = (uint32_t)IB.idx[pb]; break;
                            }
                            slot = (slot + 1) & mask;
                        }
                    }
                }
            }
            #pragma omp parallel for schedule(dynamic, 64) reduction(+:ncd, fnd) if(npairs[1] > 65536)
            for (long u = nm[0]; u < nm[0] + nm[1]; u++) {
                long gc = mp[2 * u], gd = mp[2 * u + 1];
                for (long pc = IC.start[gc]; pc < IC.start[gc + 1]; pc++) {
                    uint64_t xc = hT - hs[2][pc];
                    for (long pd = ID.start[gd]; pd < ID.start[gd + 1]; pd++) {
                        ncd++;
                        uint64_t x = xc - hs[3][pd]; if (x == ~0ULL) x--;
                        uint64_t slot = (x * 0x9E3779B97F4A7C15ULL) & mask;
                        while (hk[slot] != ~0ULL) {
                            if (hk[slot] == x) {
                                long a = ha[slot], b = hb[slot], c = IC.idx[pc], d = ID.idx[pd];
                                const int16_t *va = V[0] + (size_t)D * a, *vb = V[1] + (size_t)D * b, *vc = V[2] + (size_t)D * c, *vd = V[3] + (size_t)D * d;
                                int ok = 1; for (int q = 0; q < D; q++) if (va[q] + vb[q] + vc[q] + vd[q] != T[q]) { ok = 0; break; }
                                if (ok) {
                                    fnd++;
                                    #pragma omp critical
                                    { if (found < K) { printf("SOL %ld %ld %ld %ld\n", a, b, c, d); fflush(stdout); } found++; }
                                }
                            }
                            slot = (slot + 1) & mask;
                        }
                    }
                }
            }
            free(hk); free(ha); free(hb);
        }
        totAB += cnt; totCD += ncd;
        fprintf(stderr, "kt4: sigma %ld/%ld (", ci + 1, NSIG);
        for (int s = 0; s < NS; s++) fprintf(stderr, "%d%s", sig[s], s + 1 < NS ? "," : "");
        fprintf(stderr, "): pairs AB %ld CD %lld, found %ld (total %ld) [%.0fs]\n", cnt, ncd, fnd, found, now() - t0);
    }
    printf("TOTAL %ld %lld %lld\n", found, totAB, totCD);
    return 0;
}
