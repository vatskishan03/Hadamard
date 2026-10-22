// exact 4-list zero-sum join with linear 64-bit hashing, (A,C) x (B,D), R passes by low hash bits.
// vectors are int16, row-major, d columns. returns number of solutions found (may exceed maxout; only maxout stored).
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static inline uint64_t hv(const int16_t *v, int d, const uint64_t *coef) {
    uint64_t h = 0;
    for (int i = 0; i < d; i++) h += coef[i] * (uint64_t)(int64_t)v[i];
    return h;
}

typedef struct { uint64_t key; uint32_t a, c; } ent;
static inline int clsf(const int16_t *v, int d, const int64_t *w, int R) { int64_t a_ = 0; for (int j_ = 0; j_ < d; j_++) a_ += w[j_] * (int64_t)v[j_]; a_ %= R; if (a_ < 0) a_ += R; return (int)a_; }

long mitm4(int nA, int nB, int nC, int nD, int d,
           const int16_t *VA, const int16_t *VB, const int16_t *VC, const int16_t *VD,
           const int16_t *target, int logR, uint32_t *out, long maxout, uint64_t seed, int verbose)
{
    uint64_t coef[64];
    uint64_t s = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    for (int i = 0; i < d; i++) { s = s * 6364136223846793005ULL + 1442695040888963407ULL; coef[i] = (s >> 7) | 1ULL; }
    uint64_t *hA = malloc(sizeof(uint64_t) * nA), *hB = malloc(sizeof(uint64_t) * nB);
    uint64_t *hC = malloc(sizeof(uint64_t) * nC), *hD = malloc(sizeof(uint64_t) * nD);
    for (int i = 0; i < nA; i++) hA[i] = hv(VA + (long)i * d, d, coef);
    for (int i = 0; i < nB; i++) hB[i] = hv(VB + (long)i * d, d, coef);
    for (int i = 0; i < nC; i++) hC[i] = hv(VC + (long)i * d, d, coef);
    for (int i = 0; i < nD; i++) hD[i] = hv(VD + (long)i * d, d, coef);
    uint64_t ht = hv(target, d, coef);
    int R = logR; if (R < 1) R = 1;   // number of classes (any integer)
    int64_t wcls[64]; { uint64_t s2 = seed * 2862933555777941757ULL + 3037000493ULL; for (int i = 0; i < d; i++) { s2 = s2 * 6364136223846793005ULL + 1442695040888963407ULL; wcls[i] = 1 + (int64_t)((s2 >> 33) % (R > 1 ? (R - 1) : 1)); } }
    int *kA = malloc(sizeof(int) * nA), *kB = malloc(sizeof(int) * nB), *kC = malloc(sizeof(int) * nC), *kD = malloc(sizeof(int) * nD);
    #define CLS(v) clsf((v), d, wcls, R)
    for (int i = 0; i < nA; i++) kA[i] = CLS(VA + (long)i * d);
    for (int i = 0; i < nB; i++) kB[i] = CLS(VB + (long)i * d);
    for (int i = 0; i < nC; i++) kC[i] = CLS(VC + (long)i * d);
    for (int i = 0; i < nD; i++) kD[i] = CLS(VD + (long)i * d);
    int kT = CLS(target);
    // bucket C and D by low bits of their hashes
    uint32_t *cC = calloc(R + 1, sizeof(uint32_t)), *cD = calloc(R + 1, sizeof(uint32_t));
    for (int i = 0; i < nC; i++) cC[kC[i] + 1]++;
    for (int i = 0; i < nD; i++) cD[kD[i] + 1]++;
    for (int r = 0; r < R; r++) { cC[r + 1] += cC[r]; cD[r + 1] += cD[r]; }
    uint32_t *oC = malloc(sizeof(uint32_t) * (nC + 1)), *oD = malloc(sizeof(uint32_t) * (nD + 1));
    { uint32_t *pos = malloc(sizeof(uint32_t) * (R + 1)); memcpy(pos, cC, sizeof(uint32_t) * (R + 1));
      for (int i = 0; i < nC; i++) oC[pos[kC[i]]++] = i;
      memcpy(pos, cD, sizeof(uint32_t) * (R + 1));
      for (int i = 0; i < nD; i++) oD[pos[kD[i]]++] = i; free(pos); }
    long found = 0;
    for (int r = 0; r < R; r++) {
        // AC pairs with (hA+hC)&Rm == r
        long npairs = 0;
        for (int a = 0; a < nA; a++) { int need = ((r - kA[a]) % R + R) % R; npairs += cC[need + 1] - cC[need]; }
        if (npairs == 0) continue;
        long cap = 1; while (cap < npairs * 2) cap <<= 1;
        ent *tab = malloc(sizeof(ent) * cap);
        if (!tab) { fprintf(stderr, "oom table %ld\n", cap); return -1; }
        memset(tab, 0xff, sizeof(ent) * cap);   // key = all ones marks empty with a=c=0xffffffff
        uint64_t cm = cap - 1;
        for (int a = 0; a < nA; a++) {
            int need = ((r - kA[a]) % R + R) % R;
            for (uint32_t p = cC[need]; p < cC[need + 1]; p++) {
                uint32_t c = oC[p]; uint64_t key = hA[a] + hC[c];
                uint64_t i = (key * 0x9E3779B97F4A7C15ULL) >> 20; i &= cm;
                while (tab[i].a != 0xffffffffu) i = (i + 1) & cm;
                tab[i].key = key; tab[i].a = a; tab[i].c = c;
            }
        }
        // probe: BD pairs with (ht - hB - hD) & Rm == r
        for (int b = 0; b < nB; b++) {
            int need = (((kT - r - kB[b]) % R) + R) % R;   // kD must equal need
            for (uint32_t p = cD[need]; p < cD[need + 1]; p++) {
                uint32_t dd = oD[p]; uint64_t key = ht - hB[b] - hD[dd];
                uint64_t i = (key * 0x9E3779B97F4A7C15ULL) >> 20; i &= cm;
                while (tab[i].a != 0xffffffffu) {
                    if (tab[i].key == key) {
                        uint32_t a = tab[i].a, c = tab[i].c; int ok = 1;
                        for (int j = 0; j < d; j++)
                            if ((int)VA[(long)a * d + j] + VB[(long)b * d + j] + VC[(long)c * d + j] + VD[(long)dd * d + j] != target[j]) { ok = 0; break; }
                        if (ok) { if (found < maxout) { out[4 * found] = a; out[4 * found + 1] = b; out[4 * found + 2] = c; out[4 * found + 3] = dd; } found++; }
                    }
                    i = (i + 1) & cm;
                }
            }
        }
        free(tab);
        if (verbose) fprintf(stderr, "pass %d/%d pairs %ld found %ld\n", r + 1, R, npairs, found);
    }
    free(hA); free(hB); free(hC); free(hD); free(kA); free(kB); free(kC); free(kD); free(cC); free(cD); free(oC); free(oD);
    return found;
}
