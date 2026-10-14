// liftsolve.c : EXHAUSTIVE solver for the Legendre-sign lift of a Q-level Spence solution (Z_2^2 x Z_11 x Z_13, fibre Z_11).
// Problem: signs sA (nA bits), sB (nB bits), s = +-1 (bit=1 <=> s=-1).  Linear system  EA sA + EB sB = rhs  (E equations, small ints)  and
// full class conditions key_c(sA,sB) = T_c, key = KA_c(sA) + KB_c(sB), K_c(s) = c + lin.s + sum_{i<j} Q_ij s_i s_j  (all exact int64).
// Method: 4-list Schroeppel-Shamir-style join on the linear system with modular bucketing on an additive hash (exact, no sampling):
//   Left = A x B1, Right = B2 x B3; bucket r in [0,M): left pairs with hm(l)=r into a chained hash table keyed by a 64-bit additive hash,
//   right pairs with hm = hm(rhs)-r probe it; every hash match is re-verified on the exact E-vector, then the exact 27-class check runs.
// usage: liftsolve problem.txt [logM=12] [r_lo r_hi] [seed]   (fixbit in the problem file: restrict sA[0]=+1 (valid for sigma-symmetric targets))
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
typedef uint64_t u64; typedef int64_t i64; typedef uint32_t u32; typedef int8_t i8;
#define MAXE 12
#define MAXC 32
#define MAXN 64
static int nA, nB, E, NC, fixbit;
static int EAm[MAXE][MAXN], EBm[MAXE][MAXN]; static i64 rhs[MAXE], T[MAXC], cA[MAXC], cB[MAXC];
static int linA[MAXC][MAXN], linB[MAXC][MAXN];
static int (*QA)[MAXN][MAXN], (*QB)[MAXN][MAXN];
typedef struct { uint8_t i, j; int16_t w; } pr_t;
static pr_t *pairsB[MAXC]; static int npB[MAXC]; static pr_t *pairsA[MAXC]; static int npA[MAXC];
static int ord[MAXC];
static int32_t *keyA;  // [listA][NC]
static u64 cH[MAXE], dH[MAXE];
static u64 rng_s = 88172645463325252ULL;
static u64 rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }

typedef struct { int nbits, off; int n; int8_t *u; u64 *H; u32 *hm; u32 *start; u32 *ordx; u64 *sH; } list_t;
// build list over bits [off, off+nbits) of matrix Em (E x N); u = (v - p)/2 with p parity offset
static int par_total[MAXE];
static void build_list(list_t *L, int (*Em)[MAXN], int off, int nbits, int bitoff_in_mask, int M, int bucketed, int skipfirst) {
  // skipfirst: first bit of the group is pinned to +1 (fixbit); list enumerates the remaining bits
  L->nbits = nbits; L->off = off;
  int nfree = nbits - (skipfirst ? 1 : 0); int n = 1 << nfree; L->n = n;
  L->u = malloc((size_t)n * E); L->H = malloc(8 * (size_t)n); L->hm = malloc(4 * (size_t)n);
  int p[MAXE]; int base[MAXE];
  for (int e = 0; e < E; e++) { int s = 0; for (int j = 0; j < nbits; j++) s += Em[e][off + j]; p[e] = ((s % 2) + 2) % 2; par_total[e] += p[e]; }
  int sh = skipfirst ? 1 : 0;
  int8_t *v = malloc((size_t)n * E * 2 > 0 ? (size_t)n * E : 1); (void)v;
  int *vi = malloc(sizeof(int) * (size_t)n * E);
  for (int e = 0; e < E; e++) { int s = 0; for (int j = 0; j < nbits; j++) s += Em[e][off + j]; base[e] = s; vi[e] = s; }   // mask 0: all +1 (pinned bit +1 too)
  for (int m = 1; m < n; m++) {
    int lb = __builtin_ctz(m); int pm = m & (m - 1); int col = off + lb + sh;
    for (int e = 0; e < E; e++) vi[(size_t)m * E + e] = vi[(size_t)pm * E + e] - 2 * Em[e][col];
  }
  for (int m = 0; m < n; m++) {
    u64 H = 0, h = 0;
    for (int e = 0; e < E; e++) { int val = (vi[(size_t)m * E + e] - p[e]) / 2; L->u[(size_t)m * E + e] = (int8_t)val; H += cH[e] * (u64)(i64)val; h += dH[e] * (u64)(i64)val; }
    L->H[m] = H; L->hm[m] = (u32)(h & (u64)(M - 1));
  }
  free(vi); free(v);
  if (bucketed) {
    L->start = calloc(M + 1, 4); L->ordx = malloc(4 * (size_t)n); L->sH = malloc(8 * (size_t)n);
    for (int m = 0; m < n; m++) L->start[L->hm[m] + 1]++;
    for (int b = 0; b < M; b++) L->start[b + 1] += L->start[b];
    u32 *fill = malloc(4 * (size_t)M); memcpy(fill, L->start, 4 * (size_t)M);
    for (int m = 0; m < n; m++) { u32 pos = fill[L->hm[m]]++; L->ordx[pos] = m; L->sH[pos] = L->H[m]; }
    free(fill);
  }
}
static i64 evalB(int c, const int *s) {
  i64 v = cB[c]; for (int i = 0; i < nB; i++) v += (i64)linB[c][i] * s[i];
  pr_t *P = pairsB[c]; for (int k = 0; k < npB[c]; k++) v += (i64)P[k].w * s[P[k].i] * s[P[k].j];
  return v;
}
static i64 evalA(int c, const int *s) {
  i64 v = cA[c]; for (int i = 0; i < nA; i++) v += (i64)linA[c][i] * s[i];
  pr_t *P = pairsA[c]; for (int k = 0; k < npA[c]; k++) v += (i64)P[k].w * s[P[k].i] * s[P[k].j];
  return v;
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); if (!f) { fprintf(stderr, "no input\n"); return 1; }
  int logM = argc > 2 ? atoi(argv[2]) : 12; int M = 1 << logM;
  int rlo = argc > 3 ? atoi(argv[3]) : 0, rhi = argc > 4 ? atoi(argv[4]) : M; u64 seed = argc > 5 ? strtoull(argv[5], 0, 10) : 12345;
  rng_s ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) rnd();
  long long tmp;
  fscanf(f, "%d %d %d %d %d", &nA, &nB, &E, &NC, &fixbit);
  for (int e = 0; e < E; e++) for (int j = 0; j < nA; j++) fscanf(f, "%d", &EAm[e][j]);
  for (int e = 0; e < E; e++) for (int j = 0; j < nB; j++) fscanf(f, "%d", &EBm[e][j]);
  for (int e = 0; e < E; e++) { fscanf(f, "%lld", &tmp); rhs[e] = tmp; }
  for (int c = 0; c < NC; c++) { fscanf(f, "%lld", &tmp); T[c] = tmp; }
  for (int c = 0; c < NC; c++) { fscanf(f, "%lld", &tmp); cA[c] = tmp; }
  for (int c = 0; c < NC; c++) { fscanf(f, "%lld", &tmp); cB[c] = tmp; }
  for (int c = 0; c < NC; c++) for (int j = 0; j < nA; j++) fscanf(f, "%d", &linA[c][j]);
  for (int c = 0; c < NC; c++) for (int j = 0; j < nB; j++) fscanf(f, "%d", &linB[c][j]);
  QA = calloc((size_t)NC, sizeof(int[MAXN][MAXN])); QB = calloc((size_t)NC, sizeof(int[MAXN][MAXN]));
  for (int c = 0; c < NC; c++) for (int i = 0; i < nA; i++) for (int j = 0; j < nA; j++) fscanf(f, "%d", &QA[c][i][j]);
  for (int c = 0; c < NC; c++) for (int i = 0; i < nB; i++) for (int j = 0; j < nB; j++) fscanf(f, "%d", &QB[c][i][j]);
  fclose(f);
  for (int c = 0; c < NC; c++) {
    pairsA[c] = malloc(sizeof(pr_t) * (size_t)(nA * nA + 1)); pairsB[c] = malloc(sizeof(pr_t) * (size_t)(nB * nB + 1)); npA[c] = npB[c] = 0;
    for (int i = 0; i < nA; i++) for (int j = i + 1; j < nA; j++) if (QA[c][i][j]) { pairsA[c][npA[c]++] = (pr_t){(uint8_t)i, (uint8_t)j, (int16_t)QA[c][i][j]}; }
    for (int i = 0; i < nB; i++) for (int j = i + 1; j < nB; j++) if (QB[c][i][j]) { pairsB[c][npB[c]++] = (pr_t){(uint8_t)i, (uint8_t)j, (int16_t)QB[c][i][j]}; }
  }
  // class check order: decreasing estimated variance of the B part
  { double var[MAXC]; for (int c = 0; c < NC; c++) { double v = 0; for (int i = 0; i < nB; i++) v += (double)linB[c][i] * linB[c][i]; for (int k = 0; k < npB[c]; k++) v += (double)pairsB[c][k].w * pairsB[c][k].w; var[c] = v; ord[c] = c; }
    for (int a = 0; a < NC; a++) for (int b = a + 1; b < NC; b++) if (var[ord[b]] > var[ord[a]]) { int t = ord[a]; ord[a] = ord[b]; ord[b] = t; } }
  for (int e = 0; e < E; e++) { cH[e] = rnd() | 1; dH[e] = rnd() | 1; }
  memset(par_total, 0, sizeof par_total);
  // split bits
  int nAf = nA - (fixbit ? 1 : 0); int ntot = nAf + nB; int nR = ntot / 2; int g2 = nR / 2, g3 = nR - g2, g1 = nB - nR;
  if (g1 < 0) { fprintf(stderr, "bad split\n"); return 1; }
  list_t LA, L1, L2, L3;
  build_list(&LA, EAm, 0, nA, 0, M, 0, fixbit);
  build_list(&L1, EBm, 0, g1, 0, M, 1, 0);
  build_list(&L2, EBm, g1, g2, 0, M, 0, 0);
  build_list(&L3, EBm, g1 + g2, g3, 0, M, 1, 0);
  // rhs' = (rhs - sum parities)/2
  i64 rp[MAXE]; for (int e = 0; e < E; e++) { i64 t = rhs[e] - par_total[e]; if (t & 1) { printf("PARITY_IMPOSSIBLE\nSTATS lin 0 full 0 leftpairs 0 probes 0\n"); return 0; } rp[e] = t / 2; }
  u64 HR = 0, hR = 0; for (int e = 0; e < E; e++) { HR += cH[e] * (u64)rp[e]; hR += dH[e] * (u64)rp[e]; } hR &= (u64)(M - 1);
  // A tables of exact keys
  int nAl = LA.n; keyA = malloc(4 * (size_t)nAl * NC);
  for (int m = 0; m < nAl; m++) { int s[MAXN]; int mm = fixbit ? (m << 1) : m; for (int i = 0; i < nA; i++) s[i] = 1 - 2 * ((mm >> i) & 1);
    for (int c = 0; c < NC; c++) keyA[(size_t)m * NC + c] = (int32_t)evalA(c, s); }
  fprintf(stderr, "lists: A %d (nA=%d fix=%d) B1 %d (g=%d) B2 %d (g=%d) B3 %d (g=%d) M=%d E=%d\n", LA.n, nA, fixbit, L1.n, g1, L2.n, g2, L3.n, g3, M, E);
  // hash table
  double avg = (double)LA.n * L1.n / M; int hb = 4; while ((1 << hb) < 2 * avg) hb++; if (hb > 24) hb = 24; size_t HS = (size_t)1 << hb; int *head = malloc(4 * HS); for (size_t i = 0; i < HS; i++) head[i] = -1;
  size_t cap = (size_t)(avg * 2) + 1024; u64 *nkey = malloc(8 * cap); u32 *nia = malloc(4 * cap), *nib = malloc(4 * cap); int *nnext = malloc(4 * cap);
  u64 nlin = 0, nfull = 0, nleft = 0, nprobe = 0, nhash = 0; clock_t t0 = clock();
  for (int r = rlo; r < rhi; r++) {
    size_t nent = 0;
    for (int ia = 0; ia < LA.n; ia++) {
      u32 bk = ((u32)r - LA.hm[ia]) & (M - 1); u64 hA = LA.H[ia];
      for (u32 p = L1.start[bk]; p < L1.start[bk + 1]; p++) {
        if (nent == cap) { cap *= 2; nkey = realloc(nkey, 8 * cap); nia = realloc(nia, 4 * cap); nib = realloc(nib, 4 * cap); nnext = realloc(nnext, 4 * cap); }
        u64 key = hA + L1.sH[p]; size_t ix = (key * 0x9E3779B97F4A7C15ULL) >> (64 - hb);
        nkey[nent] = key; nia[nent] = ia; nib[nent] = L1.ordx[p]; nnext[nent] = head[ix]; head[ix] = (int)nent; nent++;
      }
    }
    nleft += nent;
    for (int i2 = 0; i2 < L2.n; i2++) {
      u32 bk = (((u32)hR - (u32)r) - L2.hm[i2]) & (M - 1); u64 hh = HR - L2.H[i2];
      for (u32 p = L3.start[bk]; p < L3.start[bk + 1]; p++) {
        u64 key = hh - L3.sH[p]; size_t ix = (key * 0x9E3779B97F4A7C15ULL) >> (64 - hb); nprobe++;
        for (int nd = head[ix]; nd >= 0; nd = nnext[nd]) if (nkey[nd] == key) {
          nhash++;
          u32 ia = nia[nd], i1 = nib[nd], i3 = L3.ordx[p];
          int ok = 1; for (int e = 0; e < E; e++) { int sum = LA.u[(size_t)ia * E + e] + L1.u[(size_t)i1 * E + e] + L2.u[(size_t)i2 * E + e] + L3.u[(size_t)i3 * E + e]; if (sum != rp[e]) { ok = 0; break; } }
          if (!ok) continue;
          nlin++;
          u64 bmask = (u64)i1 | ((u64)i2 << g1) | ((u64)i3 << (g1 + g2));
          int sB[MAXN]; for (int i = 0; i < nB; i++) sB[i] = 1 - 2 * (int)((bmask >> i) & 1);
          int good = 1; for (int k = 0; k < NC; k++) { int c = ord[k]; i64 v = evalB(c, sB) + keyA[(size_t)ia * NC + c]; if (v != T[c]) { good = 0; break; } }
          if (good) { nfull++; int mm = fixbit ? (int)(ia << 1) : (int)ia; printf("SOL %d %llu\n", mm, (unsigned long long)bmask); fflush(stdout); }
        }
      }
    }
    for (size_t k = 0; k < nent; k++) head[(nkey[k] * 0x9E3779B97F4A7C15ULL) >> (64 - hb)] = -1;
    if (((r - rlo) & 255) == 255) fprintf(stderr, "r=%d/%d lin=%llu full=%llu  %.1fs\n", r + 1, rhi, (unsigned long long)nlin, (unsigned long long)nfull, (double)(clock() - t0) / CLOCKS_PER_SEC);
  }
  printf("STATS lin %llu full %llu leftpairs %llu probes %llu hashmatches %llu fixbit %d cpu %.1f\n", (unsigned long long)nlin, (unsigned long long)nfull, (unsigned long long)nleft, (unsigned long long)nprobe, (unsigned long long)nhash, fixbit, (double)(clock() - t0) / CLOCKS_PER_SEC);
  return 0;
}
