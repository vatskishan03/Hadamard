// Hierarchical stage-1 engine: B-array with prescribed COLUMN sums (columns = Z_13 points, 4 Q-points each) is walked by single-column moves;
// every neighbour's class-key is looked up in a hash of A candidates (A with prescribed Z_13 compression).  Hits are Q-level solutions.
// input: V D ncol ; colvars (4 each) ; options per column ; MS_c (V x V) ; target ; nA ; A hashes (uint64) ; seconds seed
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#define MV 64
#define MD 16
static int V, D, ncol; static int colv[16][4]; static int nopt[16]; static int opt[16][600][4];
static long MS[MD][MV][MV], Md[MD][MV], tgt;
static uint64_t PR[16] = {0x9E3779B97F4A7C15ULL, 0xC2B2AE3D27D4EB4FULL, 0x165667B19E3779F9ULL, 0xD6E8FEB86659FD93ULL, 0xFF51AFD7ED558CCDULL, 0xC4CEB9FE1A85EC53ULL, 0x2545F4914F6CDD1DULL, 0x94D049BB133111EBULL, 0xBF58476D1CE4E5B9ULL, 0x9FB21C651E98DF25ULL, 0x1B873593CC9E2D51ULL, 0xE6546B64AB4F5C0FULL, 3, 5, 7, 11};
static uint64_t *tab; static uint64_t tmask; static long nA;
static uint64_t s[2];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t nxt(void) { const uint64_t s0 = s[0]; uint64_t s1 = s[1]; const uint64_t r = s0 + s1; s1 ^= s0; s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); s[1] = rotl(s1, 37); return r; }
static long x[MV], Vv[MD][MV], T[MD];
static int cur[16];
static inline int present(uint64_t h) { uint64_t p = h & tmask; while (tab[p]) { if (tab[p] == h) return 1; p = (p + 1) & tmask; } return 0; }
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); double secs = atof(argv[2]); uint64_t seed = atol(argv[3]); int maxhits = argc > 4 ? atoi(argv[4]) : 50;
  s[0] = seed * 0x9E3779B97F4A7C15ULL + 7; s[1] = seed ^ 0xD1B54A32D192ED03ULL; for (int i = 0; i < 20; i++) nxt();
  fscanf(f, "%d %d %d", &V, &D, &ncol);
  for (int t = 0; t < ncol; t++) { for (int k = 0; k < 4; k++) fscanf(f, "%d", &colv[t][k]); fscanf(f, "%d", &nopt[t]); for (int o = 0; o < nopt[t]; o++) for (int k = 0; k < 4; k++) fscanf(f, "%d", &opt[t][o][k]); }
  for (int c = 0; c < D; c++) for (int i = 0; i < V; i++) for (int j = 0; j < V; j++) fscanf(f, "%ld", &MS[c][i][j]);
  for (int c = 0; c < D; c++) for (int i = 0; i < V; i++) Md[c][i] = MS[c][i][i] / 2;
  fscanf(f, "%ld", &tgt); fscanf(f, "%ld", &nA);
  uint64_t sz = 1; while (sz < 2 * (uint64_t)nA) sz <<= 1; tmask = sz - 1; tab = calloc(sz, 8);
  for (long i = 0; i < nA; i++) { unsigned long long h; fscanf(f, "%llu", &h); if (h == 0) h = 1; uint64_t p = h & tmask; while (tab[p] && tab[p] != h) p = (p + 1) & tmask; tab[p] = h; }
  fclose(f);
  for (int t = 0; t < ncol; t++) { cur[t] = nxt() % nopt[t]; for (int k = 0; k < 4; k++) x[colv[t][k]] = opt[t][cur[t]][k]; }
  for (int c = 0; c < D; c++) { for (int i = 0; i < V; i++) { long a = 0; for (int j = 0; j < V; j++) a += MS[c][i][j] * x[j]; Vv[c][i] = a; } long tt = 0; for (int i = 0; i < V; i++) tt += x[i] * Vv[c][i]; T[c] = tt / 2; }
  long probes = 0, hits = 0, steps = 0; clock_t t0 = clock();
  while ((double)(clock() - t0) / CLOCKS_PER_SEC < secs) {
    for (int rep = 0; rep < 200; rep++) {
      for (int t = 0; t < ncol; t++) {
        int *cv = colv[t]; long xo[4]; for (int k = 0; k < 4; k++) xo[k] = x[cv[k]];
        long lin[MD];
        for (int c = 0; c < D; c++) { for (int k = 0; k < 4; k++) { lin[c] = 0; } }
        for (int o = 0; o < nopt[t]; o++) {
          long d[4]; for (int k = 0; k < 4; k++) d[k] = opt[t][o][k] - xo[k];
          uint64_t h = 0;
          for (int c = 0; c < D; c++) {
            long dT = 0;
            for (int k = 0; k < 4; k++) { dT += d[k] * Vv[c][cv[k]] + d[k] * d[k] * Md[c][cv[k]]; for (int l = k + 1; l < 4; l++) dT += d[k] * d[l] * MS[c][cv[k]][cv[l]]; }
            long need = tgt - (T[c] + dT);
            h += (uint64_t)need * PR[c];
          }
          if (h == 0) h = 1;
          probes++;
          if (present(h)) { hits++; if (hits <= maxhits) { printf("HIT hash %llu B", (unsigned long long)h); for (int tt = 0; tt < ncol; tt++) { for (int k = 0; k < 4; k++) printf(" %ld", tt == t ? (long)opt[t][o][k] : x[colv[tt][k]]); } printf("\n"); fflush(stdout); } }
        }
        // move: random option on a random column
        int tm = nxt() % ncol; int om = nxt() % nopt[tm]; int *cm = colv[tm]; long dm[4]; for (int k = 0; k < 4; k++) dm[k] = opt[tm][om][k] - x[cm[k]];
        for (int k = 0; k < 4; k++) { if (!dm[k]) continue; for (int c = 0; c < D; c++) { for (int j = 0; j < V; j++) Vv[c][j] += dm[k] * MS[c][j][cm[k]]; } }
        // update T from scratch-free formula: recompute T via x
        for (int k = 0; k < 4; k++) x[cm[k]] += dm[k];
        for (int c = 0; c < D; c++) { long tt = 0; for (int i = 0; i < V; i++) tt += x[i] * Vv[c][i]; T[c] = tt / 2; }
        steps++;
      }
    }
  }
  printf("DONE probes %ld hits %ld steps %ld\n", probes, hits, steps);
  return 0;
}
