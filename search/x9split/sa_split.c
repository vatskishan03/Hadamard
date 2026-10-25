// Simulated annealing for a circulant 4-split of a Goethals-Seidel quadruple (x9 plug-in).
// Input file: n, then 4 lines of n entries (+-1), then 4 lines of n slot ids (slots partition positions per sequence;
// slot ids are global), then the list of fixed (slot, layer) pairs: "nfix" followed by pairs.
// Layers 0..3 each need sum_k PAF(layer_k) = (total weight) delta, total weight per layer = n (enforced by swaps
// that preserve the per-layer weight, given an initial balanced assignment).
// usage: sa_split infile seed maxiters T0 T1 restarts [shiftmode: 0 = all d in 1..(n-1)/2, else list of reps file]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
static uint64_t rs;
static inline uint64_t rnd() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand() { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
int n, x[4][512], slotof[4][512], nslot, slotk[4096], slotsz[4096], slotpos[4096][16], fixedl[4096];
int lab[4096]; int D, dl[512]; long P[4][512];
static inline int L(int k, int t) { return lab[slotof[k][t]]; }
// apply: move slot s to layer b, updating P
static void apply(int s, int b) {
  int a = lab[s], k = slotk[s];
  if (a == b) return;
  for (int di = 0; di < D; di++) {
    int d = dl[di]; long dA = 0, dB = 0, in = 0;
    for (int q = 0; q < slotsz[s]; q++) {
      int t = slotpos[s][q], tp = (t + d) % n, tm = (t - d + n) % n;
      int lp = L(k, tp), lm = L(k, tm);
      long xp = x[k][t] * x[k][tp], xm = x[k][t] * x[k][tm];
      if (lp == a) dA += xp; if (lm == a) dA += xm;
      if (lp == b) dB += xp; if (lm == b) dB += xm;
      if (slotof[k][tp] == s) in += xp;
    }
    // pairs inside the slot: counted twice in dA (both currently a), want once; not counted in dB
    P[a][di] += -dA + in;
    P[b][di] += dB + in;
  }
  lab[s] = b;
}
static long energy() { long e = 0; for (int a = 0; a < 4; a++) for (int di = 0; di < D; di++) e += P[a][di] * P[a][di]; return e; }
static void recompute() {
  for (int a = 0; a < 4; a++) for (int di = 0; di < D; di++) {
    long s = 0; int d = dl[di];
    for (int k = 0; k < 4; k++) for (int t = 0; t < n; t++) if (L(k, t) == a && L(k, (t + d) % n) == a) s += x[k][t] * x[k][(t + d) % n];
    P[a][di] = s;
  }
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); uint64_t seed = strtoull(argv[2], 0, 10); long maxit = atol(argv[3]);
  double T0 = atof(argv[4]), T1 = atof(argv[5]); int restarts = atoi(argv[6]);
  rs = seed * 0x9E3779B97F4A7C15ull + 12345; for (int i = 0; i < 10; i++) rnd();
  if (fscanf(f, "%d", &n) != 1) return 1;
  for (int k = 0; k < 4; k++) for (int t = 0; t < n; t++) if (fscanf(f, "%d", &x[k][t]) != 1) return 1;
  nslot = 0;
  for (int k = 0; k < 4; k++) for (int t = 0; t < n; t++) { int s; if (fscanf(f, "%d", &s) != 1) return 1; slotof[k][t] = s; if (s + 1 > nslot) nslot = s + 1; }
  for (int s = 0; s < nslot; s++) { slotsz[s] = 0; fixedl[s] = -1; }
  for (int k = 0; k < 4; k++) for (int t = 0; t < n; t++) { int s = slotof[k][t]; slotk[s] = k; slotpos[s][slotsz[s]++] = t; }
  int nfix; if (fscanf(f, "%d", &nfix) != 1) return 1;
  for (int i = 0; i < nfix; i++) { int s, l; if (fscanf(f, "%d %d", &s, &l) != 2) return 1; fixedl[s] = l; }
  if (fscanf(f, "%d", &D) != 1) return 1;
  for (int i = 0; i < D; i++) if (fscanf(f, "%d", &dl[i]) != 1) return 1;
  fclose(f);
  // free slots grouped by size; balanced initial assignment: each layer gets weight n
  int freeS[4096], nf = 0; for (int s = 0; s < nslot; s++) if (fixedl[s] < 0) freeS[nf++] = s;
  long best_overall = -1;
  for (int r = 0; r < restarts; r++) {
    // random balanced assignment: fixed first, then fill layers greedily in random order by size
    long wt[4] = {0, 0, 0, 0};
    for (int s = 0; s < nslot; s++) if (fixedl[s] >= 0) { lab[s] = fixedl[s]; wt[fixedl[s]] += slotsz[s]; }
    for (int i = nf - 1; i > 0; i--) { int j = rnd() % (i + 1); int tmp = freeS[i]; freeS[i] = freeS[j]; freeS[j] = tmp; }
    for (int i = 0; i < nf; i++) { int s = freeS[i]; int best = -1;
      int st = (int)(rnd() % 4); for (int a0 = 0; a0 < 4; a0++) { int a = (a0 + st) % 4; if (wt[a] + slotsz[s] <= n) { best = a; break; } }
      if (best < 0) { fprintf(stderr, "cannot balance\n"); return 1; }
      lab[s] = best; wt[best] += slotsz[s]; }
    for (int a = 0; a < 4; a++) if (wt[a] != n) { fprintf(stderr, "unbalanced %ld\n", wt[a]); return 1; }
    recompute(); long E = energy(), bestE = E;
    for (long it = 0; it < maxit; it++) {
      double T = T0 * pow(T1 / T0, (double)it / maxit);
      int s1 = freeS[rnd() % nf], s2 = freeS[rnd() % nf];
      if (lab[s1] == lab[s2] || slotsz[s1] != slotsz[s2]) continue;
      int a = lab[s1], b = lab[s2];
      apply(s1, b); apply(s2, a);
      long E2 = energy();
      if (E2 <= E || urand() < exp((E - E2) / T)) { E = E2; if (E < bestE) bestE = E;
        if (E == 0) break; }
      else { apply(s2, b); apply(s1, a); }
    }
    printf("restart %d seed %llu final E %ld best %ld\n", r, (unsigned long long)seed, E, bestE); fflush(stdout);
    if (E == 0) {
      printf("SOLUTION");
      for (int s = 0; s < nslot; s++) printf(" %d", lab[s]);
      printf("\n"); fflush(stdout);
      return 0;
    }
  }
  return 0;
}
