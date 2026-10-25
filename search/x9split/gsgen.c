// SA for a Goethals-Seidel quadruple of order n (sum of PAFs = 4n delta). usage: gsgen n seed iters
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static uint64_t rs; static inline uint64_t rnd() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int n, x[4][512]; long P[512];
int main(int argc, char **argv) {
  n = atoi(argv[1]); rs = strtoull(argv[2], 0, 10) * 2654435761ull + 1; long it = atol(argv[3]);
  for (int tries = 0; tries < 1000; tries++) {
    for (int k = 0; k < 4; k++) for (int t = 0; t < n; t++) x[k][t] = (rnd() & 1) ? 1 : -1;
    for (int d = 1; d < n; d++) { long s = 0; for (int k = 0; k < 4; k++) for (int t = 0; t < n; t++) s += x[k][t] * x[k][(t + d) % n]; P[d] = s; }
    long E = 0; for (int d = 1; d < n; d++) E += P[d] * P[d];
    for (long i = 0; i < it && E; i++) {
      double T = 8.0 * pow(0.05 / 8.0, (double)i / it);
      int k = rnd() % 4, t = rnd() % n; long dE = 0;
      static long dP[512];
      for (int d = 1; d < n; d++) { dP[d] = -2 * x[k][t] * (x[k][(t + d) % n] + x[k][(t - d + n) % n]); }
      for (int d = 1; d < n; d++) dE += (P[d] + dP[d]) * (P[d] + dP[d]) - P[d] * P[d];
      if (dE <= 0 || (rnd() >> 11) * (1.0 / 9007199254740992.0) < exp(-dE / T)) { x[k][t] = -x[k][t]; for (int d = 1; d < n; d++) P[d] += dP[d]; E += dE; }
    }
    if (E == 0) { printf("%d\n", n); for (int k = 0; k < 4; k++) { for (int t = 0; t < n; t++) printf("%d ", x[k][t]); printf("\n"); } return 0; }
  }
  return 1;
}
