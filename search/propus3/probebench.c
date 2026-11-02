/* Microbenchmark of the join kernel of p3.c: open-addressing table of nd 64-bit keys (empty = 0, load <= 1/4,
 * linear probing, multiplicative hash), probed by key = base - 2 h_B[j] over a contiguous B group.  Reports ns/probe
 * for table sizes that occur at v = 239 (ns = 1..3 split classes).   cc -O3 -march=native probebench.c -o probebench */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
typedef unsigned long long u64;
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static u64 s = 88172645463325252ULL;
static u64 rnd(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s * 0x2545F4914F6CDD1DULL; }
int main(int argc, char **argv) {
    u64 probes = argc > 1 ? strtoull(argv[1], 0, 10) : 40000000ULL;
    u64 sizes[] = {2000, 18000, 70000, 200000, 600000, 2400000};
    for (int z = 0; z < 6; z++) {
        u64 nd = sizes[z];
        u64 size = 64; while (size < 4 * nd) size <<= 1;
        u64 mask = size - 1; int sh = 64 - __builtin_ctzll(size);
        u64 *tk = calloc(size, 8); uint32_t *ti = malloc(4 * size);
        for (u64 i = 0; i < nd; i++) { u64 key = rnd() | 1; u64 p = (key * 0x9E3779B97F4A7C15ULL) >> sh; while (tk[p]) p = (p + 1) & mask; tk[p] = key; ti[p] = (uint32_t)i; }
        u64 nb = 1 << 16; u64 *hb = malloc(8 * nb); for (u64 i = 0; i < nb; i++) hb[i] = rnd();
        u64 found = 0, done = 0; double t0 = now();
        while (done < probes) {
            u64 base = rnd();
            for (u64 j = 0; j < nb; j++) {
                u64 key = base - 2 * hb[j];
                u64 p = (key * 0x9E3779B97F4A7C15ULL) >> sh, tk_;
                while ((tk_ = tk[p]) != 0) { if (tk_ == key) found++; p = (p + 1) & mask; }
            }
            done += nb;
        }
        double dt = now() - t0;
        printf("nd %8llu  slots %8llu (%6.1f MB keys)  %.2f ns/probe  (%.3g probes/s)  found %llu\n", nd, size, size * 8 / 1e6, 1e9 * dt / done, done / dt, found);
        free(tk); free(ti); free(hb);
    }
    return 0;
}
