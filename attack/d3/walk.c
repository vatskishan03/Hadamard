/* d=3 walk search: table of all X0 class vectors (exact hash), tabu walks over X1 (H-invariant, row sum fixed)
 * evaluating every row-sum preserving swap; each neighbour with centrality score <= SLOOK is looked up.
 * usage: ./walk seconds threads slook seed     ; prints HIT lines "HIT idx0 <86 signs of X1>" */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define CK 15
static int N0, k1; static int16_t *V0; static int *sizes; static int *MF; static int DW[CK]; static double mu0[CK], sd0[CK];
static uint64_t R[CK]; static uint64_t *tkey; static uint32_t *tidx; static uint64_t tmask;
static inline uint64_t hkey(const int *v){ uint64_t h = 0; for (int k = 0; k < CK; k++) h += R[k] * (uint64_t)(int64_t)v[k]; return h; }
static inline uint64_t mix(uint64_t x){ x ^= x >> 31; x *= 0x7fb5d329728ea185ULL; x ^= x >> 27; x *= 0x81dadef4bc2dd44dULL; x ^= x >> 33; return x; }
static long lookup(const int *v){             /* returns pool index with V0 == v exactly, or -1 */
    uint64_t h = hkey(v), s = mix(h) & tmask;
    while (tkey[s]){ if (tkey[s] == (h | 1)){ uint32_t i = tidx[s]; int ok = 1; for (int k = 0; k < CK && ok; k++) ok = V0[(long)i*CK+k] == v[k]; if (ok) return i; } s = (s + 1) & tmask; }
    return -1;
}
static inline uint64_t xr(uint64_t *s){ *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17; return *s; }
int main(int argc, char **argv){
    double secs = atof(argv[1]); int nth = atoi(argv[2]); double slook = atof(argv[3]); uint64_t seed = atoll(argv[4]);
    FILE *f = fopen("pool0.bin", "rb"); if (fread(&N0, 4, 1, f) != 1) return 2;
    V0 = malloc(2L * N0 * CK); if (fread(V0, 2, (size_t)N0 * CK, f) != (size_t)N0 * CK) return 2; fclose(f);
    f = fopen("walk.bin", "rb"); int ck; if (fread(&k1, 4, 1, f) != 1 || fread(&ck, 4, 1, f) != 1 || ck != CK) return 2;
    sizes = malloc(4 * k1); MF = malloc(4L * CK * k1 * k1);
    if (fread(sizes, 4, k1, f) != (size_t)k1 || fread(MF, 4, (size_t)CK*k1*k1, f) != (size_t)CK*k1*k1 || fread(DW, 4, CK, f) != CK
        || fread(mu0, 8, CK, f) != CK || fread(sd0, 8, CK, f) != CK) return 2; fclose(f);
    uint64_t s0 = 0x9E3779B97F4A7C15ULL; for (int k = 0; k < CK; k++){ s0 = mix(s0 + k); R[k] = s0 | 1; }
    long tsz = 1; while (tsz < 3L * N0) tsz <<= 1; tmask = tsz - 1;
    tkey = calloc(tsz, 8); tidx = malloc(4 * tsz);
    for (long i = 0; i < N0; i++){ int v[CK]; for (int k = 0; k < CK; k++) v[k] = V0[i*CK+k];
        uint64_t h = hkey(v), s = mix(h) & tmask; while (tkey[s]) s = (s + 1) & tmask; tkey[s] = h | 1; tidx[s] = i; }
    fprintf(stderr, "table built: %d entries, %ld slots\n", N0, tsz);
    double w[CK]; for (int k = 0; k < CK; k++) w[k] = 1.0 / (sd0[k] * sd0[k]);
    long tot_eval = 0, tot_look = 0, tot_steps = 0; double t0 = omp_get_wtime();
    #pragma omp parallel num_threads(nth) reduction(+:tot_eval, tot_look, tot_steps)
    {
        uint64_t rs = seed * 1000003ULL + 77 + omp_get_thread_num() * 0x9E3779B97F4A7C15ULL; xr(&rs);
        int *y = malloc(4 * k1); int *G = malloc(4L * CK * k1); long *tabu = calloc(k1, 8);
        int W[CK]; int idx15[64], idx30[64], n15 = 0, n30 = 0;
        for (int o = 0; o < k1; o++){ if (sizes[o] == 15) idx15[n15++] = o; else if (sizes[o] == 30) idx30[n30++] = o; }
        long step = 0;
        while (omp_get_wtime() - t0 < secs){
            /* random start: y(size1) = -1, y(size2) = +1, half of each size class positive */
            for (int o = 0; o < k1; o++) y[o] = sizes[o] == 1 ? -1 : 1;
            int *grp[2] = {idx15, idx30}; int gn[2] = {n15, n30};
            for (int g = 0; g < 2; g++){ int m = gn[g]; int perm[64]; for (int i = 0; i < m; i++) perm[i] = grp[g][i];
                for (int i = m - 1; i > 0; i--){ int j = xr(&rs) % (i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
                for (int i = 0; i < m; i++) y[perm[i]] = i < m / 2 ? 1 : -1; }
            for (int k = 0; k < CK; k++){ long q = 0; for (int P = 0; P < k1; P++){ long gg = 0; const int *row = MF + ((long)k*k1 + P)*k1;
                    for (int O = 0; O < k1; O++) gg += row[O] * y[O]; G[k*k1+P] = (int)gg; q += gg * y[P]; } W[k] = DW[k] + (int)(q / 2); }
            memset(tabu, 0, 8 * k1); double best = 1e300; long stall = 0;
            while (stall < 4000 && omp_get_wtime() - t0 < secs){
                double bsc = 1e300; int b1 = -1, b2 = -1, bcnt = 0;
                for (int g = 0; g < 2; g++) for (int i = 0; i < gn[g]; i++){ int o1 = grp[g][i]; if (y[o1] < 0) continue;
                    for (int j = 0; j < gn[g]; j++){ int o2 = grp[g][j]; if (y[o2] > 0) continue;
                        int Wn[CK]; double sc = 0; int y1 = y[o1], y2 = y[o2];
                        for (int k = 0; k < CK; k++){ int d = -2 * (y1 * G[k*k1+o1] + y2 * G[k*k1+o2]) + 4 * y1 * y2 * MF[((long)k*k1 + o1)*k1 + o2];
                            Wn[k] = W[k] + d; double z = Wn[k] + mu0[k]; sc += z * z * w[k]; }
                        tot_eval++;
                        if (sc <= slook){ int tv[CK]; for (int k = 0; k < CK; k++) tv[k] = -Wn[k]; tot_look++;
                            long hit = lookup(tv);
                            if (hit >= 0){
                                #pragma omp critical
                                { printf("HIT %ld", hit); for (int o = 0; o < k1; o++){ int v = o == o1 || o == o2 ? -y[o] : y[o]; printf(" %d", v); } printf("\n"); fflush(stdout); }
                            } }
                        int istabu = tabu[o1] > step || tabu[o2] > step;
                        if (istabu && sc >= best) continue;
                        if (sc < bsc){ bsc = sc; b1 = o1; b2 = o2; bcnt = 1; } else if (sc == bsc && (xr(&rs) % ++bcnt) == 0){ b1 = o1; b2 = o2; }
                    } }
                if (b1 < 0) break;
                int fl[2] = {b1, b2};
                for (int t = 0; t < 2; t++){ int o = fl[t]; int yo = y[o];
                    for (int k = 0; k < CK; k++){ int *Gk = G + k*k1; const int *col = MF + (long)k*k1*k1;
                        W[k] += -2 * yo * Gk[o] + 0;   /* sequential single flips: dQ = -2 y_o g_o */
                        for (int P = 0; P < k1; P++) Gk[P] += -2 * yo * col[(long)P*k1 + o]; }
                    y[o] = -yo; }
                step++; tot_steps++;
                tabu[b1] = step + 5 + xr(&rs) % 10; tabu[b2] = step + 5 + xr(&rs) % 10;
                if (bsc < best){ best = bsc; stall = 0; } else stall++;
            }
        }
        free(y); free(G); free(tabu);
    }
    double dt = omp_get_wtime() - t0;
    fprintf(stderr, "steps %ld  evaluated %ld (%.3g/s)  lookups %ld (%.3g/s)  %.0fs\n", tot_steps, tot_eval, tot_eval / dt, tot_look, tot_look / dt, dt);
    return 0;
}
