/*
 * Four-list variant of search/match.c (extended k-tree matcher) for quadruples whose
 * four roles all have different candidate lists, e.g. skew Goethals-Seidel
 * (skew-type A, and B, C, D with distinct row sums).
 *
 * Input: four candidate files (records: uint64 mask, then C int8 class values v).
 * We look for a in A, b in B, c in C, d in D with
 *
 *     v_a + v_b + v_c + v_d = TARGET   on every class,
 *
 * where TARGET = -(n mod 4). Pairs (a, d) form the left side and (b, c) the right.
 *
 * Algorithm (Wagner; Minder-Sinclair): fix the values tau of v_a + v_d on the first
 * S classes. Pairs (a, d) with that partial sum, and pairs (b, c) with partial sum
 * TARGET - tau, are enumerated by merging lists sorted on an additive 64-bit hash
 * of the first S classes. Both pair lists are sorted on a hash of the remaining
 * classes and joined, and every hit is checked exactly. Values of tau are taken
 * in decreasing order of probability (estimated from the class histograms).
 *
 *     cc -O3 -march=native -fopenmp attack/match4.c -o match4 -lm
 *     ./match4 -s 6 -t -1 A.bin B.bin C.bin D.bin 12 > solutions.txt
 *
 * Output lines: "SOL maskA maskB maskC maskD" (hexadecimal). Progress on stderr.
 */
#include <fcntl.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAXC 128
#define TOPMAX 16

static int C, S, TARGET, RECLEN;
static uint64_t R1[MAXC], R2[MAXC];

typedef struct {
    long n;
    const unsigned char *rec;
} List;

typedef struct {   /* list sorted on the clamp hash */
    long n;
    uint64_t *key, *h2;
    uint32_t *idx;
} Index;

typedef struct {
    uint64_t h;
    uint32_t i, j;
} Pair;

static uint64_t mix(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    return x ^ (x >> 33);
}

static const int8_t *values(const List *L, long i) { return (const int8_t *)(L->rec + i * RECLEN + 8); }

static uint64_t mask(const List *L, long i)
{
    uint64_t m;
    memcpy(&m, L->rec + i * RECLEN, 8);
    return m;
}

/* LSD radix sort on Pair.h, after a parallel partition on the top 12 bits. */
static void sort_pairs(Pair *a, Pair *tmp, long n)
{
    const int NB = 4096;
    int nt = omp_get_max_threads();
    long *tc = calloc((size_t)nt * NB, sizeof(long)), *start = calloc(NB + 1, sizeof(long));
#pragma omp parallel
    {
        int t = omp_get_thread_num();
        long lo = n * t / nt, hi = n * (t + 1) / nt, *c = tc + (size_t)t * NB;
        for (long i = lo; i < hi; i++)
            c[a[i].h >> 52]++;
    }
    long s = 0;
    for (int b = 0; b < NB; b++) {
        start[b] = s;
        for (int t = 0; t < nt; t++) {
            long v = tc[(size_t)t * NB + b];
            tc[(size_t)t * NB + b] = s;
            s += v;
        }
    }
    start[NB] = n;
#pragma omp parallel
    {
        int t = omp_get_thread_num();
        long lo = n * t / nt, hi = n * (t + 1) / nt, *c = tc + (size_t)t * NB;
        for (long i = lo; i < hi; i++)
            tmp[c[a[i].h >> 52]++] = a[i];
    }
#pragma omp parallel for schedule(dynamic, 8)
    for (int b = 0; b < NB; b++) {
        long lo = start[b], m = start[b + 1] - lo;
        if (m <= 1) {
            if (m == 1)
                a[lo] = tmp[lo];
            continue;
        }
        Pair *src = tmp + lo, *dst = a + lo;
        long *cnt = malloc(8193 * sizeof(long));
        for (int sh = 0; sh < 52; sh += 13) {
            memset(cnt, 0, 8193 * sizeof(long));
            for (long i = 0; i < m; i++)
                cnt[((src[i].h >> sh) & 8191) + 1]++;
            for (int k = 0; k < 8192; k++)
                cnt[k + 1] += cnt[k];
            for (long i = 0; i < m; i++)
                dst[cnt[(src[i].h >> sh) & 8191]++] = src[i];
            Pair *t = src;
            src = dst;
            dst = t;
        }
        if (src != a + lo)
            memcpy(a + lo, src, m * sizeof(Pair));
        free(cnt);
    }
    free(tc);
    free(start);
}

static void load(const char *fn, List *L)
{
    int fd = open(fn, O_RDONLY);
    struct stat st;
    if (fd < 0 || fstat(fd, &st)) {
        perror(fn);
        exit(1);
    }
    L->rec = mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    madvise((void *)L->rec, st.st_size, MADV_WILLNEED);
    L->n = st.st_size / RECLEN;
}

static void build_index(const List *L, Index *I)
{
    long n = L->n;
    Pair *a = malloc(n * sizeof(Pair)), *t = malloc(n * sizeof(Pair));
    uint64_t *h2 = malloc(n * 8);
#pragma omp parallel for schedule(static)
    for (long i = 0; i < n; i++) {
        const int8_t *v = values(L, i);
        uint64_t x = 0, y = 0;
        for (int c = 0; c < S; c++)
            x += R1[c] * (uint64_t)(int64_t)v[c];
        for (int c = S; c < C; c++)
            y += R2[c] * (uint64_t)(int64_t)v[c];
        a[i] = (Pair){x, (uint32_t)i, 0};
        h2[i] = y;
    }
    sort_pairs(a, t, n);
    free(t);
    I->n = n;
    I->key = malloc(n * 8);
    I->h2 = malloc(n * 8);
    I->idx = malloc(n * 4);
#pragma omp parallel for
    for (long i = 0; i < n; i++) {
        I->key[i] = a[i].h;
        I->idx[i] = a[i].i;
        I->h2[i] = h2[a[i].i];
    }
    free(a);
    free(h2);
}

/*
 * Pairs (x, y) with key(x) + key(y) = H (mod 2^64), by a two-pointer merge.
 * write = 0 counts; write = 1 stores Pair{h2(x) + h2(y)} (left) or
 * Pair{H2 - h2(x) - h2(y)} (right, H2 != 0).
 */
static long merge(const Index *X, const Index *Y, uint64_t H, int write, uint64_t H2, Pair *out, int nth)
{
    long *cnt = calloc(nth + 1, sizeof(long));
    for (int pass = write ? 0 : 1; pass < 2; pass++) {
#pragma omp parallel num_threads(nth)
        {
            int t = omp_get_thread_num();
            long lo = X->n * t / nth, hi = X->n * (t + 1) / nth, c = 0, w = (pass && write) ? cnt[t] : 0;
            long i = lo;
            while (i < hi) {
                /* keys at or below H need targets in [0, H]; the rest wrap around */
                int side = X->key[i] <= H;
                uint64_t tgt = H - X->key[i];
                long a = 0, b = Y->n;
                while (a < b) {
                    long mid = (a + b) >> 1;
                    if (Y->key[mid] < tgt)
                        a = mid + 1;
                    else
                        b = mid;
                }
                long j = a;
                while (i < hi && (X->key[i] <= H) == side) {
                    uint64_t kx = X->key[i];
                    tgt = H - kx;
                    while (j > 0 && Y->key[j - 1] >= tgt)
                        j--;
                    while (j < Y->n && Y->key[j] < tgt)
                        j++;
                    long j2 = j;
                    while (j2 < Y->n && Y->key[j2] == tgt)
                        j2++;
                    long i2 = i;
                    while (i2 < hi && X->key[i2] == kx)
                        i2++;
                    if (j2 > j) {
                        if (pass == 0 || !write)
                            c += (j2 - j) * (i2 - i);
                        else
                            for (long u = i; u < i2; u++)
                                for (long v = j; v < j2; v++) {
                                    uint64_t s = X->h2[u] + Y->h2[v];
                                    out[w++] = (Pair){H2 ? H2 - s : s, X->idx[u], Y->idx[v]};
                                }
                    }
                    i = i2;
                }
            }
            if (pass == 0 || !write)
                cnt[t + 1] = c;
        }
        if (pass == 0)
            for (int u = 0; u < nth; u++)
                cnt[u + 1] += cnt[u];
        if (!write)
            break;
    }
    long total = 0;
    if (write)
        total = cnt[nth];
    else
        for (int u = 1; u <= nth; u++)
            total += cnt[u];
    free(cnt);
    return total;
}

static double *g_weight;
static int by_weight(const void *a, const void *b)
{
    double x = g_weight[*(const long *)a], y = g_weight[*(const long *)b];
    return x < y ? 1 : x > y ? -1 : 0;
}

static void histogram(const List *L, double h[][256])
{
    long step = L->n > 20000000 ? L->n / 20000000 : 1, m = 0;
    for (long i = 0; i < L->n; i += step, m++) {
        const int8_t *v = values(L, i);
        for (int c = 0; c < C; c++)
            h[c][v[c] + 128] += 1;
    }
    for (int c = 0; c < C; c++)
        for (int x = 0; x < 256; x++)
            h[c][x] /= m;
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "usage: %s [-s clamp] [-t target] [-n max_tau] [-f first_tau] [-j threads]\n"
            "          [-k max_solutions] [-m values_per_class] A.bin B.bin C.bin D.bin CLASSES\n",
            prog);
    exit(2);
}

int main(int argc, char **argv)
{
    long max_tau = 1000000, first_tau = 0, max_sol = 16;
    int nth = omp_get_max_threads(), top = 6, opt;
    S = 8;
    TARGET = -3;
    while ((opt = getopt(argc, argv, "s:t:n:f:j:k:m:")) != -1) {
        switch (opt) {
        case 's': S = atoi(optarg); break;
        case 't': TARGET = atoi(optarg); break;
        case 'n': max_tau = atol(optarg); break;
        case 'f': first_tau = atol(optarg); break;
        case 'j': nth = atoi(optarg); break;
        case 'k': max_sol = atol(optarg); break;
        case 'm': top = atoi(optarg); break;
        default: usage(argv[0]);
        }
    }
    if (argc - optind != 5)
        usage(argv[0]);
    C = atoi(argv[optind + 4]);
    if (C > MAXC || S < 1 || S >= C || top < 1 || top > TOPMAX)
        usage(argv[0]);
    RECLEN = 8 + C;
    omp_set_num_threads(nth);
    for (int c = 0; c < C; c++) {
        R1[c] = mix(0x1234567ULL + c * 7919ULL) | 1;
        R2[c] = mix(0xabcdef1ULL + c * 104729ULL) | 1;
    }

    double t0 = omp_get_wtime();
    const char *fn[4] = {argv[optind], argv[optind + 1], argv[optind + 2], argv[optind + 3]};
    List Ls[4];
    Index Is[4];
    for (int k = 0; k < 4; k++) {
        int same = -1;
        for (int q = 0; q < k; q++)
            if (!strcmp(fn[q], fn[k]))
                same = q;
        if (same >= 0)
            Ls[k] = Ls[same], Is[k] = Is[same];
        else {
            load(fn[k], &Ls[k]);
            build_index(&Ls[k], &Is[k]);
        }
    }
    List A = Ls[0], B = Ls[1], Cl = Ls[2], D = Ls[3];
    Index IA = Is[0], IB = Is[1], IC = Is[2], ID = Is[3];
    fprintf(stderr, "lists A=%ld B=%ld C=%ld D=%ld indexed in %.1fs\n", A.n, B.n, Cl.n, D.n, omp_get_wtime() - t0);

    /* candidate partial sums per clamped class, ranked by P(left = t) P(right = TARGET - t) */
    static double hA[MAXC][256], hB[MAXC][256], hC[MAXC][256], hD[MAXC][256];
    histogram(&A, hA);
    histogram(&B, hB);
    histogram(&Cl, hC);
    histogram(&D, hD);
    int cand[MAXC][TOPMAX];
    double cw[MAXC][TOPMAX], pl[MAXC][TOPMAX], pr[MAXC][TOPMAX], logp_rest = 0;
    for (int c = 0; c < C; c++) {
        static double L2[512], R2s[512];
        memset(L2, 0, sizeof L2);
        memset(R2s, 0, sizeof R2s);
        for (int x = 0; x < 256; x++)
            for (int y = 0; y < 256; y++) {
                L2[x + y] += hA[c][x] * hD[c][y];
                R2s[x + y] += hB[c][x] * hC[c][y];
            }
        double p4 = 0;
        for (int t = -256; t < 256; t++) {
            int u = TARGET - t;
            if (u >= -256 && u < 256)
                p4 += L2[t + 256] * R2s[u + 256];
        }
        if (c >= S) {
            logp_rest += log(p4);
            continue;
        }
        for (int q = 0; q < top; q++)
            cw[c][q] = -1, cand[c][q] = 0;
        for (int t = -256; t < 256; t++) {
            int u = TARGET - t;
            if (u < -256 || u >= 256)
                continue;
            double w = L2[t + 256] * R2s[u + 256];
            for (int q = 0; q < top; q++)
                if (w > cw[c][q]) {
                    for (int z = top - 1; z > q; z--) {
                        cw[c][z] = cw[c][z - 1];
                        cand[c][z] = cand[c][z - 1];
                        pl[c][z] = pl[c][z - 1];
                        pr[c][z] = pr[c][z - 1];
                    }
                    cw[c][q] = w;
                    cand[c][q] = t;
                    pl[c][q] = L2[t + 256];
                    pr[c][q] = R2s[u + 256];
                    break;
                }
        }
    }
    long ncomb = 1;
    for (int c = 0; c < S; c++)
        ncomb *= top;
    double *weight = malloc(ncomb * sizeof(double));
    long *order = malloc(ncomb * sizeof(long));
    for (long z = 0; z < ncomb; z++) {
        long r = z;
        double w = 1;
        for (int c = 0; c < S; c++, r /= top)
            w *= cw[c][r % top];
        weight[z] = w;
        order[z] = z;
    }
    g_weight = weight;
    qsort(order, ncomb, sizeof(long), by_weight);
    {
        long r = order[0];
        double el = 1, er = 1;
        for (int c = 0; c < S; c++, r /= top) {
            el *= pl[c][r % top];
            er *= pr[c][r % top];
        }
        /* one class is determined by the others (row-sum identity): undo its factor */
        double dep = exp(-logp_rest / (C - S));
        fprintf(stderr, "best tau: expected pairs %.3g / %.3g, expected solutions %.3g per tau\n",
                el * A.n * D.n, er * B.n * Cl.n, el * A.n * D.n * er * B.n * Cl.n * exp(logp_rest) * dep);
    }

    uint64_t H2 = 0;
    for (int c = S; c < C; c++)
        H2 += R2[c] * (uint64_t)(int64_t)TARGET;
    Pair *left = NULL, *right = NULL, *tmp = NULL;
    long cap = 0, nsol = 0;
    for (long ti = first_tau; ti < max_tau && ti < ncomb && nsol < max_sol; ti++) {
        long r = order[ti];
        int tau[MAXC];
        for (int c = 0; c < S; c++, r /= top)
            tau[c] = cand[c][r % top];
        uint64_t HL = 0, HR = 0;
        for (int c = 0; c < S; c++) {
            HL += R1[c] * (uint64_t)(int64_t)tau[c];
            HR += R1[c] * (uint64_t)(int64_t)(TARGET - tau[c]);
        }
        long nl = merge(&IA, &ID, HL, 0, 0, NULL, nth), nr = merge(&IB, &IC, HR, 0, 0, NULL, nth);
        long need = nl > nr ? nl : nr;
        if (need > cap) {
            cap = need + need / 8 + 1024;
            free(left);
            free(right);
            free(tmp);
            left = malloc(cap * sizeof(Pair));
            right = malloc(cap * sizeof(Pair));
            tmp = malloc(cap * sizeof(Pair));
            if (!left || !right || !tmp) {
                fprintf(stderr, "out of memory (%ld pairs)\n", cap);
                return 1;
            }
        }
        merge(&IA, &ID, HL, 1, 0, left, nth);
        merge(&IB, &IC, HR, 1, H2, right, nth);
        sort_pairs(left, tmp, nl);
        sort_pairs(right, tmp, nr);
        long i = 0, j = 0, found = 0;
        while (i < nl && j < nr) {
            if (left[i].h < right[j].h)
                i++;
            else if (left[i].h > right[j].h)
                j++;
            else {
                uint64_t h = left[i].h;
                long i2 = i, j2 = j;
                while (i2 < nl && left[i2].h == h)
                    i2++;
                while (j2 < nr && right[j2].h == h)
                    j2++;
                for (long x = i; x < i2; x++)
                    for (long y = j; y < j2; y++) {
                        uint32_t a = left[x].i, d = left[x].j, b = right[y].i, c2 = right[y].j;
                        const int8_t *va = values(&A, a), *vd = values(&D, d), *vb = values(&B, b), *vc = values(&Cl, c2);
                        int ok = 1;
                        for (int c = 0; c < C && ok; c++)
                            ok = va[c] + vb[c] + vc[c] + vd[c] == TARGET;
                        if (ok) {
                            printf("SOL %016llx %016llx %016llx %016llx\n", (unsigned long long)mask(&A, a),
                                   (unsigned long long)mask(&B, b), (unsigned long long)mask(&Cl, c2),
                                   (unsigned long long)mask(&D, d));
                            fflush(stdout);
                            found++;
                            nsol++;
                        }
                    }
                i = i2;
                j = j2;
            }
        }
        fprintf(stderr, "tau #%ld weight %.3g pairs %ld/%ld solutions %ld (total %ld) %.1fs\n", ti, weight[order[ti]],
                nl, nr, found, nsol, omp_get_wtime() - t0);
    }
    fprintf(stderr, "done: %ld solutions in %.1fs\n", nsol, omp_get_wtime() - t0);
    return 0;
}
