/*
 * Enlarge a candidate list with multiplier images.
 *
 * If x is constant on the cosets of H = <g> in Z_p^*, so is x_u(t) = x(ut) for every
 * unit u, and P_{x_u}(s) = P_x(us). The class values of x_u are therefore a
 * permutation of those of x, and descent-selected candidates stay selected. We
 * apply one u from each coset of <H, -1> (u = -1 only reverses x and leaves the
 * class values unchanged), so each input record yields up to C records.
 *
 *     cc -O2 search/expand.c -o expand
 *     ./expand p g IMAGES in.bin out.bin     (IMAGES <= C, including the identity; "-" for stdout)
 *
 * Records use the format of generate.py; run dedupe on the output.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    if (argc != 6) {
        fprintf(stderr, "usage: %s p g IMAGES in.bin out.bin\n", argv[0]);
        return 2;
    }
    long p = atol(argv[1]), g = atol(argv[2]);
    int images = atoi(argv[3]);

    /* coset (orbit) and class numbering exactly as in generate.py */
    int *orb = malloc(p * sizeof(int)), *cls = malloc(p * sizeof(int));
    long *orep = malloc(p * sizeof(long)), *crep = malloc(p * sizeof(long));
    int no = 0, nc = 0;
    for (long j = 0; j < p; j++)
        orb[j] = cls[j] = -1;
    for (long j = 0; j < p; j++) {
        if (orb[j] >= 0)
            continue;
        orep[no] = j;
        for (long x = j; orb[x] < 0; x = x * g % p)
            orb[x] = no;
        no++;
    }
    for (long s = 1; s < p; s++) {
        if (cls[s] >= 0)
            continue;
        crep[nc] = s;
        for (long x = s; cls[x] < 0; x = x * g % p)
            cls[x] = cls[p - x] = nc;
        nc++;
    }
    if (no - 1 > 64 || images < 1 || images > nc) {
        fprintf(stderr, "unsupported parameters (cosets %d, classes %d)\n", no - 1, nc);
        return 2;
    }
    long len = 8 + nc;

    /* one multiplier per coset of <H, -1>: the class representatives themselves */
    uint64_t (*bits)[8][256] = calloc(images, sizeof *bits);
    int *vperm = malloc((size_t)images * nc * sizeof(int));
    for (int k = 0; k < images; k++) {
        long u = crep[k];   /* crep[0] = 1 is the identity */
        int perm[64];
        for (int o = 1; o < no; o++)        /* x_u negative on coset o  <=>  x negative on coset of u r_o */
            perm[o - 1] = orb[u * orep[o] % p] - 1;
        for (int byte = 0; byte < 8; byte++)
            for (int val = 0; val < 256; val++) {
                uint64_t m = 0;
                for (int b = 0; b < 8; b++) {
                    int src = byte * 8 + b;  /* source bit in the input mask */
                    if (!((val >> b) & 1) || src >= no - 1)
                        continue;
                    for (int o = 0; o < no - 1; o++)
                        if (perm[o] == src)
                            m |= 1ULL << o;
                }
                bits[k][byte][val] = m;
            }
        for (int c = 0; c < nc; c++)
            vperm[k * nc + c] = cls[u * crep[c] % p];
    }

    FILE *in = fopen(argv[4], "rb"), *out = strcmp(argv[5], "-") ? fopen(argv[5], "wb") : stdout;
    if (!in || !out) {
        perror("open");
        return 1;
    }
    unsigned char *rec = malloc(len), *img = malloc((size_t)images * len);
    long n = 0;
    while (fread(rec, len, 1, in) == 1) {
        uint64_t mask;
        memcpy(&mask, rec, 8);
        for (int k = 0; k < images; k++) {
            unsigned char *r = img + (size_t)k * len;
            uint64_t m = 0;
            for (int byte = 0; byte < 8; byte++)
                m |= bits[k][byte][(mask >> (8 * byte)) & 255];
            memcpy(r, &m, 8);
            for (int c = 0; c < nc; c++)
                r[8 + c] = rec[8 + vperm[k * nc + c]];
        }
        fwrite(img, len, images, out);
        n++;
    }
    fclose(out);
    fprintf(stderr, "%ld records in, %ld out (%d images each)\n", n, n * images, images);
    return 0;
}
