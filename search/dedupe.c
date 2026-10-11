/*
 * Remove duplicate candidates (equal masks) from a file written by generate.py.
 *
 *     cc -O2 search/dedupe.c -o dedupe
 *     ./dedupe CLASSES in.bin out.bin
 *
 * CLASSES is the number of autocorrelation classes C; records are 8 + C bytes.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint64_t mask;
    uint64_t index;
} Key;

static int by_mask(const void *a, const void *b)
{
    uint64_t x = ((const Key *)a)->mask, y = ((const Key *)b)->mask;
    return x < y ? -1 : x > y;
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr, "usage: %s CLASSES in.bin out.bin\n", argv[0]);
        return 2;
    }
    long len = 8 + atol(argv[1]);
    FILE *in = fopen(argv[2], "rb");
    if (!in) {
        perror(argv[2]);
        return 1;
    }
    fseek(in, 0, SEEK_END);
    long size = ftell(in);
    fseek(in, 0, SEEK_SET);
    long n = size / len;
    unsigned char *buf = malloc(size);
    Key *keys = malloc(n * sizeof(Key));
    if (!buf || !keys || fread(buf, 1, size, in) != (size_t)size) {
        fprintf(stderr, "cannot read %s\n", argv[2]);
        return 1;
    }
    fclose(in);

    for (long i = 0; i < n; i++) {
        memcpy(&keys[i].mask, buf + i * len, 8);
        keys[i].index = i;
    }
    qsort(keys, n, sizeof(Key), by_mask);

    FILE *out = fopen(argv[3], "wb");
    long kept = 0;
    for (long i = 0; i < n; i++) {
        if (i && keys[i].mask == keys[i - 1].mask)
            continue;
        fwrite(buf + keys[i].index * len, len, 1, out);
        kept++;
    }
    fclose(out);
    fprintf(stderr, "%ld records, %ld distinct\n", n, kept);
    return 0;
}
