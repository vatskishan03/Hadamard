/*
 * Independent exact check of a Hadamard matrix given as CSV (one row per line,
 * entries 1 / -1 separated by commas). Reads a file, or standard input if no
 * argument is given, so compressed matrices can be piped in:
 *
 *     cc -O2 -fopenmp verify/check_matrix.c -o check_matrix
 *     gzip -dc 2092/H2092.csv.gz | ./check_matrix
 *
 * Every pair of rows is compared with integer arithmetic; no floating point and
 * no external libraries are involved. Exit status 0 means H H^T = nI.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    FILE *f = argc > 1 ? fopen(argv[1], "r") : stdin;
    if (!f) {
        perror(argv[1]);
        return 2;
    }
    size_t cap = (size_t)1 << 26;
    char *line = malloc(cap);
    signed char *H = NULL;
    long n = 0, rows = 0;

    while (fgets(line, (int)cap, f)) {
        long cnt = 1;
        for (char *p = line; *p; p++)
            cnt += *p == ',';
        if (!n) {
            n = cnt;
            H = malloc((size_t)n * n);
        }
        if (cnt != n || rows >= n) {
            fprintf(stderr, "row %ld: expected %ld entries\n", rows, n);
            return 1;
        }
        char *p = line;
        for (long j = 0; j < n; j++) {
            long v = strtol(p, &p, 10);
            if (v != 1 && v != -1) {
                fprintf(stderr, "row %ld, column %ld: entry %ld is not +-1\n", rows, j, v);
                return 1;
            }
            H[rows * n + j] = (signed char)v;
            if (*p == ',')
                p++;
        }
        rows++;
    }
    if (rows != n) {
        fprintf(stderr, "matrix is %ld x %ld, not square\n", rows, n);
        return 1;
    }

    long bad = 0;
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : bad)
    for (long i = 0; i < n; i++) {
        const signed char *a = H + i * n;
        for (long k = i; k < n; k++) {
            const signed char *b = H + k * n;
            int s = 0;
            for (long j = 0; j < n; j++)
                s += a[j] * b[j];
            if (s != (i == k ? n : 0))
                bad++;
        }
    }
    printf("n = %ld, %ld row pairs, %ld violations: %s\n", n, n * (n + 1) / 2, bad,
           bad ? "NOT a Hadamard matrix" : "H H^T = nI");
    return bad != 0;
}
