#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

/* Benchmark naive integer matrix multiplication.
 * No external libraries needed.
 * Usage: matmul_bench <n>
 */

int main(int argc, char *argv[]) {
  if (argc < 2) { fprintf(stderr, "Usage: %s <n>\n", argv[0]); return 1; }
  int n = atoi(argv[1]);
  struct timeval t0, t1;

  long *A = calloc(n * n, sizeof(long));
  long *C = calloc(n * n, sizeof(long));
  if (!A || !C) { fprintf(stderr, "malloc failed\n"); return 1; }

  /* Build A[i,j] = i^j mod n */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      long v = 1;
      for (int k = 0; k < j; k++) v = (v * i) % n;
      A[i * n + j] = v;
    }

  /* Timed: 3 trials of A*A */
  for (int t = 0; t < 3; t++) {
    gettimeofday(&t0, NULL);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        long s = 0;
        for (int k = 0; k < n; k++)
          s += A[i * n + k] * A[k * n + j];
        C[i * n + j] = s;
      }
    gettimeofday(&t1, NULL);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) * 1e-6;
    printf("%d,%d,%.6f\n", n, t + 1, elapsed);
  }

  free(A);
  free(C);
  return 0;
}
