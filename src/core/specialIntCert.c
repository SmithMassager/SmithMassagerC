#include "specialIntCert.h"
#include "fastIntCert.h"

int specialIntCert(fmpz_mat_t ret, fmpz_t s, fmpz_mat_t A, fmpz_mat_t B, int n, int m, int r) {
  fmpz_mat_t BB, AA, res;
  fmpz_mat_init(BB, n+r, m);
  fmpz_mat_init(AA, n+r, n+r);
  fmpz_mat_init(res, n+r, m);
  fmpz_mat_zero(BB);

  // BB[1..n, 1..m] = B[1..n, 1..m]
  // AA[1..n, 1..n] = A[1..n, 1..n]
  // AA[1..n, -r..-1] := A[1..n, -r..-1]
  for (int i = 0; i < n; ++i) {
    _fmpz_vec_set(fmpz_mat_entry(BB, i, 0), fmpz_mat_entry(B, i, 0), B->c);
    _fmpz_vec_set(fmpz_mat_entry(AA, i, 0), fmpz_mat_entry(A, i, 0), n);
    _fmpz_vec_set(fmpz_mat_entry(AA, i, 0) + AA->c - r, fmpz_mat_entry(A, i, 0) + A->c - r, r);
  }
  // AA[-r..-1, 1..n] := A[-r..-1, 1..n]
  // AA[-r..-1, -r..-1] := A[-r..-1, -r..-1]
  for (int i = 1; i <= r; ++i) {
    _fmpz_vec_set(fmpz_mat_entry(AA, n+r-i, 0), fmpz_mat_entry(A, 2*n-i, 0), n);
    _fmpz_vec_set(fmpz_mat_entry(AA, n+r-i, 0) + AA->c - r, fmpz_mat_entry(A, 2*n-i, 0) + A->c - r, r);
  }

  int success = fmpzFastIntCert(s, res, AA, BB, n+r, m);
  if (!success) {
    printf("Warning: specialInt is not successful\n");
    goto clear;
  }

  // ret[1..n, 1..m] := res[1..n, 1..m]
  for (int i = 0; i < n; ++i) {
    _fmpz_vec_set(fmpz_mat_entry(ret, i, 0), fmpz_mat_entry(res, i, 0), m);
  }
  // ret[-r..-1, 1..m] := res[-r..-1, 1..m]
  for (int i = 1; i <= r; ++i) {
    _fmpz_vec_set(fmpz_mat_entry(ret, 2*n-i, 0), fmpz_mat_entry(res, n+r-i, 0), m);
  }

clear:
  fmpz_mat_clear(BB);
  fmpz_mat_clear(AA);
  fmpz_mat_clear(res);

  return success;
}
