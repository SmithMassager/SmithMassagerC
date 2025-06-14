#include "pl_matrix.h"

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#include "cblas.h"
#include "gmp.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIZ(x) ((x)->_mp_size)
#define ABSIZ(x) ABS (SIZ (x))
#define PTR(x) ((x)->_mp_d)
#define ALLOC(x) ((x)->_mp_alloc)
#define max(x, y) ((x) >= (y) ? (x) : (y))
#define min(x, y) ((x) <= (y) ? (x) : (y))
#define _ONES(x) ((x) == 64-1 ? 0xffffffffffffffff : (1ul<<(x+1))-1)
#define ONES(l, r) ((l == 0) ? _ONES(r) : _ONES(r) & ~_ONES(l-1))

// Set out to represent [-l, -r]-th bits from in's binary reprsentation.
void extractBits(fmpz_t out, const fmpz_t in, long l, long r) {
  if (l > r) return;
  if (!COEFF_IS_MPZ(*in)) {
    fmpz_set_si(out, ((ONES(l, r) & (fmpz_get_si(in) * fmpz_sgn(in))) >> l) * fmpz_sgn(in));
    return;
  }

  long numBits = r - l + 1;
  long numLimbs = (numBits + FLINT_BITS - 1)/FLINT_BITS;
  int bitsWritten = 0;
  int outIdx = 0;
  int outOffset = 0;
  mpz_t outMpz;
  mpz_init2(outMpz, numBits);
  SIZ(outMpz) = numLimbs * fmpz_sgn(in);
  mp_limb_t *inLimbs = PTR(COEFF_TO_PTR(*in)), *outLimbs = PTR(outMpz);
  for (int i = 0; i < ALLOC(outMpz); ++i) { outLimbs[i] = 0; }

  while (bitsWritten < numBits) {
    int curBit = l + bitsWritten;
    int inIdx = curBit/FLINT_BITS;
    int inOffset = curBit % FLINT_BITS;

    int inBitsLeft = FLINT_BITS - inOffset;
    int outBitsLeft = FLINT_BITS - outOffset;
    int bitsRemain = numBits - bitsWritten;
    int bitsToCopy = min(inBitsLeft, min(outBitsLeft, bitsRemain));

    long extracted = (inLimbs[inIdx] >> inOffset) & ONES(0, bitsToCopy-1);
    outLimbs[outIdx] |= (extracted << outOffset);
    bitsWritten += bitsToCopy;
    outOffset += bitsToCopy;

    if (bitsWritten == numBits) break;

    if (outOffset == FLINT_BITS) {
      ++outIdx;
      outOffset = 0;
    }
  }
  fmpz_set_mpz(out, outMpz);
  mpz_clear(outMpz);
}

void reverse(char *str) {
  long slen = strlen(str);
  for (long l = 0; l < slen / 2; l++) {
      char temp = str[l];
      str[l] = str[slen - l - 1];
      str[slen - l - 1] = temp;
  }
}

void fmpzMatrix_maxColEntry(fmpz_t ret, fmpz_mat_t const A, long col)
{
  fmpz_set_ui(ret, 0);
  for (long i = 0; i < fmpz_mat_nrows(A); ++i) {
    if (fmpz_cmpabs(fmpz_mat_entry(A, i, col), ret) > 0) {
      fmpz_abs(ret, fmpz_mat_entry(A, i, col));
    }
  }
}

void fmpzMatrix_maxRowEntry(fmpz_t ret, fmpz_mat_t const A, long row)
{
  fmpz_set_ui(ret, 0);
  for (long i = 0; i < fmpz_mat_ncols(A); ++i) {
    if (fmpz_cmpabs(fmpz_mat_entry(A, row, i), ret) > 0) {
      fmpz_abs(ret, fmpz_mat_entry(A, row, i));
    }
  }
}

long calculateBase(fmpz_mat_t const src) {
  long ans = 0, exp;
  fmpz_t mx;
  fmpz_init(mx);
  for (long i = 0; i < fmpz_mat_ncols(src); ++i) {
    fmpzMatrix_maxColEntry(mx, src, i);
    fmpz_get_d_2exp(&exp, mx);
    ans += exp;
  }
  fmpz_clear(mx);
  return ans/fmpz_mat_ncols(src);
}

plMatrix_t * colPLMatrix_initFromFmpzMatrix(fmpz_mat_t const src, long X) {
  long origRows = fmpz_mat_nrows(src);
  long origCols = fmpz_mat_ncols(src);
  long *eparams = malloc((origCols+1) * sizeof(long));
  plMatrix_t *M = malloc(sizeof(plMatrix_t));
  M->origRows = origRows; M->origCols = origCols; M->eparams = eparams; M->X = X;

  fmpz_t mx, P;
  fmpz_init(mx);
  fmpz_init(P);

  eparams[0] = 0;
  for (long i = 0; i < origCols; ++i) {
    eparams[i+1] = eparams[i];
    fmpzMatrix_maxColEntry(mx, src, i);
    fmpz_add_ui(mx, mx, 1);
    fmpz_set_si(P, 1);
    while(fmpz_cmpabs(P, mx) < 0) {
      fmpz_mul_2exp(P, P, X);
      ++eparams[i+1];
    }
  }

  long newCols = eparams[origCols];
  fmpz_mat_t matrix;
  fmpz_mat_init(matrix, origRows, newCols);
  M->matrix[0] = matrix[0];
  fmpz *entry, *tar;
  long sz;
  long i, j, k, l, len, bitcnts;

  for (i = 0; i < origRows; ++i) {
    for (j = 0; j < origCols; ++j) {
      len = eparams[j+1] - eparams[j];
      if (len == 0) continue;
      else {
        entry = fmpz_mat_entry(src, i, j);
        sz = fmpz_size(entry);
        for (k = 0; k < len; ++k) {
          tar = fmpz_mat_entry(matrix, i, eparams[j] + k);
          if ((k*X)/FLINT_BITS < sz) {
            extractBits(tar, entry, k*X, min((k+1)*X-1, sz*FLINT_BITS-1));
          } else {
            fmpz_set_ui(entry, 0);
          }
        }
      }
    }
  }
  fmpz_clear(mx);
  fmpz_clear(P);

  return M;
}


plMatrix_t * rowPLMatrix_initFromFmpzMatrix(fmpz_mat_t const src, long X) {
  long origRows = fmpz_mat_nrows(src);
  long origCols = fmpz_mat_ncols(src);
  long *eparams = malloc((origRows+1) * sizeof(long));
  plMatrix_t *M = malloc(sizeof(plMatrix_t));
  M->origRows = origRows; M->origCols = origCols; M->eparams = eparams;

  fmpz_t mx, P;
  fmpz_init(mx);
  fmpz_init(P);

  eparams[0] = 0;
  for (long i = 0; i < origRows; ++i) {
    eparams[i+1] = eparams[i];
    fmpzMatrix_maxRowEntry(mx, src, i);
    fmpz_add_ui(mx, mx, 1);
    fmpz_set_si(P, 1);
    while(fmpz_cmpabs(P, mx) < 0) {
      fmpz_mul_2exp(P, P, X);
      ++eparams[i+1];
    }
  }

  long newRows = eparams[origRows];
  fmpz_mat_t matrix;
  fmpz_mat_init(matrix, newRows, origCols);
  M->matrix[0] = matrix[0];
  fmpz *entry, *tar;
  long sz;
  long i, j, k, l, len, bitcnts;

  for (i = 0; i < origRows; ++i) {
    len = eparams[i+1] - eparams[i];
    if (len == 0) { continue;
    } else {
      for (j = 0; j < origCols; ++j) {
	entry = fmpz_mat_entry(src, i, j);
	sz = fmpz_size(entry);
	for (k = 0; k < len; ++k) {
	  tar = fmpz_mat_entry(matrix, eparams[i]+k, j);
	  if ((k * X)/FLINT_BITS < sz) {
	    extractBits(tar, entry, k*X, min((k+1)*X-1, sz*FLINT_BITS-1));
	  } else {
	    fmpz_set_ui(entry, 0);
	  }
	}
      }
    }
  }
  fmpz_clear(mx);
  fmpz_clear(P);

  return M;

}

void plMatrix_fini(plMatrix_t *M) {
  //if (M != NULL) fmpz_mat_clear(M->matrix);
  fmpz_mat_clear(M->matrix);
  if (M->eparams != NULL) free(M->eparams);
  M->origRows = 0;
  M->origCols = 0;
  M->X = 0;
  //M->matrix = NULL;
  M->eparams = NULL;
  free(M);
}

void plMatrix_print(FILE *f, plMatrix_t const * A) {
  fmpz_mat_fprint(f, A->matrix);
}

void MultiplyViaColPL(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const B) {
  cmodMulviaColPL(dst, A, B, NULL);
}

void MultiplyViaPL(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const B) {
  cmodMulviaPL(dst, A, B, NULL);
}

void cmodMulviaColPL(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const B, fmpz_mat_struct const *F) {
  assert(fmpz_mat_nrows(dst) == fmpz_mat_nrows(A) && fmpz_mat_ncols(dst) == fmpz_mat_ncols(B) && fmpz_mat_ncols(A) == fmpz_mat_nrows(B));
  fmpz_mat_zero(dst);
  if (F != NULL) { cmod(B, B, F); }
  //time_t start, end;
  struct timeval start, end;
  gettimeofday(&start,0);
  long X = max(calculateBase(B), calculateBase(A));
  printf("Base chose to be %d\n", X);
  if (X == 0) { fmpz_mat_zero(dst); return; }
  printf("Starting to init Apl\n");
  plMatrix_t *Apl = colPLMatrix_initFromFmpzMatrix(A, X);

  if (Apl->eparams[Apl->origCols] == 0) { fmpz_mat_zero(dst); return; }
  gettimeofday(&end,0);
  printf("Apl %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
  printf("Dimension of Apl %d \n", Apl->eparams[Apl->origCols]);

  //time(&start);
  gettimeofday(&start,0);
  fmpz_mat_t Be; 
  fmpz_mat_init(Be, Apl->eparams[Apl->origCols], fmpz_mat_ncols(B));
  long offset = 0;
  for (long i = 0; i < Apl->origCols; ++i) {
    long len = Apl->eparams[i+1] - Apl->eparams[i];
    if (len == 0) continue;
    for (long j = 0; j < fmpz_mat_ncols(B); ++j) {
      fmpz_set(fmpz_mat_entry(Be, offset, j), fmpz_mat_entry(B, i, j));
    }
    ++offset;
    for (long k = 1; k < len; ++k) {
      for (long j = 0; j < fmpz_mat_ncols(B); ++j) {
	fmpz_mul_2exp(fmpz_mat_entry(Be, offset, j), fmpz_mat_entry(Be, offset-1, j), X);
	if (F != NULL) fmpz_mod(fmpz_mat_entry(Be, offset, j), fmpz_mat_entry(Be, offset, j), fmpz_mat_entry(F, j, j));
      }
      ++offset;
    }
  }
  //time(&end);
  //printf("Be %ld\n", (end - start));
  gettimeofday(&end,0);
  printf("Be %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));

  gettimeofday(&start,0);
  //time(&start);
  plMatrix_t *Bpl = colPLMatrix_initFromFmpzMatrix(Be, X);
  // TODO: uncomment this
  fmpz_mat_clear(Be);
  gettimeofday(&end,0);
  printf("Bpl %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
  gettimeofday(&start,0);
  fmpz_mat_t Cl;
  fmpz_mat_init(Cl, fmpz_mat_nrows(A), Bpl->eparams[fmpz_mat_ncols(B)]);
  printf("Dimension of Apl, Bpl, Cl: %d x %d, %d x %d, %d x %d \n", fmpz_mat_nrows(Apl->matrix), fmpz_mat_ncols(Apl->matrix), fmpz_mat_nrows(Bpl->matrix), fmpz_mat_ncols(Bpl->matrix), fmpz_mat_nrows(Cl), fmpz_mat_ncols(Cl));
  printf("Starting to multiply matrices\n");
  fmpz_mat_mul(Cl, Apl->matrix, Bpl->matrix);
  printf("Finished multipliying matrices\n");
  long *bplEparams = Bpl->eparams;
  Bpl->eparams = NULL;
  printf("starting to deconstruct plMatrices\n");
  plMatrix_fini(Apl);
  plMatrix_fini(Bpl);
  printf("finished to deconstruct plMatrices\n");
  //time(&end);
  //printf("Cl %ld\n", (end - start));
  gettimeofday(&end,0);
  printf("Cl %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
  //time(&start);
  gettimeofday(&start, 0);
  printf("Starting to colCompress\n");
  fflush(stdout);
  colCompress(dst, Cl, X, bplEparams, F, true);
  //time(&end);
  //printf("colCompress %ld\n", (end - start));
  gettimeofday(&end,0);
  printf("colCompress %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));

  free(bplEparams);
  fmpz_mat_clear(Cl);
}

void cmodMulviaPL(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const B, fmpz_mat_struct const *F) {
  fmpz_mat_zero(dst);
  long X = max(calculateBase(B), calculateBase(A));
  if (X == 0) { fmpz_mat_zero(dst); return; }
  plMatrix_t *Apl = rowPLMatrix_initFromFmpzMatrix(A, X);
  if (Apl->eparams[Apl->origRows] == 0) { fmpz_mat_zero(dst); return; }
  plMatrix_t *Bpl = colPLMatrix_initFromFmpzMatrix(B, X);
  if (Bpl->eparams[Bpl->origCols] == 0) { fmpz_mat_zero(dst); return; }

  fmpz_mat_t Cl;
  fmpz_mat_init(Cl, Apl->eparams[Apl->origRows], Bpl->eparams[Bpl->origCols]);
  fmpz_mat_mul(Cl, Apl->matrix, Bpl->matrix);
  fmpz_mat_t C;
  fmpz_mat_init(C, Apl->eparams[Apl->origRows], Bpl->origCols);
  fmpz_mat_zero(C);
  colCompress(C, Cl, X, Bpl->eparams, F, true);
  rowCompress(dst, C, X, Apl->eparams, F, true);

  plMatrix_fini(Apl);
  plMatrix_fini(Bpl);
  fmpz_mat_clear(Cl);
}

void colCompress(fmpz_mat_t dst, fmpz_mat_t const src, long X, long *eparams, fmpz_mat_struct const *F, bool cmod) {
  long rows = fmpz_mat_nrows(dst);
  long cols = fmpz_mat_ncols(dst);
  fmpz_mat_zero(dst);
  fmpz_t val, tmp;
  fmpz_init(val);
  fmpz_init(tmp);

  long offset = 0;
  for (long j = 0; j < cols; ++j) {
    long len = eparams[j+1] - eparams[j];
    if (len == 0) continue;
    fmpz_set_ui(val, 1);
    for (long k = 0; k < len; ++k) {
      for (long i = 0; i < rows; ++i) {
	fmpz *entry = fmpz_mat_entry(dst, i, j);
	fmpz_mul(tmp, val, fmpz_mat_entry(src, i, offset));
	fmpz_add(entry, entry, tmp);
      }
      fmpz_mul_2exp(val, val, X);
      ++offset;
    }
    if (F != NULL) {
      for (long i = 0; i < rows; ++i) {
	fmpz *entry = fmpz_mat_entry(dst, i, j);
	if (cmod) fmpz_mod(entry, entry, fmpz_mat_entry(F, j, j));
	else fmpz_mod(entry, entry, fmpz_mat_entry(F, i, i));
      }
    }
  }

  fmpz_clear(val);
  fmpz_clear(tmp);
}

void rowCompress(fmpz_mat_t dst, fmpz_mat_t const src, long X, long *eparams, fmpz_mat_struct const *F, bool cmod) {
  long rows = fmpz_mat_nrows(dst);
  long cols = fmpz_mat_ncols(dst);
  fmpz_mat_zero(dst);
  fmpz_t val, tmp;
  fmpz_init(val);
  fmpz_init(tmp);

  long offset = 0;
  for (long i = 0; i < rows; ++i) {
    long len = eparams[i+1] - eparams[i];
    if (len == 0) continue;
    fmpz_set_ui(val, 1);
    for (long k = 0; k < len; ++k) {
      for (long j = 0; j < cols; ++j) {
	fmpz *entry = fmpz_mat_entry(dst, i, j);
	fmpz_mul(tmp, val, fmpz_mat_entry(src, offset, j));
	fmpz_add(entry, entry, tmp);
      }
      fmpz_mul_2exp(val, val, X);
      ++offset;
    }
    if (F != NULL) {
      for (long j = 0; j < cols; ++j) {
	fmpz *entry = fmpz_mat_entry(dst, i, j);
	if (cmod) fmpz_mod(entry, entry, fmpz_mat_entry(F, j, j));
	else fmpz_mod(entry, entry, fmpz_mat_entry(F, i, i));
      }
    }
  }

  fmpz_clear(val);
  fmpz_clear(tmp);
}

void cmod(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const C) {
  for (long i = 0; i < fmpz_mat_nrows(A); ++i) {
    for (long j = 0; j < fmpz_mat_ncols(A); ++j) {
      if (fmpz_is_zero(fmpz_mat_entry(C, j, j))) { fmpz_set_ui(fmpz_mat_entry(dst, i, j), 0); }
      else { fmpz_mod(fmpz_mat_entry(dst, i, j), fmpz_mat_entry(A, i, j), fmpz_mat_entry(C, j, j)); }
    }
  }
}

void rmod(fmpz_mat_t dst, fmpz_mat_t const A, fmpz_mat_t const C) {
  for (long i = 0; i < fmpz_mat_nrows(A); ++i) {
    for (long j = 0; j < fmpz_mat_ncols(A); ++j) {
      if (fmpz_is_zero(fmpz_mat_entry(C, i, i))) { fmpz_set_ui(fmpz_mat_entry(dst, i, j), 0); }
      else { fmpz_mod(fmpz_mat_entry(dst, i, j), fmpz_mat_entry(A, i, j), fmpz_mat_entry(C, i, i)); }
    }
  }
}
