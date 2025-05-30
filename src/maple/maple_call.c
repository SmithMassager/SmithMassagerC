#include "maple_call.h"
#include "fmpz.h"
#include "fmpz_mat.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mpz_matrix.h"
#include "highorder.h"
#include "iherm.h"
//#include "linsys.h"
//#include "horsolve.h"
#include "imlsolve.h"
#include "timer.h"
#include "unicert.h"
#include "lift.h"
#include "pl_matrix.h"

/*static void whattype( MKernelVector kv, ALGEB s )
{
  if( IsMapleAssignedName(kv, s) )        MaplePrintf(kv, "AssignedName\n");
  if( IsMapleComplexNumeric(kv, s) )      MaplePrintf(kv, "ComplexNumeric\n");
  if( IsMapleComplex64(kv, s) )           MaplePrintf(kv, "Complex64\n");
  if( IsMapleNumeric(kv, s) )             MaplePrintf(kv, "Numeric\n");
  if( IsMapleFloat64(kv, s) )             MaplePrintf(kv, "Float64\n");
  if( IsMapleInteger(kv, s) )             MaplePrintf(kv, "Integer\n");
  if( IsMapleInteger8(kv, s) )            MaplePrintf(kv, "Integer8\n");
  if( IsMapleInteger16(kv, s) )           MaplePrintf(kv, "Integer16\n");
  if( IsMapleInteger32(kv, s) )           MaplePrintf(kv, "Integer32\n");
  if( IsMapleInteger64(kv, s) )           MaplePrintf(kv, "Integer64\n");
  if( IsMapleList(kv, s) )                MaplePrintf(kv, "List\n");
  if( IsMapleExpressionSequence(kv, s) )  MaplePrintf(kv, "ExprSeq\n");
  if( IsMapleName(kv, s) )                MaplePrintf(kv, "Name\n");
  if( IsMapleNULL(kv, s) )                MaplePrintf(kv, "Null\n");
  if( IsMaplePointer(kv, s) )             MaplePrintf(kv, "Pointer\n");
  if( IsMaplePointerNULL(kv, s) )         MaplePrintf(kv, "PointerNull\n");
  if( IsMapleProcedure(kv, s) )           MaplePrintf(kv, "Proc\n");
  if( IsMapleRTable(kv, s) )              MaplePrintf(kv, "RTable\n");
  if( IsMapleSet(kv, s) )                 MaplePrintf(kv, "Set\n");
  if( IsMapleStop(kv, s) )                MaplePrintf(kv, "Stop\n");
  if( IsMapleString(kv, s) )              MaplePrintf(kv, "String\n");
  if( IsMapleTable(kv, s) )               MaplePrintf(kv, "Table\n");
  if( IsMapleUnassignedName(kv, s) )      MaplePrintf(kv, "Unassignedname\n");
  if( IsMapleUnnamedZero(kv, s) )         MaplePrintf(kv, "UnnamedZero\n");
}*/

static mpzMatrix_t * mpzMatrix_fromRTable( MKernelVector kv, ALGEB rt )
{
  long src_idx, dst_idx = 0, nrows, ncols;
  mpzMatrix_t * A;
  ALGEB * A_mpl;
  mpz_srcptr tmp;
  RTableSettings rts;

  assert(IsMapleRTable(kv, rt));

  RTableGetSettings(kv, &rts, rt);
  assert(rts.subtype == RTABLE_MATRIX);
  assert(rts.num_dimensions == 2);

  nrows = RTableUpperBound(kv, rt, 1);
  ncols = RTableUpperBound(kv, rt, 2);
  assert(nrows*ncols == RTableNumElements(kv, rt));

  A_mpl = (ALGEB*)RTableDataBlock(kv, rt);

  A = mpzMatrix_init(nrows, ncols);

  for(src_idx = 0; src_idx < nrows*ncols; ++src_idx) {
      tmp = MapleToGMPInteger(kv, A_mpl[src_idx]);
      switch (rts.order) {
        case RTABLE_C:
          dst_idx = src_idx;
          break;
        case RTABLE_FORTRAN:
          dst_idx = (src_idx/nrows)+(src_idx%nrows)*ncols;
          break;
        default:
          assert(0);
          break;
      }
      mpz_set(A->data[dst_idx], tmp);
  }
  return A;
}

static ALGEB mpzMatrix_toRTable( MKernelVector kv, mpzMatrix_t const * A )
{
  long i;
  RTableSettings rts;
  M_INT bounds[4];
  ALGEB * dst;
  mpz_t const * src;
  ALGEB rt;

  RTableGetDefaults(kv, &rts);
  rts.order = RTABLE_C;
  rts.num_dimensions = 2;
  rts.subtype = RTABLE_MATRIX;
  rts.data_type = RTABLE_DAG;

  bounds[0] = 1;
  //bounds[1] = A->m;
  bounds[1] = A->nrows;
  bounds[2] = 1;
  //bounds[3] = A->n;
  bounds[3] = A->ncols;

  rt = RTableCreate(kv,&rts, NULL, bounds);
  dst = (ALGEB*)RTableDataBlock(kv, rt);
  src = mpzMatrix_constData(A);

  for (i = 0; i < A->ncols * A->nrows; ++i) {
    dst[i] = GMPIntegerToMaple(kv, (mpz_ptr)src[i]);
  }
  return rt;
}

static void fmpzMat_fromRTable(fmpz_mat_t Af, MKernelVector kv, ALGEB rt) {
  mpzMatrix_t *A;
  A = mpzMatrix_fromRTable(kv, rt);
  fmpz_mat_init(Af, A->nrows, A->ncols);
  struct timeval start, end;
  gettimeofday(&start, 0);
  for (long i = 0; i < A->nrows; ++i) {
    for (long j = 0; j < A->ncols; ++j) {
      fmpz_set_mpz(fmpz_mat_entry(Af, i, j), mmget(A, i, j));
    }
  }
  gettimeofday(&end, 0);
  printf("Time to convert from mpzMatrix to fmpzMatrix %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
  mpzMatrix_fini(A);
}

static ALGEB fmpzMat_toRTable(MKernelVector kv, fmpz_mat_t Af) {
  mpzMatrix_t *A;
  ALGEB rt;
  A = mpzMatrix_init(fmpz_mat_nrows(Af), fmpz_mat_ncols(Af));
  struct timeval start, end;
  gettimeofday(&start, 0);

  for (long i = 0; i < A->nrows; ++i) {
    for (long j = 0; j < A->ncols; ++j) {
      fmpz_get_mpz(mmget(A, i, j), fmpz_mat_entry(Af, i, j));
    }
  }
  gettimeofday(&end, 0);
  printf("Time to convert from fmpzMatrix to mpzMatrix %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));

  rt =  mpzMatrix_toRTable(kv, A);
  mpzMatrix_fini(A);
  return rt;
}

static openblasResidue_t * openblasResidue_fromRTable( MKernelVector kv, ALGEB rt, double p)
{
  long src_idx, dst_idx = 0, nrows, ncols;
  openblasResidue_t * A;
  ALGEB * A_mpl;
  double tmp;
  RTableSettings rts;

  assert(IsMapleRTable(kv, rt));

  RTableGetSettings(kv, &rts, rt);
  assert(rts.subtype == RTABLE_MATRIX);
  assert(rts.num_dimensions == 2);

  nrows = RTableUpperBound(kv, rt, 1);
  ncols = RTableUpperBound(kv, rt, 2);
  assert(nrows*ncols == RTableNumElements(kv, rt));

  A_mpl = (ALGEB*)RTableDataBlock(kv, rt);

  A = openblas_init(nrows, ncols, p);

  for(src_idx = 0; src_idx < nrows*ncols; ++src_idx) {
      tmp = MapleToFloat64(kv, A_mpl[src_idx]);
      switch (rts.order) {
        case RTABLE_C:
          dst_idx = src_idx;
          break;
        case RTABLE_FORTRAN:
          dst_idx = (src_idx/nrows)+(src_idx%nrows)*ncols;
          break;
        default:
          assert(0);
          break;
      }
      openblas_setEntry(A, dst_idx, (long)tmp);
  }
  return A;
}

static ALGEB openblasResidue_toRTable( MKernelVector kv, openblasResidue_t const * A )
{
  long i;
  RTableSettings rts;
  M_INT bounds[4];
  ALGEB * dst;
  ALGEB rt;

  RTableGetDefaults(kv, &rts);
  rts.order = RTABLE_C;
  rts.num_dimensions = 2;
  rts.subtype = RTABLE_MATRIX;
  rts.data_type = RTABLE_DAG;

  bounds[0] = 1;
  //bounds[1] = A->m;
  bounds[1] = A->nrows;
  bounds[2] = 1;
  //bounds[3] = A->n;
  bounds[3] = A->ncols;

  rt = RTableCreate(kv,&rts, NULL, bounds);
  dst = (ALGEB*)RTableDataBlock(kv, rt);

  for (i = 0; i < A->ncols * A->nrows; ++i) {
    dst[i] = ToMapleInteger64(kv, openblas_getEntry(A, i));
  }
  return rt;
}

static void * alloc_func(size_t sz) { return malloc(sz); }
static void * realloc_func(void * ptr, size_t old_sz, size_t new_sz) { (void)old_sz; return realloc(ptr, new_sz); }
static void free_func(void * ptr, size_t sz) { (void)sz; free(ptr); }

ALGEB M_DECL highOrderResidue_maple( MKernelVector kv, ALGEB args )
{
  int argc;
  mpzMatrix_t * A;
  ALGEB A_rt;
  mpzMatrix_t * R;
  ALGEB R_rt;

  argc = MapleNumArgs(kv,args);
  if (argc != 1 || !IsMapleRTable(kv, (ALGEB)args[1])) {
    MapleRaiseError(kv,"one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }

  R = highOrderResidue(A);
  R_rt = mpzMatrix_toRTable(kv, R);

  mpzMatrix_fini(A);
  mpzMatrix_fini(R);

  MaplePopGMPAllocators(kv);

  return R_rt;
}

ALGEB M_DECL iherm_maple( MKernelVector kv, ALGEB args )
{
  int argc;
  mpzMatrix_t * A;
  ALGEB A_rt;
  mpzMatrix_t * H;
  ALGEB H_rt;

  argc = MapleNumArgs(kv,args);
  if (argc != 1 || !IsMapleRTable(kv, (ALGEB)args[1])) {
    MapleRaiseError(kv,"one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }

  //H = pkMatrix_toFull(myHermite(A));
  H = hermiteRect(A);

  H_rt = mpzMatrix_toRTable(kv, H);

  mpzMatrix_fini(A);
  mpzMatrix_fini(H);

  MaplePopGMPAllocators(kv);

  return H_rt;
}

// (A, B, s)
// Calculates A^(-1)sB and if it's integral return A^(-1)sB mods s
// A is nxn, B is nxm.
ALGEB M_DECL fastIntCert_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  mpzMatrix_t *A, *B, *R, *Res, *sBR;
  rnsMatrix_t *Ap, *Bp, *Cp;
  basis_t *P;
  ALGEB A_rt, B_rt, s_m;
  mpz_ptr s, h;
  mpz_t thres, BMX, tmp, bound, AMX;
  int nprimes, k;
  long *primes, n;

  argc = MapleNumArgs(kv,args);
  if (argc != 3 
      || !IsMapleRTable(kv, (ALGEB)args[1]) 
      || !IsMapleRTable(kv, (ALGEB)args[2]) 
      || !IsMapleInteger(kv, args[3])) {
    MapleRaiseError(kv,"one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];
  s_m = (ALGEB)args[3];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  B = mpzMatrix_fromRTable(kv, B_rt);
  Res = mpzMatrix_init(B->nrows, B->ncols);
  sBR = mpzMatrix_init(B->nrows, B->ncols);
  s = MapleToGMPInteger(kv, s_m);
  n = A->nrows;
  mpz_init_set_ui(thres, 1);
  mpz_init(BMX);
  mpz_init(AMX);
  mpz_init_set_ui(bound, 1);
  mpz_init(tmp);
  if (!A || !B) { return NULL; }

  // Claculate thres >= 1.2 * s * n^2 * ||A|| * ||B||
  mpz_mul(thres, thres, s);
  mpzMatrix_max(AMX, A);
  mpz_mul(thres, thres, AMX);
  mpzMatrix_max(BMX, B);
  mpz_mul(thres, thres, BMX);
  mpz_mul_ui(thres, thres, 12);
  mpz_mul_ui(thres, thres, n);
  mpz_mul_ui(thres, thres, n);
  mpz_cdiv_q_ui(thres, thres, 10);

  //R = highOrderResidue(A);
  mpzMatrix_t * rslt;
  lift_info_t * info;
  // Bound = 2 s n^(n/2) ||A||^(n-1) ||B||
  mpz_pow_ui(tmp, AMX, n-1);
  mpz_ui_pow_ui(bound, n, n/2+1);
  mpz_mul(bound, bound, BMX);
  mpz_mul_ui(bound, bound, 2);
  mpz_mul(bound, bound, s);
  mpz_mul(bound, bound, tmp);
  //info = horBaseFromBound(A, NULL, bound);
  info = horBase(A, NULL);
  R = mpzMatrix_init(n, n);
  mpzMatrix_reconstruct(R, liftInfo_R(info));
  h = liftInfo_modulus(info);
  k = liftInfo_iter(info);

  // Calculate sRB
  //mpzMatrix_gemm(Res, R, B);
  mpzMatrix_rnsGemm(sBR, R, B);
  mpzMatrix_scale(sBR, s);
  int returnFalse;
  mpz_set_ui(tmp, 10); mpz_pow_ui(tmp, tmp, 60);
  // Calculate A^(-1)sRB and store it into Res 
  if (mpz_cmp(s, tmp) <= 0) {
    primes = genCoPrimes(pickStartModulus(n), thres, &nprimes, s);
    P = basis_init(primes, nprimes);  
    Ap = rnsMatrix_init(n, n, P);
    Bp = rnsMatrix_init(B->nrows, B->ncols, P);
    Cp = rnsMatrix_init(B->nrows, B->ncols, P);
    rnsMatrix_fromMpzMatrix(Ap, A);
    rnsMatrix_fromMpzMatrix(Bp, sBR);
    rnsMatrix_inverse(Ap);
    rnsMatrix_gemm(Cp, Ap, Bp);
    //rnsMatrix_print(Ap);
    //rnsMatrix_print(Bp);
    //rnsMatrix_print(Cp);
    mpzMatrix_reconstruct(Res, Cp);
    //mpzMatrix_print(stdout, Res);
    mpzMatrix_max(tmp, Res);
    // Check the answer A^(-1)sRB is integral iff ||rem(A^(-1)sRB, thres)||  < 0.6*s*n*||B||
    mpz_mul(BMX, BMX, s);
    mpz_mul_ui(BMX, BMX, 6);
    mpz_mul_ui(BMX, BMX, n);
    mpz_cdiv_q_ui(BMX, BMX, 10);
    returnFalse = (mpz_cmp(tmp, BMX) >= 0);
    if (!returnFalse)  {
      // Calculate Res * h^k mod s
      mpzMatrix_modp(Res, s);
      mpz_powm_ui(h, h, k, s);
      mpzMatrix_scale(Res, h);
      mpzMatrix_modp(Res, s);
      B_rt = mpzMatrix_toRTable(kv, Res);
    }
    rnsMatrix_fini(Ap);
    rnsMatrix_fini(Bp);
    rnsMatrix_fini(Cp);
    basis_fini(P);
    free(primes);
  } else {
    imlSolve(Res, tmp, A, sBR);
    returnFalse = (mpz_cmp_ui(tmp, 1) != 0);
    if (!returnFalse) { mpzMatrix_modp(Res, s); B_rt = mpzMatrix_toRTable(kv, Res); }
  }

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  mpzMatrix_fini(R);
  mpzMatrix_fini(Res);
  mpzMatrix_fini(sBR);
  mpz_clear(thres);
  mpz_clear(BMX);
  mpz_clear(tmp);
  MaplePopGMPAllocators(kv);
  finiLift(info);

  if (returnFalse) return ToMapleBoolean(kv, 0);
  else return ToMapleExpressionSequence(kv, 1, B_rt);
  //return ToMapleExpressionSequence(kv, 1, B_rt);

}

ALGEB M_DECL unicert_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  mpzMatrix_t * A;
  ALGEB A_rt;
  int rslt;

  argc = MapleNumArgs(kv,args);
  if (argc != 1 || !IsMapleRTable(kv, (ALGEB)args[1])) {
    MapleRaiseError(kv,"one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }

  rslt = uniCert(A);
  mpzMatrix_fini(A);

  MaplePopGMPAllocators(kv);

  return ToMapleBoolean(kv, rslt);
}

ALGEB M_DECL horSolveIML_maple(MKernelVector kv, ALGEB args)
{
  int argc, Rzero;
  ALGEB A_rt, b_rt, y_rt, d_rt, Rzero_maple;
  mpzMatrix_t * A, * b;
  mpz_t d;
  mpzMatrix_t * y;

  argc = MapleNumArgs(kv,args);
  if (argc != 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  b_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  b = mpzMatrix_fromRTable(kv, b_rt);
  if (!b) { return NULL; }
  mpz_init(d);

  y = horSolveIML(d, &Rzero, A, b);

  y_rt = mpzMatrix_toRTable(kv, y);
  d_rt = GMPIntegerToMaple(kv, d);
  Rzero_maple = ToMapleBoolean(kv, Rzero);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  mpz_clear(d);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 3, Rzero_maple, y_rt, d_rt);
}

ALGEB M_DECL horSolveSpinv_maple(MKernelVector kv, ALGEB args)
{
  int argc, Rzero;
  ALGEB A_rt, b_rt, y_rt, d_rt, Rzero_maple;
  mpzMatrix_t * A, * b;
  mpz_t d;
  mpzMatrix_t * y;

  argc = MapleNumArgs(kv,args);
  if (argc != 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  b_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  b = mpzMatrix_fromRTable(kv, b_rt);
  if (!b) { return NULL; }

  mpz_init(d);

  y = horSolveSpinv(d, &Rzero, A, b);

  y_rt = mpzMatrix_toRTable(kv, y);
  d_rt = GMPIntegerToMaple(kv, d);
  Rzero_maple = ToMapleBoolean(kv, Rzero);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  mpz_clear(d);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 3, Rzero_maple, y_rt, d_rt);
}

ALGEB M_DECL imlSolve_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, b_rt, y_rt, d_rt;
  mpzMatrix_t * A, * b;
  mpz_t d;
  mpzMatrix_t * y;

  argc = MapleNumArgs(kv,args);
  if (argc != 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  b_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  b = mpzMatrix_fromRTable(kv, b_rt);
  if (!b) { return NULL; }
  mpz_init(d);
  y = mpzMatrix_init(b->nrows, b->ncols);

  imlSolve(y, d, A, b);
  fflush(stdout);

  y_rt = mpzMatrix_toRTable(kv, y);
  d_rt = GMPIntegerToMaple(kv, d);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  mpz_clear(d);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 2, y_rt, d_rt);
}

ALGEB M_DECL maple_passing(MKernelVector kv, ALGEB args)
{
  int argc;
  mpzMatrix_t * A;
  mpzMatrix_t * B;
  ALGEB A_rt, B_rt;

  argc = MapleNumArgs(kv,args);
  if (argc != 1 || !IsMapleRTable(kv, (ALGEB)args[1])) {
    MapleRaiseError(kv, "one argument expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  B = mpzMatrix_init(A->nrows, A->ncols);

  TIMER("mpz_set", mpzMatrix_set(B, A));

  B_rt = mpzMatrix_toRTable(kv, B);

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);

  MaplePopGMPAllocators(kv);

  return B_rt;
}

// A, args[1] reprsenting a matrix with appropiate size.
// b, args[2] reprsenting a vector/matrix with appropiate size.
// d, args[3] optinal integer.
// s, args[4] optional modulus.
// Returns dAb mods s
ALGEB M_DECL imlMultiply_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, b_rt, y_rt, d_int, mod_int;
  mpzMatrix_t * A, * b;
  mpzMatrix_t * y;
  mpz_ptr mod, d;

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  b_rt = (ALGEB)args[2];
  if (argc >= 3 && IsMapleInteger(kv, args[3])) { d_int = (ALGEB)args[3]; }
  if (argc >= 4 && IsMapleInteger(kv, args[4])) { mod_int = (ALGEB)args[4]; }

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  b = mpzMatrix_fromRTable(kv, b_rt);
  if (!b) { return NULL; }
  y = mpzMatrix_init(A->nrows, b->ncols);

  if (argc >= 3) { d = MapleToGMPInteger(kv, d_int); }
  if (argc >= 4) { mod = MapleToGMPInteger(kv, mod_int); }
  if (argc >= 4) { mpzMatrix_mods(A, mod); mpzMatrix_mods(b, mod); }

  mpzMatrix_rnsGemm(y, A, b);
  if (argc >= 3) { mpzMatrix_scale(y, d); }
  if (argc >= 4) { mpzMatrix_mods(y, mod); }
  
  y_rt = mpzMatrix_toRTable(kv, y);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  if (argc >= 3) mpz_clear(mod);
  if (argc >= 4) mpz_clear(d);
  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 1, y_rt);
}

ALGEB M_DECL colplMultiply_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, B_rt, y_rt;
  mpzMatrix_t *A, *B, *y;

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  B = mpzMatrix_fromRTable(kv, B_rt);
  if (!B) { return NULL; }
  y = mpzMatrix_init(A->nrows, B->ncols);

  MultiplyViaColPL(y, A, B);
  
  y_rt = mpzMatrix_toRTable(kv, y);

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  mpzMatrix_fini(y);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 1, y_rt);
}

ALGEB M_DECL plMultiply_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, B_rt, y_rt;
  mpzMatrix_t *A, *B, *y;

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  B = mpzMatrix_fromRTable(kv, B_rt);
  if (!B) { return NULL; }
  y = mpzMatrix_init(A->nrows, B->ncols);

  MultiplyViaPL(y, A, B);
  
  y_rt = mpzMatrix_toRTable(kv, y);

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  mpzMatrix_fini(y);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 1, y_rt);
}

ALGEB M_DECL cmodMulviaColPL_maple(MKernelVector kv, ALGEB args)
{
  struct timeval start, end;
  gettimeofday(&start, 0);
  int argc;
  ALGEB A_rt, B_rt, F_rt, y_rt;
  fmpz_mat_t A, B, y;
  fmpz_mat_struct *F = malloc(sizeof(fmpz_mat_struct));

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least three arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];
  if (argc >= 3) F_rt = (ALGEB)args[3];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  printf("starting to convert from rtable to fmpz_mat\n");
  fmpzMat_fromRTable(A, kv, A_rt);
  if (!A) { return NULL; }
  fmpzMat_fromRTable(B, kv, B_rt);
  if (!B) { return NULL; }
  if (argc >= 3) {
    fmpzMat_fromRTable(F, kv, F_rt);
    if (!F) { return NULL; }
  } else {
    F = NULL;
  }
  fmpz_mat_init(y, fmpz_mat_nrows(A), fmpz_mat_ncols(B));

  gettimeofday(&end, 0);
  printf("Time to convert from Rtable to mpzMatrix %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));

  gettimeofday(&start, 0);
  cmodMulviaColPL(y, A, B, F);
  gettimeofday(&end, 0);
  printf("Time to call cmodMulviaColPL %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
 
  gettimeofday(&start, 0);
  y_rt = fmpzMat_toRTable(kv, y);

  fmpz_mat_clear(A);
  fmpz_mat_clear(B);
  fmpz_mat_clear(F);
  free(F);
  fmpz_mat_clear(y);

  MaplePopGMPAllocators(kv);

  gettimeofday(&end, 0);
  printf("Time deconstruct and return back to maple %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
  return ToMapleExpressionSequence(kv, 1, y_rt);
}

ALGEB M_DECL cmodMulviaPL_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, B_rt, F_rt, y_rt;
  mpzMatrix_t *A, *B, *F, *y;

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];
  if (argc >= 3) F_rt = (ALGEB)args[3];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  B = mpzMatrix_fromRTable(kv, B_rt);
  if (!B) { return NULL; }
  if (argc >= 3) {
    F = mpzMatrix_fromRTable(kv, F_rt);
    if (!F) { return NULL; }
  } else {
    F = NULL;
  }

  y = mpzMatrix_init(A->nrows, B->ncols);

  cmodMulviaPL(y, A, B, F);
  
  y_rt = mpzMatrix_toRTable(kv, y);

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  mpzMatrix_fini(F);
  mpzMatrix_fini(y);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 1, y_rt);
}

ALGEB M_DECL cmodMulviaPLmpz_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, B_rt, F_rt, y_rt;
  mpzMatrix_t *A, *B, *F, *y;

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];
  if (argc >= 3) F_rt = (ALGEB)args[3];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  B = mpzMatrix_fromRTable(kv, B_rt);
  if (!B) { return NULL; }
  if (argc >= 3) {
    F = mpzMatrix_fromRTable(kv, F_rt);
    if (!F) { return NULL; }
  } else {
    F = NULL;
  }

  y = mpzMatrix_init(A->nrows, B->ncols);

  cmodMulviaColPLMpz(y, A, B, F);
  
  y_rt = mpzMatrix_toRTable(kv, y);

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  mpzMatrix_fini(F);
  mpzMatrix_fini(y);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 1, y_rt);
}

// (A, B, l, mod) Return mods(A + lB, mod)
// Require l to be non-neg.
ALGEB M_DECL add_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, B_rt, C_rt, l_int, mod_int;
  mpzMatrix_t * A, * B;
  mpz_ptr mod, l;

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];
  if (argc >= 3 && IsMapleInteger(kv, args[3])) { l_int = (ALGEB)args[3]; }
  if (argc >= 4 && IsMapleInteger(kv, args[4])) { mod_int = (ALGEB)args[4]; }

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  A = mpzMatrix_fromRTable(kv, A_rt);
  if (!A) { return NULL; }
  B = mpzMatrix_fromRTable(kv, B_rt);
  if (!B) { return NULL; }
  if (argc >= 3) { l = MapleToGMPInteger(kv, l_int); }
  if (argc >= 4) { mod = MapleToGMPInteger(kv, mod_int); }

  if (argc >= 3) {  mpzMatrix_scale(B, l); }
  if (argc >= 4) {  mpzMatrix_mods(A, mod); mpzMatrix_mods(B, mod); }

  mpzMatrix_add(A, B);
  if (argc >= 4) { mpzMatrix_mods(A, mod); }
 
  C_rt = mpzMatrix_toRTable(kv, A);

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
 
  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 1, C_rt);
}

ALGEB M_DECL fmpzMult_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, B_rt, C_rt, l_int, mod_int;
  fmpz_mat_t A, B, C;

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  fmpzMat_fromRTable(A, kv, A_rt);
  fmpzMat_fromRTable(B, kv, B_rt);
  fmpz_mat_init(C, fmpz_mat_nrows(A), fmpz_mat_ncols(B));

  fmpz_mat_mul(C, A, B);
 
  C_rt = fmpzMat_toRTable(kv, C);

  fmpz_mat_clear(A);
  fmpz_mat_clear(B);
  fmpz_mat_clear(C);
 
  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 1, C_rt);
}
