#include "fmpq_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "maple_call.h"
#include "nmod.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

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
#include "rns_conversion.h"
#include "reconstruct.h"

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

void mapleToFmpzInteger(MKernelVector kv, fmpz_t dst, ALGEB val) {
  mpz_ptr tmp;
  tmp = MapleToGMPInteger(kv, val);
  fmpz_set_mpz(dst, tmp);
  mpz_clear(tmp);
}

ALGEB FMPZIntegerToMaple(MKernelVector kv, fmpz_t val) {
  mpz_t tmp;
  ALGEB ret;
  mpz_init(tmp);

  fmpz_get_mpz(tmp, val);
  ret = GMPIntegerToMaple(kv, (mpz_ptr)tmp);
  mpz_clear(tmp);

  return ret;
}

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

static void fmpzMat_fromMpzMat(fmpz_mat_t Af, mpzMatrix_t *A) {
  for (long i = 0; i < A->nrows; ++i) {
    for (long j = 0; j < A->ncols; ++j) {
      fmpz_set_mpz(fmpz_mat_entry(Af, i, j), mmget(A, i, j));
    }
  }
}

static void fmpzMat_fromRTable(fmpz_mat_t Af, MKernelVector kv, ALGEB rt) {
  mpzMatrix_t *A;
  A = mpzMatrix_fromRTable(kv, rt);
  fmpz_mat_init(Af, A->nrows, A->ncols);
  struct timeval start, end;
  //gettimeofday(&start, 0);
  fmpzMat_fromMpzMat(Af, A);
  //gettimeofday(&end, 0);
  //printf("Time to convert from mpzMatrix to fmpzMatrix %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
  mpzMatrix_fini(A);
}

static mpzMatrix_t* fmpzMat_toMpzMat(fmpz_mat_t Af) {
  mpzMatrix_t *A;
  A = mpzMatrix_init(fmpz_mat_nrows(Af), fmpz_mat_ncols(Af));
  for (long i = 0; i < A->nrows; ++i) {
    for (long j = 0; j < A->ncols; ++j) {
      fmpz_get_mpz(mmget(A, i, j), fmpz_mat_entry(Af, i, j));
    }
  }
  return A;
}

static ALGEB fmpzMat_toRTable(MKernelVector kv, fmpz_mat_t Af) {
  mpzMatrix_t *A;
  ALGEB rt;
  struct timeval start, end;
  //gettimeofday(&start, 0);
  A = fmpzMat_toMpzMat(Af);
  //gettimeofday(&end, 0);
  //printf("Time to convert from fmpzMatrix to mpzMatrix %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));

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

unsigned long fmpz_mat_find_good_prime_and_solve(nmod_mat_t Xmod,
		                 nmod_mat_t Amod, nmod_mat_t Bmod,
                const fmpz_mat_t A, const fmpz_mat_t B, const fmpz_t det_bound)
{
    ulong p;
    fmpz_t tested;

    p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
    fmpz_init(tested);
    fmpz_one(tested);

    while (1)
    {
        p = n_nextprime(p, 0);
        nmod_mat_set_mod(Xmod, p);
        nmod_mat_set_mod(Amod, p);
        nmod_mat_set_mod(Bmod, p);
        fmpz_mat_get_nmod_mat(Amod, A);
        fmpz_mat_get_nmod_mat(Bmod, B);
        if (nmod_mat_solve(Xmod, Amod, Bmod))
            break;
        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
        {
            p = 0;
            break;
        }
    }

    fmpz_clear(tested);
    return p;
}

void
_fmpz_mat_solve_multi_mod(fmpz_mat_t Ret,
                        const fmpz_mat_t A, const fmpz_mat_t B,
                     nmod_mat_t Xmod, nmod_mat_t Amod, nmod_mat_t Bmod,
		                   ulong p, const fmpz_t N, const fmpz_t D)
{
    fmpz_t bound, pprod;
    fmpz_mat_t x;
    fmpq_mat_t AX, X;
    slong i, n, nexti, cols;
    int stabilised; /* has CRT stabilised */

    n = A->r;
    cols = B->c;

    fmpz_init(bound);
    fmpz_init(pprod);

    fmpq_mat_init(AX, B->r, B->c);
    fmpq_mat_init(X, n, cols);
    fmpz_mat_init(x, n, cols);

    /* Compute bound for the needed modulus. TODO: if one of N and D
       is much smaller than the other, we could use a tighter bound (i.e. 2ND).
       This would require the ability to forward N and D to the
       CRT routine.
     */
    if (fmpz_cmpabs(N, D) < 0)
        fmpz_mul(bound, D, D);
    else
        fmpz_mul(bound, N, N);
    fmpz_mul_ui(bound, bound, UWORD(2));  /* signs */
    fmpz_set(bound, D);

    fmpz_set_ui(pprod, p);
    fmpz_mat_set_nmod_mat(x, Xmod);

    i = 1; /* working with i primes */
    nexti = 1; /* when to do next termination test */

    while (fmpz_cmp(pprod, bound) <= 0)
    {
        stabilised = i == nexti;
        stabilised = 0;
        if (stabilised) /* set next termination test iteration */
            nexti = (slong)(i*1.4) + 1;

        /* full matrix stabilisation check */
        if (stabilised)
        {
            stabilised = fmpq_mat_set_fmpz_mat_mod_fmpz(X, x, pprod);

	        if (stabilised)
            {
                if (_fmpq_mat_check_solution_fmpz_mat(X, A, B))
                    goto multi_mod_done;
            }
        }
        i++;

        while (1)
        {
           p = n_nextprime(p, 1);

           nmod_mat_set_mod(Xmod, p);
           nmod_mat_set_mod(Amod, p);
           nmod_mat_set_mod(Bmod, p);
           fmpz_mat_get_nmod_mat(Amod, A);
           fmpz_mat_get_nmod_mat(Bmod, B);
           if (nmod_mat_solve(Xmod, Amod, Bmod))
              break;
        }

        fmpz_mat_CRT_ui(x, x, pprod, Xmod, 0);
        fmpz_mul_ui(pprod, pprod, p);
    }

multi_mod_done:
    fmpz_mat_set(Ret, x);
    fmpz_mat_scalar_smod(Ret, Ret, pprod);

    fmpz_clear(bound);
    fmpz_clear(pprod);

    fmpq_mat_clear(AX);
    fmpz_mat_clear(x);

}

int
fmpz_mat_solve_fmpz_mat_multi_mod_bound(fmpz_mat_t X,
                        const fmpz_mat_t A, const fmpz_mat_t B, const fmpz_t bound)
{
    nmod_mat_t Xmod, Amod, Bmod;
    fmpz_t N, D;
    ulong p;

    if (!fmpz_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_mat_solve_fmpz_mat_multi_mod). Non-square system matrix.\n");
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
        return 1;

    fmpz_init(N);
    fmpz_init(D);
    fmpz_mat_solve_bound(N, D, A, B);

    nmod_mat_init(Amod, A->r, A->c, 1);
    nmod_mat_init(Bmod, B->r, B->c, 1);
    nmod_mat_init(Xmod, B->r, B->c, 1);

    //printf("calling fmpz_mat_find_good_prime\n");
    p = fmpz_mat_find_good_prime_and_solve(Xmod, Amod, Bmod, A, B, D);
    //printf("calling _solve_multi\n");
    if (p != 0)
        _fmpz_mat_solve_multi_mod(X, A, B, Xmod, Amod, Bmod, p, bound, bound);

    nmod_mat_clear(Xmod);
    nmod_mat_clear(Bmod);
    nmod_mat_clear(Amod);
    fmpz_clear(N);
    fmpz_clear(D);

    return p != 0;
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
  fmpz_mat_t Af, sBRf, Resf;
  fmpz_t tmpf, denf;
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
  fmpz_init(tmpf);
  fmpz_init(denf);

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

  fmpz_mat_init(Af, A->nrows, A->ncols);
  fmpz_mat_init(sBRf, sBR->nrows, sBR->ncols);
  fmpz_mat_init(Resf, Res->nrows, Res->ncols);
  fmpzMat_fromMpzMat(Af, A);
  fmpzMat_fromMpzMat(sBRf, sBR);

  // Calculate A^(-1)sRB and store it into Res 
  //if (false) {
  if (mpz_cmp(s, tmp) <= 0) {
    ////printf("calling reconstruct\n");
    //fmpz_t thressqrt;

    //fmpz_init(thressqrt);
    //fmpz_set_mpz(thressqrt, thres);
    ////fmpz_sqrt(thressqrt, thressqrt);
    ////fmpz_add_si(thressqrt, thressqrt, 1);
    ////printf("calling fmpz_mat_solve\n");
    //fmpz_mat_solve_fmpz_mat_multi_mod_bound(Resf, Af, sBRf, thressqrt);
    ////printf("converting Q \n");
    ////fmpq_mat_get_fmpz_mat_matwise(Resf, denf, Q);
    ////printf("finished converting Q \n");
    ////fmpz_mat_scalar_mul_fmpz(Resf, Resf, denf);
    ////printf("converting to mpzMat\n");
    //Res = fmpzMat_toMpzMat(Resf);

    ////printf("clearing fmp\n");
    //fmpz_clear(thressqrt);
    ////printf("finished clearing fmp\n");

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

    // Check the answer A^(-1)sRB is integral iff ||rem(A^(-1)sRB, thres)||  < 0.6*s*n*||B||
    mpzMatrix_max(tmp, Res);
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
    //printf("calling imlSolve\n");
    fmpz_mat_solve_dixon_den(Resf, denf, Af, sBRf);
    //fmpz_mat_solve(Resf, denf, Af, sBRf);
    //fmpz_mat_solve_multi_mod_den(Resf, denf, Af, sBRf);
    Res = fmpzMat_toMpzMat(Resf);
    fmpz_get_mpz(tmp, denf);
    //imlSolve(Res, tmp, A, sBR);
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
  fmpz_mat_clear(Af);
  fmpz_mat_clear(sBRf);
  fmpz_mat_clear(Resf);
  fmpz_clear(tmpf);
  fmpz_clear(denf);

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

  y_rt = mpzMatrix_toRTable(kv, y);
  d_rt = GMPIntegerToMaple(kv, d);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  mpz_clear(d);

  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 2, y_rt, d_rt);
}

ALGEB M_DECL dixonSolve_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, b_rt, y_rt, d_rt;
  fmpz_mat_t A, b, y;
  fmpz_t d, gcd;

  argc = MapleNumArgs(kv,args);
  if (argc != 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  b_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  fmpzMat_fromRTable(A, kv, A_rt);
  if (!A) { return NULL; }
  fmpzMat_fromRTable(b, kv, b_rt);
  if (!b) { return NULL; }
  fmpz_init(d);
  fmpz_init(gcd);
  fmpz_mat_init(y, fmpz_mat_nrows(b), fmpz_mat_ncols(b));

  //printf("Calling fmpz_mat_solve_dixon_den\n");
  fmpz_mat_solve_dixon_den(y, d, A, b);
  fmpz_mat_content(gcd, y);
  fmpz_gcd(gcd, gcd, d);
  fmpz_mat_scalar_divexact_fmpz(y, y, gcd);
  fmpz_divexact(d, d, gcd);

  y_rt = fmpzMat_toRTable(kv, y);
  //printf("calling fmpz to maple\n");
  d_rt = FMPZIntegerToMaple(kv, d);
  //printf("finish fmpz to maple\n");

  fmpz_mat_clear(A);
  fmpz_mat_clear(b);
  fmpz_mat_clear(y);
  fmpz_clear(d);
  fmpz_clear(gcd);

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
  if (argc >= 4) { 
    mod = MapleToGMPInteger(kv, mod_int);
    mpzMatrix_mods(A, mod);
    mpzMatrix_mods(b, mod);
  }

  mpzMatrix_rnsGemm(y, A, b);
  if (argc >= 3) { mpzMatrix_scale(y, d); }
  if (argc >= 4) { mpzMatrix_mods(y, mod); }
  
  y_rt = mpzMatrix_toRTable(kv, y);

  mpzMatrix_fini(A);
  mpzMatrix_fini(b);
  mpzMatrix_fini(y);
  if (argc >= 3) mpz_clear(d);
  if (argc >= 4) mpz_clear(mod);
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

  //printf("starting to convert from rtable to fmpz_mat\n");
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

  //gettimeofday(&end, 0);
  //printf("Time to convert from Rtable to mpzMatrix %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));

  //gettimeofday(&start, 0);
  cmodMulviaColPL(y, A, B, F);
  //gettimeofday(&end, 0);
  //printf("Time to call cmodMulviaColPL %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
 
  //gettimeofday(&start, 0);
  y_rt = fmpzMat_toRTable(kv, y);

  fmpz_mat_clear(A);
  fmpz_mat_clear(B);
  fmpz_mat_clear(F);
  free(F);
  fmpz_mat_clear(y);

  MaplePopGMPAllocators(kv);

  //gettimeofday(&end, 0);
  //printf("Time deconstruct and return back to maple %.3f\n", (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6));
  return ToMapleExpressionSequence(kv, 1, y_rt);
}

ALGEB M_DECL cmodMulviaPL_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, B_rt, F_rt, y_rt;
  fmpz_mat_t A, B, y;
  fmpz_mat_struct *F = malloc(sizeof(fmpz_mat_struct));

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

  cmodMulviaPL(y, A, B, F);
  
  y_rt = fmpzMat_toRTable(kv, y);

  fmpz_mat_clear(A);
  fmpz_mat_clear(B);
  fmpz_mat_clear(F);
  free(F);
  fmpz_mat_clear(y);

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

// A, args[1] reprsenting a matrix with appropiate size.
// b, args[2] reprsenting a vector/matrix with appropiate size.
// d, args[3] optinal integer.
// s, args[4] optional modulus.
// Returns dAb mods s
ALGEB M_DECL fmpzMult_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB A_rt, B_rt, C_rt, d_int, mod_int;
  fmpz_mat_t A, B, C;
  fmpz_t d, mod;

  argc = MapleNumArgs(kv,args);
  if (argc < 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "at least two arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  B_rt = (ALGEB)args[2];
  if (argc >= 3 && IsMapleInteger(kv, args[3])) { d_int = (ALGEB)args[3]; }
  if (argc >= 4 && IsMapleInteger(kv, args[4])) { mod_int = (ALGEB)args[4]; }


  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);

  fmpzMat_fromRTable(A, kv, A_rt);
  fmpzMat_fromRTable(B, kv, B_rt);
  fmpz_mat_init(C, fmpz_mat_nrows(A), fmpz_mat_ncols(B));

  if (argc >= 3) { mapleToFmpzInteger(kv, d, d_int); }
  if (argc >= 4) { mapleToFmpzInteger(kv, mod, mod_int); }
  if (argc >= 4) { fmpz_mat_scalar_smod(A, A, mod); fmpz_mat_scalar_smod(B, B, mod); }

  fmpz_mat_mul(C, A, B);
  if (argc >= 3) { fmpz_mat_scalar_mul_fmpz(C, C, d); }
  if (argc >= 4) { fmpz_mat_scalar_smod(C, C, mod); }
 
 
  C_rt = fmpzMat_toRTable(kv, C);

  fmpz_mat_clear(A);
  fmpz_mat_clear(B);
  fmpz_mat_clear(C);
  fmpz_clear(mod);
  fmpz_clear(d);
 
  MaplePopGMPAllocators(kv);

  return ToMapleExpressionSequence(kv, 1, C_rt);
}
