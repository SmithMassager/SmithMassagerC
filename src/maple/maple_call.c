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
#include "fastIntCert.h"
#include "fmpz_conv.h"
#include "specialIntCert.h"
#include "indexMassager.h"
#include "largestInvariantFactor.h"
#include "smithMassager.h"

#define min(x, y) ((x < y) ? x : y)

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
  mpz_ptr tmp; // Tmp should not be freed, since it is managed by MKernelVector
  tmp = MapleToGMPInteger(kv, val);
  fmpz_set_mpz(dst, tmp);
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

static mpzMatrix_t *mpzMatrix_fromRTable( MKernelVector kv, ALGEB rt )
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
  fmpzMat_fromMpzMat(Af, A);
  mpzMatrix_fini(A);
}

static ALGEB fmpzMat_toRTable(MKernelVector kv, fmpz_mat_t Af) {
  mpzMatrix_t *A;
  ALGEB rt;
  struct timeval start, end;
  A = fmpzMat_toMpzMat(Af);
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
  mpzMatrix_t *A, *B, *Res;
  ALGEB A_rt, B_rt, s_m;
  mpz_ptr s;

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
  s = MapleToGMPInteger(kv, s_m);
  Res = mpzMatrix_init(B->nrows, B->ncols);

  int success = fastIntCert(s, Res, A, B, A->nrows, B->ncols);

  if (success) {
    B_rt = mpzMatrix_toRTable(kv, Res);
  }

  mpzMatrix_fini(A);
  mpzMatrix_fini(B);
  MaplePopGMPAllocators(kv);
  mpzMatrix_fini(Res);

  if (!success) return ToMapleBoolean(kv, 0);
  else return ToMapleExpressionSequence(kv, 1, B_rt);
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

// Assume P is n x (r+k)
ALGEB M_DECL computeProjBasis_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB P_rt, n_rt, r_rt, s_rt, k_rt, S_rt, U_rt, M_rt, T_rt;
  fmpz_mat_t P, S, U, M, T;
  int n, r, k;
  fmpz_t s;

  argc = MapleNumArgs(kv,args);
  if (argc != 5 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleInteger(kv, (ALGEB)args[2])
                || !IsMapleInteger(kv, (ALGEB)args[3])
                || !IsMapleInteger(kv, (ALGEB)args[4])
                || !IsMapleInteger(kv, (ALGEB)args[5])) {
    MapleRaiseError(kv, "exactly five arguments expected");
    return NULL;
  }

  P_rt = (ALGEB)args[1];
  n_rt = (ALGEB)args[2];
  r_rt = (ALGEB)args[3];
  k_rt = (ALGEB)args[4];
  s_rt = (ALGEB)args[5];
  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);
  printf("initializing all the vars from RTable\n");

  n = MapleToInteger32(kv, n_rt);
  r = MapleToInteger32(kv, r_rt);
  k = MapleToInteger32(kv, k_rt);
  fmpz_init(s);
  mapleToFmpzInteger(kv, s, s_rt);
  fmpz_mat_init(P, n, r+k);
  fmpzMat_fromRTable(P, kv, P_rt);
  fmpz_mat_init(S, r, 1);
  fmpz_mat_init(U, min(r, n), n);
  fmpz_mat_init(M, n, min(r, n));
  fmpz_mat_init(T, min(r, n), min(r, n));
  //printf("finish initializing all the vars from RTable\n");
  
  //printf("calling computeProjBasis\n");
  //printf("%d\n", s);
  //fflush(stdout);
  int success = computeProjBasis(S, U, M, T, P, n, r, s, k);
  //printf("finish calling computeProjBasis\n");
  //fflush(stdout);
  //fmpz_mat_neg(M, M);

  S_rt = fmpzMat_toRTable(kv, S);
  U_rt = fmpzMat_toRTable(kv, U);
  M_rt = fmpzMat_toRTable(kv, M);
  T_rt = fmpzMat_toRTable(kv, T);

  fmpz_mat_clear(P);
  fmpz_mat_clear(S);
  fmpz_mat_clear(U);
  fmpz_mat_clear(M);
  fmpz_mat_clear(T);
  fmpz_clear(s);
  MaplePopGMPAllocators(kv);

  if (success) {
    return ToMapleExpressionSequence(kv, 4, S_rt, U_rt, M_rt, T_rt);
  } else {
    return ToMapleBoolean(kv, 0);
  }
}

ALGEB M_DECL specialIntCert_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB s_rt, A_rt, B_rt, n_rt, m_rt, Res_rt, r_rt;
  fmpz_mat_t A, B, Res;
  int n, m, r;
  fmpz_t s;

  //printf("calling specialIntCert\n");
  argc = MapleNumArgs(kv,args);
  if (argc != 8 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleRTable(kv, (ALGEB)args[2])
                || !IsMapleRTable(kv, (ALGEB)args[3])
                || !IsMapleInteger(kv, (ALGEB)args[4])
                || !IsMapleInteger(kv, (ALGEB)args[5])
                || !IsMapleInteger(kv, (ALGEB)args[6])
                || !IsMapleInteger(kv, (ALGEB)args[7])
                || !IsMapleInteger(kv, (ALGEB)args[8])) {
    MapleRaiseError(kv, "exactly six arguments expected");
    return NULL;
  }

  s_rt = (ALGEB)args[1];
  A_rt = (ALGEB)args[2];
  B_rt = (ALGEB)args[3];
  n_rt = (ALGEB)args[4];
  m_rt = (ALGEB)args[5];
  r_rt = (ALGEB)args[6];
  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);
  n = MapleToInteger32(kv, n_rt);
  m = MapleToInteger32(kv, m_rt);
  r = MapleToInteger32(kv, r_rt);

  mapleToFmpzInteger(kv, s, s_rt);
  fmpz_mat_init(A, 2*n, 2*n);
  fmpzMat_fromRTable(A, kv, A_rt);
  fmpz_mat_init(B, 2*n, m);
  fmpzMat_fromRTable(B, kv, B_rt);
  fmpz_mat_init(Res, 2*n, m);
  
  int success = specialIntCert(Res, s, A, B, n, m, r);

  Res_rt = fmpzMat_toRTable(kv, Res);

  fmpz_mat_clear(A);
  fmpz_mat_clear(B);
  fmpz_mat_clear(Res);
  fmpz_clear(s);
  MaplePopGMPAllocators(kv);

  if (!success) {
    return ToMapleBoolean(kv, 0);
  } else {
    return ToMapleExpressionSequence(kv, 1, Res_rt);
  }
}

ALGEB M_DECL indexMassager_maple(MKernelVector kv, ALGEB args)
{
  int argc;
  ALGEB B_rt, n_rt, m_rt, r_rt, eps_rt, Q_rt, U_rt, M_rt, S_rt, T_rt, k_rt, s_rt;
  fmpz_mat_t B, Q, U, M, S, T;
  int n, m, r, k;
  fmpz_t s;
  double eps;

  //printf("calling specialIntCert\n");
  argc = MapleNumArgs(kv,args);
  if (argc != 8 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleInteger(kv, (ALGEB)args[2])
                || !IsMapleInteger(kv, (ALGEB)args[3])
                || !IsMapleInteger(kv, (ALGEB)args[4])
                || !IsMapleInteger(kv, (ALGEB)args[5])
                //|| !IsMapleFloat64(kv, (ALGEB)args[6])
                //|| !(IsMapleRTable(kv, (ALGEB)args[7]) || IsMapleBoolean(kv, (ALGEB)args[7]))
                || !IsMapleInteger(kv, (ALGEB)args[8])) {
    MapleRaiseError(kv, "exactly eight arguments expected");
    return NULL;
  }

  B_rt = (ALGEB)args[1];
  n_rt = (ALGEB)args[2];
  m_rt = (ALGEB)args[3];
  r_rt = (ALGEB)args[4];
  s_rt = (ALGEB)args[5];
  eps_rt = (ALGEB)args[6];
  Q_rt = (ALGEB)args[7];
  k_rt = (ALGEB)args[8];
  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);
  n = MapleToInteger32(kv, n_rt);
  m = MapleToInteger32(kv, m_rt);
  r = MapleToInteger32(kv, r_rt);
  k = MapleToInteger32(kv, k_rt);
  //eps = MapleToFloat64(kv, eps_rt);
  if (IsMapleRTable(kv, (ALGEB)args[7])) {
    fmpz_mat_init(Q, n, r+k);
    fmpzMat_fromRTable(Q, kv, Q_rt);
  } else {
    fmpz_mat_init(Q, 0, 0);
  }
  fmpz_init(s);
  mapleToFmpzInteger(kv, s, s_rt);
  fmpz_mat_init(B, 2*n, 2*n);
  fmpzMat_fromRTable(B, kv, B_rt);
  fmpz_mat_init(U, r, n);
  fmpz_mat_init(M, n, r);
  fmpz_mat_init(T, r, r);
  fmpz_mat_init(S, r, 1);
  
  //printf("Calling indexMassage\n");
  int success = indexMassager(S, U, M, T, B, n, m, r, s, k, Q);

  U_rt = fmpzMat_toRTable(kv, U);
  M_rt = fmpzMat_toRTable(kv, M);
  S_rt = fmpzMat_toRTable(kv, S);
  T_rt = fmpzMat_toRTable(kv, T);

  fmpz_mat_clear(Q);
  fmpz_mat_clear(B);
  fmpz_mat_clear(S);
  fmpz_mat_clear(M);
  fmpz_mat_clear(U);
  fmpz_mat_clear(T);
  fmpz_clear(s);
  MaplePopGMPAllocators(kv);

  //printf("Returning from indexMassager_maple: %d \n", success);
  if (!success) {
    return ToMapleBoolean(kv, 0);
  } else {
    return ToMapleExpressionSequence(kv, 4, U_rt, M_rt, T_rt, S_rt);
  }
}

ALGEB M_DECL largestInvariantFactor_maple(MKernelVector kv, ALGEB args)
{
  int argc, n, startDimension;
  ALGEB A_rt, n_rt, startDimension_rt, X_rt, s_rt;
  fmpz_mat_t A, X;
  fmpz_t s;

  argc = MapleNumArgs(kv,args);
  if (argc != 3 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleInteger(kv, (ALGEB)args[2])
                || !IsMapleInteger(kv, (ALGEB)args[3])) {
    MapleRaiseError(kv, "exactly eight arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  n_rt = (ALGEB)args[2];
  startDimension_rt = (ALGEB)args[3];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);
  startDimension = MapleToInteger32(kv, startDimension_rt);
  n = MapleToInteger32(kv, n_rt);

  fmpz_init(s);
  fmpz_mat_init(A, n, n);
  fmpzMat_fromRTable(A, kv, A_rt);
  fmpz_mat_init(X, n, startDimension);
  
  int success = largestInvariantFactor(s, X, A, startDimension);

  X_rt = fmpzMat_toRTable(kv, X);
  s_rt = FMPZIntegerToMaple(kv, s);

  fmpz_mat_clear(A);
  fmpz_mat_clear(X);
  fmpz_clear(s);
  MaplePopGMPAllocators(kv);

  if (!success) {
    return ToMapleBoolean(kv, 0);
  } else {
    return ToMapleExpressionSequence(kv, 2, s_rt, X_rt);
  }
}

ALGEB M_DECL smithMassager_maple(MKernelVector kv, ALGEB args)
{
  int argc, n;
  ALGEB A_rt, n_rt, U_rt, M_rt, S_rt, T_rt;
  fmpz_mat_t A, U, M, S, T;

  argc = MapleNumArgs(kv,args);
  if (argc != 2 || !IsMapleRTable(kv, (ALGEB)args[1])
                || !IsMapleInteger(kv, (ALGEB)args[2])) {
    MapleRaiseError(kv, "exactly eight arguments expected");
    return NULL;
  }

  A_rt = (ALGEB)args[1];
  n_rt = (ALGEB)args[2];

  MaplePushGMPAllocators(kv, alloc_func, realloc_func, free_func);
  n = MapleToInteger32(kv, n_rt);

  fmpz_mat_init(A, n, n);
  fmpzMat_fromRTable(A, kv, A_rt);
  fmpz_mat_init(U, n, n);
  fmpz_mat_init(M, n, n);
  fmpz_mat_init(T, n, n);
  fmpz_mat_init(S, n, 1);
  
  printf("calling smithMassager\n");
  int success = smithMassager(U, M, T, S, A);

  U_rt = fmpzMat_toRTable(kv, U);
  M_rt = fmpzMat_toRTable(kv, M);
  T_rt = fmpzMat_toRTable(kv, T);
  S_rt = fmpzMat_toRTable(kv, S);

  fmpz_mat_clear(A);
  fmpz_mat_clear(U);
  fmpz_mat_clear(M);
  fmpz_mat_clear(S);
  fmpz_mat_clear(T);
  MaplePopGMPAllocators(kv);

  if (!success) {
    return ToMapleBoolean(kv, 0);
  } else {
    return ToMapleExpressionSequence(kv, 4, U_rt, M_rt, T_rt, S_rt);
  }
}
