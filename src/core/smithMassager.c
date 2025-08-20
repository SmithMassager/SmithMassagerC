#include <assert.h>
#include <math.h>

#include "fmpz_mat.h"
#include "fmpz.h"
#include "flint.h"

#define min(a, b) (((a) < (b)) ? (a) : (b))

int smithMassager(fmpz_mat_t U, fmpz_mat_t M, fmpz_mat_t T, fmpz_mat_t S, fmpz_mat_t A) {
  assert(A->r == A->c);
  fmpz_t s, thres, p, maxA, maxU, maxUp, tmp;
  fmpz_mat_t B;
  fmpz_mat_struct *Q;
  int i, sumR, n, startDim, k, oldsbits, success, r;

  i = sumR = r = 0;
  n = A->r;
  k = startDim = 5;
  success = 1;
  Q = malloc(sizeof(fmpz_mat_struct));

  printf("SmithMassager: Initiating variables\n");
  fmpz_mat_init(B, 2*n, 2*n);
  printf("SmithMassager: Initiating Q\n");
  fmpz_mat_init(Q, n, startDim+k);
  printf("SmithMassager: Finish initiating Q\n");
  fmpz_init(s);
  fmpz_init(thres);
  fmpz_init(p);
  fmpz_init(maxA);
  fmpz_init(maxU);
  fmpz_init(maxUp);
  fmpz_init(tmp);

  fmpz_mat_cpy(B, 0, 0, n, n, A, 0, 0, n, n);
  printf("SmithMassager: Finish initiating variables\n");

  printf("SmithMassager: Calling largestInv\n");
  largestInvariantFactor(s, Q, A, startDim+k);
  printf("SmithMassager: Finish calling largestInv\n");

  while (sumR < n) {
    fmpz_mat_t Up, Mp, Tp, Sp, SpInv, MpSpInv, Update, TpSpInv, UpMp, NegUpM, Msub, Usub;
    if (i == 0) {
      r = min(startDim, n-sumR);
    } else {
      int sbits = fmpz_bits(s);
      r = min((oldsbits - sbits)*i, n-sumR);
    }
    oldsbits = fmpz_bits(s);

    printf("SmithMassager: iteration %d\n", i);
    printf("SmithMassager: Inner Initiating n, r: %d %d\n", n, r);
    printf("%d %d %d %d\n", B->r - sumR - r, B->c - sumR - r, B->r - sumR, B->c - sumR);
    fmpz_mat_init(Up, r, n);
    printf(".1\n");
    fmpz_mat_init(Mp, n, r);
    printf(".2\n");
    fmpz_mat_init(Tp, r, r);
    printf(".3\n");
    fmpz_mat_init(Sp, r, 1);
    printf(".4\n");
    fmpz_mat_init(SpInv, r, 1);
    printf(".5\n");
    fmpz_mat_init(MpSpInv, n, r);
    printf(".6\n");
    fmpz_mat_init(Update, n+sumR+r, n);
    printf(".7\n");
    fmpz_mat_init(TpSpInv, r, r);
    printf(".8\n");
    fmpz_mat_window_init(UpMp, B, B->r - sumR - r, B->c - sumR - r, B->r - sumR, B->c - sumR);
    printf(".9\n");
    fmpz_mat_init(NegUpM, r, sumR);
    printf(".10\n");
    fmpz_mat_window_init(Msub, M, 0, M->c - sumR, M->r, M->c);
    printf(".11\n");
    fmpz_mat_window_init(Usub, U, U->c - sumR, 0, U->r, U->c);
    printf("SmithMassager: Inner Finish Initiating\n");

    printf("SmithMassager: Calling index massager\n");
    success = indexMassager(Sp, Up, Mp, Tp, B, n, sumR, r, s, 0, Q, k);
    printf("SmithMassager: Finish calling index massager and returned %d \n", success);

    if (i == 0) {
      fmpz_mat_clear(Q);
      free(Q);
      Q = NULL;
    }

    if (!success) goto cleanInner;

    // Apply the (m,r)-index massager on B and update (U, M, T, S) to become
    // (0, r+m)-index massager for diag(A, I_n).
    // Pre:
    // B looks like:
    // [A       A M S^-1   ]
    // [   I               ]
    // [U    (T+ U M) S^-1 ]
    //
    // i.e. B has a (0, sumR)-index smith massager (U, M, T, S).
    //
    // Post:
    // B looks like:
    // [A         A Mp Sp^-1               A M S^-1        ]
    // [   I                                               ]
    // [Up    (Up Mp + Tp) Sp^-1                           ]
    // [U          U Mp Sp^-1          (T + U M)S^-1       ]

    // claculate, thres = max(abs(A), abs(U), abs(Up)) * n * 2 + 1
    printf("SmithMassager: Applying index massager\n");
    printf("SmithMassager: Finding max A, U, Up\n");
    fmpz_set_ui(thres, 0);
    printf("SmithMassager: Calling mat_max\n");
    fmpz_mat_max(maxUp, Up);
    printf("SmithMassager: Finish calling mat_max\n");
    fmpz_max(thres, thres, maxA);
    fmpz_max(maxU, maxU, maxUp);
    fmpz_max(thres, thres, maxU);
    printf("SmithMassager: Finish Finding max A, U, Up\n");

    fmpz_set_ui(p, 1709);
    while (1) { 
      fmpz_gcd(tmp, p, fmpz_mat_entry(Sp, r-1, 0));
      if (fmpz_cmp_ui(tmp, 1) == 0) {
	break;
      }
      fmpz_nextprime(p, p, 1);
    }
    while (fmpz_cmp(p, thres) < 0) { fmpz_mul(p, p, p); }
    printf("SmithMassager: Finish finding p\n");
    
    printf("SmithMassager: 1\n");
    printf("%d %d %d\n", B->r, sumR, r);
    printf("%d %d %d %d %d %d %d %d\n", B->r - sumR - r, 0, B->r - sumR, n, 0, 0, Up->r, Up->c);
    fmpz_mat_cpy(B, B->r - sumR - r, 0, B->r - sumR, n, Up, 0, 0, Up->r, Up->c);
    printf("SmithMassager: 2\n");
    fmpz_mat_modDiagInv(SpInv, p, Sp);
    printf("SmithMassager: 3\n");
    fmpz_mat_modDiagMul(MpSpInv, p, Mp, SpInv);
    printf("SmithMassager: 4\n");
    printf("n, sumR, r: %d %d %d\n", n, sumR, r);
    fmpz_mat_concat_vertical3(Update, A, Up, Usub);
    printf("SmithMassager: 5\n");
    fmpz_mat_mul(Update, Update, MpSpInv);
    printf("SmithMassager: 6\n");
    fmpz_mat_modDiagMul(TpSpInv, p, Tp, SpInv);
    printf("SmithMassager: 7\n");
    fmpz_mat_mods(Update, Update, p);
    printf("SmithMassager: 8\n");
    fmpz_mat_cpy(B, 0, B->c - sumR - r, n, B->c - sumR, Update, 0, 0, n, r);
    printf("SmithMassager: 9\n");
    fmpz_mat_cpy(B, B->r - sumR, B->c - sumR - r, B->r, B->c - sumR, Update, n + r, 0, Update->r, r);
    printf("SmithMassager: 10\n");
    fmpz_mat_cpy(UpMp, 0, 0, UpMp->r, UpMp->c, Update, n, 0, n+r, Update->c);
    printf("SmithMassager: 11\n");
    fmpz_mat_add(UpMp, UpMp, TpSpInv);
    printf("SmithMassager: 12\n");
    fmpz_mat_mods(UpMp, UpMp, p);
    printf("SmithMassager: 13\n");
    fmpz_mat_cpy(B, B->r - sumR - r, B->c - sumR - r, B->r - sumR, B->c - sumR, UpMp, 0, 0, UpMp->r, UpMp->c);

    printf("SmithMassager: Updating U, M, S, T\n");
    fmpz_mat_cpy(U, U->r - sumR - r, 0, U->r - sumR, U->c, Up, 0, 0, Up->r, Up->c);
    fmpz_mat_cpy(M, 0, M->c - sumR - r, M->r, M->c - sumR, Mp, 0, 0, Mp->r, Mp->c);
    fmpz_mat_cpy(S, S->r - sumR - r, 0, S->r - sumR, 1, Sp, 0, 0, Sp->r, Sp->c);
    fmpz_mat_cpy(T, T->r - sumR - r, T->r - sumR - r, T->r - sumR, T->r - sumR, Tp, 0, 0, Tp->r, Tp->c);
    if (sumR != 0) {
      fmpz_mat_mul(NegUpM, Up, Msub);
      fmpz_mat_neg(NegUpM, NegUpM);
      fmpz_mat_cpy(T, T->r - sumR - r, T->c - sumR , T->r - sumR, T->c, NegUpM, 0, 0, r, sumR);
    }
    printf("SmithMassager: Finish Applying index massager\n");

    fmpz_set(s, fmpz_mat_entry(Sp, 0, 0));
    sumR += r;
    ++i;
cleanInner:
    fmpz_mat_clear(Up);
    fmpz_mat_clear(Mp);
    fmpz_mat_clear(Tp);
    fmpz_mat_clear(Sp);
    fmpz_mat_clear(SpInv);
    fmpz_mat_clear(MpSpInv);
    fmpz_mat_clear(Update);
    fmpz_mat_clear(TpSpInv);
    fmpz_mat_window_clear(UpMp);
    fmpz_mat_clear(NegUpM);
    fmpz_mat_window_clear(Msub);
    fmpz_mat_window_clear(Usub);

    if (!success) goto cleanOuter;
    printf("SmithMassager: Finish cleaning going into next loop\n");
  }

  printf("SmithMassager: Building Unicert\n");
// Now need todo unicertificate.
  int l, m;
  l = -1, r = n-1;
  while (l < r) {
    m = ceil((double)(l+r)/2);
    if (fmpz_cmp_ui(fmpz_mat_entry(S, m, 0), 1) > 0) {
      r = m-1;
    } else {
      l = m;
    }
  }

  if (l != -1) {
    fmpz_mat_t E;
    fmpz_mat_init(E, 2*n-l, 2*n - l);
    fmpz_mat_cpy(E, 0, 0, n, n, B, 0, 0, n, n);
    fmpz_mat_cpy(E, n, 0, 2*n-l, n, B, n+l, 0, 2*n, n);
    fmpz_mat_cpy(E, 0, n, n, 2*n-l, B, 0, n+l, n, 2*n);
    fmpz_mat_cpy(E, n, n, 2*n-l, 2*n-l, B, n+l, n+l, 2*n, 2*n);

    success = fmpz_uniCert(E);

    fmpz_mat_clear(E);
  } else {
    success = fmpz_uniCert(B);
  }

cleanOuter:
  fmpz_mat_clear(B);
  if (Q) { fmpz_mat_clear(Q); free(Q); }
  fmpz_clear(s);
  fmpz_clear(thres);
  fmpz_clear(p);
  fmpz_clear(maxA);
  fmpz_clear(maxU);
  fmpz_clear(maxUp);
  fmpz_clear(tmp);

  return success;
}
