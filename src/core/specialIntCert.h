#pragma once

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_vec.h"

// specialIntCert(ret, s, A, B, n, m, r)
//
// Input:
//   A, 2n x 2n matrix over Z, nonsingular.
//   s, positive integer.
//   B, 2n x m matrix over Z(s).
//   n, m, r, positive integer.
//
// Pre:
//   NOTE: Empty space means zero entries.
//
//   A = [ nxn    nxr ]
//       [            ]
//       [            ]
//       [ rxn    rxr ]
//
//   B = [ n x m ]
//       [       ]
//
// Output:
//   If sA^-1B is integral then return modp(sA^-1B, s) otherwise return false.
//
int specialIntCert(fmpz_mat_t ret, fmpz_t s, fmpz_mat_t A, fmpz_mat_t B, int n, int m, int r);
