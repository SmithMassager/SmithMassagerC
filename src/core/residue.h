#pragma once

/*#ifdef OPENBLAS*/

#include "openblas_residue.h"

typedef openblasResidue_t residue_t;

#define residue_init openblas_init
#define residue_fini openblas_fini

#define residue_getEntry openblas_getEntry
#define residue_setEntry openblas_setEntry


#define residue_copy openblas_copy
#define residue_fromMpzMatrix openblas_fromMpzMatrix
#define residue_identity openblas_identity
#define residue_zero openblas_zero

#define residue_isZero openblas_isZero
#define residue_print openblas_print

#define residue_add openblas_add
#define residue_determinant openblas_determinant
#define residue_inverse openblas_inverse
#define residue_gemm openblas_gemm
#define residue_gemm_helper openblas_gemm_helper
#define residue_mods openblas_mods
#define residue_quadLift openblas_quadLift
#define residue_scale openblas_scale
#define residue_set openblas_set
#define residue_get openblas_get

/*#endif*/
