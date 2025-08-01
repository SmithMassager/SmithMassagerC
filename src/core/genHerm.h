#pragma once

#include "fmpz_mat.h"

// Assume A is p x p where p is a prime and generates A such that it has intersting HNF
void interestHerm(fmpz_mat_t A);
