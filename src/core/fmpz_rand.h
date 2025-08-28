#pragma once

#include "flint.h"

extern int initialized;
extern flint_rand_t rand;

void fmpz_rand_init();
void fmpz_rand_clear();
