#include "fmpz_mat.h"
#include "smithMassager.h"
#include "genHerm.h"

#include "timer.h"

int my_getnbr(char *str)
{
  int result;
  int puiss;

  result = 0;
  puiss = 1;
  while (('-' == (*str)) || ((*str) == '+'))
  {
      if (*str == '-')
        puiss = puiss * -1;
      str++;
  }
  while ((*str >= '0') && (*str <= '9'))
  {
      result = (result * 10) + ((*str) - '0');
      str++;
  }
  return (result * puiss);
}

int main(int argc, char *argv[]) {
  fmpz_mat_t A, U, M, S, T;
  int n = my_getnbr(argv[1]);
  fmpz_mat_init(A, n, n);
  fmpz_mat_init(U, n, n);
  fmpz_mat_init(M, n, n);
  fmpz_mat_init(T, n, n);
  fmpz_mat_init(S, n, 1);
  interestHerm(A);

  int success;
  REAL_TIMER("smithMassager", success = smithMassager(U, M, T, S, A));
  printf("success: %d\n", success);

  fmpz_mat_clear(A);
  fmpz_mat_clear(U);
  fmpz_mat_clear(M);
  fmpz_mat_clear(T);
  fmpz_mat_clear(S);
}
