#include "fmpz_mat.h"
#include "smithMassager.h"
#include "genHerm.h"
#include "timer.h"

#include <string.h>

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
  if (argc == 2 && strcmp(argv[1], "-h") == 0) {
        printf("Help: This program processes matrices based on a given prime number n.\n");
        printf("Usage: %s <prime_number>\n", argv[0]);
        printf("Example: %s 5\n", argv[0]);
        return 0;
    }

    if (argc != 2) {
        fprintf(stderr, "Error: Please specify a prime number as the argument.\n");
        fprintf(stderr, "Usage: %s <prime_number>\n", argv[0]);
        fprintf(stderr, "Use '%s -h' for help.\n", argv[0]);
        return 1;
    }

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
  printf("U, M, T, S: \n");
  fmpz_mat_print_pretty(U);
  printf("\n");
  fmpz_mat_print_pretty(M);
  printf("\n");
  fmpz_mat_print_pretty(T);
  printf("\n");
  fmpz_mat_print_pretty(S);
  printf("\n");
  printf("success: %d\n", success);

  fmpz_mat_clear(A);
  fmpz_mat_clear(U);
  fmpz_mat_clear(M);
  fmpz_mat_clear(T);
  fmpz_mat_clear(S);
}
