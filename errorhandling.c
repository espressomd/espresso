#include <mpi.h>
#include <string.h>
#include "utils.h"
#include "errorhandling.h"

/******************* exported variables **********************/
/** buffer for error messages during the integration process. */
char *error_msg;
int n_error_msg = 0;

/******************* exported functions **********************/

char *runtime_error(int errlen)
{
  /* the true length of the string will be in general shorter than n_error_msg,
     at least if numbers are involved */
  int curend = error_msg ? strlen(error_msg) : 0;
  n_error_msg = curend + errlen + 1;
 
  error_msg = realloc(error_msg, n_error_msg);
  return error_msg + curend;
}

int check_runtime_errors()
{
  int n_all_error_msg;
  MPI_Allreduce(&n_error_msg, &n_all_error_msg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return n_all_error_msg;
}
