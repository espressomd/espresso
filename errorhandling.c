// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file errorhandling.c
    Implementation of \ref errorhandling.h "errorhandling.h".
*/
#include <mpi.h>
#include <string.h>
#include "utils.h"
#include "errorhandling.h"

/******************* exported variables **********************/
/** buffer for error messages during the integration process. NULL if no errors occured. */
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
