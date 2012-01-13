/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
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
  MPI_Allreduce(&n_error_msg, &n_all_error_msg, 1, MPI_INT, MPI_SUM, comm_cart);
  return n_all_error_msg;
}
