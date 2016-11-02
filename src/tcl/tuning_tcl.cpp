/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file tuning_tcl.cpp
 *
 *  Implements the callback for the timings global variable used during
 *  tuning of e.g. P3M or mmm1d.
 */
#include "parser.hpp"
#include "tuning.hpp"

int tclcallback_timings(Tcl_Interp *interp, void *data)
{
  if (*(int *)data <= 0)
    timing_samples = 0;
  else 
    timing_samples = *(int *)data;
  return TCL_OK;
}

int tclcommand_time_integration(ClientData data, Tcl_Interp *interp, int argc, char *argv[]) {
  char buffer[10+TCL_DOUBLE_SPACE];
  double t;
  int n = 1;
  if(argc > 2) {
    Tcl_AppendResult(interp, "time_integration expects zero or one argument.", (char *)NULL);
    return TCL_ERROR;
  }
  if(argc == 2) {
    if(!(ARG1_IS_I(n) && (n > 1))) {
      return TCL_ERROR;
    }
  }

  t = time_force_calc(n);

  sprintf(buffer, "%lf", t);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

int tclcommand_tune_skin(ClientData data, Tcl_Interp *interp, int argc, char *argv[]) {
  if(argc != 5) {
    puts("usage:");
    return TCL_ERROR;
  }

  double min, max, tol;
  int steps;
  if(!(ARG_IS_D(1, min) && ARG_IS_D(2, max) && ARG_IS_D(3, tol) && ARG_IS_I(4, steps))) {
    puts("usage:");
    return TCL_ERROR;
  }

    tune_skin(min, max, tol, steps);
  return TCL_OK;
}







