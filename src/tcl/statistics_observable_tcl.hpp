/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  
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
#ifndef STATISTICS_OBSERVABLE_TCL_H
#define STATISTICS_OBSERVABLE_TCL_H

#include "config.hpp"
#include <tcl.h>
#include "statistics_observable.hpp"

int tclcommand_observable(ClientData data, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_observable_print_formatted(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs, double* values);
int parse_id_list(Tcl_Interp* interp, int argc, char** argv, int* change, IntList** ids ); 

int observable_calc_tclcommand(observable* self);

typedef struct {
  Tcl_Interp* interp;
  int n_A;
  char* command;
} Observable_Tclcommand_Arg_Container;

#endif
