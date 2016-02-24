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

#include "twist_stack_tcl.hpp"
#include "twist_stack.hpp"

#ifdef TWIST_STACK

static void twist_stack_usage(Tcl_Interp *interp) { 
  puts("usage: twist_stack { rm epsilon a0 a1 a2 a3 a4 a5 a6 a7 b1 b2 b3 b4 b5 b6 b7 }");
}

int tclcommand_inter_parse_twist_stack(Tcl_Interp *interp, int bond_type, int argc, char **argv) {
  DoubleList params;
  
  init_doublelist(&params);

  argc--;
  argv++;

  if(!ARG0_IS_DOUBLELIST(params)) {
    twist_stack_usage(interp);
    return ES_ERROR;
  }

  if(params.n != 17) {
    puts("Wrong number of parameters");
    twist_stack_usage(interp);
    return ES_ERROR;
  }

  twist_stack_set_params(bond_type, &params);

  return ES_OK; 
}

#endif
