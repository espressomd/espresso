/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

#include "utils.hpp"
#include "parser.hpp"
#include "cg_dna_tcl.hpp"

int tclprint_to_result_cg_dnaIA(Tcl_Interp *interp, Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, params->p.fene.k, buffer);
  Tcl_AppendResult(interp, "FENE ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.fene.drmax, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  Tcl_PrintDouble(interp, params->p.fene.r0, buffer);
  Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
  return (TCL_OK);
}

void cg_dna_basepair_usage(Tcl_Interp *interp) { 
  puts("usage: cg_dna_basepair_usage { r0 alpha E0 kd E01 E02 sigma1 sigma2 theta01 theta02 }");
}

int tclcommand_inter_parse_cg_dna_basepair(Tcl_Interp *interp, int bond_type, int argc, char **argv) {   
  DoubleList params;
  
  init_doublelist(&params);

  argc--;
  argv++;

  if(!ARG0_IS_DOUBLELIST(params)) {
    cg_dna_basepair_usage(interp);
    return ES_ERROR;
  }

  if(params.n != 10) {
    puts("Wrong number of parameters");
    cg_dna_basepair_usage(interp);
    return ES_ERROR;
  }

  cg_dna_basepair_set_params(bond_type, &params);

  return ES_OK; 
}

int tclcommand_inter_parse_cg_dna_stacking(Tcl_Interp *interp, int bond_type, int argc, char **argv) { return 0; }
int tclcommand_inter_parse_cg_dna_backbone(Tcl_Interp *interp, int bond_type, int argc, char **argv) { return 0;}

