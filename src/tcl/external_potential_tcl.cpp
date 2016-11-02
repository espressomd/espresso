/*
  Copyright (C) 2014,2015,2016 The ESPResSo project
  
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


#include "external_potential_tcl.hpp"

int tclcommand_external_potential_tabulated(Tcl_Interp* interp, int argc, char **argv, ExternalPotential* e);

int tclcommand_external_potential(ClientData _data, Tcl_Interp *interp,
    int argc, char **argv) {
  ExternalPotential* e;
  int error = generate_external_potential(&e);
  if (error == ES_ERROR)
    return TCL_ERROR;

  if (argc<2) {
    Tcl_AppendResult(interp, "Usage: external_potential <tabulate|rod|...>\n" , (char *)NULL);
    return TCL_ERROR;
  }
  if (ARG1_IS_S("tabulated")) {
    return tclcommand_external_potential_tabulated(interp, argc-2, argv+2, e);
  }
  Tcl_AppendResult(interp, "Usage: external_potential <tabulate|rod|...>\n" , (char *)NULL);
  return TCL_ERROR;

}



int tclcommand_external_potential_tabulated(Tcl_Interp* interp, int argc, char **argv, ExternalPotential* e) 
{
  char* filename =0;

  DoubleList scalelist;

  init_doublelist(&scalelist);

  while (argc>0) {
    if (ARG0_IS_S("file") ) {
      if (argc>1) {
        filename = argv[1];
        argc-=2;
        argv+=2;
      } else {
        Tcl_AppendResult(interp, "Usage: external_potential file <filename>\n" , (char *)NULL);
        return TCL_ERROR;
      }
    } else if (ARG0_IS_S("scale")) {
      if (argc>1  && ARG_IS_DOUBLELIST(1, scalelist)) {
        argc-=2;
        argv+=2;
      } else {
        Tcl_AppendResult(interp, "Usage: external_potential tabulated scale <float>\n" , (char *)NULL);
        return TCL_ERROR;
      }
    } else {
      Tcl_AppendResult(interp, "Unknown argument to external_potential: " , argv[0], "\n", (char *)NULL);
      return TCL_ERROR;
    }
  }
  if (filename == 0) {
    Tcl_AppendResult(interp, "No filename given to external_potential tabulated\n" , (char *)NULL);
    return TCL_ERROR;
  }
  if (external_potential_tabulated_init(n_external_potentials-1, filename, scalelist.n, scalelist.e)==ES_ERROR) {
    Tcl_AppendResult(interp, "Error creating external potential\n" , (char *)NULL);
    return TCL_ERROR;
  }
  return gather_runtime_errors(interp, TCL_OK);
}




