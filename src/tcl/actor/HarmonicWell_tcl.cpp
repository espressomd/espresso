/*
  Copyright (C) 2012,2013 The ESPResSo project
  
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

#include "actor/HarmonicWell_tcl.hpp"

#ifdef CUDA

#include "forces.hpp"
#include "actor/HarmonicWell.hpp"

int tclcommand_HarmonicWell(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  DoubleList dl;

  init_doublelist(&dl);

  if(!ARG1_IS_DOUBLELIST(dl)) {
    puts("Expected double list");
    return TCL_ERROR;
  }

  if(dl.n != 4) {
    puts("Wrong # of args");
    for(int i = 0; i < dl.n; i++)
      printf("%d %e\n", i, dl.e[i]);

    return TCL_ERROR;
  }

  // printf("x %e %e %e, k %e\n", dl.e[0], dl.e[1],dl.e[2],dl.e[3]);

  addHarmonicWell(dl.e[0], dl.e[1], dl.e[2], dl.e[3]);

  return TCL_OK;
}


#endif
