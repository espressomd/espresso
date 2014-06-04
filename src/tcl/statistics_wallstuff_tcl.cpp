/*
  Copyright (C) 2010,2011,2013 The ESPResSo project
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
/** \file statistics_wallstuff_tcl.cpp
    This is the place for analysis (so far...).
*/
#include "statistics_wallstuff_tcl.hpp"

#include "parser.hpp"
#include "statistics_wallstuff.hpp"
#include "statistics.hpp"

int tclcommand_analyze_wallstuff(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze wallmsd -xy|-z <min> <max>' */
  /******************************************************************************/
  char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 2];
  DoubleList g;
  int job, bin;
  double rmin, rmax,rclocal;
  int rbins, boxes;
  enum { BINS, MX, MYZ, RDFYZ,BONDYZ,SCALE,SCALE2, PRINT };

  if (argc < 2) {
    Tcl_AppendResult(interp, "expected: analyze wallstuff -bins <binboundaries> | -myz <bin> |-mx <bin> | -rdfyz <bin> <rmin> <rmax> <rdfbins> | -bondyz | -scale | -scale2 | -print",
		     (char *)NULL);
    return TCL_ERROR;
  }

  // 1. what do we do?
  if (ARG0_IS_S("-bins")) {
    job = BINS;
  }
  else if (ARG0_IS_S("-mx") && argc == 2) {
    job = MX;
  }
  else if (ARG0_IS_S("-myz") && argc == 2) {
    job = MYZ;
  }
  else if (ARG0_IS_S("-rdfyz") && argc == 5) {
    job = RDFYZ;
  }
  else if (ARG0_IS_S("-bondyz") && argc == 5) {
    job = BONDYZ;
  }
  else if (ARG0_IS_S("-scale") && argc == 4) {
    job = SCALE;
  }
  else if (ARG0_IS_S("-scale2") && argc == 4) {
    job = SCALE2;
  }
  else if (ARG0_IS_S("-print") && argc == 2) {
    job = PRINT;
  }
  else {
    Tcl_AppendResult(interp, ": analyze wallstuff -bins|-myz|-mx|-rdfyz|-bondyz|-scale|-scale2...", (char *)NULL);
    return TCL_ERROR;
  }
  
  // 2. parameters
  // 2. a) 1. parameter, bin or boundaries
  switch (job) {
  case BINS:
    realloc_doublelist(&wallstuff_boundaries, wallstuff_boundaries.n = 0);
    if (!ARG_IS_DOUBLELIST(1, wallstuff_boundaries)) {
      return TCL_ERROR;
    }
    if (wallstuff_boundaries.n < 2) {
      return (TCL_ERROR);
    }
    break;
  case MX:
  case MYZ:
  case RDFYZ:
  case BONDYZ:
  case SCALE:
  case SCALE2:
  case PRINT:
    if (!ARG_IS_I(1, bin)) {
      return (TCL_ERROR);
    }
    if (bin < 0 || bin >= wallstuff_boundaries.n-1) {
      return (TCL_ERROR);
    }
    break;
  }

  // 2. other parameters, only for rdf
  switch (job) {
  case RDFYZ:
    if (!ARG_IS_D(2, rmin)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_D(3, rmax)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_I(4, rbins)) {
      return (TCL_ERROR);
    }
    break;
  case BONDYZ:
    if (!ARG_IS_D(2, rclocal)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_D(3, rmax)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_I(4, rbins)) {
      return (TCL_ERROR);
    }
    break;
  case SCALE:
    if (!ARG_IS_I(2, boxes)) {
      return (TCL_ERROR);
    }
  case SCALE2:
    if (!ARG_IS_I(2, boxes)) {
      return (TCL_ERROR);
    }
    if (!ARG_IS_D(3, rclocal)) {
      return (TCL_ERROR);
    }
    break;
  
  }

  // result double list
  init_doublelist(&g);

  // check that data is there
  switch (job) {
  case BINS:
  case RDFYZ:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  case BONDYZ:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  case SCALE:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  case SCALE2:
    // these cases use partCfg
    updatePartCfg(WITHOUT_BONDS);
    break;
  case MX:
  case MYZ:
    // these cases use the positions array
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! Use 'analyze append' to save some!",
		       (char *)NULL);
      return TCL_ERROR;
    }
    break;
  }

  // finally, do what is necessary
  switch (job) {
  case BINS:
    wall_sort_particles();
    break;
  case MX:
    realloc_doublelist(&g, g.n = n_configs);
    calc_wallmsdx(g.e, bin);    
    break;
  case MYZ:
    realloc_doublelist(&g, g.n = n_configs);
    calc_wallmsdyz(g.e, bin);    
    break;
  case RDFYZ:
    realloc_doublelist(&g, g.n = rbins);
    calc_wallrdfyz(g.e, bin, rmin, rmax, rbins);
    break;
  case BONDYZ:
    realloc_doublelist(&g, g.n = rbins+2);
    calc_wallbondyz(g.e, bin, rclocal, rmax, rbins);
    break;
  case SCALE:
    realloc_doublelist(&g, g.n = 3*pow(4,boxes));
    calc_scaling (g.e,bin, boxes, rclocal);
    break;
  case SCALE2:
    realloc_doublelist(&g, g.n = 14*pow(4,boxes));
    calc_scaling2 (g.e,bin, boxes, rclocal);
    break;
  case PRINT:
    // just write out what wall_sort_particles has put into
    // this bin
    for (int i = 1; i < wallstuff_part_in_bin[bin].n; i++) { 
      sprintf(buffer," %d",wallstuff_part_in_bin[bin].e[i]);
      Tcl_AppendResult(interp, buffer, (char *)NULL); 
    }
    break;
  }

  // print out double results, if any
  if (g.n) {
    sprintf(buffer,"%f",g.e[0]);
    Tcl_AppendResult(interp, buffer, (char *)NULL); 

    for (int i = 1; i < g.n; i++) { 
      sprintf(buffer," %f",g.e[i]);
      Tcl_AppendResult(interp, buffer, (char *)NULL); 
    }
    realloc_doublelist(&g, g.n = 0);
  }

  return (TCL_OK);
}
