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
/** \file energy_tcl.cpp
 *
 *  Implements the analyze energy command.
 */
#include "utils.hpp"
#include "parser.hpp"
#include "energy.hpp"

/**********************************************************************
 *                                 parser
 **********************************************************************/

static void tclcommand_analyze_print_all(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  double value;
  int i, j;

  value = total_energy.data.e[0];
  for (i = 1; i < total_energy.data.n; i++)
    value += total_energy.data.e[i];

  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, "{ energy ", buffer, " } ", (char *)NULL);

  Tcl_PrintDouble(interp, total_energy.data.e[0], buffer);
  Tcl_AppendResult(interp, "{ kinetic ", buffer, " } ", (char *)NULL);

  for(i=0;i<n_bonded_ia;i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      sprintf(buffer, "%d ", i);
      Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
      Tcl_PrintDouble(interp, *obsstat_bonded(&total_energy, i), buffer);
      Tcl_AppendResult(interp,
		       get_name_of_bonded_ia(bonded_ia_params[i].type),
		       " ", buffer, " } ", (char *) NULL);
    }
  }

  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, *obsstat_nonbonded(&total_energy, i, j), buffer);
	Tcl_AppendResult(interp, "nonbonded ", buffer, " } ", (char *)NULL);	    
      }
    }

#if defined(ELECTROSTATICS) || defined(DIPOLES)
  if(
#if defined(ELECTROSTATICS) && defined(DIPOLES) 
     coulomb.method != COULOMB_NONE || coulomb.Dmethod != DIPOLAR_NONE
#elif defined(ELECTROSTATICS)
     coulomb.method != COULOMB_NONE
#elif defined(DIPOLES)     
     coulomb.Dmethod != DIPOLAR_NONE
#endif
     ) {
    /* total Coulomb energy */
    value = 0;
    for (i = 0; i < total_energy.n_coulomb; i++)
      value += total_energy.coulomb[i];
    for (i = 0; i < total_energy.n_dipolar; i++)
      value += total_energy.dipolar[i];
    Tcl_PrintDouble(interp, value, buffer);
    
#if defined(ELECTROSTATICS) && defined(DIPOLES) 
    Tcl_AppendResult(interp, "{ coulomb+magdipoles ", buffer, (char *)NULL);  
#elif defined(ELECTROSTATICS)
    Tcl_AppendResult(interp, "{ coulomb ", buffer, (char *)NULL);
#elif defined(DIPOLES)
    Tcl_AppendResult(interp, "{ magdipoles ", buffer, (char *)NULL);  
#endif

    /* if it is split up, then print the split up parts */
#ifdef ELECTROSTATICS
    if (total_energy.n_coulomb > 1) {
      for (i = 0; i < total_energy.n_coulomb; i++) {
	Tcl_PrintDouble(interp, total_energy.coulomb[i], buffer);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
      }
    }
#endif

#ifdef DIPOLES
    if (total_energy.n_dipolar > 1) {
      for (i = 0; i < total_energy.n_dipolar; i++) {
 	Tcl_PrintDouble(interp, total_energy.dipolar[i], buffer);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
      }
    }
#endif
    Tcl_AppendResult(interp, " }", (char *)NULL);
  }
#endif
}

/************************************************************/

int tclcommand_analyze_parse_and_print_energy(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze energy [{ fene <type_num> | harmonic <type_num> | subt_lj_harm <type_num> | subt_lj_fene <type_num> | subt_lj <type_num> | lj <type1> <type2> | ljcos <type1> <type2> | ljcos2 <type1> <type2> | gb <type1> <type2> | coulomb | kinetic | total }]' */
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  int i, j;
  double value;
  value = 0.0;
  if (n_part == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  if (total_energy.init_status == 0) {
    init_energies(&total_energy);
    master_energy_calc();
  }

  if (argc == 0)
    tclcommand_analyze_print_all(interp);
  else {

    if      (ARG0_IS_S("kinetic"))
      value = total_energy.data.e[0];
    else if (ARG0_IS_S("bonded") ||
	     ARG0_IS_S("fene") ||
	     ARG0_IS_S("subt_lj_harm") ||
	     ARG0_IS_S("subt_lj_fene") ||
	     ARG0_IS_S("subt_lj") ||
	     ARG0_IS_S("harmonic") ||
	     ARG0_IS_S("endangledist")) {
      if(argc<2 || ! ARG1_IS_I(i)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze energy bonded <type_num>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_bonded_ia) {
	Tcl_AppendResult(interp, "bond type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_bonded(&total_energy, i);
    }
    else if (ARG0_IS_S("nonbonded") ||
	     ARG0_IS_S("lj") ||
	     ARG0_IS_S("buckingham") ||
	     ARG0_IS_S("lj-cos") ||
             ARG0_IS_S("lj-cos2") ||
	     ARG0_IS_S("gb") ||
	     ARG0_IS_S("tabulated")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze energy nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded(&total_energy, i, j);
    }
 
    else if( ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      value = 0;
      for (i = 0; i < total_energy.n_coulomb; i++)
	value += total_energy.coulomb[i];
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see myconfig.hpp)\n", (char *)NULL);
#endif
    }    
    else if( ARG0_IS_S("magnetic")) {
#ifdef DIPOLES
      value = 0;
      for (i = 0; i < total_energy.n_dipolar; i++)
	value += total_energy.dipolar[i];
#else
      Tcl_AppendResult(interp, "DIPOLES not compiled (see myconfig.hpp)\n", (char *)NULL);
#endif
    }
    
    else if (ARG0_IS_S("total")) {
      value = total_energy.data.e[0];
      for (i = 1; i < total_energy.data.n; i++)
	value += total_energy.data.e[i];
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of: analyze energy",
		       (char *)NULL);
      return (TCL_ERROR);
    }
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }

  return (TCL_OK);
}
