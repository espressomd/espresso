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
/** \file pressure.c
    Implementation of \ref pressure.h "pressure.h".
*/
#include "pressure.hpp"
#include "parser.hpp"

/************************************************************/
/* callbacks for setmd                                      */
/************************************************************/

int tclcallback_npt_piston(Tcl_Interp *interp, void *_data) {
  double data = *(double *)_data;
  if (data < 0.0) { Tcl_AppendResult(interp, "the piston's mass must be positive.", (char *) NULL); return (TCL_ERROR); }
  nptiso.piston = data;
  mpi_bcast_parameter(FIELD_NPTISO_PISTON);
  return (TCL_OK);
}

int tclcallback_p_ext(Tcl_Interp *interp, void *_data) {
  double data = *(double *)_data;
  nptiso.p_ext = data;
  mpi_bcast_parameter(FIELD_NPTISO_PEXT);
  return (TCL_OK);
}

int tclcallback_npt_p_diff(Tcl_Interp *interp, void *_data) {
  double data = *(double *)_data;
  nptiso.p_diff = data;
  mpi_bcast_parameter(FIELD_NPTISO_PDIFF);
  return (TCL_OK);
}


/****************************************************************************************
 *                                 parser
 ****************************************************************************************/

static void tclcommand_analyze_print_pressure_all(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  double value;
  int i, j;
  value = 0.0;

  value = total_pressure.data.e[0];
  for (i = 1; i < total_pressure.data.n; i++)
    value += total_pressure.data.e[i];

  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, "{ pressure ", buffer, " } ", (char *)NULL);

  Tcl_PrintDouble(interp, total_pressure.data.e[0], buffer);
  Tcl_AppendResult(interp, "{ ideal ", buffer, " } ", (char *)NULL);

  for(i=0;i<n_bonded_ia;i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      sprintf(buffer, "%d ", i);
      Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
      Tcl_PrintDouble(interp, *obsstat_bonded(&total_pressure, i), buffer);
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
	Tcl_PrintDouble(interp, *obsstat_nonbonded(&total_pressure, i, j), buffer);
	Tcl_AppendResult(interp, "nonbonded ", buffer, " } ", (char *)NULL);
      }
    }
  
/* In case we need intra- and inter- nonbonded (nb) contribution of total pressure  */
  value = 0.0;
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      value += *obsstat_nonbonded_intra(&total_pressure_non_bonded, i, j);
    }
  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, " { total_nb_intra ", (char *)NULL);
  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_AppendResult(interp, "} ", (char *)NULL);
 
  value = 0.0;
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      value += *obsstat_nonbonded_inter(&total_pressure_non_bonded, i, j);
    }
  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, " { total_nb_inter ", (char *)NULL);
  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_AppendResult(interp, "} ", (char *)NULL);
  
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, *obsstat_nonbonded_intra(&total_pressure_non_bonded, i, j), buffer);
	Tcl_AppendResult(interp, "nb_intra ", buffer, " } ", (char *)NULL);
      }
    }
  
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, *obsstat_nonbonded_inter(&total_pressure_non_bonded, i, j), buffer);
	Tcl_AppendResult(interp, "nb_inter ", buffer, " } ", (char *)NULL);
      }
    } 
  
#if  defined(ELECTROSTATICS) || defined (DIPOLES)
  if( 
#ifdef ELECTROSTATICS
      coulomb.method != COULOMB_NONE
#else
      0
#endif
      ||
#ifdef DIPOLES
      coulomb.Dmethod != DIPOLAR_NONE
#else
      0
#endif     
   ) {
    /* total Coulomb pressure */
    value = total_pressure.coulomb[0];
    for (i = 1; i < total_pressure.n_coulomb; i++)
      value += total_pressure.coulomb[i];
    for (i = 0; i < total_pressure.n_dipolar; i++)
      value += total_pressure.dipolar[i];
    Tcl_PrintDouble(interp, value, buffer);
#if  defined(ELECTROSTATICS) && defined (DIPOLES)
    Tcl_AppendResult(interp, "{ coulomb+magdipoles ", buffer, (char *)NULL);
#else
#if defined(ELECTROSTATICS)
    Tcl_AppendResult(interp, "{ coulomb ", buffer, (char *)NULL);
#endif	
#if defined(DIPOLES)
    Tcl_AppendResult(interp, "{ magdipoles ", buffer, (char *)NULL);
#endif		
#endif

    /* if it is split up, then print the split up parts */
    if (total_pressure.n_coulomb > 1) {
      for (i = 0; i < total_pressure.n_coulomb; i++) {
	Tcl_PrintDouble(interp, total_pressure.coulomb[i], buffer);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
      }
    }
    Tcl_AppendResult(interp, " } ", (char *)NULL);
  }
#endif
#ifdef VIRTUAL_SITES_RELATIVE
    buffer[0] = 0;
    Tcl_AppendResult(interp, "{ vs_relative ", buffer, (char *)NULL);
    Tcl_PrintDouble(interp, total_pressure.vs_relative[0], buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    Tcl_AppendResult(interp, " }", (char *)NULL);
#endif
}

/************************************************************/

int tclcommand_analyze_parse_and_print_pressure(Tcl_Interp *interp, int v_comp, int argc, char **argv)
{
  /* 'analyze pressure [{ bond <type_num> | nonbonded <type1> <type2> | coulomb | ideal | total }]' */
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  int i, j;
  double value, p_vel[3];
  value = 0.0;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  /* if desired (v_comp==1) replace ideal component with instantaneous one */
  if (total_pressure.init_status != 1+v_comp ) {
    init_virials(&total_pressure);
    init_p_tensor(&total_p_tensor);
    
    init_virials_non_bonded(&total_pressure_non_bonded);
    init_p_tensor_non_bonded(&total_p_tensor_non_bonded);
    
    if(v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) && !(nptiso.invalidate_p_vel)) {
      if (total_pressure.init_status == 0)
	master_pressure_calc(0);
      total_pressure.data.e[0] = 0.0;
      MPI_Reduce(nptiso.p_vel, p_vel, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      for(i=0; i<3; i++)
	if(nptiso.geometry & nptiso.nptgeom_dir[i])
	  total_pressure.data.e[0] += p_vel[i];
      total_pressure.data.e[0] /= (nptiso.dimension*nptiso.volume);
      total_pressure.init_status = 1+v_comp;   }
    else
      master_pressure_calc(v_comp);
  }

  if (argc == 0)
    tclcommand_analyze_print_pressure_all(interp);
  else {

    if      (ARG0_IS_S("ideal"))
      value = total_pressure.data.e[0];
    else if (ARG0_IS_S("bonded") ||
	     ARG0_IS_S("fene") ||
	     ARG0_IS_S("subt_lj_harm") ||
	     ARG0_IS_S("subt_lj_fene") ||
	     ARG0_IS_S("subt_lj") ||
	     ARG0_IS_S("harmonic") ||
	     ARG0_IS_S("endangledist")) {
      if(argc<2 || ! ARG1_IS_I(i)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure bonded <type_num>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_bonded_ia) {
	Tcl_AppendResult(interp, "bond type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_bonded(&total_pressure, i);
    }
    else if (ARG0_IS_S("nonbonded") ||
	     ARG0_IS_S("lj") ||
	     ARG0_IS_S("buckingham") ||
             ARG0_IS_S("morse") ||
             ARG0_IS_S("soft-sphere") ||
	     ARG0_IS_S("lj-cos") ||
	     ARG0_IS_S("lj-cos2") ||
	     ARG0_IS_S("tabulated") ||
	     ARG0_IS_S("gb")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded(&total_pressure, i, j);
    }
    else if( ARG0_IS_S("tot_nb_intra") ||
	     ARG0_IS_S("tot_nonbonded_intra")) {
      value = 0.0;
      for (i = 0; i < n_particle_types; i++)
        for (j = i; j < n_particle_types; j++)
        value += *obsstat_nonbonded_intra(&total_pressure_non_bonded, i, j);
    }
    else if( ARG0_IS_S("tot_nb_inter") ||
	     ARG0_IS_S("tot_nonbonded_inter")) {
      value = 0.0;
      for (i = 0; i < n_particle_types; i++)
        for (j = i; j < n_particle_types; j++)
        value += *obsstat_nonbonded_inter(&total_pressure_non_bonded, i, j);
    }
    else if( ARG0_IS_S("nb_intra") ||
	     ARG0_IS_S("nonbonded_intra")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nb_intra <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded_intra(&total_pressure_non_bonded, i, j);
    }   
    else if( ARG0_IS_S("nb_inter") ||
	     ARG0_IS_S("nonbonded_inter")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nb_inter <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded_inter(&total_pressure_non_bonded, i, j);
    }
    else if( ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      value = 0;
      for (i = 0; i < total_pressure.n_coulomb; i++)
	value += total_pressure.coulomb[i];
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if( ARG0_IS_S("dipolar")) {
#ifdef DIPOLES
      value = 0;
      for (i = total_pressure.n_coulomb-1; i < total_pressure.n_coulomb; i++)  /*when DLC will be installed this has to be changed */
        value += total_pressure.coulomb[i];
#else
      Tcl_AppendResult(interp, "DIPOLES not compiled (see config.h)\n", (char *)NULL);
#endif
    }
#ifdef VIRTUAL_SITES_RELATIVE
    else if (ARG0_IS_S("vs_relative")) {
      value =total_pressure.vs_relative[0];
    }
#endif
    else if (ARG0_IS_S("total")) {
      value = total_pressure.data.e[0];
      for (i = 1; i < total_pressure.data.n; i++) {
	value += total_pressure.data.e[i];
      }
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of: analyze pressure",
		       (char *)NULL);
      return (TCL_ERROR);
    }
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }

  return (TCL_OK);
}

int tclcommand_analyze_parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "analyze p_IK no longer exists", (char *)NULL);
  return (TCL_ERROR);
}
 
int tclcommand_analyze_parse_bins(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "analyze bins no longer exists", (char *)NULL);
  return (TCL_ERROR);
}

/****************************************************************************************
 *                                 parser
 ****************************************************************************************/
/* Parser routines used in the "analyze stress_tensor" command */ 
static void tclcommand_analyze_print_stress_tensor_all(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  double value, tvalue[9];
  int i, j, k;
  value = 0.0;

  Tcl_AppendResult(interp, "{ pressure ", (char *)NULL);
  for(j=0; j<9; j++) {
    value = total_p_tensor.data.e[j];
    for (i = 1; i < total_p_tensor.data.n/9; i++) value += total_p_tensor.data.e[9*i + j];
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);

  Tcl_AppendResult(interp, "{ ideal ", (char *)NULL);
  for(j=0; j<9; j++) {
    Tcl_PrintDouble(interp, total_p_tensor.data.e[j], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);

  for(i=0;i<n_bonded_ia;i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      sprintf(buffer, "%d ", i);
      Tcl_AppendResult(interp, "{ ", buffer, get_name_of_bonded_ia(bonded_ia_params[i].type)," ", (char *)NULL);
      for(j=0; j<9; j++) {
	Tcl_PrintDouble(interp, obsstat_bonded(&total_p_tensor, i)[j], buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      }
      Tcl_AppendResult(interp, "} ", (char *)NULL);
    }
  }

  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, "nonbonded ", (char *)NULL);
	for(k=0; k<9; k++) {
	  Tcl_PrintDouble(interp, obsstat_nonbonded(&total_p_tensor, i, j)[k], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, "} ", (char *)NULL);
      }
    }

  Tcl_AppendResult(interp, " { total_nb_intra ", (char *)NULL);
  for(k=0; k<9; k++) {
    tvalue[k] = 0.0;
    for (i = 0; i < n_particle_types; i++)
      for (j = i; j < n_particle_types; j++) {
        tvalue[k] += obsstat_nonbonded_intra(&total_p_tensor_non_bonded, i, j)[k];
      }
    Tcl_PrintDouble(interp, tvalue[k], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
  
  Tcl_AppendResult(interp, " { total_nb_inter ", (char *)NULL);
  for(k=0; k<9; k++) {
    tvalue[k] = 0.0;
    for (i = 0; i < n_particle_types; i++)
      for (j = i; j < n_particle_types; j++) {
        tvalue[k] += obsstat_nonbonded_inter(&total_p_tensor_non_bonded, i, j)[k];
      }
    Tcl_PrintDouble(interp, tvalue[k], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
 
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, "nb_intra_tensor ", (char *)NULL);
	for(k=0; k<9; k++) {
	  Tcl_PrintDouble(interp, obsstat_nonbonded_intra(&total_p_tensor_non_bonded, i, j)[k], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, "} ", (char *)NULL);
      }
    }
  
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, "nb_inter_tensor ", (char *)NULL);
	for(k=0; k<9; k++) {
	  Tcl_PrintDouble(interp, obsstat_nonbonded_inter(&total_p_tensor_non_bonded, i, j)[k], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, "} ", (char *)NULL);
      }
    }

#ifdef ELECTROSTATICS
  if(coulomb.method != COULOMB_NONE) {
    Tcl_AppendResult(interp, "{ coulomb ", (char *)NULL);
    for(j=0; j<9; j++) {
      sprintf(buffer, " %lf ", total_p_tensor.coulomb[j]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    Tcl_AppendResult(interp, "} ", (char *)NULL);
  }
#endif

#ifdef DIPOLES
  if(coulomb.Dmethod != DIPOLAR_NONE) {
    fprintf(stderr,"tensor magnetostatics, something should go here, file pressure.c ... \n");
  }  
#endif

#ifdef VIRTUAL_SITES_RELATIVE
    buffer[0] = 0;
    Tcl_AppendResult(interp, "{ vs_relative ", buffer, (char *)NULL);
    for (j=0;j<9;j++) {
      sprintf(buffer, "%g", total_p_tensor.vs_relative[j]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " ", (char *)NULL);
    }
    Tcl_AppendResult(interp, " }", (char *)NULL);
#endif

  
}

/************************************************************/
int tclcommand_analyze_parse_and_print_stress_tensor(Tcl_Interp *interp, int v_comp, int argc, char **argv)
{
  /* 'analyze stress_tensor [{ bond <type_num> | nonbonded <type1> <type2> | coulomb | ideal | total }]' */
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  int i, j, k;
  double p_vel[3], tvalue[9];
  for(j=0; j<9; j++)  tvalue[j] = 0.0;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  /* if desired (v_comp==1) replace ideal component with instantaneous one */
   if (total_pressure.init_status != 1+v_comp ) {
    init_virials(&total_pressure);
    init_p_tensor(&total_p_tensor);

    init_virials_non_bonded(&total_pressure_non_bonded);
    init_p_tensor_non_bonded(&total_p_tensor_non_bonded);

    if(v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) && !(nptiso.invalidate_p_vel)) {
      if (total_pressure.init_status == 0)
	master_pressure_calc(0);
      p_tensor.data.e[0] = 0.0;
      MPI_Reduce(nptiso.p_vel, p_vel, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      for(i=0; i<3; i++)
	if(nptiso.geometry & nptiso.nptgeom_dir[i])
	  p_tensor.data.e[0] += p_vel[i];
      p_tensor.data.e[0] /= (nptiso.dimension*nptiso.volume);
      total_pressure.init_status = 1+v_comp;   }
    else
      master_pressure_calc(v_comp);
  }

  if (argc == 0)
    tclcommand_analyze_print_stress_tensor_all(interp);
  else {

    if      (ARG0_IS_S("ideal")) {
      for(j=0; j<9; j++)  tvalue[j] = total_p_tensor.data.e[j];
    }
    else if (ARG0_IS_S("bonded") ||
	     ARG0_IS_S("fene") ||
	     ARG0_IS_S("subt_lj_harm") ||
	     ARG0_IS_S("subt_lj_fene") ||
	     ARG0_IS_S("subt_lj") ||
	     ARG0_IS_S("harmonic") ||
	     ARG0_IS_S("endangledist")) {
      if(argc<2 || ! ARG1_IS_I(i)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure bonded <type_num>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_bonded_ia) {
	Tcl_AppendResult(interp, "bond type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      for(k=0; k<9; k++) tvalue[k] = obsstat_bonded(&total_p_tensor, i)[k];
    }
    else if (ARG0_IS_S("nonbonded") ||
	     ARG0_IS_S("lj") ||
	     ARG0_IS_S("buckingham") ||
             ARG0_IS_S("morse") ||
             ARG0_IS_S("soft-sphere") ||
	     ARG0_IS_S("lj-cos") ||
	     ARG0_IS_S("lj-cos2") ||
	     ARG0_IS_S("tabulated") ||
	     ARG0_IS_S("gb")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      for(k=0; k<9; k++) tvalue[k] = obsstat_nonbonded(&total_p_tensor, i, j)[k];
    }
    else if( ARG0_IS_S("tot_nb_intra")) {
      for(k=0; k<9; k++) {
        for (i = 0; i < n_particle_types; i++)
          for (j = i; j < n_particle_types; j++) {
            tvalue[k] += obsstat_nonbonded_intra(&total_p_tensor_non_bonded, i, j)[k];
          }
      }
    }
    else if( ARG0_IS_S("tot_nb_inter")) {
      for(k=0; k<9; k++) {
        for (i = 0; i < n_particle_types; i++)
          for (j = i; j < n_particle_types; j++) {
            tvalue[k] += obsstat_nonbonded_inter(&total_p_tensor_non_bonded, i, j)[k];
          }
      }
    }
    else if( ARG0_IS_S("nb_intra")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze stress tensor nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      for(k=0; k<9; k++) tvalue[k] = obsstat_nonbonded_intra(&total_p_tensor_non_bonded, i, j)[k];
    }
    else if( ARG0_IS_S("nb_inter")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze stress tensor nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      for(k=0; k<9; k++) tvalue[k] = obsstat_nonbonded_inter(&total_p_tensor_non_bonded, i, j)[k];
    }
    else if( ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      for(j=0; j<9; j++) tvalue[j] = total_p_tensor.coulomb[j];
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if( ARG0_IS_S("dipolar")) {
#ifdef DIPOLES
      /* for(j=0; j<9; j++) tvalue[j] = total_p_tensor.coulomb[j];*/
      fprintf(stderr," stress tensor, magnetostatics, something should go here, file pressure.c ");
#else
      Tcl_AppendResult(interp, "DIPOLES not compiled (see config.h)\n", (char *)NULL);
#endif
    }
#ifdef VIRTUAL_SITES_RELATIVE
    else if (ARG0_IS_S("VS_RELATIVE")) {
      for(j=0; j<9; j++) tvalue[j] = total_p_tensor.vs_relative[j];
    }
#endif
    else if (ARG0_IS_S("total")) {
      for(j=0; j<9; j++) {
        tvalue[j] = total_p_tensor.data.e[j];
        for (i = 1; i < total_p_tensor.data.n/9; i++) tvalue[j] += total_p_tensor.data.e[9*i + j];
     }
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of: analyze pressure",
		       (char *)NULL);
      return (TCL_ERROR);
    }

    Tcl_AppendResult(interp, *argv, (char *)NULL);
    Tcl_AppendResult(interp, " ", (char *)NULL);
    for(j=0; j<9; j++) {
      Tcl_PrintDouble(interp, tvalue[j], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    }
  }

  return (TCL_OK);
}

int tclcommand_analyze_parse_local_stress_tensor(Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[TCL_DOUBLE_SPACE];
  int periodic[3];
  double range_start[3];
  double range[3];
  int bins[3];
  int i,j,k,l;
  DoubleList *TensorInBin;
  PTENSOR_TRACE(fprintf(stderr,"%d: Running tclcommand_analyze_parse_local_stress_tensor\n",this_node));
  /* 'analyze stress profile ' */
  if (argc != 12) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "local_stress_tensor requires 12 inputs: x_periodic, y_periodic, z_periodic, x_range_start, y_range_start, z_range_start, x_range, y_range, z_range, x_bins, y_bins, z_bins", (char *)NULL);
    return(TCL_ERROR);
  }
  const char *usage = "usage: analyse local_stress_tensor <x_periodic> <y_periodic> <z_periodic> <x_range_start> <y_range_start> <z_range_start> <x_range> <y_range> <z_range> <x_bins> <y_bins> <z_bins>";
 
  for (i=0;i<3;i++) {
    if ( !ARG0_IS_I(periodic[i]) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,usage, (char *)NULL);
      return (TCL_ERROR);
    } else {
      argc -= 1;
      argv += 1;
    }
  }
  for (i=0;i<3;i++) {
    if ( !ARG0_IS_D(range_start[i]) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,usage, (char *)NULL);
      return (TCL_ERROR);
    } else {
      argc -= 1;
      argv += 1;
    }
  }
  for (i=0;i<3;i++) {
    if ( !ARG0_IS_D(range[i]) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,usage, (char *)NULL);
      return (TCL_ERROR);
    } else {
      argc -= 1;
      argv += 1;
    }
  }
  for (i=0;i<3;i++) {
    if ( !ARG0_IS_I(bins[i]) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,usage, (char *)NULL);
      return (TCL_ERROR);
    } else {
      argc -= 1;
      argv += 1;
    }
  }

  /* Allocate a doublelist of bins to keep track of stress profile */
  TensorInBin = (DoubleList *)malloc(bins[0]*bins[1]*bins[2]*sizeof(DoubleList));
  if ( TensorInBin ) {
  /* Initialize the stress profile */
    for ( i = 0 ; i < bins[0]*bins[1]*bins[2]; i++ ) {
      init_doublelist(&TensorInBin[i]);
      alloc_doublelist(&TensorInBin[i],9);
      for ( j = 0 ; j < 9 ; j++ ) {
	TensorInBin[i].e[j] = 0.0;
      }
    }
  } else {
    Tcl_AppendResult(interp,"could not allocate memory for local_stress_tensor", (char *)NULL);
    return (TCL_ERROR);
  }

  mpi_local_stress_tensor(TensorInBin, bins, periodic,range_start, range);
  PTENSOR_TRACE(fprintf(stderr,"%d: tclcommand_analyze_parse_local_stress_tensor: finished mpi_local_stress_tensor \n",this_node));

  /* Write stress profile to Tcl export */
  Tcl_AppendResult(interp, "{ LocalStressTensor } ", (char *)NULL);
  for ( i = 0 ; i < bins[0] ; i++) {
    for ( j = 0 ; j < bins[1] ; j++) {
      for ( k = 0 ; k < bins[2] ; k++) {
	Tcl_AppendResult(interp, " { ", (char *)NULL);
	sprintf(buffer," { %d %d %d } ",i,j,k);
	Tcl_AppendResult(interp,buffer, (char *)NULL);
	Tcl_AppendResult(interp, " { ", (char *)NULL);
	for ( l = 0 ; l < 9 ; l++) {
	  Tcl_PrintDouble(interp,TensorInBin[i*bins[1]*bins[2]+j*bins[2]+k].e[l],buffer);
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
	  Tcl_AppendResult(interp, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, " } ", (char *)NULL);
	Tcl_AppendResult(interp, " } ", (char *)NULL);
      }
    }
  }
  
  /* Free memory */
  for ( i = 0 ; i < bins[0]*bins[1]*bins[2] ; i++ ) {
    realloc_doublelist(&TensorInBin[i],0);
  }
  free(TensorInBin);
  return TCL_OK;
}
