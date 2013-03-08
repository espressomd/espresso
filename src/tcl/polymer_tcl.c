/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file polymer.c
    This file contains everything needed to create a start-up configuration
    of (partially charged) polymer chains with counterions and salt molecules,
    assigning velocities to the particles and crosslinking the polymers if necessary.
 
    The corresponding header file is polymer.h.
 
    Created:       27.02.2003 by BAM
       Based upon 'polymer.tcl' by BAM (20.02.2003).
*/

#include "polymer.h"
#include "communication.h"
#include "interaction_data.h"
#include "parser.h"
#include "grid.h"
#include "constraint.h"

/************************************************************* 
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/
int tclcommand_polymer (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int N_P, MPC; 
  double bond_length; 
  int part_id = 0; 
  double *posed = NULL; 
  double *posed2 = NULL;
  /* mode==0 equals "SAW", mode==1 equals "RW", mode==2 equals "PSAW" */
  int mode = 1; 
  double shield = 1.0; 
  int tmp_try,max_try = 30000;       
  double val_cM = 0.0; 
  int cM_dist = 1, type_nM = 0, type_cM = 1, type_bond = 0;
  double angle = -1.0, angle2 = -1.0;
  int constr = 0;
  char buffer[128 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;
#ifdef CONSTRAINTS
  int j;
#endif

  if (argc < 4) { Tcl_AppendResult(interp, "Wrong # of args! Usage: polymer <N_P> <MPC> <bond_length> [start <n> | pos <x> <y> <z> | mode | charge | distance | types | bond | angle | constraints]", (char *)NULL); return (TCL_ERROR); }
  if (!ARG_IS_I(1, N_P)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Number of polymers must be integer (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
 }
  else {
    if(N_P < 0) {
      Tcl_AppendResult(interp, "Number of polymers must be positive (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (!ARG_IS_I(2, MPC)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Number of monomers must be integer (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR); }
  else {
    if(MPC < 2) {
      Tcl_AppendResult(interp, "Number of monomers must be greater than 1 (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (!ARG_IS_D(3, bond_length)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Bondlength must be double (got: ", argv[3],")!", (char *)NULL); return (TCL_ERROR); }
  else {
    if(bond_length < 0.0) {
      Tcl_AppendResult(interp, "Bondlength  must be positive (got: ", argv[3],")!", (char *)NULL); return (TCL_ERROR); }
  }
  for (i=4; i < argc; i++) {
    /* [start <part_id>] */
    if (ARG_IS_S(i, "start")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, part_id)) {	
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "Index of polymer chain's first monomer must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if (part_id < 0) {
	    Tcl_AppendResult(interp, "Index of polymer chain's first monomer must be positive (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else { Tcl_AppendResult(interp, "Not enough arguments for start", (char *)NULL); return (TCL_ERROR); }
    }
    /* [pos <x> <y> <z>] */
    else if (ARG_IS_S(i, "pos")) {
      if (i+3 < argc) { 
	posed = malloc(3*sizeof(double));
	if (!(ARG_IS_D(i+1, posed[0]) && ARG_IS_D(i+2, posed[1]) && ARG_IS_D(i+3, posed[2]))) {
	  Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "The first start monomers position must be double (got: ",argv[i+1],",",argv[i+2],",",argv[i+3],")!", (char *)NULL);
	  return (TCL_ERROR); } else { i+=3; } }
      else { Tcl_AppendResult(interp, "The first start monomers position must be 3D!", (char *)NULL); return (TCL_ERROR); }
    }
    /* [mode { SAW | RW | PSAW } [<shield> [max_try]]] */
    else if (ARG_IS_S(i, "mode")) {
      if (i+1 < argc) {
	if (ARG_IS_S(i+1, "SAW") || ARG_IS_S(i+1, "PSAW")) {
	  if (ARG_IS_S(i+1, "SAW")) mode = 0; else mode = 2;
	  if ((i+2 >= argc) || !ARG_IS_D(i+2, shield)) { Tcl_ResetResult(interp); i++; }
	  else {
	    if (shield < 0) { Tcl_AppendResult(interp, "The SAW-shield must be positive (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
	    if ((i+3 >= argc) || !ARG_IS_I(i+3, max_try)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } } }
	else if (ARG_IS_S(i+1, "RW")) {
	  mode = 1;
	  if ((i+2 >= argc) || !ARG_IS_D(i+2, shield)) { Tcl_ResetResult(interp); i++; }
	      else if ((i+3 >= argc) || !ARG_IS_I(i+3, max_try)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } }
	else {
	  Tcl_AppendResult(interp, "The mode you specified does not exist (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      }
      else {Tcl_AppendResult(interp, "Not enough arguments for mode", (char *)NULL); return (TCL_ERROR); }
    }
    /* [charge <val_cM>] */
    else if (ARG_IS_S(i, "charge")) {
      if(i+1 < argc) {
	if (!ARG_IS_D(i+1, val_cM)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The charge of the chain's monomers must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else { i++; }
      }
      else { Tcl_AppendResult(interp, "Not enough arguments for charge", (char *)NULL); return (TCL_ERROR); }
    }
    /* [distance <cM_dist>] */
    else if (ARG_IS_S(i, "distance")) {
      if(i+1 <argc) {
	if (!ARG_IS_I(i+1, cM_dist)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The distance between two charged monomers' indices must be integer (got: ",argv[i+1],")!", (char *)NULL); 
	  return (TCL_ERROR); }
	else {
	  if(cM_dist < 0) { Tcl_AppendResult(interp, "The charge of the chain's monomers  must be positive (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
	  else { i++; }
	}
      }
      else { Tcl_AppendResult(interp, "Not enough arguments for distance", (char *)NULL); return (TCL_ERROR); }
    }
    /* [types <type_nM> [<type_cM>]] */
    else if (ARG_IS_S(i, "types")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, type_nM)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The type-# of neutral monomers must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((i+2 >= argc) || !ARG_IS_I(i+2, type_cM)) { Tcl_ResetResult(interp); i++; } else { i+=2; } }
      }
      else {Tcl_AppendResult(interp, "Not enough arguments for types", (char *)NULL); return (TCL_ERROR); }
    }
    /* [bond <type_bond>] */
    else if (ARG_IS_S(i, "bond") || ARG_IS_S(i, "FENE")) {
      if (i+1 < argc) {
	if (!ARG_IS_I(i+1, type_bond)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The type-# of the bond-interaction must be integer (got: ", argv[i+1],")!", (char *)NULL); 
	  return (TCL_ERROR); 
	}
	else {
	  i++;
	  
	}
      } else {
	Tcl_AppendResult(interp, "Not enough arguments for bond", 
			 (char *)NULL); 
	return (TCL_ERROR); 
      }
    }
    /* [angle <angle> [\<angle2\>]] */
    else if (ARG_IS_S(i, "angle")) {
      Tcl_AppendResult(interp, "Warning: The angle definition has been changed. See RELEASE-NOTES.\n", (char *)NULL);
      if (i+1 < argc) {
	if (ARG_IS_D(i+1, angle)) {
	  if (angle < 0.0) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "The angle phi must be positive (got: ",argv[i+1],")!", (char *)NULL); 
	    return (TCL_ERROR);
	  }
	  while (angle >= 2.0*PI) angle -= 2.0*PI;
	  i++;
    	  if (i+1 < argc) {
	    if (ARG_IS_D(i+1, angle2)) {
	      if (angle2 < 0.0) {
		Tcl_ResetResult(interp);
		Tcl_AppendResult(interp, "The angle theta must be positive (got: ",argv[i+1],")!", (char *)NULL); 
		return (TCL_ERROR); 
	      }
	      while(angle2 >= 2.0*PI) angle2 -= 2.0*PI;
	      i++;
	      if (i+3 < argc) {
		posed2=malloc(3*sizeof(double));
		if (ARG_IS_D(i+1, posed2[0]) && ARG_IS_D(i+2, posed2[1]) && ARG_IS_D(i+3, posed2[2])) {
		  i+=3; 
		} else {
		  free(posed2); 
		}
	      }
	    }
	  }
	}
      }
      else {Tcl_AppendResult(interp, "Not enough arguments for angle", (char *)NULL); return (TCL_ERROR); }
    }
    /* [constraints] */
    else if (ARG_IS_S(i, "constraints")) {
#ifndef CONSTRAINTS
      Tcl_AppendResult(interp, "Constraints are not compiled in!", (char *)NULL); return (TCL_ERROR); }
#else
      constr=1;
      tmp_try=0;
      for(j=0;j<n_constraints;j++){
	if(constraints[j].type==CONSTRAINT_MAZE || constraints[j].type==CONSTRAINT_PORE || constraints[j].type==CONSTRAINT_PLATE || constraints[j].type==CONSTRAINT_RHOMBOID)
	  tmp_try++;
      }
      if (tmp_try>0) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "Warning: Only constraints of type WALL/SPHERE/CYLINDER are respected!", (char *)NULL); return (TCL_ERROR); }
      else { i++; }
    }
#endif
    /* Default */
  else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
}

if(type_bond<0 || type_bond>=n_bonded_ia){
  Tcl_AppendResult(interp, "Please define a bonded interaction before setting up polymers!", (char *)NULL);
  return(TCL_ERROR);
 }

  if (fabs(val_cM) < 1e-10) { val_cM = 0.0; type_cM = type_nM; }

  POLY_TRACE(if (posed!=NULL) {if (posed2!=NULL) printf("int N_P %d, int MPC %d, double bond_length %f, int part_id %d, double posed (%f,%f,%f), int mode %d, double shield %f, int max_try %d, double val_cM %f, int cM_dist %d, int type_nM %d, int type_cM %d, int type_bond %d, double angle %f, double angle2 %f, double posed (%f,%f,%f), int constraints %d\n", N_P, MPC, bond_length, part_id, posed[0],posed[1],posed[2], mode, shield, max_try, val_cM, cM_dist, type_nM, type_cM, type_bond, angle,angle2, posed2[0], posed2[1], posed2[2], constr); else printf("int N_P %d, int MPC %d, double bond_length %f, int part_id %d, double posed (%f,%f,%f), int mode %d, double shield %f, int max_try %d, double val_cM %f, int cM_dist %d, int type_nM %d, int type_cM %d, int type_bond %d, double angle %f, double angle2 %f, double posed2 NULL, constraints %d\n", N_P, MPC, bond_length, part_id, posed[0],posed[1],posed[2], mode, shield, max_try, val_cM, cM_dist, type_nM, type_cM, type_bond,angle,angle2,constr);} else {if (posed2!=NULL) printf("int N_P %d, int MPC %d, double bond_length %f, int part_id %d, double posed NULL, int mode %d, double shield %f, int max_try %d, double val_cM %f, int cM_dist %d, int type_nM %d, int type_cM %d, int type_bond %d, double angle %f, double angle2 %f, double posed2 (%f,%f,%f), int constraints %d\n", N_P, MPC, bond_length, part_id, mode, shield, max_try, val_cM, cM_dist, type_nM, type_cM, type_bond, angle, angle2,posed2[0],posed2[1],posed2[2],constr); else printf("int N_P %d, int MPC %d, double bond_length %f, int part_id %d, double posed NULL, int mode %d, double shield %f, int max_try %d, double val_cM %f, int cM_dist %d, int type_nM %d, int type_cM %d, int type_bond %d, double angle %f, double angle2 %f, double posed2 NULL, int constraints %d\n", N_P, MPC, bond_length, part_id, mode, shield, max_try, val_cM, cM_dist, type_nM, type_cM, type_bond, angle, angle2, constr);});

  tmp_try = polymerC(N_P, MPC, bond_length, part_id, posed, mode, shield, max_try, val_cM, cM_dist, type_nM, type_cM, type_bond, angle, angle2, posed2, constr);
  if (tmp_try == -1) {
    sprintf(buffer, "Failed to find a suitable place for the start-monomer for %d times!\nUse option 'mode { SAW | RW | PSAW } <shield> <max_try>' to increase this limit...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -2) {
    sprintf(buffer, "Failed to place current polymer chain in the simulation box for %d times!\nUse option 'mode { SAW | RW | PSAW } <shield> <max_try>' to increase this limit...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -3) {
    sprintf(buffer, "Failed upon creating one of the monomers in Espresso!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try == -4) {
    sprintf(buffer, "Failed upon removing one of the monomers in Espresso while trying to reset current chain!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return gather_runtime_errors(interp, tmp_try);
}

int tclcommand_counterions (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int N_CI; int part_id = n_total_particles; 
  int mode = 0; double shield = 0.0; int tmp_try,max_try = 30000;                             /* mode==0 equals "SAW", mode==1 equals "RW" */
  double val_CI = -1.0; int type_CI = 2;
  char buffer[128 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 2) { Tcl_AppendResult(interp, "Wrong # of args! Usage: counterions <N_CI> [options]", (char *)NULL); return (TCL_ERROR); }
  if (!ARG_IS_I(1, N_CI)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Number of conterions must be integer (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
 }
  else {
    if(N_CI < 0) {
      Tcl_AppendResult(interp, "Number of counterions must be positive (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR); }
  }
  for (i=2; i < argc; i++) {
    /* [start <part_id>] */
    if (ARG_IS_S(i, "start")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, part_id)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "Index of first counterion must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if (part_id < 0) {
	    Tcl_AppendResult(interp, "Index of first counterion must be positive (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for start!", (char *)NULL); return (TCL_ERROR); }
    }
    /* [mode { SAW | RW } [<shield> [max_try]]] */
    else if (ARG_IS_S(i, "mode")) {
      if(i+1 < argc) {
	if (ARG_IS_S(i+1, "SAW")) {
	  mode = 0;
	  if ((i+2 >= argc ) || !ARG_IS_D(i+2, shield)) { Tcl_ResetResult(interp); i++; }
	  else {
	    if (shield < 0) { Tcl_AppendResult(interp, "The SAW-shield must be positive (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
	    if ((i+3) >= argc || !ARG_IS_I(i+3, max_try)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } }
        }
	else if (ARG_IS_S(i+1, "RW")) {
	  mode = 1;
	  if ((i+2 >= argc) || !ARG_IS_D(i+2, shield)) { Tcl_ResetResult(interp); i++; }
	  else if ((i+3 >= argc) || !ARG_IS_I(i+3, max_try)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } }
        else {
	  Tcl_AppendResult(interp, "The mode you specified does not exist (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      }
      else {
        Tcl_AppendResult(interp, "Not enough arguments for mode!", (char *)NULL); return (TCL_ERROR); }
    }
    /* [charge <val_CI>] */
    else if (ARG_IS_S(i, "charge")) {
      if (!ARG_IS_D(i+1, val_CI)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The charge of the counterions must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else { i++; }
    }
    /* [type <type_CI>] */
    else if (ARG_IS_S(i, "type")) {
      if (!ARG_IS_I(i+1, type_CI)) { 
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The type-# of the counterions must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else { i++; }
    }
    /* default */
    else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }

  POLY_TRACE(printf("int N_CI %d, int part_id %d, int mode %d, double shield %f, int max_try %d, double val_CI %f, int type_CI %d\n", N_CI, part_id, mode, shield, max_try, val_CI, type_CI));

  tmp_try = counterionsC(N_CI, part_id, mode, shield, max_try, val_CI, type_CI);
  if (tmp_try == -1) {
    sprintf(buffer, "Failed to place current counterion in the simulation box for %d times!\nAborting...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -3) {
    sprintf(buffer, "Failed upon creating one of the monomers in Espresso!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
  else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return gather_runtime_errors(interp, tmp_try);
}

int tclcommand_salt (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int N_pS, N_nS; int part_id = n_total_particles; 
  int mode = 0; double shield = 0.0; int tmp_try,max_try = 30000;                             /* mode==0 equals "SAW", mode==1 equals "RW" */
  double val_pS = 1.0, val_nS = -1.0; int type_pS = 3, type_nS = 4;
  double rad=0.;
  char buffer[128 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 3) { Tcl_AppendResult(interp, "Wrong # of args! Usage: salt <N_pS> <N_nS> [options]", (char *)NULL); return (TCL_ERROR); }
  if (!ARG_IS_I(1, N_pS)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Number of positive salt-ions must be integer (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
  }
  else {
    if(N_pS < 0) {
      Tcl_AppendResult(interp, "Number of positive salt-ions must be positive (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (!ARG_IS_I(2, N_nS)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Number of negative salt-ions must be integer (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR);
  }
  else {
    if(N_nS < 0) {
      Tcl_AppendResult(interp, "Number of negative salt-ions must be positive (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR); }
  }
  for (i=3; i < argc; i++) {
    /* [start <part_id>] */
    if (ARG_IS_S(i, "start")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, part_id)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "Index of first salt ion must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if (part_id < 0) {
	    Tcl_AppendResult(interp, "Index of first salt ion must be positive (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for start!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [mode { SAW | RW } [<shield> [max_try]]] */
    else if (ARG_IS_S(i, "mode")) {
      if(i+1 < argc) {
	if (ARG_IS_S(i+1, "SAW")) {
	  mode = 0;
	  if ((i+2 >= argc) || !ARG_IS_D(i+2, shield)) { Tcl_ResetResult(interp); i++; }
	  else {
	    if (shield < 0) { Tcl_AppendResult(interp, "The SAW-shield must be positive (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
	    if ((i+3 >= argc) || !ARG_IS_I(i+3, max_try)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } } }
	else if (ARG_IS_S(i+1, "RW")) {
	  mode = 1;
	  if ((i+2 >= argc) || !ARG_IS_D(i+2, shield)) { Tcl_ResetResult(interp); i++; }
	  else if ((i+3 >= argc) || !ARG_IS_I(i+3, max_try)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } }
	else {
	  Tcl_AppendResult(interp, "The mode you specified does not exist (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for mode!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [charges <val_pS> [val_nS]] */
    else if (ARG_IS_S(i, "charges")) {
      if(i+1 < argc) {
	if (!ARG_IS_D(i+1, val_pS)) { 
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The charge of positive salt ions must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((i+2 >= argc) || !ARG_IS_D(i+2, val_nS)) { Tcl_ResetResult(interp); val_nS = -1.*val_pS; i++; } 
	  else { i+=2; } }
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for charges!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [types <type_pS> [<type_nS>]] */
    else if (ARG_IS_S(i, "types")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, type_pS)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The type-# of positive salt ions must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((i+2 >= argc) || !ARG_IS_I(i+2, type_nS)) { Tcl_ResetResult(interp); type_nS = type_pS; i++; } 
	  else { i+=2; } }
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for types!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [rad <rad> ] */
    else if (ARG_IS_S(i, "rad")) {
      if(i+1 < argc) {
	if ((!ARG_IS_D(i+1, rad)) || rad < 0.)  {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The radius for the cell model must be positive double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else { i+=2; }
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for rad!",(char *)NULL); return (TCL_ERROR); }
    }
    /* default */
  else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }
  
  POLY_TRACE(printf("int N_pS %d, int N_nS %d, int part_id %d, int mode %d, double shield %f, int max_try %d, double val_pS %f, double val_nS %f, int type_pS %d, int type_nS %d, double rad %f\n", N_pS, N_nS, part_id, mode, shield, max_try, val_pS, val_nS, type_pS, type_nS, rad));

  tmp_try = saltC(N_pS, N_nS, part_id, mode, shield, max_try, val_pS, val_nS, type_pS, type_nS, rad);
  if (tmp_try == -1) {
    sprintf(buffer, "Failed to place current positive salt ion in the simulation box for %d times!\nAborting...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -2) {
    sprintf(buffer, "Failed to place current negative salt ion in the simulation box for %d times!\nAborting...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -3) {
    sprintf(buffer, "Failed upon creating one of the monomers in Espresso!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
  else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return gather_runtime_errors(interp, tmp_try);
}

int tclcommand_velocities (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  double v_max; int part_id = 0, N_T = n_total_particles;
  double tmp_try;
  char buffer[128 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 2) { Tcl_AppendResult(interp, "Wrong # of args! Usage: velocities <v_max> [options]", (char *)NULL); return (TCL_ERROR); }
  if (!ARG_IS_D(1, v_max)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Maximum velocity must be double (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
  }
  else {
    if(v_max < 0) {
      Tcl_AppendResult(interp, "Maximum velocity must be positive (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR); }
  }
  for (i=2; i < argc; i++) {
    /* [start <part_id>] */
    if (ARG_IS_S(i, "start")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, part_id)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "Index of first particle must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((part_id < 0) || (part_id>=n_total_particles)) {
	    sprintf(buffer,"Index of first particle must be in [0,%d[ (got: ", n_total_particles);
	    Tcl_AppendResult(interp, buffer, argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for start!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [count <N_T>] */
    else if (ARG_IS_S(i, "count")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, N_T)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The amount of particles to be set must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((N_T < 0) || (part_id+N_T > n_total_particles)) {
	    sprintf(buffer,"The amount of particles to be set must be in [0,%d] (got: ",n_total_particles-part_id);
	    Tcl_AppendResult(interp, buffer, argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for count!",(char *)NULL); return (TCL_ERROR); }

    }
    /* default */
    else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (part_id+N_T > n_total_particles) N_T = n_total_particles - part_id;

  POLY_TRACE(printf("double v_max %f, int part_id %d, int N_T %d\n", v_max, part_id, N_T));

  tmp_try = velocitiesC(v_max, part_id, N_T);
  sprintf(buffer, "%f", tmp_try); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  return gather_runtime_errors(interp, TCL_OK);
}

int tclcommand_maxwell_velocities (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int part_id = 0, N_T = n_total_particles;
  double tmp_try;
  char buffer[128 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 1) { Tcl_AppendResult(interp, "Wrong # of args! Usage: maxwell_velocities [options]", (char *)NULL); return (TCL_ERROR); }
  for (i=1; i < argc; i++) {
    /* [start <part_id>] */
    if (ARG_IS_S(i, "start")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, part_id)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "Index of first particle must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((part_id < 0) || (part_id>=n_total_particles)) {
	    sprintf(buffer,"Index of first particle must be in [0,%d[ (got: ", n_total_particles);
	    Tcl_AppendResult(interp, buffer, argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for start!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [count <N_T>] */
    else if (ARG_IS_S(i, "count")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, N_T)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "The amount of particles to be set must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((N_T < 0) || (part_id+N_T > n_total_particles)) {
	    sprintf(buffer,"The amount of particles to be set must be in [0,%d] (got: ",n_total_particles-part_id);
	    Tcl_AppendResult(interp, buffer, argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
    }
    /* default */
    else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (part_id+N_T > n_total_particles) N_T = n_total_particles - part_id;

  POLY_TRACE(printf("int part_id %d, int N_T %d\n", part_id, N_T));

  tmp_try = maxwell_velocitiesC(part_id, N_T);
  sprintf(buffer, "%f", tmp_try); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  return gather_runtime_errors(interp, TCL_OK);
}

int tclcommand_crosslink (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int N_P, MPC; int part_id=0;
  double r_catch=1.9; int link_dist=2, chain_dist, type_bond=0, tmp_try,max_try=30000; 
  char buffer[128 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 3) { Tcl_AppendResult(interp, "Wrong # of args! Usage: crosslink <N_P> <MPC> [options]", (char *)NULL); return (TCL_ERROR); }
  if (!ARG_IS_I(1, N_P)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Number of polymer chains must be integer (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
  }
  else {
    if(N_P <= 1) {
      Tcl_AppendResult(interp, "Need at least 2 Polymers to crosslink (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (!ARG_IS_I(2, MPC)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Monomers per chain must be integer (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR);
  }
  else {
    if(MPC <= 1) {
      Tcl_AppendResult(interp, "Polymers must consist of at least 2 monomers per chain (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR); }
  }
  chain_dist = MPC;
  for (i=3; i < argc; i++) {
    /* [start <part_id>] */
    if (ARG_IS_S(i, "start")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, part_id)) {	
	  Tcl_AppendResult(interp, "Index of first particle must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((part_id < 0) || (part_id > n_total_particles - N_P*MPC)) {
	    sprintf(buffer,"Index of first particle must be in [0,%d] (got: ", n_total_particles - N_P*MPC);
	    Tcl_AppendResult(interp, buffer, argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for start!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [catch <r_catch>] */
    else if (ARG_IS_S(i, "catch")) {
      if(i+1 < argc) {
	if (!ARG_IS_D(i+1, r_catch)) {	
	  Tcl_AppendResult(interp, "Catching radius must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((r_catch < 0.) || (r_catch > dmax(dmax(box_l[0],box_l[1]),box_l[2]) )) {
	    sprintf(buffer,"Catching radius must be in [0,%f] (got: ", dmax(dmax(box_l[0],box_l[1]),box_l[2]) );
	    Tcl_AppendResult(interp, buffer, argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for catch!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [distLink <link_dist>] */
    else if (ARG_IS_S(i, "distLink")) {
      if (!ARG_IS_I(i+1, link_dist)) {	
	Tcl_AppendResult(interp, "Minimum distance between bonds must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if ((link_dist < 0) || (link_dist > MPC-1)) {
	  sprintf(buffer,"Minimum distance between bonds must be in [0,%d] (got: ",MPC-1);
	  Tcl_AppendResult(interp, buffer,argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
      i++;
    }
    /* [distChain <chain_dist>] */
    else if (ARG_IS_S(i, "distChain")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, chain_dist)) {	
	  Tcl_AppendResult(interp, "Minimum distance between partners must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
	else {
	  if ((chain_dist < 0) || (chain_dist > MPC)) {
	    sprintf(buffer,"Minimum distance between partners must be in [0,%d] (got: ",MPC);
	    Tcl_AppendResult(interp, buffer,argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      }
      else {
	Tcl_AppendResult(interp, "Not enough arguments for distChain!",(char *)NULL); return (TCL_ERROR); }
    }
    /* [bond <type_bond>] */
    else if (ARG_IS_S(i, "bond") || ARG_IS_S(i, "FENE")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, type_bond)) { 
	  Tcl_AppendResult(interp, "The type-# of the bind-interaction must be integer (got: ",
			   argv[i+1],")!", (char *)NULL); 
	  return (TCL_ERROR); 
	} else { i++; }
      } else {
	Tcl_AppendResult(interp, "Not enough arguments for bond!",(char *)NULL); 
	return (TCL_ERROR); 
      }
    }
    /* [trials <max_try>] */
    else if (ARG_IS_S(i, "trials")) {
      if(i+1 < argc) {
	if (!ARG_IS_I(i+1, max_try)) {	
	  Tcl_AppendResult(interp, "Amount of retries must be integer (got: ",argv[i+1],")!", 
			   (char *)NULL); 
	  return (TCL_ERROR); 
	} else {
	  if (max_try < 0) {
	    sprintf(buffer,"Amount of retries must be positive (got: ");
	    Tcl_AppendResult(interp, buffer,argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
	i++;
      } else {
	Tcl_AppendResult(interp, "Not enough arguments for trials!",(char *)NULL); 
	return (TCL_ERROR); 
      }
    }
    /* default */
    else { 
      Tcl_AppendResult(interp, 
		       "The parameters you supplied do not seem to be valid (stuck at: ",
		       argv[i],")!", (char *)NULL); 
      return (TCL_ERROR); 
    }
  }

  POLY_TRACE(printf("int N_P %d, int MPC %d, int part_id %d, double r_catch %f, int link_dist %d, int chain_dist %d, int type_bond %d, int max_try %d\n", N_P, MPC, part_id, r_catch, link_dist, chain_dist, type_bond, max_try));

  tmp_try = crosslinkC(N_P, MPC, part_id, r_catch, link_dist, chain_dist, type_bond, max_try);
  if (tmp_try == -1) {
    sprintf(buffer, "Failed to crosslink current system for %d times!\nAborting...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -2) {
    sprintf(buffer, "An error occured while crosslinking the system!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try == -3) {
    sprintf(buffer, "Failed upon submitting a bond to Espresso!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
  else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return gather_runtime_errors(interp, tmp_try);
}

int tclcommand_diamond (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  double a, bond_length; int MPC, N_CI = 0; double val_nodes = 0.0, val_cM = 0.0, val_CI = 0.0; int cM_dist = 1; int nonet = 0;
  char buffer[128 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i, tmp_try;

  if (argc < 4) { Tcl_AppendResult(interp, "Wrong # of args! Usage: diamond <a> <bond_length> <MPC> [options]", (char *)NULL); return (TCL_ERROR); }
  if (!ARG_IS_D(1, a)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Unit cell spacing must be double (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
  }
  if (!ARG_IS_D(2, bond_length)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Bond length must be double (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR);
  }
  else {
    if(bond_length < 0) {
      Tcl_AppendResult(interp, "Bond length must be positive (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (!ARG_IS_I(3, MPC)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Monomers per chain must be integer (got: ", argv[3],")!", (char *)NULL); return (TCL_ERROR);
  }
  else {
    if(MPC < 0) {
      Tcl_AppendResult(interp, "Monomers per chain must be positive (got: ", argv[3],")!", (char *)NULL); return (TCL_ERROR); }
  }
  for (i=4; i < argc; i++) {
    /* [counterions <N_CI>] */
    if (ARG_IS_S(i, "counterions")) {
      if (!ARG_IS_I(i+1, N_CI)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The number of counterions must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else { i++; }
    }
    /* [charges <val_nodes> <val_cM> <val_CI>] */
    else if (ARG_IS_S(i, "charges")) {
      if (i+3 >= argc) { Tcl_AppendResult(interp, "Wrong # of args! Usage: charges <val_nodes> <val_cM> <val_CI>!", (char *)NULL); return (TCL_ERROR); }
      if (!ARG_IS_D(i+1, val_nodes)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The charge of the nodes must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      if (!ARG_IS_D(i+2, val_cM)) { 
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The charge of the monomers must be double (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
      if (!ARG_IS_D(i+3, val_CI)) { 
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The charge of the counterions must be double (got: ",argv[i+3],")!", (char *)NULL); return (TCL_ERROR); }
      i+=3;
    }
    /* [distance <cM_dist>] */
    else if (ARG_IS_S(i, "distance")) {
      if (!ARG_IS_I(i+1, cM_dist)) { 
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The distance between two charged monomers' indices must be integer (got: ",argv[i+1],")!", (char *)NULL); 
	return (TCL_ERROR); }
      else { i++; }
    }
    /* [nonet] */
    else if (ARG_IS_S(i, "nonet")) {
      nonet = 1;
    }
    /* default */
    else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }

if(0 == n_bonded_ia){
  Tcl_AppendResult(interp, "Please define a bonded interaction before setting up polymers!", (char *)NULL);
  return(TCL_ERROR);
 }
 else{
   if(bonded_ia_params[0].num > 1){
     Tcl_AppendResult(interp, "Only two-body interacions are allowed with this command!", (char *)NULL);
     return(TCL_ERROR);
   }
 }

  POLY_TRACE(printf("double a %f, bond_length %f, int MPC %d, N_CI %d, double val_nodes %f, val_cM %f, val_CI %f, int cM_dist %d, nonet %d\n", a, bond_length, MPC, N_CI, val_nodes, val_cM, val_CI, cM_dist,nonet));

  tmp_try = diamondC(a, bond_length, MPC, N_CI, val_nodes, val_cM, val_CI, cM_dist, nonet);
  if (tmp_try == -3) {
    sprintf(buffer, "Failed upon creating one of the monomers in Espresso!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
  else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return gather_runtime_errors(interp, tmp_try);
}

int tclcommand_icosaeder (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  double a; int MPC, N_CI = 0; double val_cM = 0.0, val_CI = 0.0; int cM_dist = 1; 
  char buffer[128 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i, tmp_try;

  if (argc < 3) { Tcl_AppendResult(interp, "Wrong # of args! Usage: icosaeder <a> <MPC> [options]", (char *)NULL); return (TCL_ERROR); }
  if (!ARG_IS_D(1, a)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "a must be double (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
  }
  if (!ARG_IS_I(2, MPC)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Monomers per chain must be integer (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR);
  }
  else {
    if(MPC < 1) {
      Tcl_AppendResult(interp, "Monomers per chain must be positive (got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR); }
  }
  for (i=3; i < argc; i++) {
    /* [counterions <N_CI>] */
    if (ARG_IS_S(i, "counterions")) {
      if (!ARG_IS_I(i+1, N_CI)) { 
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The number of counterions must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else { i++; }
    }
    /* [charges <val_cM> <val_CI>] */
    else if (ARG_IS_S(i, "charges")) {
      if (i+2 >= argc) { Tcl_AppendResult(interp, "Wrong # of args! Usage: charges <val_cM> <val_CI>!", (char *)NULL); return (TCL_ERROR); }
      if (!ARG_IS_D(i+1, val_cM)) { 
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The charge of the monomers must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      if (!ARG_IS_D(i+2, val_CI)) { 
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The charge of the counterions must be double (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
      i+=2;
    }
    /* [distance <cM_dist>] */
    else if (ARG_IS_S(i, "distance")) {
      if (!ARG_IS_I(i+1, cM_dist)) { 
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "The distance between two charged monomers' indices must be integer (got: ",argv[i+1],")!", (char *)NULL); 
	return (TCL_ERROR); }
      else { i++; }
    }
    /* default */
    else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }

if(0 == n_bonded_ia){
  Tcl_AppendResult(interp, "Please define a bonded interaction before setting up polymers!", (char *)NULL);
  return(TCL_ERROR);
 }
 else{
   if(bonded_ia_params[0].num > 1){
     Tcl_AppendResult(interp, "Only two-body interacions are allowed with this command!", (char *)NULL);
     return(TCL_ERROR);
   }
 }

  POLY_TRACE(printf("double a %f, int MPC %d, N_CI %d, double val_cM %f, val_CI %f, int cM_dist %d\n", a, MPC, N_CI, val_cM, val_CI, cM_dist));

  tmp_try = icosaederC(a, MPC, N_CI, val_cM, val_CI, cM_dist);
  if (tmp_try == -3) {
    sprintf(buffer, "Failed upon creating one of the monomers in Espresso!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try == -2) {
    sprintf(buffer, "Failed upon creating a bond between edges around the vertices!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try == -1) {
    sprintf(buffer, "Failed upon connecting loose edges around vertices with chains along the middle third!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
  else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return gather_runtime_errors(interp, tmp_try);
}

