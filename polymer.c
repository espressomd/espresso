/** \file polymer.c
    This file contains everything needed to create a start-up configuration
    of (partially charged) polymer chains with counterions and salt molecules,
    assigning velocities to the particles and crosslinking the polymers if necessary.
 
    The corresponding header file is \ref polymer.h "polymer.h".
 
    Created:       27.02.2003 by BAM
       Based upon 'polymer.tcl' by BAM (20.02.2003).
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "polymer.h"
#include "grid.h"
#include "communication.h"
#include "random.h"
#include "debug.h"
#include "utils.h"




/************************************************************* 
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/



int mindist3(int part_id, double r_catch, int *ids) {
  /** C implementation of 'mindist <part_id> <r_catch>',
      which returns the size of an array <ids> of indices of particles which are 
      less than <r_catch> away from the position of the particle <part_id>. */
  Particle *partCfgMD;
  double dx,dy,dz;
  int i, caught=0;
  ids = (int *)malloc(n_total_particles*sizeof(int));

  partCfgMD = malloc(n_total_particles*sizeof(Particle));
  mpi_get_particles(partCfgMD); 
  for (i=0; i<n_total_particles; i++) {
    if (i != part_id) {
      dx = partCfgMD[part_id].r.p[0] - partCfgMD[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
      dy = partCfgMD[part_id].r.p[1] - partCfgMD[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
      dz = partCfgMD[part_id].r.p[2] - partCfgMD[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
      if (sqrt(SQR(dx)+SQR(dy)+SQR(dz)) < r_catch) ids[++caught]=partCfgMD[i].r.identity;
    }
  }
  ids = (int *)realloc(ids,caught*sizeof(int));
  free(partCfgMD); 
  return (caught);
}



double mindist4(double pos[3]) {
  /** C implementation of 'mindist <posx> <posy> <posz>',
      which returns the minimum distance of all current particles
      to position (<posx>, <posy>, <posz>) as a double.
      If it fails, return value equals -1. */
  Particle *partCfgMD;
  double mindist=30000.0, dx,dy,dz;
  int i;

  if (n_total_particles ==0) return (dmin(dmin(box_l[0],box_l[1]),box_l[2]));
  partCfgMD = malloc(n_total_particles*sizeof(Particle));
  mpi_get_particles(partCfgMD); 
  for (i=0; i<n_total_particles; i++) {
    dx = pos[0] - partCfgMD[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
    dy = pos[1] - partCfgMD[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
    dz = pos[2] - partCfgMD[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
    mindist = dmin(mindist, SQR(dx)+SQR(dy)+SQR(dz));
  }
  free(partCfgMD); 
  if (mindist<30000.0) return (sqrt(mindist)); else return (-1.0);
}



int collision(double pos[3], double shield) {
  /** Checks whether a particle at coordinates (<posx>, <posy>, <posz>) collides
      with any other particle due to a minimum image distance smaller than <shield>. 
      Returns '1' if there is a collision, '0' otherwise. */
  if (mindist4(pos) > shield) return (0); else return (1);
}



int polymer (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  /** Implementation of the tcl-command
      polymer <N_P> <MPC> <bond_length> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] 
                                        [charge <val_cM>] [distance <cM_dist>] [types <type_nM> [<type_cM>]] [FENE <type_FENE>]
      Creates some polymer chains within the simulation box,
      and returns how often the attempt to place a monomer failed in the worst case.
      Parameters:  <N_P>         = how many polymers to create
                   <MPC>         = monomers per chain
                   <bond_length> = length of the bonds between two monomers
                   <box_length>  = length of the simulation box
                   <part_id>     = particle number of the start monomer (defaults to '0')
		   <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
		   <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0')
		   <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
		   <val_cM>      = valency of charged monomers (defaults to '0.0')
		   <cM_dist>     = distance between two charged monomers' indices (defaults to '1')
		   <type_{n|c}P> = type number of {neutral|charged} monomers to be used with "part" (default to '0' and '1')
		   <type_FENE>   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0') 
      If <val_cM> < 1e-10, the charge is assumed to be zero, and <type_cM> = <type_nM>.                                                    */
  int N_P, MPC; double bond_length; int part_id = 0; 
  int mode = 0; double shield = 0.0; int tmp_try,max_try = 30000;                             /* mode==0 equals "SAW", mode==1 equals "RW" */
  double val_cM = 0.0; int cM_dist = 1, type_nM = 0, type_cM = 1, type_FENE = 0;
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 4) { Tcl_AppendResult(interp, "Wrong # of args! Usage: polymer <N_P> <MPC> <bond_length> [options]", (char *)NULL); return (TCL_ERROR); }
  N_P = atoi(argv[1]); MPC = atoi(argv[2]); bond_length = atof(argv[3]);
  for (i=4; i < argc; i++) {
    /* [start <part_id>] */
    if (!strncmp(argv[i], "start", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &part_id) == TCL_ERROR) {	
	Tcl_AppendResult(interp, "Index of polymer chain's first monomer must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if (part_id < 0) {
	  Tcl_AppendResult(interp, "Index of polymer chain's first monomer must be positive (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
      i++;
    }
    /* [mode { SAW | RW } [<shield> [max_try]]] */
    else if (!strncmp(argv[i], "mode", strlen(argv[i]))) {
      if (!strncmp(argv[i+1], "SAW", strlen(argv[i+1]))) {
	mode = 0;
	if ((i+2 >= argc) || (Tcl_GetDouble(interp, argv[i+2], &shield) == TCL_ERROR)) { Tcl_ResetResult(interp); i++; }
	else {
	  if (shield < 0) { Tcl_AppendResult(interp, "The SAW-shield must be positive (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
	  if ((i+3 >= argc) || (Tcl_GetInt(interp, argv[i+3], &max_try) == TCL_ERROR)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } } }
      else if (!strncmp(argv[i+1], "RW", strlen(argv[i+1]))) {
	mode = 1; i++; }
      else {
	Tcl_AppendResult(interp, "The mode you specified does not exist (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
    }
    /* [charge <val_cM>] */
    else if (!strncmp(argv[i], "charge", strlen(argv[i]))) {
      if (Tcl_GetDouble(interp, argv[i+1], &val_cM) == TCL_ERROR) { 
	Tcl_AppendResult(interp, "The charge of the chain's monomers must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else { i++; }
    }
    /* [distance <cM_dist>] */
    else if (!strncmp(argv[i], "distance", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &cM_dist) == TCL_ERROR) { 
	Tcl_AppendResult(interp, "The distance between two charged monomers' indices must be integer (got: ",argv[i+1],")!", (char *)NULL); 
	return (TCL_ERROR); }
      else { i++; }
    }
    /* [types <type_nM> [<type_cM>]] */
    else if (!strncmp(argv[i], "types", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &type_nM) == TCL_ERROR) { 
	Tcl_AppendResult(interp, "The type-# of neutral monomers must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if ((i+2 >= argc) || (Tcl_GetInt(interp, argv[i+2], &type_cM) == TCL_ERROR)) { Tcl_ResetResult(interp); i++; } else { i+=2; } }
    }
    /* [FENE <type_FENE>] */
    else if (!strncmp(argv[i], "FENE", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &type_FENE) == TCL_ERROR) { 
	Tcl_AppendResult(interp, "The type-# of the FENE-interaction must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else { i++; }
    }
    /* default */
    else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (val_cM < 1e-10) { val_cM = 0.0; type_cM = type_nM; }

  POLY_TRACE(printf("int N_P %d, int MPC %d, double bond_length %f, int part_id %d, int mode %d, double shield %f, int max_try %d, double val_cM %f, int cM_dist %d, int type_nM %d, int type_cM %d, int type_FENE %d\n", N_P, MPC, bond_length, part_id, mode, shield, max_try, val_cM, cM_dist, type_nM, type_cM, type_FENE));

  tmp_try = polymerC(N_P, MPC, bond_length, part_id, mode, shield, max_try, val_cM, cM_dist, type_nM, type_cM, type_FENE);
  if (tmp_try == -1) {
    sprintf(buffer, "Failed to find a suitable place for the start-monomer for %d times!\nAborting...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -2) {
    sprintf(buffer, "Failed to place current polymer chain in the simulation box for %d times!\nAborting...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -3) {
    sprintf(buffer, "Failed upon creating one of the monomers in tcl_md!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
  else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (tmp_try);
}



int polymerC(int N_P, int MPC, double bond_length, int part_id, int mode, double shield, int max_try, 
	     double val_cM, int cM_dist, int type_nM, int type_cM, int type_FENE) {
  /** C implementation of 'polymer <N_P> <MPC> <bond_length> [options]',
      which returns how often the attempt to place a monomer failed in the worst case. */
  int p,n, cnt1,cnt2,max_cnt, bond[2];
  double theta,phi,pos[3],poz[3];

  cnt1 = cnt2 = max_cnt = 0;
  for (p=0; p<N_P; p++) {
    for (cnt2=0; cnt2<max_try; cnt2++) {
      /* place start monomer */
      for (cnt1=0; cnt1<max_try; cnt1++) {
	pos[0]=box_l[0]*d_random();
	pos[1]=box_l[1]*d_random();
	pos[2]=box_l[2]*d_random();
	if ((mode!=0) || (collision(pos, shield)==0)) break;
	POLY_TRACE(printf("s"); fflush(NULL));
      }
      if (cnt1 >= max_try) return (-1);
      if (place_particle(part_id, pos)==TCL_ERROR) return (-3);
      if (set_particle_q(part_id, val_cM)==TCL_ERROR) return (-3);
      if (set_particle_type(part_id, type_cM)==TCL_ERROR) return (-3);
      part_id++; max_cnt=imax(cnt1, max_cnt);
      POLY_TRACE(printf("S"); fflush(NULL));
      
      /* place remaining monomers */
      for (n=1; n<MPC; n++) {
	poz[0]=pos[0]; poz[1]=pos[1]; poz[2]=pos[2];
	for (cnt1=0; cnt1<max_try; cnt1++) {
	  theta  =     PI*d_random();
	  phi    = 2.0*PI*d_random();
	  pos[0] = poz[0] + bond_length*sin(theta)*cos(phi);
	  pos[1] = poz[1] + bond_length*sin(theta)*sin(phi);
	  pos[2] = poz[2] + bond_length*cos(theta);
	  if ((mode!=0) || (collision(pos, shield)==0)) break;
	  POLY_TRACE(printf("m"); fflush(NULL));
	}
	if (cnt1 >= max_try) {
	  fprintf(stderr,"Warning! Attempt #%d to build polymer %d failed after %d unsuccessful trials to place monomer %d!\n",cnt2+1,p,cnt1,n);
	  fprintf(stderr,"         Retrying by re-setting the start-monomer of current chain...\n");
	  part_id = part_id - n; n=0; break;
	}
	bond[0] = type_FENE; bond[1] = n-1;
	if (place_particle(part_id, pos)==TCL_ERROR) return (-3);
	if (set_particle_q(part_id, ((n % cM_dist==0) ? val_cM : 0.0) )==TCL_ERROR) return (-3);
	if (set_particle_type(part_id, ((n % cM_dist==0) ? type_cM : type_nM) )==TCL_ERROR) return (-3);
	if (change_particle_bond(part_id, bond, 0)==TCL_ERROR) return (-3);
	part_id++; max_cnt=imax(cnt1, max_cnt);
	POLY_TRACE(printf("M"); fflush(NULL));
      }
      if ((mode!=0) || (n>0)) break;
    }
    POLY_TRACE(printf(" %d/%d->%d \n",cnt1,cnt2,max_cnt));
    if (cnt2 >= max_try) return(-2); else max_cnt = imax(max_cnt,imax(cnt1,cnt2));
  }
  return(imax(max_cnt,cnt2));
}



int counterions (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  /** Implementation of the tcl-command
      counterions <N_CI> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] [charge <val_CI>] [type <type_CI>]
      Creates <N_CI> counterions of charge <val_CI> within the simulation box,
      and returns how often the attempt to place a particle failed in the worst case.
      Parameters:  <N_CI>        = number of counterions to create
                   <part_id>     = particle number of the first counterion (defaults to 'n_total_particles')
		   <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
		   <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0')
		   <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
		   <val_CI>      = valency of the counterions (defaults to '-1.0')
		   <type_CI>     = type number of the counterions to be used with "part" (default to '2') */
  int N_CI; int part_id = n_total_particles; 
  int mode = 0; double shield = 0.0; int tmp_try,max_try = 30000;                             /* mode==0 equals "SAW", mode==1 equals "RW" */
  double val_CI = -1.0; int type_CI = 2;
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 2) { Tcl_AppendResult(interp, "Wrong # of args! Usage: counterions <N_CI> [options]", (char *)NULL); return (TCL_ERROR); }
  N_CI = atoi(argv[1]);
  for (i=2; i < argc; i++) {
    /* [start <part_id>] */
    if (!strncmp(argv[i], "start", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &part_id) == TCL_ERROR) {	
	Tcl_AppendResult(interp, "Index of first counterion must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if (part_id < 0) {
	  Tcl_AppendResult(interp, "Index of first counterion must be positive (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
      i++;
    }
    /* [mode { SAW | RW } [<shield> [max_try]]] */
    else if (!strncmp(argv[i], "mode", strlen(argv[i]))) {
      if (!strncmp(argv[i+1], "SAW", strlen(argv[i+1]))) {
	mode = 0;
	if ((i+2 >= argc) || (Tcl_GetDouble(interp, argv[i+2], &shield) == TCL_ERROR)) { Tcl_ResetResult(interp); i++; }
	else {
	  if (shield < 0) { Tcl_AppendResult(interp, "The SAW-shield must be positive (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
	  if ((i+3 >= argc) || (Tcl_GetInt(interp, argv[i+3], &max_try) == TCL_ERROR)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } } }
      else if (!strncmp(argv[i+1], "RW", strlen(argv[i+1]))) {
	mode = 1; i++; }
      else {
	Tcl_AppendResult(interp, "The mode you specified does not exist (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
    }
    /* [charge <val_CI>] */
    else if (!strncmp(argv[i], "charge", strlen(argv[i]))) {
      if (Tcl_GetDouble(interp, argv[i+1], &val_CI) == TCL_ERROR) { 
	Tcl_AppendResult(interp, "The charge of the counterions must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else { i++; }
    }
    /* [type <type_CI>] */
    else if (!strncmp(argv[i], "type", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &type_CI) == TCL_ERROR) { 
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
    sprintf(buffer, "Failed upon creating one of the monomers in tcl_md!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
  else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (tmp_try);
}



int counterionsC(int N_CI, int part_id, int mode, double shield, int max_try, double val_CI, int type_CI) {
  /** C implementation of 'counterions <N_CI> [options]',
      which returns how often the attempt to place a counterion failed in the worst case. */
  int n, cnt1,max_cnt;
  double pos[3];

  cnt1 = max_cnt = 0;
  for (n=0; n<N_CI; n++) {
    for (cnt1=0; cnt1<max_try; cnt1++) {
      pos[0]=box_l[0]*d_random();
      pos[1]=box_l[1]*d_random();
      pos[2]=box_l[2]*d_random();
      if ((mode!=0) || (collision(pos, shield)==0)) break;
      POLY_TRACE(printf("c"); fflush(NULL));
    }
    if (cnt1 >= max_try) return (-1);
    if (place_particle(part_id, pos)==TCL_ERROR) return (-3);
    if (set_particle_q(part_id, val_CI)==TCL_ERROR) return (-3);
    if (set_particle_type(part_id, type_CI)==TCL_ERROR) return (-3);
    part_id++; max_cnt=imax(cnt1, max_cnt);
    POLY_TRACE(printf("C"); fflush(NULL));
  }
  POLY_TRACE(printf(" %d->%d \n",cnt1,max_cnt));
  if (cnt1 >= max_try) return(-1); else return(imax(max_cnt,cnt1));
}



int salt (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  /** Implementation of the tcl-command
      salt <N_pS> <N_nS> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] [charges <val_pS> [<val_nS>]] [types <type_pS> [<type_nS>]]
      Creates <N_pS> positively and <N_nS> negatively charged salt ions of charge <val_pS> and <val_nS> within the simulation box,
      and returns how often the attempt to place a particle failed in the worst case.
      Parameters:  <N_pS>/<N_nS> = number of salt ions to create
                   <box_length>  = length of the simulation box
		   <part_id>     = particle number of the first salt ion (defaults to 'n_total_particles')
		   <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
		   <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0')
		   <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
		   <val_{p|n}S>  = valencies of the salt ions (default to '1' and '-1', respectively); if <val_nS> is not given, <val_nS> = -1*<val_pS>
		   <type_{p|n}S> = type numbers to be used with "part" (default to '3' and '4'); if <type_nS> is not given, <type_nS> = <type_pS> is assumed. */
  int N_pS, N_nS; int part_id = n_total_particles; 
  int mode = 0; double shield = 0.0; int tmp_try,max_try = 30000;                             /* mode==0 equals "SAW", mode==1 equals "RW" */
  double val_pS = 1.0, val_nS = -1.0; int type_pS = 3, type_nS = 4;
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 3) { Tcl_AppendResult(interp, "Wrong # of args! Usage: salt <N_pS> <N_nS> [options]", (char *)NULL); return (TCL_ERROR); }
  N_pS = atoi(argv[1]); N_nS = atoi(argv[2]);
  for (i=3; i < argc; i++) {
    /* [start <part_id>] */
    if (!strncmp(argv[i], "start", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &part_id) == TCL_ERROR) {	
	Tcl_AppendResult(interp, "Index of first salt ion must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if (part_id < 0) {
	  Tcl_AppendResult(interp, "Index of first salt ion must be positive (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
      i++;
    }
    /* [mode { SAW | RW } [<shield> [max_try]]] */
    else if (!strncmp(argv[i], "mode", strlen(argv[i]))) {
      if (!strncmp(argv[i+1], "SAW", strlen(argv[i+1]))) {
	mode = 0;
	if ((i+2 >= argc) || (Tcl_GetDouble(interp, argv[i+2], &shield) == TCL_ERROR)) { Tcl_ResetResult(interp); i++; }
	else {
	  if (shield < 0) { Tcl_AppendResult(interp, "The SAW-shield must be positive (got: ",argv[i+2],")!", (char *)NULL); return (TCL_ERROR); }
	  if ((i+3 >= argc) || (Tcl_GetInt(interp, argv[i+3], &max_try) == TCL_ERROR)) { Tcl_ResetResult(interp); i+=2; } else { i+=3; } } }
      else if (!strncmp(argv[i+1], "RW", strlen(argv[i+1]))) {
	mode = 1; i++; }
      else {
	Tcl_AppendResult(interp, "The mode you specified does not exist (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
    }
    /* [charges <val_pS> [val_nS]] */
    else if (!strncmp(argv[i], "charges", strlen(argv[i]))) {
      if (Tcl_GetDouble(interp, argv[i+1], &val_pS) == TCL_ERROR) { 
	Tcl_AppendResult(interp, "The charge of positive salt ions must be double (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if ((i+2 >= argc) || (Tcl_GetDouble(interp, argv[i+2], &val_nS) == TCL_ERROR)) { Tcl_ResetResult(interp); val_nS = -1.*val_pS; i++; } 
	else { i+=2; } }
    }
    /* [types <type_pS> [<type_nS>]] */
    else if (!strncmp(argv[i], "types", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &type_pS) == TCL_ERROR) { 
	Tcl_AppendResult(interp, "The type-# of positive salt ions must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if ((i+2 >= argc) || (Tcl_GetInt(interp, argv[i+2], &type_nS) == TCL_ERROR)) { Tcl_ResetResult(interp); type_nS = type_pS; i++; } 
	else { i+=2; } }
    }
    /* default */
    else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }

  POLY_TRACE(printf("int N_pS %d, int N_nS %d, int part_id %d, int mode %d, double shield %f, int max_try %d, double val_pS %f, double val_nS %f, int type_pS %d, int type_nS %d\n", N_pS, N_nS, part_id, mode, shield, max_try, val_pS, val_nS, type_pS, type_nS));

  tmp_try = saltC(N_pS, N_nS, part_id, mode, shield, max_try, val_pS, val_nS, type_pS, type_nS);
  if (tmp_try == -1) {
    sprintf(buffer, "Failed to place current positive salt ion in the simulation box for %d times!\nAborting...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -2) {
    sprintf(buffer, "Failed to place current negative salt ion in the simulation box for %d times!\nAborting...\n",max_try); tmp_try = TCL_ERROR; }
  else if (tmp_try == -3) {
    sprintf(buffer, "Failed upon creating one of the monomers in tcl_md!\nAborting...\n"); tmp_try = TCL_ERROR; }
  else if (tmp_try >= 0) {
    sprintf(buffer, "%d", tmp_try); tmp_try = TCL_OK; }
  else {
    sprintf(buffer, "Unknown error %d occured!\nAborting...\n",tmp_try); tmp_try = TCL_ERROR; }
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (tmp_try);
}



int saltC(int N_pS, int N_nS, int part_id, int mode, double shield, int max_try, double val_pS, double val_nS, int type_pS, int type_nS) {
  /** C implementation of 'salt <N_pS> <N_nS> [options]',
      which returns how often the attempt to place a salt ion failed in the worst case. */
  int n, cnt1,max_cnt;
  double pos[3];

  cnt1 = max_cnt = 0;

  /* Place positive salt ions */
  for (n=0; n<N_pS; n++) {
    for (cnt1=0; cnt1<max_try; cnt1++) {
      pos[0]=box_l[0]*d_random();
      pos[1]=box_l[1]*d_random();
      pos[2]=box_l[2]*d_random();
      if ((mode!=0) || (collision(pos, shield)==0)) break;
      POLY_TRACE(printf("p"); fflush(NULL));
    }
    if (cnt1 >= max_try) return (-1);
    if (place_particle(part_id, pos)==TCL_ERROR) return (-3);
    if (set_particle_q(part_id, val_pS)==TCL_ERROR) return (-3);
    if (set_particle_type(part_id, type_pS)==TCL_ERROR) return (-3);
    part_id++; max_cnt=imax(cnt1, max_cnt);
    POLY_TRACE(printf("P"); fflush(NULL));
  }
  POLY_TRACE(printf(" %d->%d \n",cnt1,max_cnt));
  if (cnt1 >= max_try) return(-1);

  /* Place negative salt ions */
  for (n=0; n<N_nS; n++) {
    for (cnt1=0; cnt1<max_try; cnt1++) {
      pos[0]=box_l[0]*d_random();
      pos[1]=box_l[1]*d_random();
      pos[2]=box_l[2]*d_random();
      if ((mode!=0) || (collision(pos, shield)==0)) break;
      POLY_TRACE(printf("n"); fflush(NULL));
    }
    if (cnt1 >= max_try) return (-1);
    if (place_particle(part_id, pos)==TCL_ERROR) return (-3);
    if (set_particle_q(part_id, val_nS)==TCL_ERROR) return (-3);
    if (set_particle_type(part_id, type_nS)==TCL_ERROR) return (-3);
    part_id++; max_cnt=imax(cnt1, max_cnt);
    POLY_TRACE(printf("N"); fflush(NULL));
  }
  POLY_TRACE(printf(" %d->%d \n",cnt1,max_cnt));
  if (cnt1 >= max_try) return(-2); else return(imax(max_cnt,cnt1));
}



int velocities (ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  /** Implementation of the tcl-command
      velocities <v_max> [start <part_id>] [count <N_T>]
      Sets the velocities of <N_T> particles to a random value [-vmax,vmax],
      and returns the averaged velocity when done.
      Parameters:  <v_max>       = maximum velocity to be used
		   <part_id>     = particle number of the first of the <N_T> particles (defaults to '0') 
		   <N_T>         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id') */
  double v_max; int part_id = 0, N_T = n_total_particles;
  double tmp_try;
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i;

  if (argc < 2) { Tcl_AppendResult(interp, "Wrong # of args! Usage: velocities <v_max> [options]", (char *)NULL); return (TCL_ERROR); }
  v_max = atof(argv[1]);
  for (i=2; i < argc; i++) {
    /* [start <part_id>] */
    if (!strncmp(argv[i], "start", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &part_id) == TCL_ERROR) {	
	Tcl_AppendResult(interp, "Index of first particle must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if ((part_id < 0) || (part_id>=n_total_particles)) {
	  sprintf(buffer,"Index of first particle must be in [0,%d[ (got: ", n_total_particles);
	  Tcl_AppendResult(interp, buffer ,argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
      i++;
    }
    /* [count <N_T>] */
    else if (!strncmp(argv[i], "count", strlen(argv[i]))) {
      if (Tcl_GetInt(interp, argv[i+1], &N_T) == TCL_ERROR) {	
	Tcl_AppendResult(interp, "The amount of particles to be set must be integer (got: ",argv[i+1],")!", (char *)NULL); return (TCL_ERROR); }
      else {
	if ((N_T < 0) || (part_id+N_T > n_total_particles)) {
	  sprintf(buffer,"The amount of particles to be set must be in [0,%d] (got: ",n_total_particles-part_id);
	  Tcl_AppendResult(interp, buffer,argv[i+1],")!", (char *)NULL); return (TCL_ERROR); } }
      i++;
    }
    /* default */
    else { Tcl_AppendResult(interp, "The parameters you supplied do not seem to be valid (stuck at: ",argv[i],")!", (char *)NULL); return (TCL_ERROR); }
  }
  if (part_id+N_T > n_total_particles) N_T = n_total_particles - part_id;

  POLY_TRACE(printf("double v_max %f, int part_id %d, int N_T %d\n", v_max, part_id, N_T));

  tmp_try = velocitiesC(v_max, part_id, N_T);
  sprintf(buffer, "%f", tmp_try); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  return (TCL_OK);
}



double velocitiesC(double v_max, int part_id, int N_T) {
  /** C implementation of 'velocities <v_max> [options]',
      which returns the averaged velocity assigned. */
  double r,theta,phi, v[3], v_av[3];
  int i;

  v_av[0] = v_av[1] = v_av[2] = 0.0;
  for (i=part_id; i < part_id+N_T; i++) {
    r     = v_max*d_random();
    theta =    PI*d_random();
    phi   = 2.*PI*d_random();
    v[0] = r*sin(theta)*cos(phi);  v_av[0]+=v[0];
    v[1] = r*sin(theta)*sin(phi);  v_av[1]+=v[1];
    v[2] = r*cos(theta);           v_av[2]+=v[2];
    if (set_particle_v(i, v)==TCL_ERROR) {
      fprintf(stderr, "Failed upon setting one of the velocities in tcl_md (current average: %f)!\n",sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2]))); 
      fprintf(stderr, "Aborting...\n"); errexit();
    }
  }
  return (sqrt(SQR(v_av[0]) + SQR(v_av[1]) + SQR(v_av[2])));
}



/*
proc crosslink { N_P MPC  {p1 NA} {p2 NA} {p3 NA} {p4 NA} {p5 NA} {p6 NA} } {
# crosslink <N_P> <MPC> [catch <r_catch>] [distance <link_dist>] [trials <max_try>]
# Evaluates the current configuration and connects each chain's end to a random monomer of another chain at most <r_catch> away,
# if the next crosslink from there is at least <link_dist> monomers away;
# returns how many ends have been successfully linked.
# Parameters:  <N_P>         = number of polymer chains
#              <MPC>         = monomers per chain
#              <r_catch>     = maximum length of a crosslink (defaults to '1.9')
#              <link_dist>   = minimum distance between the indices of two crosslinked monomers (defaults to '2')
#              <max_try>     = how often crosslinks should be removed if they are too close to other links (defaults to '30000')
    set param [list $p1 $p2 $p3 $p4 $p5 $p6]
    set r_catch 1.9; set link_dist 2; set max_try 30000
    for {set i 0} {$i < 6} {incr i} {
	switch [lindex $param $i] {
	    "catch" { incr i; set r_catch [lindex $param $i] 
		      if { ![string is double $r_catch] } { puts "Maximum length of a crosslink must be double (got: $r_catch)!\nAborting...\n"; exit } }
	    "distance" { incr i; set link_dist [lindex $param $i]
		      if { ![string is integer $link_dist] } { 
			  puts "The distance between two crosslinked monomers' indices must be integer (got: $link_dist)!\nAborting...\n"; exit } }
	    "trials" { incr i; set max_try [lindex $param $i]
		      if { ![string is integer $max_try] } { puts "The maximum \# of trials must be integer (got: $max_try)!\nAborting...\n"; exit } }
	    default { if { [lindex $param $i]!="NA" } {
		          puts "The parameter set you supplied ($param) does not seem to be valid (stuck at: [lindex $param $i])!\nAborting...\n"; exit }
	            }
	}
    }
    set bonds [countBonds [part]]
    for { set i 0 } { $i<$N_P } { incr i } {
	lappend ends [lrange [lindex $bonds [expr $i*$MPC]] 1 end]
	lappend link [mindist [expr $i*$MPC] $r_catch]
	lappend ends [lrange [lindex $bonds [expr ($i+1)*$MPC-1]] 1 end]
	lappend link [mindist [expr ($i+1)*$MPC-1] $r_catch]
    }
puts "$ends"
puts "$link"
    for { set i 0 } { $i<$N_P } { incr i } {
	for { set k 0 } { $k<2 } { incr k } {
	    set tmp_link1 [lindex $link [expr 2*$i+$k]]
	    set tmp_link2 ""
puts -nonewline "$i,$k: $tmp_link1 & $tmp_link2 => "; flush stdout
	    for { set j 0 } { $j<[llength $tmp_link1] } { incr j } {
		set tmp_link [lindex $tmp_link1 $j]
puts -nonewline "$tmp_link "
		if { ($tmp_link<[expr $i*$MPC] || $tmp_link>[expr ($i+1)*$MPC]) && $tmp_link<[expr $N_P*$MPC] } { lappend $tmp_link2 $tmp_link }
	    }
puts "$tmp_link2"
	    set link [lreplace $link [expr 2*$i+$k] [expr 2*$i+$k] $tmp_link2]
	}
    }
puts "$link"
    
}
*/
