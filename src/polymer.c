/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "polymer.h"
#include "grid.h"
#include "communication.h"
#include "interaction_data.h"
#include "random.h"
#include "parser.h"
#include "integrate.h"
#include "constraint.h"




/************************************************************* 
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/



int mindist3(int part_id, double r_catch, int *ids) {
  Particle *partCfgMD;
  double dx,dy,dz;
  int i, me, caught=0;

  partCfgMD = malloc(n_total_particles*sizeof(Particle));
  mpi_get_particles(partCfgMD, NULL);
  me = -1; /* Since 'mpi_get_particles' returns the particles unsorted, it's most likely that 'partCfgMD[i].p.identity != i'
	      --> prevent that! */
  for(i=0; i<n_total_particles; i++) if (partCfgMD[i].p.identity == part_id) me = i; 
  if (me == -1) {
    char *errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{049 failed to find desired particle %d} ",part_id);
    return 0;
  }
  for (i=0; i<n_total_particles; i++) {
    if (i != me) {
      dx = partCfgMD[me].r.p[0] - partCfgMD[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
      dy = partCfgMD[me].r.p[1] - partCfgMD[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
      dz = partCfgMD[me].r.p[2] - partCfgMD[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
      if (sqrt(SQR(dx)+SQR(dy)+SQR(dz)) < r_catch) ids[caught++]=partCfgMD[i].p.identity;
    }
  }
  free(partCfgMD); 
  return (caught);
}



double mindist4(double pos[3]) {
  Particle *partCfgMD;
  double mindist=30000.0, dx,dy,dz;
  int i;

  if (n_total_particles ==0) return (dmin(dmin(box_l[0],box_l[1]),box_l[2]));
  partCfgMD = malloc(n_total_particles*sizeof(Particle));
  mpi_get_particles(partCfgMD, NULL); 
  for (i=0; i<n_total_particles; i++) {
    dx = pos[0] - partCfgMD[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
    dy = pos[1] - partCfgMD[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
    dz = pos[2] - partCfgMD[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
    mindist = dmin(mindist, SQR(dx)+SQR(dy)+SQR(dz));
  }
  free(partCfgMD); 
  if (mindist<30000.0)
    return (sqrt(mindist));
  return (-1.0);
}

double buf_mindist4(double pos[3], int n_add, double *add) {
  double mindist=30000.0, dx,dy,dz;
  int i;

  if (n_add == 0) return (dmin(dmin(box_l[0],box_l[1]),box_l[2]));
  for (i=0; i<n_add; i++) {
    dx = pos[0] - add[3*i + 0];   dx -= dround(dx/box_l[0])*box_l[0];
    dy = pos[1] - add[3*i + 1];   dy -= dround(dy/box_l[1])*box_l[1];
    dz = pos[2] - add[3*i + 2];   dz -= dround(dz/box_l[2])*box_l[2];
    mindist = dmin(mindist, SQR(dx)+SQR(dy)+SQR(dz));
  }
  if (mindist<30000.0) return (sqrt(mindist));
  return (-1.0);
}



int collision(double pos[3], double shield, int n_add, double *add) {
  if (mindist4(pos) > shield && buf_mindist4(pos, n_add, add) > shield) return (0);
  return (1);
}



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
  return mpi_gather_runtime_errors(interp, tmp_try);
}

#ifdef CONSTRAINTS

int constraint_collision(double *p1, double *p2){
  Particle part1,part2;
  double d1,d2,v[3];
  Constraint *c;
  int i;
  double folded_pos1[3];
  double folded_pos2[3];
  int img[3];

  memcpy(folded_pos1, p1, 3*sizeof(double));
  fold_position(folded_pos1, img);

  memcpy(folded_pos2, p2, 3*sizeof(double));
  fold_position(folded_pos2, img);

  for(i=0;i<n_constraints;i++){
    c=&constraints[i];
    switch(c->type){
    case CONSTRAINT_WAL:
      calculate_wall_dist(&part1,folded_pos1,&part1,&c->c.wal,&d1,v);
      calculate_wall_dist(&part2,folded_pos2,&part2,&c->c.wal,&d2,v);
      if(d1*d2<=0.0)
	return 1;
      break;
    case CONSTRAINT_SPH:
      calculate_sphere_dist(&part1,folded_pos1,&part1,&c->c.sph,&d1,v);
      calculate_sphere_dist(&part2,folded_pos2,&part2,&c->c.sph,&d2,v);
      if(d1*d2<0.0)
	return 1;
      break;
    case CONSTRAINT_CYL:
      calculate_cylinder_dist(&part1,folded_pos1,&part1,&c->c.cyl,&d1,v);
      calculate_cylinder_dist(&part2,folded_pos2,&part2,&c->c.cyl,&d2,v);
      if(d1*d2<0.0)
	return 1;
      break;
    case CONSTRAINT_MAZE:
    case CONSTRAINT_PORE:
    case CONSTRAINT_PLATE:
    case CONSTRAINT_RHOMBOID:
      break;
    }
  }
  return 0;
}

#endif

int polymerC(int N_P, int MPC, double bond_length, int part_id, double *posed, 
	     int mode, double shield, int max_try, double val_cM, int cM_dist, 
	     int type_nM, int type_cM, int type_bond, 
	     double angle, double angle2, double *posed2, int constr) {
  int p,n, cnt1,cnt2,max_cnt, bond_size, *bond, i;
  double phi,zz,rr;
  double *poly;
  double pos[3];
  double poz[3];
  double poy[3] = {0, 0, 0};
  double pox[3] = {0, 0, 0};
  double a[3] = {0, 0, 0};
  double b[3],c[3]={0., 0., 0.},d[3];
  double absc;
  poly = malloc(3*MPC*sizeof(double));

  bond_size = bonded_ia_params[type_bond].num;
  bond = malloc(sizeof(int) * (bond_size + 1));
  bond[0] = type_bond;

  cnt1 = cnt2 = max_cnt = 0;
  for (p=0; p < N_P; p++) {
    for (cnt2=0; cnt2 < max_try; cnt2++) {
      /* place start monomer */
      if (posed!=NULL) {
	/* if position of 1st monomer is given */
	if (p > 0) {
	  free(posed);
	  posed=NULL;
	} else {
	  pos[0]=posed[0];
	  pos[1]=posed[1];
	  pos[2]=posed[2];
	}
      } else {
	/* randomly set position */
	for (cnt1=0; cnt1<max_try; cnt1++) {
	  pos[0]=box_l[0]*d_random();
	  pos[1]=box_l[1]*d_random();
	  pos[2]=box_l[2]*d_random();
	  if ((mode==1) || (collision(pos, shield, 0, NULL)==0)) break;
	  POLY_TRACE(printf("s"); fflush(NULL));
	}
	if (cnt1 >= max_try) { free(poly); return (-1); }
      }
      poly[0] = pos[0]; poly[1] = pos[1]; poly[2] = pos[2];
      max_cnt=imax(cnt1, max_cnt);
      POLY_TRACE(printf("S"); fflush(NULL));
      //POLY_TRACE(/* printf("placed Monomer 0 at (%f,%f,%f)\n",pos[0],pos[1],pos[2]) */);

      poz[0]=pos[0]; poz[1]=pos[1]; poz[2]=pos[2];

      /* place 2nd monomer */
      n=1;
      if (posed2 != NULL && posed != NULL && angle2 > -1.0) {
	/* if position of 2nd monomer is given */
	pos[0]=posed2[0];
	pos[1]=posed2[1];
	pos[2]=posed2[2];
	/* calculate preceding monomer so that bond_length is correct */
	absc=sqrt(SQR(pos[0]-poz[0])+SQR(pos[1]-poz[1])+SQR(pos[2]-poz[2]));
	poz[0]=pos[0]+(poz[0]-pos[0])*bond_length/absc;
	poz[1]=pos[1]+(poz[1]-pos[1])*bond_length/absc;
	poz[2]=pos[2]+(poz[2]-pos[2])*bond_length/absc;
	//POLY_TRACE(/* printf("virtually shifted position of first monomer to (%f,%f,%f)\n",poz[0],poz[1],poz[2]) */);
      } else {
	/* randomly place 2nd monomer */
	for (cnt1=0; cnt1<max_try; cnt1++) {
	  zz     = (2.0*d_random()-1.0)*bond_length;
          rr     = sqrt(SQR(bond_length)-SQR(zz));
	  phi    = 2.0*PI*d_random();
	  pos[0] = poz[0]+rr*cos(phi);
	  pos[1] = poz[1]+rr*sin(phi);
	  pos[2] = poz[2]+zz;
#ifdef CONSTRAINTS
	  if(constr==0 || constraint_collision(pos,poly+3*(n-1))==0){
#endif

	    if (mode==1 || collision(pos, shield, n, poly)==0) break;
	    if (mode==0) { cnt1 = -1; break; }
#ifdef CONSTRAINTS
	  }
#endif
	  POLY_TRACE(printf("m"); fflush(NULL));
	}
	if (cnt1 >= max_try) {
	  fprintf(stderr,"\nWarning! Attempt #%d to build polymer %d failed while placing monomer 2!\n",cnt2+1,p);
	  fprintf(stderr,"         Retrying by re-setting the start-monomer of current chain...\n");
	}
	if (cnt1 == -1 || cnt1 >= max_try) {
	  continue; /* continue the main loop */
	}
      }
      if(posed2!=NULL && p>0) {
	free(posed2);
	posed2=NULL;
      }
      poly[3*n] = pos[0]; poly[3*n+1] = pos[1]; poly[3*n+2] = pos[2];
      max_cnt=imax(cnt1, max_cnt);
      POLY_TRACE(printf("M"); fflush(NULL));
      //POLY_TRACE(/* printf("placed Monomer 1 at (%f,%f,%f)\n",pos[0],pos[1],pos[2]) */);
      
      /* place remaining monomers */
      for (n=2; n<MPC; n++) { 
	if (angle2 > -1.0) {
	  if(n==2) { /* if the 2nd angle is set, construct preceding monomer 
		       with resulting plane perpendicular on the xy-plane */
	    poy[0]=2*poz[0]-pos[0];
	    poy[1]=2*poz[1]-pos[1];
	    if(pos[2]==poz[2])
	      poy[2]=poz[2]+1;
	    else
	      poy[2]=poz[2];
	  } else { 
	    /* save 3rd last monomer */
	    pox[0]=poy[0]; pox[1]=poy[1]; pox[2]=poy[2]; 
	  }
	}
	if (angle > -1.0) { 
	  /* save one but last monomer */
	  poy[0]=poz[0]; poy[1]=poz[1]; poy[2]=poz[2]; 
	}
	/* save last monomer */
	poz[0]=pos[0]; poz[1]=pos[1]; poz[2]=pos[2];

	if(angle > -1.0){
	  a[0]=poy[0]-poz[0];
	  a[1]=poy[1]-poz[1];
	  a[2]=poy[2]-poz[2];

	  b[0]=pox[0]-poy[0];
	  b[1]=pox[1]-poy[1];
	  b[2]=pox[2]-poy[2];

	  vector_product(a,b,c);	  
	}

	for (cnt1=0; cnt1<max_try; cnt1++) {
	  if(angle > -1.0) {
	    if (sqrlen(c) < ROUND_ERROR_PREC) {
	      fprintf(stderr, "WARNING: rotation axis is 0,0,0, check the angles given to the polymer command\n");
	      c[0] = 1; c[1] = 0; c[2] = 0;
	    }
	    if(angle2 > -1.0 && n>2) {
	      vec_rotate(a,angle2,c,d);
	    } else {
	      phi = 2.0*PI*d_random();
	      vec_rotate(a,phi,c,d);
	    }

	    vec_rotate(d,angle,a,b);

	    pos[0] = poz[0] + b[0];
	    pos[1] = poz[1] + b[1];
	    pos[2] = poz[2] + b[2];

	  } else {
            zz     = (2.0*d_random()-1.0)*bond_length;
            rr     = sqrt(SQR(bond_length)-SQR(zz));
            phi    = 2.0*PI*d_random();
            pos[0] = poz[0]+rr*cos(phi);
            pos[1] = poz[1]+rr*sin(phi);
            pos[2] = poz[2]+zz;
	  }
	  
	  //POLY_TRACE(/* printf("a=(%f,%f,%f) absa=%f M=(%f,%f,%f) c=(%f,%f,%f) absMc=%f a*c=%f)\n",a[0],a[1],a[2],sqrt(SQR(a[0])+SQR(a[1])+SQR(a[2])),M[0],M[1],M[2],c[0],c[1],c[2],sqrt(SQR(M[0]+c[0])+SQR(M[1]+c[1])+SQR(M[2]+c[2])),a[0]*c[0]+a[1]*c[1]+a[2]*c[2]) */);
	  //POLY_TRACE(/* printf("placed Monomer %d at (%f,%f,%f)\n",n,pos[0],pos[1],pos[2]) */);

#ifdef CONSTRAINTS
	  if(constr==0 || constraint_collision(pos,poly+3*(n-1))==0){
#endif
	    if (mode==1 || collision(pos, shield, n, poly)==0) break;
	    if (mode==0) { cnt1 = -2; break; }
#ifdef CONSTRAINTS
	  }
#endif
	  POLY_TRACE(printf("m"); fflush(NULL));
	}
	if (cnt1 >= max_try) {
	  fprintf(stderr,"\nWarning! Attempt #%d to build polymer %d failed after %d unsuccessful trials to place monomer %d!\n",cnt2+1,p,cnt1,n);
	  fprintf(stderr,"         Retrying by re-setting the start-monomer of current chain...\n");
	}
	if (cnt1 == -2 || cnt1 >= max_try) {
	  n=0; break;
	}
	poly[3*n] = pos[0]; poly[3*n+1] = pos[1]; poly[3*n+2] = pos[2];
	max_cnt=imax(cnt1, max_cnt);
	POLY_TRACE(printf("M"); fflush(NULL));
      }
      if (n>0) break;
    } /* cnt2 */
    POLY_TRACE(printf(" %d/%d->%d \n",cnt1,cnt2,max_cnt));
    if (cnt2 >= max_try) { free(poly); return(-2); } else max_cnt = imax(max_cnt,imax(cnt1,cnt2));

    /* actually creating current polymer in ESPResSo */
    for (n=0; n<MPC; n++) {
      
      pos[0] = poly[3*n]; pos[1] = poly[3*n+1]; pos[2] = poly[3*n+2];
      if (place_particle(part_id, pos)==TCL_ERROR ||
	  (set_particle_q(part_id, ((n % cM_dist==0) ? val_cM : 0.0) )==TCL_ERROR) ||
	  (set_particle_type(part_id, ((n % cM_dist==0) ? type_cM : type_nM) )==TCL_ERROR))
	{ free(poly); return (-3); }
      
      if(n>=bond_size){
	bond[1] = part_id - bond_size;
	for(i=2;i<=bond_size;i++){
	  bond[i] = part_id - bond_size + i;
	}
	if(change_particle_bond(part_id-bond_size+1, bond, 0)==TCL_ERROR)
	  { free(poly); return (-3); }
      }
      part_id++;
      //POLY_TRACE(/* printf("placed Monomer %d at (%f,%f,%f)\n",n,pos[0],pos[1],pos[2]) */);
    }
  }
  free(poly);
  return(imax(max_cnt,cnt2));
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
  return mpi_gather_runtime_errors(interp, tmp_try);
}



int counterionsC(int N_CI, int part_id, int mode, double shield, int max_try, double val_CI, int type_CI) {
  int n, cnt1,max_cnt;
  double pos[3];

  cnt1 = max_cnt = 0;
  for (n=0; n<N_CI; n++) {
    for (cnt1=0; cnt1<max_try; cnt1++) {
      pos[0]=box_l[0]*d_random();
      pos[1]=box_l[1]*d_random();
      pos[2]=box_l[2]*d_random();
      if ((mode!=0) || (collision(pos, shield, 0, NULL)==0)) break;
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
  if (cnt1 >= max_try) return(-1);
  return(imax(max_cnt,cnt1));
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
  return mpi_gather_runtime_errors(interp, tmp_try);
}



int saltC(int N_pS, int N_nS, int part_id, int mode, double shield, int max_try, double val_pS, double val_nS, int type_pS, int type_nS, double rad) {
  int n, cnt1,max_cnt;
  double pos[3], dis2;

  cnt1 = max_cnt = 0;

  /* Place positive salt ions */
  for (n=0; n<N_pS; n++) {
    for (cnt1=0; cnt1<max_try; cnt1++) {
      if (rad > 0.) {
        pos[0]=rad*(2.*d_random()-1.);
        pos[1]=rad*(2.*d_random()-1.);
        pos[2]=rad*(2.*d_random()-1.);
        dis2 = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
        pos[0] += box_l[0]*0.5;
        pos[1] += box_l[1]*0.5;
        pos[2] += box_l[2]*0.5;
        if (((mode!=0) || (collision(pos, shield, 0, NULL)==0)) && (dis2 < (rad * rad))) break;
      } else {
        pos[0]=box_l[0]*d_random();
        pos[1]=box_l[1]*d_random();
        pos[2]=box_l[2]*d_random();
        if ((mode!=0) || (collision(pos, shield, 0, NULL)==0)) break;
      }
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
      if (rad > 0.) {
        pos[0]=rad*(2.*d_random()-1.);
        pos[1]=rad*(2.*d_random()-1.);
        pos[2]=rad*(2.*d_random()-1.);
        dis2 = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
        pos[0] += box_l[0]*0.5;
        pos[1] += box_l[1]*0.5;
        pos[2] += box_l[2]*0.5;
        if (((mode!=0) || (collision(pos, shield, 0, NULL)==0)) && (dis2 < (rad * rad))) break;
      } else {
        pos[0]=box_l[0]*d_random();
        pos[1]=box_l[1]*d_random();
        pos[2]=box_l[2]*d_random();
        if ((mode!=0) || (collision(pos, shield, 0, NULL)==0)) break;
      }
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
  if (cnt1 >= max_try) return(-2);
  return(imax(max_cnt,cnt1));
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
  return mpi_gather_runtime_errors(interp, TCL_OK);
}



double velocitiesC(double v_max, int part_id, int N_T) {
  double v[3], v_av[3];
  int i;

  v_av[0] = v_av[1] = v_av[2] = 0.0;
  for (i=part_id; i < part_id+N_T; i++) {
    do {
      v[0] = v_max * 2.*(d_random()-.5);
      v[1] = v_max * 2.*(d_random()-.5);
      v[2] = v_max * 2.*(d_random()-.5);
    } while ( sqrt(SQR(v[0])+SQR(v[1])+SQR(v[2])) > v_max);
    v_av[0]+=v[0]; v_av[1]+=v[1]; v_av[2]+=v[2];
    if (set_particle_v(i, v)==TCL_ERROR) {
      fprintf(stderr, "INTERNAL ERROR: failed upon setting one of the velocities in Espresso (current average: %f)!\n",
	      sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2]))); 
      fprintf(stderr, "Aborting...\n"); errexit();
    }
  }
  return ( sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2])) );
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
  return mpi_gather_runtime_errors(interp, TCL_OK);
}

double maxwell_velocitiesC(int part_id, int N_T) {
  double v[3], v_av[3],uniran[2];
  int i;
  int flag=1;
  uniran[0]=d_random();
  uniran[1]=d_random();
  v_av[0] = v_av[1] = v_av[2] = 0.0;
  for (i=part_id; i < part_id+N_T; i++) {
    if(flag == 1 ) {
      v[0] = pow((-2. * log(uniran[0])),0.5) * cos (2. * PI * uniran[1]) * time_step;
      v[1] = pow((-2. * log(uniran[1])),0.5) * sin (2. * PI * uniran[0]) * time_step;
      uniran[0]=d_random();
      uniran[1]=d_random();
      v[2] = pow((-2. * log(uniran[0])),0.5) * cos (2. * PI * uniran[1]) * time_step;
      flag = 0;
    } else {
      v[0] = pow((-2. * log(uniran[1])),0.5) * sin (2. * PI * uniran[0]) * time_step;
      uniran[0]=d_random();
      uniran[1]=d_random();
      v[1] = pow((-2. * log(uniran[0])),0.5) * cos (2. * PI * uniran[1]) * time_step;
      v[2] = pow((-2. * log(uniran[1])),0.5) * sin (2. * PI * uniran[0]) * time_step;
      flag = 1;      
    }
    //printf("%f \n %f \n %f \n",v[0],v[1],v[2]);
    v_av[0]+=v[0]; v_av[1]+=v[1]; v_av[2]+=v[2];
    if (set_particle_v(i, v)==TCL_ERROR) {
      fprintf(stderr, "INTERNAL ERROR: failed upon setting one of the velocities in Espresso (current average: %f)!\n",sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2]))); 
      fprintf(stderr, "Aborting...\n"); errexit();
    }
  }
  return ( sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2])) );
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
	if (ARG_IS_I(i+1, type_bond)) { 
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
  return mpi_gather_runtime_errors(interp, tmp_try);
}



int collectBonds(int mode, int part_id, int N_P, int MPC, int type_bond, int **bond_out, int ***bonds_out) {
  int i,j,k,ii,size, *bond=NULL, **bonds=NULL;

  /* Get particle and bonding informations. */
  IntList *bl;
  Particle *prt, *sorted;
  bl  = malloc(1*sizeof(IntList));
  prt = malloc(n_total_particles*sizeof(Particle));
  mpi_get_particles(prt, bl); 

  /* Sort the received informations. */
  sorted = malloc(n_total_particles*sizeof(Particle));
  for(i = 0; i < n_total_particles; i++)
    memcpy(&sorted[prt[i].p.identity], &prt[i], sizeof(Particle));
  free(prt);
  prt = sorted;
  
  if (mode == 1) {
    /* Find all the bonds leading to and from the ending monomers of the chains. */
    bond  = malloc(2*N_P*sizeof(int));      bonds   = malloc(2*N_P*sizeof(int *));
    for (i=0; i < 2*N_P; i++) { bond[i]=0;  bonds[i]= malloc(1*sizeof(int)); }
    for (k=part_id; k < N_P*MPC + part_id; k++) {
      i=0;
      while(i < prt[k].bl.n) {
	size = bonded_ia_params[prt[k].bl.e[i]].num;
	if (prt[k].bl.e[i++] == type_bond) {
	  for(j=0; j<size; j++) {
	    if ((prt[k].p.identity % MPC == 0) || ( (prt[k].p.identity+1) % MPC == 0)) {
	      ii = prt[k].p.identity%MPC ? 2*(prt[k].p.identity+1)/MPC-1 : 2*prt[k].p.identity/MPC;
	      bonds[i] = realloc(bonds[i], (bond[i]+1)*sizeof(int));
	      bonds[ii][bond[ii]++] = prt[k].bl.e[i];
	    }
	    else if ((prt[k].bl.e[i] % MPC == 0) || ( (prt[k].bl.e[i]+1) % MPC == 0)) {
	      ii = prt[k].bl.e[i]%MPC ? 2*(prt[k].bl.e[i]+1)/MPC-1 : 2*prt[k].bl.e[i]/MPC;
	      bonds[i] = realloc(bonds[i], (bond[i]+1)*sizeof(int));
	      bonds[ii][bond[ii]++] = prt[k].p.identity;
	    }
	    i++;
	  }
	}
	else i += size;
      }
    }
    POLY_TRACE(for (i=0; i < 2*N_P; i++) {
      printf("(%d) %d:\t",i,i%2 ? (i+1)*MPC/2-1 : i*MPC/2); if(bond[i]>0) for(j=0;j<bond[i];j++) printf("%d ",bonds[i][j]); printf("\t=%d\n",bond[i]);
    });
  }
  else if (mode == 2) {
    /* Find all the bonds leading to and from each monomer. */
    bond  = malloc(N_P*MPC*sizeof(int));                bonds   = malloc(N_P*MPC*sizeof(int *));
    for (i=0; i < N_P*MPC + part_id; i++) { bond[i]=0;  bonds[i]= malloc(1*sizeof(int)); }
    for (k=part_id; k < N_P*MPC + part_id; k++) {
      i=0;
      while(i < prt[k].bl.n) {
	size = bonded_ia_params[prt[k].bl.e[i]].num;
	if (prt[k].bl.e[i++] == type_bond) {
	  for(j=0; j<size; j++) {
	    ii = prt[k].bl.e[i];
	    bonds[k] = (int *) realloc(bonds[k], (bond[k]+1)*sizeof(int));
	    bonds[k][bond[k]++] = ii;
	    bonds[ii] = (int *) realloc(bonds[ii], (bond[ii]+1)*sizeof(int));
	    bonds[ii][bond[ii]++] = k;
	    i++;
	  }
	}
	else i += size;
      }
    }
    POLY_TRACE(for (i=0; i < N_P*MPC + part_id; i++) { 
      printf("%d:\t",i); if(bond[i]>0) for(j=0;j<bond[i];j++) printf("%d ",bonds[i][j]); printf("\t=%d\n",bond[i]); 
    });
  }
  else {
    fprintf(stderr, "Unknown mode %d requested!\nAborting...\n",mode); fflush(NULL); return(-2);
  }
  free(prt); realloc_intlist(bl, 0);
  *bond_out  = bond;
  *bonds_out = bonds;
  return(0);
}



int crosslinkC(int N_P, int MPC, int part_id, double r_catch, int link_dist, int chain_dist, int type_bond, int max_try) {
  int i,j,k,ii,size, bondN[2], *bond, **bonds, *link, **links, *cross, crossL;

  /* Find all the bonds leading to and from each monomer. */
  if (collectBonds(2, part_id, N_P, MPC, type_bond, &bond, &bonds)) return(-2);
  POLY_TRACE(for (i=0; i < N_P*MPC + part_id; i++) { 
    printf("%d:\t",i); if(bond[i]>0) for(j=0;j<bond[i];j++) printf("%d ",bonds[i][j]); printf("\t=%d\n",bond[i]); 
  });
  
  /* Find all possible binding partners in the neighbourhood of the unconnected ending monomers. */
  link  = malloc(2*N_P*sizeof(int));       
  links = malloc(2*N_P*sizeof(int *));
  for (i=0; i < N_P; i++) {
    for (k=0; k<2; k++) {
      if (bond[i*MPC+k*(MPC-1)] == 1) {
	links[2*i+k] = malloc(n_total_particles*sizeof(int));
	link[2*i+k] = mindist3(i*MPC+k*(MPC-1)+part_id, r_catch, links[2*i+k]);
	links[2*i+k] = realloc(links[2*i+k],link[2*i+k]*sizeof(int));
      }
      else if (bond[i*MPC+k*(MPC-1)] == 2) link[2*i+k] = -1;  /* Note that links[2*i+k] will not be malloc()ed now (taken care of at end)!!! */
      else { fprintf(stderr,"Runaway end-monomer %d detected (has %d bonds)!\nAborting...\n", i*N_P+k*(MPC-1)+part_id, bond[i*MPC+k*(MPC-1)]); 
             fflush(NULL); return(-2); }
      POLY_TRACE(printf("%d: ",i*MPC+k*(MPC-1)+part_id); 
		 for (j=0; j<link[2*i+k]; j++) printf("%d ",links[2*i+k][j]); printf("\t=%d\n",link[2*i+k]); fflush(NULL) );
    }
  }

  /* Throw out all the monomers which are ends, which are too close to the ending monomers on the same chain, or which are no monomers at all. */
  for (i=0; i < N_P; i++) {
    for (k=0; k<2; k++) {
      size = 0;  ii = i*MPC + k*(MPC-1) + part_id;
      if (link[2*i+k] >= 0) {
	for (j=0; j < link[2*i+k]; j++) {                     /* only monomers && ((same chain, but sufficiently far away) || (different chain)) */
	  if ( (links[2*i+k][j] < N_P*MPC+part_id) && ( ((abs(links[2*i+k][j] - ii) > chain_dist) || (abs(links[2*i+k][j]-i*MPC) > (1.*MPC))) ) )
	    if ((links[2*i+k][j] % MPC != 0) && ((links[2*i+k][j]+1) % MPC != 0)) links[2*i+k][size++] = links[2*i+k][j];    /* no ends accepted */
	}
	link[2*i+k]  = size; 
	links[2*i+k] = realloc(links[2*i+k],link[2*i+k]*sizeof(int));
      }
      POLY_TRACE(printf("%d: ",ii); for (j=0; j<link[2*i+k]; j++) printf("%d ",links[2*i+k][j]); printf("\t=%d\n",link[2*i+k]); fflush(NULL) );
    }
  }

  /* Randomly choose a partner (if not available -> '-1') for each polymer chain's end if it's not already been crosslinked (-> '-2'). */
  cross = malloc(2*N_P*sizeof(int)); crossL = 0;
  for (i=0; i < 2*N_P; i++) 
    if (link[i] > 0) { cross[i] = links[i][(int)dround(d_random()*(link[i]-1))]; crossL++; }  else { cross[i] = -1+link[i]; crossL -= link[i]; }
  POLY_TRACE(for (i=0; i < 2*N_P; i++) printf("%d -> %d \t", i%2 ? (i+1)*MPC/2-1 : i*MPC/2, cross[i]); printf("=> %d\n",crossL); fflush(NULL) );

  /* Remove partners (-> '-3') if they are less than link_dist apart and retry. */
  k = 0; ii = 1;
  while ((k < max_try) && (ii > 0)) {
    POLY_TRACE(printf("Check #%d: ",k));
    for (i=0; i < 2*N_P; i++) {
      if (cross[i] >= 0) {
	for (j=0; j < 2*N_P; j++) {       /* In the neighbourhood of each partner shall be no future crosslinks (preventing stiffness). */
	  if ((j != i) && (cross[j] >=0) && (abs(cross[j]-cross[i]) < link_dist)) {
	    cross[i] = -3; cross[j] = -3; crossL -= 2; POLY_TRACE(printf("%d->%d! ",i,j)); break;
	  }
	}
	if (cross[i] == -3) continue;     /* Partners shall not be too close to the chain's ends (because these will be crosslinked at some point). */
	if ((cross[i] % MPC < link_dist) || (cross[i] % MPC >= MPC-link_dist)) {
	  cross[i] = -3; crossL--; POLY_TRACE(printf("%d->end! ",i)); }      
	else {                            /* In the neighbourhood of each partner there shall be no other crosslinks (preventing stiffness). */
	  for (j = cross[i]-link_dist+1; j < cross[i]+link_dist-1; j++) {
	    if ((j % MPC == 0) || ((j+1) % MPC == 0)) size = 1; else size = 2;
	    if ((bond[j] > size) && (j - floor(i/2.)*MPC < MPC)) {
	      cross[i] = -3; crossL--; POLY_TRACE(printf("%d->link! ",i)); break; 
	    }
	  }
	}
      }
    }   POLY_TRACE(printf("complete => %d CL left; ",crossL));
    if (k == max_try-1) break; else ii = 0;  /* Get out if max_try is about to be reached, preventing dangling unchecked bond suggestions. */
    if (crossL < 2*N_P) {
      for (i=0; i < 2*N_P; i++) {         /* If crosslinks violated the rules & had to be removed, create new ones now. */
	if (cross[i] == -3) {
	  ii++;
	  if (link[i] > 0) { cross[i] = links[i][(int)dround(d_random()*(link[i]-1))]; crossL++; } else { return(-2); }
	}
      }
    }   POLY_TRACE(printf("+ %d new = %d CL.\n",ii,crossL));
    if (ii > 0) k++;
  }
  POLY_TRACE(for (i=0; i < 2*N_P; i++) printf("%d -> %d \t", i%2 ? (i+1)*MPC/2-1 : i*MPC/2, cross[i]); printf("=> %d\n",crossL); fflush(NULL) );

  /* Submit all lawful partners as new bonds to Espresso (observing that bonds are stored with the higher-ID particle only). */
  if (k >= max_try) return(-1);

  {
    size = 0;
    for (i=0; i < N_P; i++) {
      if (cross[2*i] >= 0 ) {
	bondN[0] = type_bond; bondN[1] = i*MPC + part_id; size++;
	if (change_particle_bond(cross[2*i], bondN, 0)==TCL_ERROR) return (-3);
      }
      if (cross[2*i+1] >= 0) {
	bondN[0] = type_bond; bondN[1] = cross[2*i+1]; size++;
	if (change_particle_bond(i*MPC+(MPC-1) + part_id, bondN, 0)==TCL_ERROR) return (-3);
      }
      free(bonds[2*i]);    if (link[2*i]   >= 0) free(links[2*i]);    /* else crash(); because links[2*i]   has never been malloc()ed then */
      free(bonds[2*i+1]);  if (link[2*i+1] >= 0) free(links[2*i+1]);  /* else crash(); because links[2*i+1] has never been malloc()ed then */
    }
    free(bond); free(bonds); free(link); free(links); free(cross);
    POLY_TRACE(printf("Created %d new bonds; now %d ends are crosslinked!\n", size, crossL));
    return(crossL); 
  }
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
  return mpi_gather_runtime_errors(interp, tmp_try);
}


int diamondC(double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, int nonet) {
  int i,j,k, part_id, bond[2], type_bond=0,type_node=0,type_cM=1,type_nM=1, type_CI=2;
  double pos[3], off = bond_length/sqrt(3);
  double dnodes[8][3]  = {{0,0,0}, {1,1,1}, {2,2,0}, {0,2,2}, {2,0,2}, {3,3,1}, {1,3,3}, {3,1,3}};
  int    dchain[16][5] = {{0,1, +1,+1,+1}, {1,2, +1,+1,-1}, {1,3, -1,+1,+1}, {1,4, +1,-1,+1},
			  {2,5, +1,+1,+1}, {3,6, +1,+1,+1}, {4,7, +1,+1,+1}, {5,0, +1,+1,-1},
			  {5,3, +1,-1,+1}, {5,4, -1,+1,+1}, {6,0, -1,+1,+1}, {6,2, +1,-1,+1},
			  {6,4, +1,+1,-1}, {7,0, +1,-1,+1}, {7,2, -1,+1,+1}, {7,3, +1,+1,-1}};

  part_id = 0;
  /* place 8 tetra-functional nodes */
  for(i=0; i<8; i++) {
    for(j=0; j<3; j++) { 
      dnodes[i][j] *= a/4.; pos[j] = dnodes[i][j]; 
    }
    if (place_particle(part_id, pos)==TCL_ERROR) return (-3);
    if (set_particle_q(part_id, val_nodes)==TCL_ERROR) return (-3);
    if (set_particle_type(part_id, type_node)==TCL_ERROR) return (-3);
    part_id++;
  }

  /* place intermediate monomers on chains connecting the nodes */
  for(i=0; i<2*8; i++) {
    for(k=1; k<=MPC; k++) {
      for(j=0; j<3; j++) pos[j] = dnodes[dchain[i][0]][j] + k*dchain[i][2+j]*off;
      if (place_particle(part_id, pos)==TCL_ERROR) return (-3);
      if (set_particle_q(part_id, (k % cM_dist==0) ? val_cM : 0.0)==TCL_ERROR) return (-3);
      if (set_particle_type(part_id, (k % cM_dist==0) ? type_cM : type_nM)==TCL_ERROR) return (-3);
      bond[0] = type_bond; 
      if(k==1) { 
	if(nonet!=1) { bond[1] = dchain[i][0]; if (change_particle_bond(part_id, bond, 0)==TCL_ERROR) return (-3); } }
      else { 
	bond[1] = part_id-1; if (change_particle_bond(part_id, bond, 0)==TCL_ERROR) return (-3); }
      if((k==MPC)&&(nonet!=1)) { 
	bond[1] = dchain[i][1];
	if (change_particle_bond(part_id, bond, 0)==TCL_ERROR) return (-3);
      }
      part_id++;
    }
  }

  /* place counterions (if any) */
  if(N_CI > 0) counterionsC(N_CI, part_id, 1, 0.0, 30000, val_CI, type_CI);
  
  return(0);
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
  return mpi_gather_runtime_errors(interp, tmp_try);
}


int icosaederC(double ico_a, int MPC, int N_CI, double val_cM, double val_CI, int cM_dist) {
  int i,j,k,l, part_id, bond[2], type_bond=0,type_cM=0,type_nM=1, type_CI=2;
  double pos[3],pos_shift[3], vec[3],e_vec[3],vec_l, bond_length=(2*ico_a/3.)/(1.*MPC);
  double ico_g=ico_a*(1+sqrt(5))/2.0, shift=0.0;
  double ico_coord[12][3] = {{0,+ico_a,+ico_g}, {0,+ico_a,-ico_g}, {0,-ico_a,+ico_g}, {0,-ico_a,-ico_g},
			     {+ico_a,+ico_g,0}, {+ico_a,-ico_g,0}, {-ico_a,+ico_g,0}, {-ico_a,-ico_g,0},
			     {+ico_g,0,+ico_a}, {-ico_g,0,+ico_a}, {+ico_g,0,-ico_a}, {-ico_g,0,-ico_a}};
  int    ico_NN[12][5]    = {{2,8,4,6, 9}, {3,10,4,6,11}, {0,8,5,7, 9}, {1,10,5,7,11}, {0,6,1,10,8}, {2,7,3,10,8},
			     {0,4,1,11,9}, {2,5,3,11, 9}, {0,2,5,10,4}, {0,2,7,11, 6}, {1,3,5,8, 4}, {1,3,7,9, 6}};
  int    ico_ind[12][10];

  /* make sure that the edges in ico_NN are sorted such that NearestNeighbours are next to each other */
  /* int    ico_NN[12][5]    = {{2,4,6,8, 9}, {3,4,6,10,11}, {0,5,7,8, 9}, {1,5,7,10,11}, {0,1,6,8,10}, {2,3,7,8,10}, 
			     {0,1,4,9,11}, {2,3,5, 9,11}, {0,2,4,5,10}, {0,2,6, 7,11}, {1,3,4,5, 8}, {1,3,6,7, 9}};
  for(i=0; i<12; i++) {
    printf("%d: { ",i); 
    for(j=0; j<5; j++) printf("%d ",ico_NN[i][j]);
    printf("} -> ");
    for(j=0; j<5; j++) 
      for(l=0; l<5; l++) 
	if(j!=l) 
	  for(k=0; k<5; k++) 
	    if(ico_NN[ico_NN[i][j]][k]==ico_NN[i][l]) printf("%d = %d (@%d)  ",ico_NN[i][j],ico_NN[i][l],k);
    printf("\n");
  } */

  /* shift coordinates to not be centered around zero but rather be positive */
  if(ico_a > ico_g) shift=ico_a; else shift=ico_g;

  /* create fulleren & soccer-ball */
  part_id = 0;
  for(i=0; i<12; i++) {
    for(j=0; j<5; j++) {
      /* place chains along the 5 edges around each of the 12 icosaeder's vertices */
      if(j < 4) for(l=0; l<3; l++) vec[l] = (ico_coord[ico_NN[i][j+1]][l] - ico_coord[ico_NN[i][j]][l])/3.;
      else      for(l=0; l<3; l++) vec[l] = (ico_coord[ico_NN[i][0]][l]   - ico_coord[ico_NN[i][4]][l])/3.;
      vec_l = sqrt(SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]));
      for(l=0; l<3; l++) e_vec[l] = vec[l]/vec_l;

      ico_ind[i][j] = part_id; bond_length = vec_l/(1.*MPC);
      for(l=0; l<3; l++) pos[l] = ico_coord[i][l] + (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l])/3.;
      for(k=0; k<MPC; k++) {
	for(l=0; l<3; l++) pos_shift[l] = pos[l] + shift;
	if (place_particle(part_id, pos_shift)==TCL_ERROR) return (-3);
	if (set_particle_q(part_id, val_cM)==TCL_ERROR) return (-3);
	if (set_particle_type(part_id, type_cM)==TCL_ERROR) return (-3);
	bond[0] = type_bond;
	if (k > 0) {
	  bond[1] = part_id-1; if (change_particle_bond(part_id, bond, 0)==TCL_ERROR) return (-3); 
	}
	part_id++;
	for(l=0; l<3; l++) pos[l] += bond_length*e_vec[l];
      }

      /* place chains along the 5 edges on the middle third of the connection between two NN vertices */
      if(i < ico_NN[i][j]) {
	for(l=0; l<3; l++) vec[l] = (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l])/3.;
	vec_l = sqrt(SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]));
	for(l=0; l<3; l++) e_vec[l] = vec[l]/vec_l;
	
	ico_ind[i][j+5] = part_id; bond_length = vec_l/(1.*MPC);
	for(l=0; l<3; l++) pos[l] = ico_coord[i][l] + (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l])/3. + bond_length*e_vec[l];
	for(k=1; k<MPC; k++) {
	  for(l=0; l<3; l++) pos_shift[l] = pos[l] + shift;
	  if (place_particle(part_id, pos_shift)==TCL_ERROR) return (-3);
	  if (set_particle_q(part_id, 0.0)==TCL_ERROR) return (-3);
	  if (set_particle_type(part_id, type_nM)==TCL_ERROR) return (-3);
	  bond[0] = type_bond;
	  if (k > 1) {
	    bond[1] = part_id-1; if (change_particle_bond(part_id, bond, 0)==TCL_ERROR) return (-3); }
	  else {
	    bond[1] = ico_ind[i][j]; if (change_particle_bond(part_id, bond, 0)==TCL_ERROR) return (-3); }
	  part_id++;
	  for(l=0; l<3; l++) pos[l] += bond_length*e_vec[l];
	}
      } 
    }

    for(j=0; j<5; j++) {
      /* add bonds between the edges around the vertices */
      bond[0] = type_bond;
      //      if(j>0) bond[1] = ico_ind[i][j-1] + (MPC-1); else bond[1] = ico_ind[i][4] + (MPC-1);
      if(j>0) bond[1] = ico_ind[i][j-1] + (MPC-1); else if(MPC>0) bond[1] = ico_ind[i][4] + (MPC-1); else bond[1] = ico_ind[i][4];
      if (change_particle_bond(ico_ind[i][j], bond, 0)==TCL_ERROR) return (-2);

      /* connect loose edges around vertices with chains along the middle third already created earlier */
      if(i > ico_NN[i][j]) {
	bond[0] = type_bond;
	for(l=0; l<5; l++) if(ico_NN[ico_NN[i][j]][l] == i) break;
	if(l==5) {
	  fprintf(stderr, "INTERNAL ERROR: Couldn't find my neighbouring edge upon creating the icosaeder!\n");
	  errexit();
	}
	bond[1] = ico_ind[ico_NN[i][j]][l+5] + (MPC-2);
	if (change_particle_bond(ico_ind[i][j], bond, 0)==TCL_ERROR) return (-1);
      }
    }
  }

  /* place counterions (if any) */
  if(N_CI > 0) counterionsC(N_CI, part_id, 1, 0.0, 30000, val_CI, type_CI);
  
  return(0);
}

