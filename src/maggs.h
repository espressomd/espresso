// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file maggs.h   
 *  Calculation of the electrostatic forces using Maggs method.
 */

#include "utils.h"
#include "integrate.h"
#include "ghosts.h"

#ifndef MAGGS_H
#define MAGGS_H

#ifdef ELECTROSTATICS


#define SPACE_DIM 3

/* Directions, and a macro to give the opposite direction */
/*  These must go from 0 to 5 because they will be used to index an
 array. */
/* Also define NDIRS = number of directions */
#define XPLUS 0
#define YPLUS 1
#define ZPLUS 2
#define ZMINUS 3
#define YMINUS 4
#define XMINUS 5

#define NOWHERE -1  /* not a direction */
#define OPP_DIR(dir)	(5-(dir))	/* Opposite direction */
#define NDIRS 6				/* number of directions */

/* Useful macros for looping over 3D */
#define FOR3D(dir) for(dir=0; dir<SPACE_DIM; dir++)

#define FORALL_INNER_SITES(i,j,k) \
for(i=lparam.in_ld[0];i<lparam.in_ur[0];i++) \
for(j=lparam.in_ld[1];j<lparam.in_ur[1];j++) \
for(k=lparam.in_ld[2];k<lparam.in_ur[2];k++) 

#define FORALL_SITES(i,j,k) \
for(i=0;i<lparam.dim[0];i++) \
for(j=0;j<lparam.dim[1];j++) \
for(k=0;k<lparam.dim[2];k++) 


typedef struct {
	double f_mass;
	double invsqrt_f_mass;
	double prefactor;
	double pref2;
	double bjerrum;
	int    mesh;
	double inva;
	double a;
} MAGGS_struct;
extern MAGGS_struct maggs;

void Maggs_init();
int tclcommand_inter_coulomb_parse_maggs(Tcl_Interp * interp, int argc, char ** argv);
int maggs_set_parameters(Tcl_Interp *interp, double bjerrum, double f_mass, int mesh);

void maggs_propagate_B_field(double dt);
void maggs_calc_forces();

double maggs_electric_energy();
double maggs_magnetic_energy();

int maggs_count_charged_particles();

void Maggs_exit();
int tclprint_to_result_Maggs(Tcl_Interp *interp);


#endif
#endif