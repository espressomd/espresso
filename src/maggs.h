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

/** \file maggs.h
 *  Maxwell Equations Molecular Dynamics (MEMD) method for electrostatic
 *  interactions.
 *
 *  We use a local update scheme to propagate artificial B-fields on a
 *  lattice in the system. In principal, the algorithm simulates full
 *  electrodynamics, but with a tunable speed of light.
 *
 *  The method is very usable for large particle numbers or highly
 *  parallel architectures, since it is local and scales linearly.
 *  It is not suited for high-precision calculation of forces, since
 *  the simple interpolation scheme produces errors in the order of
 *  10^-5 in the force.
 *
 *  The chosen mesh should roughly be of the size of the particles.
 *
 *  Further reading on the algorithm:
 *  <ul>
 *  <li> I. Pasichnyk and B. Dunweg, Coulomb interaction via local dynamics: a molecular-dynamics algorithm. J. Phys: Condens. Matter, 16 ,p. 1399-4020, (2004).
 *  </ul>
 *  
 */


/* protect header file: */
#ifndef MAGGS_H
#define MAGGS_H

/* maggs structure. Contains global system information. */
typedef struct {
	double f_mass;         /* = 1/c^2    speed of light parameter. */
	double invsqrt_f_mass; /* inverse of square root of f_mass. */
	double prefactor;      /* prefactor to convert field to force. */
	double pref2;          /* prefactor / (4*pi) */
	double bjerrum;        /* bjerrum length of the system. */
	int    mesh;           /* mesh size in one dimension */
	double inva;           /* = 1/a = mesh / box_length */
	double a;              /* size of mesh cube */
} MAGGS_struct;
extern MAGGS_struct maggs;


/**********************************************/
/* Functions to be called from external code: */
/**********************************************/

/* initialize, parse command and set parameters: */
void Maggs_init(); /* called from: initialize.c */
int tclcommand_inter_coulomb_parse_maggs(Tcl_Interp * interp, int argc, char ** argv);
int maggs_set_parameters(Tcl_Interp *interp, double bjerrum, double f_mass, int mesh);

/* for integration routine. */
void maggs_propagate_B_field(double dt); /* called from integrate.c (twice, with dt/2) */
void maggs_calc_forces(); /* called from forces.c */

/* get electric and magnetic energy. Called from energy.c */
double maggs_electric_energy();
double maggs_magnetic_energy();

/* Count the number of charges in the whole system. */
/* Called from communication.c */
int maggs_count_charged_particles();

/* Clean up, print results: */
void Maggs_exit();
int tclprint_to_result_Maggs(Tcl_Interp *interp);


/*************************/
/* Internal definitions: */
/*************************/

/* Define numbers for directions and dimensions: */
#define SPACE_DIM 3                 /* number of dimensions */
#define NOWHERE -1                  /* not a direction */
#define NDIRS 6				        /* number of directions */
#define XPLUS 0                     /* add for specific direction */
#define YPLUS 1
#define ZPLUS 2
#define ZMINUS 3
#define YMINUS 4
#define XMINUS 5
#define OPP_DIR(dir)	(5-(dir))	/* Opposite direction */


/* Three often used macros for looping over 3D */
#define FOR3D(dir) for(dir=0; dir<SPACE_DIM; dir++)

#define FORALL_INNER_SITES(i,j,k) \
for(i=lparams.inner_left_down[0];i<lparams.inner_up_right[0];i++) \
for(j=lparams.inner_left_down[1];j<lparams.inner_up_right[1];j++) \
for(k=lparams.inner_left_down[2];k<lparams.inner_up_right[2];k++) 

#define FORALL_SITES(i,j,k) \
for(i=0;i<lparams.dim[0];i++) \
for(j=0;j<lparams.dim[1];j++) \
for(k=0;k<lparams.dim[2];k++) 


/* from ifndef MAGGS_H */
#endif
