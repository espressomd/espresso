// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef MAGGS_H 
#define MAGGS_H
/** \file maggs.h   
 *  Calculation of the electrostatic forces using Maggs method.
 */

#include "utils.h"
#include "integrate.h"
#include "ghosts.h"

#ifdef ELECTROSTATICS

#define SPACE_DIM 3
#define LINEAR_INTERPOLATION

/** Not an index value for arrays. */
#define NOT_AN_INDEX -99999999

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

#define PNN   8                      /*number nearest neighbor sites */

/* Useful macros */
/* Macros for looping over 3D */
#define FOR3D(dir) for(dir=0; dir<SPACE_DIM; dir++)

#define FORALL_INNER_SITES(i,j,k) \
 for(i=lparam.in_ld[0];i<lparam.in_ur[0];i++) \
  for(j=lparam.in_ld[1];j<lparam.in_ur[1];j++) \
   for(k=lparam.in_ld[2];k<lparam.in_ur[2];k++) 

#define FORALL_SITES(i,j,k) \
 for(i=0;i<lparam.dim[0];i++) \
  for(j=0;j<lparam.dim[1];j++) \
   for(k=0;k<lparam.dim[2];k++) 

#define LOOP_CUBE_VERTICES(i,j,k) \
 for(i=0;i<2;i++) \
  for(j=0;j<2;j++) \
   for(k=0;k<2;k++) 

/* "field offset" and "field pointer" */
/* used when fields are arguments to subroutines */
/* Usage:  fo = F_OFFSET( field, lattice[0] ), where "field" is the name of a field
  in lattice.
     address = F_PT( &site , fo ), where &site is the address of the
  site and fo is a field_offset.  Usually, the result will have to be
  cast to a pointer to the appropriate type. (It is naturally a char *).
*/
#define F_OFFSET(a, site) \
  ((int)(((char *)&(site. a ))-((char *)&(site)) ))
#define F_PT( site , fo )  ((char *)( site ) + (fo)) 

/********************************************************************/

typedef int      t_ivector [SPACE_DIM];
typedef double   t_dvector [SPACE_DIM];
typedef int      t_dirs    [NDIRS];

typedef struct {
  double f_mass;
  double invsqrt_f_mass;
  double prefactor;
  double pref2;
  double bjerrum;
  int    mesh;
  double fric_gamma;
  double inva;
  double a;
  int    yukawa;
  double kappa;
  double r_cut;
} MAGGS_struct;

extern MAGGS_struct maggs;

int set_maggs_params(Tcl_Interp *interp, double bjerrum, 
		     double f_mass, int mesh, double gamma, 
		     int yukawa, double mue, double r_cut);

/** Initialize all structures, parameters and arrays needed for the 
 *  Maggs algorithm.
 */
void Maggs_init();

/** Calculate number of charged particles.
 */
int maggs_count_charged_particles();

/** Update B field.
 */
void propagate_B_field(double dt);

/** \name Exported Functions */
/************************************************************/
/*@{*/

double maggs_magnetic_energy();
/** Calculate electric energy of the field.
 */
double maggs_electric_energy();
/** Functions for integration of additional fields
 */
void maggs_propagate_psi_vel(double dt);
void maggs_propagate_psi_vel_pos(double dt);
void maggs_calc_psi_forces();
void maggs_calc_e_forces();

void check_gauss_law();
void check_charge_conserv();

double calc_gauss_res();

/** implementation of analyze energy */
int parse_and_print_gauss_res(Tcl_Interp *interp, int argc, char **argv);

/// print the maggs parameters to the interpreters result
int printMaggsToResult(Tcl_Interp *interp);

/// parse the basic maggs parameters
int inter_parse_maggs(Tcl_Interp * interp, int argc, char ** argv);

///
int Maggs_sanity_checks();

/** Calculate real space contribution of Yukawa pair forces. */
MDINLINE void add_maggs_yukawa_pair_force(Particle *p1, Particle *p2,
					  double *d,double dist2,double dist,double force[3])
{
  int j;
  double fac, temp;

  if(dist < maggs.r_cut) {
    
    temp = maggs.kappa*dist;
    fac = maggs.pref2 * p1->p.q * p2->p.q * 
      exp(-temp) * (1. + temp) / (dist2*dist);

    FOR3D(j) force[j] += fac * d[j];

    MAGGS_TRACE(
      fprintf(stderr, "%d: Yukawa: pair (%d-%d) dist=%.3f, q1=%.1f, q2=%.1f, force+-: (%.3e,%.3e,%.3e)\n",
	      this_node,p1->p.identity,p2->p.identity,dist, p1->p.q, p2->p.q, fac*d[0],fac*d[1],fac*d[2]);
      );
  }
}

/*@}*/
#endif

#endif
