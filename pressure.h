// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef PRESSURE_H
#define PRESSURE_H


#include "config.h"
#include "integrate.h"
#include "statistics.h"
#include "thermostat.h"
#include "communication.h"

/************************************************
 * data types
 ************************************************/

/** Structure to hold all variables related to the isotropic NpT-integration scheme. */
typedef struct { 
  /** mass of a virtual piston representing the shaken box */
  double piston;
  /** inverse of \ref piston */
  double inv_piston;
  /** isotropic volume.  Note that we use the term volume throughout
      although for a 2d or 1d system we mean Area and Length
      respectively */
  double volume;

  /** desired pressure to which the algorithm strives to */
  double p_ext;
  /** instantaneous pressure the system currently has */
  double p_inst;
  /** instantaneous pressure averaged over current integration cycle */
  double p_inst_av;
  /** difference between \ref p_ext and \ref p_inst */
  double p_diff;
  /** virial (short-range) components of \ref p_inst */
  double p_vir[3];
  /** ideal gas components of \ref p_inst, derived from the velocities */
  double p_vel[3];
  /** flag which indicates if \ref p_vel may (0) or may not (1) be used
      in offline pressure calculations such as 'analyze p_inst' */ 
  int invalidate_p_vel;
  /** geometry information for the npt integrator.  Holds the vector
      <dir, dir ,dir> where a positive value for dir indicates that
      box movement is allowed in that direction. To check whether a
      given direction is turned on use bitwise comparison with \ref
      nptgeom_dir */
  int geometry;
  /** bitwise comparison values corresponding to different directions*/
  int nptgeom_dir[3];
  /** The number of dimensions in which npt boxlength motion is coupled to particles */
  int dimension;
  /** Set this flag if you want all box dimensions to be identical. Needed for electrostatics.  
      If the value of \ref dimension is less than 3 then box length motion in one or more
      directions will be decoupled from the particle motion */
  int cubic_box;
  /** An index to one of the non_constant dimensions. handy if you just want the variable box_l */
  int non_const_dim;
} nptiso_struct;
extern nptiso_struct nptiso;

/** Allowable values for nptiso.geometry*/
#define NPTGEOM_XDIR 1
#define NPTGEOM_YDIR 2
#define NPTGEOM_ZDIR 4

/* include the potential files */
#include "p3m.h"
#include "lj.h"
#include "ljcos.h"
#include "tab.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "subt_lj_harm.h"
#include "subt_lj_fene.h"
#include "subt_lj.h"
#include "angle.h"
#include "dihedral.h"
#include "debye_hueckel.h"
#include "mmm1d.h"


/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat virials, total_pressure;
///
extern Observable_stat p_tensor;
///
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Callback for setting \ref nptiso_struct::piston */
int piston_callback(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref nptiso_struct::p_ext */
int p_ext_callback(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref nptiso_struct::p_diff */
int p_diff_callback(Tcl_Interp *interp, void *_data);

/** Calculates the pressure in the system from a virial expansion using the terms from \ref calculate_verlet_virials or \ref nsq_calculate_virials dependeing on the used cell system.<BR>
    @param result here all the data is stored
    @param v_comp flag which enables (1) compensation of the velocities required 
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT.
*/
void pressure_calc(double *result, int v_comp);

/** Calculate non bonded energies between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
MDINLINE void add_non_bonded_pair_virials(Particle *p1, Particle *p2, double d[3],
					  double dist, double dist2)
{
  IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
  double F1[3], F2[3];
  int i;
#ifdef ELECTROSTATICS
  double ret; 
#endif
#ifdef NPT
  double p_vir[3];
#endif

  for (i = 0; i < 3; i++) {
    F1[i] = p1->f.f[i];
    p1->f.f[i] = 0;
    F2[i] = p2->f.f[i];
    p2->f.f[i] = 0;
#ifdef NPT
    p_vir[i] = nptiso.p_vir[i];
#endif
  }

  /* lennard jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(p1,p2,ia_params,d,dist);
#endif
  /* lennard jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1,p2,ia_params,d,dist);
#endif
  /* tabulated */
#ifdef TABULATED
  add_tabulated_pair_force(p1,p2,ia_params,d,dist);
#endif
 
#ifdef ROTATION  
  /* Gay-Berne */
  add_gb_pair_force(p1,p2,ia_params,d,dist);
#endif
  *obsstat_nonbonded(&virials, p1->p.type, p2->p.type) += d[0]*p1->f.f[0] + d[1]*p1->f.f[1] + d[2]*p1->f.f[2];

#ifdef ELECTROSTATICS
  /* real space coulomb */
  if (coulomb.method != COULOMB_NONE) {
    for (i = 0; i < 3; i++) {
      p1->f.f[i] = 0;
      p2->f.f[i] = 0;
    }

    switch (coulomb.method) {
    case COULOMB_P3M:
      ret = p3m_coulomb_pair_energy(p1,p2,d,dist2,dist);
      break;
    case COULOMB_DH:
      ret = dh_coulomb_pair_energy(p1,p2,dist);
      break;
    case COULOMB_MMM1D:
      ret = mmm1d_coulomb_pair_energy(p1,p2,d, dist2,dist);
      break;
    default:
      ret = 0;
    }
    virials.coulomb[0] += ret;
  }
#endif
  
  for (i = 0; i < 3; i++) {
    p1->f.f[i] = F1[i];
    p2->f.f[i] = F2[i];
#ifdef NPT
    nptiso.p_vir[i] = p_vir[i];
#endif
  }
}

/** Calculate bonded virials for one particle.
    For performance reasons the force routines add their values directly to the particles.
    So here we do some tricks to get the value out without changing the forces.
    @param p1 particle for which to calculate virials
*/
MDINLINE void add_bonded_virials(Particle *p1)
{
  char *errtxt;

  int i,j, type_num;
  double ret, F1[3], F2[3], F3[3], d[3];
  Particle *p2;
#ifdef NPT
  double p_vir[3];
#endif

  for (j = 0; j < 3; j++) {
    F1[j] = p1->f.f[j];
#ifdef NPT
    p_vir[j] = nptiso.p_vir[j];
#endif
  }
  
  i=0;
  while(i<p1->bl.n) {
    p2 = local_particles[p1->bl.e[i+1]];
    if (!p2) {
      errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
      sprintf(errtxt,"{bond broken between particles %d and %d (particles not stored on the same node)} ",
	      p1->p.identity, p1->bl.e[i+1]);
      return;
    }

    for (j = 0; j < 3; j++) {
      p1->f.f[j] = 0;

      F2[j] = p2->f.f[j];      
      p2->f.f[j] = 0;
    }

    get_mi_vector(d, p1->r.p, p2->r.p);

    type_num = p1->bl.e[i];
    switch(bonded_ia_params[type_num].type) {
    case BONDED_IA_FENE:
      add_fene_pair_force(p1,p2,type_num);
      ret = d[0]*p1->f.f[0] + d[1]*p1->f.f[1] + d[2]*p1->f.f[2];	    
      i+=2; break;
    case BONDED_IA_HARMONIC:
      add_harmonic_pair_force(p1,p2,type_num);
      ret = d[0]*p1->f.f[0] + d[1]*p1->f.f[1] + d[2]*p1->f.f[2];
      i+=2; break;
#ifdef LENNARD_JONES
    case BONDED_IA_SUBT_LJ_HARM:
      add_subt_lj_harm_pair_force(p1,p2,type_num);
      ret = d[0]*p1->f.f[0] + d[1]*p1->f.f[1] + d[2]*p1->f.f[2]; 
      i+=2; break;
    case BONDED_IA_SUBT_LJ_FENE:
      add_subt_lj_fene_pair_force(p1,p2,type_num);
      ret = d[0]*p1->f.f[0] + d[1]*p1->f.f[1] + d[2]*p1->f.f[2];
      i+=2; break;
    case BONDED_IA_SUBT_LJ:
      add_subt_lj_pair_force(p1,p2,type_num);
      ret = d[0]*p1->f.f[0] + d[1]*p1->f.f[1] + d[2]*p1->f.f[2];
      i+=2; break;
#endif
    case BONDED_IA_ANGLE: {
      Particle *p3 = local_particles[p1->bl.e[i+2]];
      if (!p3) {
	errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	sprintf(errtxt,"{bond broken between particles %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i+2]);
	return;
      }
      for (j = 0; j < 3; j++) {
	F3[j] = p3->f.f[j];
	p3->f.f[j] = 0;
      }
      add_angle_force(p1,p2,p3,type_num);
      ret = -d[0]*p2->f.f[0] - d[1]*p2->f.f[1] - d[2]*p2->f.f[2];
      get_mi_vector(d, p1->r.p, p3->r.p);
      ret += -d[0]*p3->f.f[0] - d[1]*p3->f.f[1] - d[2]*p3->f.f[2];

      for (j = 0; j < 3; j++)
	p3->f.f[j] = F3[j];

      i+=3; break;
    }
    case BONDED_IA_DIHEDRAL:
      /* since it is not clear at the moment how to handle a four body interaction here, I skip it */
      ret = 0.0;
      i+=4; break;
   default :
      fprintf(stderr,"add_bonded_virials: WARNING: Bond type %d  of atom %d unknown\n",bonded_ia_params[type_num].type,p1->p.identity);
      ret = 0;
      i = p1->bl.n;
      break;
    }
    *obsstat_bonded(&virials, type_num) += ret;

    for (j = 0; j < 3; j++)
      p2->f.f[j] = F2[j];
  }
  
  for (j = 0; j < 3; j++) {
    p1->f.f[j] = F1[j];
#ifdef NPT
    nptiso.p_vir[j] = p_vir[j];
#endif
  }
}

/** Calculate kinetic pressure (aka energy) for one particle.
    @param p1 particle for which to calculate pressure
    @param v_comp flag which enables (1) compensation of the velocities required 
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT.
*/
MDINLINE void add_kinetic_virials(Particle *p1,int v_comp)
{
  /* kinetic energy */
  if(v_comp)
    virials.data.e[0] += SQR(p1->m.v[0] - p1->f.f[0]) + SQR(p1->m.v[1] - p1->f.f[1]) + SQR(p1->m.v[2] - p1->f.f[2]);
  else
    virials.data.e[0] += SQR(p1->m.v[0]) + SQR(p1->m.v[1]) + SQR(p1->m.v[2]);
#ifdef ROTATION
  virials.data.e[0] += (SQR(p1->m.omega[0]) + SQR(p1->m.omega[1]) + SQR(p1->m.omega[2]))*SQR(time_step);
#endif
}

/** implementation of 'analyze pressure'
    @param v_comp flag which enables (1) compensation of the velocities required 
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT. */
int parse_and_print_pressure(Tcl_Interp *interp, int argc, char **argv, int v_comp);

/** Implementation of 'analyze bins' */
int parse_bins(Tcl_Interp *interp, int argc, char **argv);

/** implementation of 'analyze p_IK1' */
int parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
