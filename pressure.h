// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file pressure.h
    Pressure calculation. Really similar to \ref energy.h "energy.h".
*/

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
  /** inverse of piston */
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
      If the value of dimension is less than 3 then box length motion in one or more
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
#include "buckingham.h"
#include "soft_sphere.h"
#include "ljcos.h"
#include "tab.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
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
  double force[3] = {0, 0, 0};
#ifdef ELECTROSTATICS
  double ret; 
#endif
#ifdef ROTATION
  double t1[3], t2[3]; /* dummies */
#endif

  /* lennard jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(p1,p2,ia_params,d,dist, force);
#endif
  /* buckingham*/
#ifdef BUCKINGHAM
  add_buck_pair_force(p1,p2,ia_params,d,dist,force);
#endif
 /*soft-sphere potential*/
#ifdef SOFT_SPHERE
  add_soft_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* lennard jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* tabulated */
#ifdef TABULATED
  add_tabulated_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* Gay-Berne */
#ifdef ROTATION
  add_gb_pair_force(p1,p2,ia_params,d,dist,force,t1,t2);
#endif

  *obsstat_nonbonded(&virials, p1->p.type, p2->p.type) += d[0]*force[0] + d[1]*force[1] + d[2]*force[2];

#ifdef ELECTROSTATICS
  /* real space coulomb */
  if (coulomb.method != COULOMB_NONE) {
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
}

/** Calculate bonded virials for one particle.
    For performance reasons the force routines add their values directly to the particles.
    So here we do some tricks to get the value out without changing the forces.
    @param p1 particle for which to calculate virials
*/
MDINLINE void add_bonded_virials(Particle *p1)
{
  double dx[3], force[3];
  char *errtxt;
  Particle *p2;
  Bonded_ia_parameters *iaparams;

  int i, type_num, type;
  
  i = 0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];
    type = bonded_ia_params[type_num].type;

    /* fetch particle 2 */
    p2 = local_particles[p1->bl.e[i++]];
    if (!p2) {
      errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
      sprintf(errtxt,"{bond broken between particles %d and %d (particles not stored on the same node)} ",
	      p1->p.identity, p1->bl.e[i-1]);
      return;
    }

    get_mi_vector(dx, p1->r.p, p2->r.p);

    switch(type) {
    case BONDED_IA_FENE:
      calc_fene_pair_force(p1,p2,iaparams,dx,force);
      break;
    case BONDED_IA_HARMONIC:
      calc_harmonic_pair_force(p1,p2,iaparams,dx,force);
      break;
#ifdef LENNARD_JONES
    case BONDED_IA_SUBT_LJ:
      calc_subt_lj_pair_force(p1,p2,iaparams,dx,force);
      break;
#endif
      /* since it is not clear at the moment how to handle a many body interaction here, I skip it */
    case BONDED_IA_ANGLE:
      i++; force[0] = force[1] = force[2] = 0; break;
    case BONDED_IA_DIHEDRAL:
      i+=2; force[0] = force[1] = force[2] = 0; break;
#ifdef BOND_CONSTRAINT
    case BONDED_IA_RIGID_BOND:
      i +=2; force[0] = force[1] = force[2] = 0; break;
#endif
    default :
      fprintf(stderr,"add_bonded_virials: WARNING: Bond type %d of atom %d unhandled\n",bonded_ia_params[type_num].type,p1->p.identity);
      force[0] = force[1] = force[2] = 0;
      break;
    }
    *obsstat_bonded(&virials, type_num) += dx[0]*force[0] + dx[1]*force[1] + dx[2]*force[2];
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
    virials.data.e[0] += (SQR(p1->m.v[0] - p1->f.f[0]) + SQR(p1->m.v[1] - p1->f.f[1]) + SQR(p1->m.v[2] - p1->f.f[2]))*PMASS(*p1);
  else
    virials.data.e[0] += (SQR(p1->m.v[0]) + SQR(p1->m.v[1]) + SQR(p1->m.v[2]))*PMASS(*p1);

#ifdef ROTATION
  virials.data.e[0] += (SQR(p1->m.omega[0]) + SQR(p1->m.omega[1]) + SQR(p1->m.omega[2]))*SQR(time_step);
#endif
}

/** implementation of 'analyze pressure'
    @param interp Tcl interpreter
    @param argc   arguments
    @param argv   arguments
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
