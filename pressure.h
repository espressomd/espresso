// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef PRESSURE_H
#define PRESSURE_H


#include "config.h"
#include "integrate.h"
#include "statistics.h"
#include "thermostat.h"
#include "communication.h"

/** FOR EACH A LINE AND A COMMENT! */
extern double piston, inv_piston, NpT_volume;
extern double p_ext, p_inst, p_diff, p_vir, p_vel;

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
#include "angle.h"
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

/** Callback for setting \ref piston */
int piston_callback(Tcl_Interp *interp, void *_data);
/** Callback for setting \ref p_ext */
int p_ext_callback(Tcl_Interp *interp, void *_data);

/** Calculates the pressure in the system from a virial expansion using the terms from \ref calculate_verlet_virials or \ref nsq_calculate_virials dependeing on the used cell system.<BR>
    @param result here all the data is stored
*/
void pressure_calc(double *result);

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
  double ret, F1[3], F2[3];
  int i;

  for (i = 0; i < 3; i++) {
    F1[i] = p1->f.f[i];
    p1->f.f[i] = 0;
    F2[i] = p2->f.f[i];
    p2->f.f[i] = 0;
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
  if (coulomb.bjerrum != 0.0) {
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
  }
}

/** Calculate bonded virials for one particle.
    For performance reasons the force routines add their values directly to the particles.
    So here we do some tricks to get the value out without changing the forces.
    @param p1 particle for which to calculate virials
*/
MDINLINE void add_bonded_virials(Particle *p1)
{
  int i,j, type_num;
  double ret, F1[3], F2[3], F3[3], d[3];
  Particle *p2, *p3;

  for (i = 0; i < 3; i++)
    F1[i] = p1->f.f[i];
  
  i=0;
  while(i<p1->bl.n) {
    p2 = checked_particle_ptr(p1->bl.e[i+1]);
    
    for (j = 0; j < 3; j++) {
      p1->f.f[i] = 0;

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
    case BONDED_IA_SUBT_LJ_HARM:
      add_subt_lj_harm_pair_force(p1,p2,type_num);
      ret = d[0]*p1->f.f[0] + d[1]*p1->f.f[1] + d[2]*p1->f.f[2]; 
      i+=2; break;
    case BONDED_IA_SUBT_LJ_FENE:
      add_subt_lj_fene_pair_force(p1,p2,type_num);
      ret = d[0]*p1->f.f[0] + d[1]*p1->f.f[1] + d[2]*p1->f.f[2];
      i+=2; break;
   case BONDED_IA_ANGLE:
      p3 = checked_particle_ptr(p1->bl.e[i+2]);
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
  
  for (j = 0; j < 3; j++)
    p1->f.f[j] = F1[j];
}

/** Calculate kinetic pressure (aka energy) for one particle.
    @param p1 particle for which to calculate pressure
*/
MDINLINE void add_kinetic_virials(Particle *p1)
{
  /* kinetic energy */
  virials.data.e[0] += SQR(p1->m.v[0]) + SQR(p1->m.v[1]) + SQR(p1->m.v[2]);
#ifdef ROTATION
  virials.data.e[0] += (SQR(p1->m.omega[0]) + SQR(p1->m.omega[1]) + SQR(p1->m.omega[2]))*SQR(time_step);
#endif
}

/** implementation of 'analyze pressure' */
int parse_and_print_pressure(Tcl_Interp *interp, int argc, char **argv);

/** Implementation of 'analyze bins' */
int parse_bins(Tcl_Interp *interp, int argc, char **argv);

/** implementation of 'analyze p_IK1' */
int parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
