// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef FORCES_H
#define FORCES_H
/** \file forces.h Force calculation. 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  \todo Preprocessor switches for all forces (Default: everything is turned on).
 *  \todo Implement more flexible thermostat, e.g. which thermostat to use.
 *
 *  For more information see \ref forces.c "forces.c".
 */
#include <tcl.h>
#include "config.h"
#include "thermostat.h"
#ifdef NPT
#include "pressure.h"
#endif
#include "communication.h"

/* include the force files */
#include "p3m.h"
#include "lj.h"
#include "tab.h"
#include "ljcos.h"
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
#include "mmm2d.h"
#include "constraint.h"
#include "comforce.h"
#include "comfixed.h"

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Calculate forces.
 *
 *  A short list, what the function is doing:
 *  <ol>
 *  <li> Initialize forces with: \ref friction_thermo_langevin (ghost forces with zero).
 *  <li> Calculate \ref tcl_bonded "bonded interaction" forces:<br>
 *       Loop all local particles (not the ghosts). 
 *       <ul>
 *       <li> FENE
 *       <li> ANGLE (cos bend potential)
 *       </ul>
 *  <li> Calculate \ref tcl_non_bonded "non-bonded short range interaction" forces:<br>
 *       Loop all \ref IA_Neighbor::vList "verlet lists" of all \ref #cells.
 *       <ul>
 *       <li> Lennard-Jones.
 *       <li> Real space part: Coulomb.
 *       <li> Ramp.
 *       </ul>
 *  <li> Calculate long range interaction forces:<br>
         Uses \ref P3M_calc_kspace_forces.
 *  </ol>
 */
void force_calc();

/** Calculate non bonded forces between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
MDINLINE void add_non_bonded_pair_force(Particle *p1, Particle *p2, 
					double d[3], double dist, double dist2)
{
  IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node, p1->p.identity, p2->p.identity, dist));

#ifdef DPD
  /* DPD thermostat forces */
  add_dpd_thermo_pair_force(p1,p2,d,dist);
#endif

#ifdef TABULATED
  /* tabulated */
  add_tabulated_pair_force(p1,p2,ia_params,d,dist);
#endif

  /* lennard jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(p1,p2,ia_params,d,dist);
#endif

  /* lennard jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1,p2,ia_params,d,dist);
#endif

#ifdef ROTATION
  /* Gay-Berne */
  add_gb_pair_force(p1,p2,ia_params,d,dist);
#endif

#ifdef ELECTROSTATICS
  /* real space coulomb */
  switch (coulomb.method) {
  case COULOMB_P3M:
    add_p3m_coulomb_pair_force(p1,p2,d,dist2,dist);
    break;
  case COULOMB_DH:
    add_dh_coulomb_pair_force(p1,p2,d,dist);
    break;
  case COULOMB_MMM1D:
    add_mmm1d_coulomb_pair_force(p1,p2,d,dist2,dist);
    break;
  case COULOMB_MMM2D:
    add_mmm2d_coulomb_pair_force(p1,p2,d,dist2,dist);
    break;
  }
#endif
}

/** Calculate bonded forces for one particle.
    @param p1 particle for which to calculate forces
*/
MDINLINE void add_bonded_force(Particle *p1)
{
  char *errtxt;
  Particle *p2;
  int i, type_num;

  i=0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i];
    p2 = local_particles[p1->bl.e[i+1]];
    if (!p2) {
      errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
      sprintf(errtxt,"{bond broken between particles %d and %d (particles not stored on the same node)} ",
	      p1->p.identity, p1->bl.e[i+1]);
      return;
    }

    switch(bonded_ia_params[type_num].type) {
    case BONDED_IA_FENE:
      add_fene_pair_force(p1, p2, type_num);
      i+=2; break;
    case BONDED_IA_HARMONIC:
      add_harmonic_pair_force(p1, p2, type_num);
      i+=2; break;
    case BONDED_IA_SUBT_LJ_HARM:
      add_subt_lj_harm_pair_force(p1, p2, type_num);
      i+=2; break;
    case BONDED_IA_SUBT_LJ_FENE:
      add_subt_lj_fene_pair_force(p1, p2, type_num);
      i+=2; break;
    case BONDED_IA_SUBT_LJ:
      add_subt_lj_pair_force(p1, p2, type_num);
      i+=2; break;
    case BONDED_IA_ANGLE: {
      Particle *p3 = local_particles[p1->bl.e[i+2]];
      if (!p3) {
	errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
	sprintf(errtxt,"{bond broken between particles %d, %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i+1], p1->bl.e[i+2]);
	return;
      }

      add_angle_force(p1, p2, p3, type_num);
      i+=3; break;
    }
    case BONDED_IA_DIHEDRAL: {
      Particle *p3 = local_particles[p1->bl.e[i+2]],
	*p4        = local_particles[p1->bl.e[i+3]];
      if (!p3 || !p4) {
	errtxt = runtime_error(128 + 4*TCL_INTEGER_SPACE);
	sprintf(errtxt,"{bond broken between particles %d, %d, %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i+1], p1->bl.e[i+2], p1->bl.e[i+3]);
	return;
      }
      add_dihedral_force(p2, p1, p3, p4, type_num);
      i+=4; break;
    }
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      switch(bonded_ia_params[type_num].p.tab.type) {
      case TAB_BOND_LENGTH:
	add_tab_bond_force(p1, p2, type_num);
	i+=2; break;
      case TAB_BOND_ANGLE: {
	Particle *p3 = local_particles[p1->bl.e[i+2]];
	if (!p3) {
	  errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
	  sprintf(errtxt,"{bond broken between particles %d, %d and %d (particles not stored on the same node)} ",
		  p1->p.identity, p1->bl.e[i+1], p1->bl.e[i+2]);
	  return;
	}	
	add_tab_angle_force(p1, p2, p3, type_num);
	i+=3; break;
      }
      case TAB_BOND_DIHEDRAL: {
	Particle *p3 = local_particles[p1->bl.e[i+2]],
	  *p4        = local_particles[p1->bl.e[i+3]];
	if (!p3 || !p4) {
	  errtxt = runtime_error(128 + 4*TCL_INTEGER_SPACE);
	  sprintf(errtxt,"{bond broken between particles %d, %d, %d and %d (particles not stored on the same node)} ",
		  p1->p.identity, p1->bl.e[i+1], p1->bl.e[i+2], p1->bl.e[i+3]);
	  return;
	}
	add_tab_dihedral_force(p2, p1, p3, p4, type_num);
	i+=4; break;
      }
      default:
	errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
	sprintf(errtxt,"{add_bonded_force: tabulated bond type of atom %d unknown\n", p1->p.identity);
	return;
      }
      break;
#endif
    default :
      errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
      sprintf(errtxt,"{add_bonded_force: bond type of atom %d unknown\n", p1->p.identity);
      return;
    }
  }
}  

/** add force to another. This is used when collecting ghost forces. */
MDINLINE void add_force(ParticleForce *F_to, ParticleForce *F_add)
{
  int i;
  for (i = 0; i < 3; i++)
    F_to->f[i] += F_add->f[i];
#ifdef ROTATION
  for (i = 0; i < 3; i++)
    F_to->torque[i] += F_add->torque[i];
#endif
}

/*@}*/

#endif
