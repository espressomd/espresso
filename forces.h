// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef FORCES_H
#define FORCES_H
/** \file forces.h Force calculation. 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  \todo Preprocessor switches for all forces (Default: everything is turned on).
 *  \todo Implement more flexible thermostat, e. g. which thermostat to use.
 *
 *  For more information see \ref forces.c "forces.c".
 */
#include <tcl.h>
#include "utils.h"
#include "thermostat.h"
#ifdef NPT
#include "pressure.h"
#endif
#include "communication.h"
#ifdef MOLFORCES
#include "topology.h"
#endif
/* include the force files */
#include "p3m.h"
#include "lj.h"
#include "buckingham.h"
#include "soft_sphere.h"
#include "maggs.h"
#include "tab.h"
#include "ljcos.h"
#include "ljcos2.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "subt_lj.h"
#include "angle.h"
#include "dihedral.h"
#include "debye_hueckel.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "constraint.h"
#include "comforce.h"
#include "comfixed.h"
#include "molforces.h"
#include "morse.h"
#include "elc.h"

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
 *       <li> Buckingham.
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
  double force[3] = { 0, 0, 0 };
  int j;

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node, p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* thermostat                                  */
  /***********************************************/

#ifdef DPD
  /* DPD thermostat forces */
  add_to_part_dpd_thermo_pair_force(p1,p2,d,dist);
#endif

  /***********************************************/
  /* non bonded pair potentials                  */
  /***********************************************/

#ifdef TABULATED
  /* tabulated */
  add_tabulated_pair_force(p1,p2,ia_params,d,dist,force);
#endif

  /* lennard jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(p1,p2,ia_params,d,dist,force);
#endif

  /* morse */
#ifdef MORSE
  add_morse_pair_force(p1,p2,ia_params,d,dist,force);
#endif

  /* buckingham */
#ifdef BUCKINGHAM
  add_buck_pair_force(p1,p2,ia_params,d,dist,force);
#endif

  /* soft-sphere */
#ifdef SOFT_SPHERE
  add_soft_pair_force(p1,p2,ia_params,d,dist,force);
#endif

  /* lennard jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1,p2,ia_params,d,dist,force);
#endif

  /* lennard jones cos2 */
#ifdef LJCOS2
  add_ljcos2_pair_force(p1,p2,ia_params,d,dist,force);
#endif

#ifdef ROTATION
  /* Gay-Berne */
  add_gb_pair_force(p1,p2,ia_params,d,dist,force,p1->f.torque,p2->f.torque);
#endif

  /***********************************************/
  /* short range electrostatics                  */
  /***********************************************/

#ifdef ELECTROSTATICS
  if (coulomb.method == COULOMB_DH)
    add_dh_coulomb_pair_force(p1,p2,d,dist,force);
#endif

  /*********************************************************************/
  /* everything before this contributes to the virial pressure in NpT, */
  /* but nothing afterwards                                            */
  /*********************************************************************/
#ifdef NPT
  for (j = 0; j < 3; j++)
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[j] += force[j] * d[j];
#endif

  /***********************************************/
  /* long range electrostatics                   */
  /***********************************************/

#ifdef ELECTROSTATICS

  /* real space coulomb */
  switch (coulomb.method) {
  case COULOMB_ELC_P3M: {
    add_p3m_coulomb_pair_force(p1->p.q*p2->p.q,d,dist2,dist,force); 
    
    // forces from the virtual charges
    // they go directly onto the particles, since they are not pairwise forces
    if (elc_params.dielectric_contrast_on)
      ELC_P3M_dielectric_layers_force_contribution(p1, p2, p1->f.f, p2->f.f);
    break;
  }
  case COULOMB_P3M: {
#ifdef NPT
    double eng = add_p3m_coulomb_pair_force(p1->p.q*p2->p.q,d,dist2,dist,force);
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += eng;
#else
    add_p3m_coulomb_pair_force(p1->p.q*p2->p.q,d,dist2,dist,force); 
#endif
    break;
  }
  case COULOMB_EWALD: {
#ifdef NPT
    double eng = add_ewald_coulomb_pair_force(p1,p2,d,dist2,dist,force);
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += eng;
#else
    add_ewald_coulomb_pair_force(p1,p2,d,dist2,dist,force);
#endif
    break;
  }
  case COULOMB_MMM1D:
    add_mmm1d_coulomb_pair_force(p1,p2,d,dist2,dist,force);
    break;
  case COULOMB_MMM2D:
    add_mmm2d_coulomb_pair_force(p1->p.q*p2->p.q,d,dist2,dist,force);
    break;
  case COULOMB_MAGGS:
    if(maggs.yukawa == 1)
      add_maggs_yukawa_pair_force(p1,p2,d,dist2,dist,force);
    break;
  case COULOMB_NONE:
    break;
  }

#endif

  /***********************************************/
  /* add total nonbonded forces to particle      */
  /***********************************************/

  for (j = 0; j < 3; j++) { 
    p1->f.f[j] += force[j];
    p2->f.f[j] -= force[j];
  }
}

/** Calculate bonded forces for one particle.
    @param p1 particle for which to calculate forces
*/
MDINLINE void add_bonded_force(Particle *p1)
{
  double dx[3], force[3], force2[3], force3[3];
  char *errtxt;
  Particle *p2, *p3 = NULL, *p4 = NULL;
  Bonded_ia_parameters *iaparams;
  int i, j, type_num, type, n_partners, bond_broken;

  i = 0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];
    type = iaparams->type;
    n_partners = iaparams->num;

    /* fetch particle 2, which is always needed */
    p2 = local_particles[p1->bl.e[i++]];
    if (!p2) {
      errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{078 bond broken between particles %d and %d (particles not stored on the same node)} ",
	      p1->p.identity, p1->bl.e[i-1]);
      return;
    }

    /* fetch particle 3 eventually */
    if (n_partners >= 2) {
      p3 = local_particles[p1->bl.e[i++]];
      if (!p3) {
	errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{079 bond broken between particles %d, %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i-2], p1->bl.e[i-1]);
	return;
      }
    }

    /* fetch particle 4 eventually */
    if (n_partners >= 3) {
      p4 = local_particles[p1->bl.e[i++]];
      if (!p4) {
	errtxt = runtime_error(128 + 4*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{080 bond broken between particles %d, %d, %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i-3], p1->bl.e[i-2], p1->bl.e[i-1]);
	return;
      }
    }

    if (n_partners == 1) {
      /* because of the NPT pressure calculation for pair forces, we need the
	 1->2 distance vector here. For many body interactions this vector is not needed,
	 and the pressure calculation not yet clear. */
      get_mi_vector(dx, p1->r.p, p2->r.p);
    }

    switch (type) {
    case BONDED_IA_FENE:
      bond_broken = calc_fene_pair_force(p1, p2, iaparams, dx, force);
      break;
    case BONDED_IA_HARMONIC:
      bond_broken = calc_harmonic_pair_force(p1, p2, iaparams, dx, force);
      break;
#ifdef LENNARD_JONES
    case BONDED_IA_SUBT_LJ:
      bond_broken = calc_subt_lj_pair_force(p1, p2, iaparams, dx, force);
      break;
#endif
    case BONDED_IA_ANGLE:
      bond_broken = calc_angle_force(p1, p2, p3, iaparams, force, force2);
      break;
    case BONDED_IA_DIHEDRAL:
      bond_broken = calc_dihedral_force(p1, p2, p3, p4, iaparams, force, force2, force3);
      break;
#ifdef BOND_CONSTRAINT
    case BONDED_IA_RIGID_BOND:
      //add_rigid_bond_pair_force(p1,p2, iaparams, force, force2);
      bond_broken = 0; 
      break;
#endif
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      switch(iaparams->p.tab.type) {
      case TAB_BOND_LENGTH:
	bond_broken = calc_tab_bond_force(p1, p2, iaparams, dx, force);
	break;
      case TAB_BOND_ANGLE:
	bond_broken = calc_tab_angle_force(p1, p2, p3, iaparams, force, force2);
	break;
      case TAB_BOND_DIHEDRAL:
	bond_broken = calc_tab_dihedral_force(p1, p2, p3, p4, iaparams, force, force2, force3);
	break;
      default:
	errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{081 add_bonded_force: tabulated bond type of atom %d unknown\n", p1->p.identity);
	return;
      }
      break;
#endif
    default :
      errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{082 add_bonded_force: bond type of atom %d unknown\n", p1->p.identity);
      return;
    }

    switch (n_partners) {
    case 1:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{083 bond broken between particles %d and %d} ",
		p1->p.identity, p2->p.identity); 
	continue;
      }
      
      for (j = 0; j < 3; j++) {
	p1->f.f[j] += force[j];
	p2->f.f[j] -= force[j];
#ifdef NPT
	if(integ_switch == INTEG_METHOD_NPT_ISO)
	  nptiso.p_vir[j] += force[j] * dx[j];
#endif
      }
      break;
    case 2:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 3*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{084 bond broken between particles %d, %d and %d} ",
		p1->p.identity, p2->p.identity, p3->p.identity); 
	continue;
      }

      for (j = 0; j < 3; j++) {
	p1->f.f[j] += force[j];
	p2->f.f[j] += force2[j];
	p3->f.f[j] -= (force[j] + force2[j]);
      }
      break;
    case 3:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 4*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{085 bond broken between particles %d, %d, %d and %d} ",
		p1->p.identity, p2->p.identity, p3->p.identity, p4->p.identity); 
	continue;
      }

      for (j = 0; j < 3; j++) {
	p1->f.f[j] += force[j];
	p2->f.f[j] += force2[j];
	p3->f.f[j] += force3[j];
	p4->f.f[j] -= (force[j] + force2[j] + force3[j]);
      }
      break;
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
