/*
  Copyright (C) 2010,2012 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
#ifndef _FORCES_H
#define _FORCES_H
/** \file forces.h Force calculation. 
 *
 *  \todo Preprocessor switches for all forces (Default: everything is turned on).
 *  \todo Implement more flexible thermostat, %e.g. which thermostat to use.
 *
 *  For more information see forces.c .
 */
#include "utils.h"
#include "thermostat.h"
#ifdef MOLFORCES
#include "topology.h"
#endif
#include "npt.h"
#include "adresso.h"
#include "virtual_sites.h"
#include "metadynamics.h"

/* include the force files */
#include "p3m.h"
#include "p3m-dipolar.h"
#include "ewald.h"
#include "lj.h"
#include "ljgen.h"
#include "steppot.h"
#include "hertzian.h"
#include "bmhtf-nacl.h"
#include "buckingham.h"
#include "soft_sphere.h"
#include "hat.h"
#include "maggs.h"
#include "tab.h"
#include "overlap.h"
#include "ljcos.h"
#include "ljcos2.h"
#include "ljangle.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "subt_lj.h"
#include "angle.h"
#include "angledist.h"
#include "dihedral.h"
#include "debye_hueckel.h"
#include "endangledist.h"
#include "reaction_field.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "comforce.h"
#include "comfixed.h"
#include "molforces.h"
#include "morse.h"
#include "elc.h"
#include "iccp3m.h"
#include "collision.h" 
/* end of force files */

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Calculate forces.
 *
 *  A short list, what the function is doing:
 *  <ol>
 *  <li> Initialize forces with: \ref friction_thermo_langevin (ghost forces with zero).
 *  <li> Calculate bonded interaction forces:<br>
 *       Loop all local particles (not the ghosts). 
 *       <ul>
 *       <li> FENE
 *       <li> ANGLE (cos bend potential)
 *       </ul>
 *  <li> Calculate non-bonded short range interaction forces:<br>
 *       Loop all \ref IA_Neighbor::vList "verlet lists" of all \ref #cells.
 *       <ul>
 *       <li> Lennard-Jones.
 *       <li> Buckingham.
 *       <li> Real space part: Coulomb.
 *       <li> Ramp.
 *       </ul>
 *  <li> Calculate long range interaction forces:<br>
         Uses <a href=P3M_calc_kspace_forces> P3M_calc_kspace_forces </a>
 *  </ol>
 */
void force_calc();

/** Set forces of all ghosts to zero
*/
void init_forces_ghosts();

/** Check if forces are NAN 
*/
void check_forces();

MDINLINE void calc_non_bonded_pair_force_parts(Particle *p1, Particle *p2, IA_parameters *ia_params,double d[3],
					 double dist, double dist2, double force[3],double torgue1[3],double torgue2[3])
{
#ifdef NO_INTRA_NB
  if (p1->p.mol_id==p2->p.mol_id) return;
#endif
  /* lennard jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(p1,p2,ia_params,d,dist, force);
#endif
  /* lennard jones generic */
#ifdef LENNARD_JONES_GENERIC
  add_ljgen_pair_force(p1,p2,ia_params,d,dist, force);
#endif
  /* Directional LJ */
#ifdef LJ_ANGLE
  /* The forces are propagated within the function */
  add_ljangle_pair_force(p1,p2,ia_params,d,dist);
#endif
  /* smooth step */
#ifdef SMOOTH_STEP
  add_SmSt_pair_force(p1,p2,ia_params,d,dist,dist2, force);
#endif
  /* Hertzian force */
#ifdef HERTZIAN
  add_hertzian_pair_force(p1,p2,ia_params,d,dist,dist2, force);
#endif
  /* BMHTF NaCl */
#ifdef BMHTF_NACL
  add_BMHTF_pair_force(p1,p2,ia_params,d,dist,dist2, force);
#endif
  /* buckingham*/
#ifdef BUCKINGHAM
  add_buck_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* morse*/
#ifdef MORSE
  add_morse_pair_force(p1,p2,ia_params,d,dist,force);
#endif
 /*soft-sphere potential*/
#ifdef SOFT_SPHERE
  add_soft_pair_force(p1,p2,ia_params,d,dist,force);
#endif
 /*hat potential*/
#ifdef HAT
  add_hat_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* lennard jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* lennard jones cosine */
#ifdef LJCOS2
  add_ljcos2_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* tabulated */
#ifdef TABULATED
  add_tabulated_pair_force(p1,p2,ia_params,d,dist,force);
#endif
  /* Gay-Berne */
#ifdef GAY_BERNE
  add_gb_pair_force(p1,p2,ia_params,d,dist,force,torgue1,torgue2);
#endif
#ifdef INTER_RF
  add_interrf_pair_force(p1,p2,ia_params,d,dist, force);
#endif
#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
  add_adress_tab_pair_force(p1,p2,ia_params,d,dist,force);
#endif
#endif
}

MDINLINE void calc_non_bonded_pair_force(Particle *p1,Particle *p2,IA_parameters *ia_params,double d[3],double dist,double dist2,double force[3],double t1[3],double t2[3]){
#ifdef MOL_CUT
   //You may want to put a correction factor and correction term for smoothing function else then theta
   if (checkIfParticlesInteractViaMolCut(p1,p2,ia_params)==1)
#endif
   {
      calc_non_bonded_pair_force_parts(p1, p2, ia_params,d, dist, dist2,force,t1,t2);
   }
}

MDINLINE void calc_non_bonded_pair_force_simple(Particle *p1,Particle *p2,double d[3],double dist,double dist2,double force[3]){
   IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
   double t1[3],t2[3];
#ifdef ADRESS
  int j;
  double force_weight=adress_non_bonded_force_weight(p1,p2);
  if (force_weight<ROUND_ERROR_PREC) return;
#endif
  calc_non_bonded_pair_force(p1,p2,ia_params,d,dist,dist2,force,t1,t2);
#ifdef ADRESS
   for (j=0;j<3;j++){
      force[j]*=force_weight;
   }
#endif
}

MDINLINE void calc_non_bonded_pair_force_from_partcfg(Particle *p1,Particle *p2,IA_parameters *ia_params,double d[3],double dist,double dist2,double force[3],double t1[3],double t2[3]){
#ifdef MOL_CUT
   //You may want to put a correction factor and correction term for smoothing function else then theta
   if (checkIfParticlesInteractViaMolCut_partcfg(p1,p2,ia_params)==1)
#endif
   {
     calc_non_bonded_pair_force_parts(p1, p2, ia_params,d, dist, dist2,force,t1,t2);
   }
}

MDINLINE void calc_non_bonded_pair_force_from_partcfg_simple(Particle *p1,Particle *p2,double d[3],double dist,double dist2,double force[3]){
   IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
   double t1[3],t2[3];
   calc_non_bonded_pair_force_from_partcfg(p1,p2,ia_params,d,dist,dist2,force,t1,t2);
}

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
  double force[3] = { 0., 0., 0. };
  double torque1[3] = { 0., 0., 0. };
  double torque2[3] = { 0., 0., 0. };
  int j;
  

#ifdef COLLISION_DETECTION
  if (collision_detection_mode > 0)
     detect_collision(p1,p2);
#endif 

#ifdef ADRESS
  double tmp,force_weight=adress_non_bonded_force_weight(p1,p2);
  if (force_weight<ROUND_ERROR_PREC) return;
#endif

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node, p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* thermostat                                  */
  /***********************************************/

#ifdef DPD
  /* DPD thermostat forces */
  if ( thermo_switch & THERMO_DPD ) add_dpd_thermo_pair_force(p1,p2,d,dist,dist2);
#endif

#ifdef INTER_DPD
  if ( thermo_switch == THERMO_INTER_DPD ) add_inter_dpd_pair_force(p1,p2,ia_params,d,dist,dist2);
#endif

  /***********************************************/
  /* non bonded pair potentials                  */
  /***********************************************/

   calc_non_bonded_pair_force(p1,p2,ia_params,d,dist,dist2,force,torque1,torque2);

  /***********************************************/
  /* short range electrostatics                  */
  /***********************************************/

#ifdef ELECTROSTATICS
  if (coulomb.method == COULOMB_DH)
    add_dh_coulomb_pair_force(p1,p2,d,dist,force);
  
  if (coulomb.method == COULOMB_RF)
    add_rf_coulomb_pair_force(p1,p2,d,dist,force);
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
  double q1q2 = p1->p.q*p2->p.q;
  if (!(iccp3m_initialized && iccp3m_cfg.set_flag)) {
    switch (coulomb.method) {
  #ifdef P3M
    case COULOMB_ELC_P3M: {
      if (q1q2) {
        p3m_add_pair_force(q1q2,d,dist2,dist,force); 
      
        // forces from the virtual charges
        // they go directly onto the particles, since they are not pairwise forces
        if (elc_params.dielectric_contrast_on)
  	ELC_P3M_dielectric_layers_force_contribution(p1, p2, p1->f.f, p2->f.f);
      }
      break;
    }
    case COULOMB_P3M: {
  #ifdef NPT
      if (q1q2) {
        double eng = p3m_add_pair_force(q1q2,d,dist2,dist,force);
        if(integ_switch == INTEG_METHOD_NPT_ISO)
  	nptiso.p_vir[0] += eng;
      }
  #else
      if (q1q2) p3m_add_pair_force(q1q2,d,dist2,dist,force); 
  #endif
      break;
    }
  #endif
    case COULOMB_EWALD: {
  #ifdef NPT
      if (q1q2) {
        double eng = add_ewald_coulomb_pair_force(p1,p2,d,dist2,dist,force);
        if(integ_switch == INTEG_METHOD_NPT_ISO)
  	nptiso.p_vir[0] += eng;
      }
  #else
      if (q1q2) add_ewald_coulomb_pair_force(p1,p2,d,dist2,dist,force);
  #endif
      break;
    }
    case COULOMB_MMM1D:
      if (q1q2) add_mmm1d_coulomb_pair_force(q1q2,d,dist2,dist,force);
      break;
    case COULOMB_MMM2D:
      if (q1q2) add_mmm2d_coulomb_pair_force(q1q2,d,dist2,dist,force);
      break;
    case COULOMB_NONE:
      break;
    }
  }

#endif /*ifdef ELECTROSTATICS */


  /***********************************************/
  /* long range magnetostatics                   */
  /***********************************************/


#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case  DIPOLAR_MDLC_P3M: 
   //fall trough 
  case DIPOLAR_P3M: {
#ifdef NPT
    double eng = dp3m_add_pair_force(p1,p2,d,dist2,dist,force);
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += eng;
#else
    dp3m_add_pair_force(p1,p2,d,dist2,dist,force);
#endif
    break;
  }
#endif /*ifdef DP3M */
  }  
#endif /* ifdef DIPOLES */

  /***********************************************/
  /* add total nonbonded forces to particle      */
  /***********************************************/

  for (j = 0; j < 3; j++) {
#ifdef ADRESS
    tmp=force_weight*force[j];
    p1->f.f[j] += tmp;
    p2->f.f[j] -= tmp;
#else
    p1->f.f[j] += force[j];
    p2->f.f[j] -= force[j];
#endif
#ifdef ROTATION
    p1->f.torque[j] += torque1[j];
    p2->f.torque[j] += torque2[j];
#endif
  }
}

/** Calculate bonded forces for one particle.
    @param p1 particle for which to calculate forces
*/
MDINLINE void add_bonded_force(Particle *p1)
{
  double dx[3]     = { 0., 0., 0. };
  double force[3]  = { 0., 0., 0. };
  double force2[3] = { 0., 0., 0. };
  double force3[3] = { 0., 0., 0. };
#ifdef ROTATION
  double torque1[3] = { 0., 0., 0. };
  double torque2[3] = { 0., 0., 0. };
#endif
  char *errtxt;
  Particle *p2, *p3 = NULL, *p4 = NULL;
  Bonded_ia_parameters *iaparams;
  int i, j, type_num, type, n_partners, bond_broken;

#ifdef ADRESS
  double tmp, force_weight=1;
  //double tmp,force_weight=adress_bonded_force_weight(p1);
  //if (force_weight<ROUND_ERROR_PREC) return;
#endif

  i = 0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];
    type = iaparams->type;
    n_partners = iaparams->num;

    /* fetch particle 2, which is always needed */
    p2 = local_particles[p1->bl.e[i++]];
    if (!p2) {
      errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{078 bond broken between particles %d and %d (particles not stored on the same node)} ",
	      p1->p.identity, p1->bl.e[i-1]);
      return;
    }

    /* fetch particle 3 eventually */
    if (n_partners >= 2) {
      p3 = local_particles[p1->bl.e[i++]];
      if (!p3) {
	errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{079 bond broken between particles %d, %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i-2], p1->bl.e[i-1]);
	return;
      }
    }

    /* fetch particle 4 eventually */
    if (n_partners >= 3) {
      p4 = local_particles[p1->bl.e[i++]];
      if (!p4) {
	errtxt = runtime_error(128 + 4*ES_INTEGER_SPACE);
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
#ifdef BOND_ANGLE
    case BONDED_IA_ANGLE:
      bond_broken = calc_angle_force(p1, p2, p3, iaparams, force, force2);
      break;
#endif
#ifdef BOND_ANGLEDIST
    case BONDED_IA_ANGLEDIST:
      bond_broken = calc_angledist_force(p1, p2, p3, iaparams, force, force2);
      break;
#endif
#ifdef BOND_ENDANGLEDIST
    case BONDED_IA_ENDANGLEDIST:
      bond_broken = calc_endangledist_pair_force(p1, p2, iaparams, dx, force, force2);
      break;
#endif
    case BONDED_IA_DIHEDRAL:
      bond_broken = calc_dihedral_force(p1, p2, p3, p4, iaparams, force, force2, force3);
      break;
#ifdef BOND_CONSTRAINT
    case BONDED_IA_RIGID_BOND:
      //add_rigid_bond_pair_force(p1,p2, iaparams, force, force2);
      bond_broken = 0; 
      force[0]=force[1]=force[2]=0.0;
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
	errtxt = runtime_error(128 + ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{081 add_bonded_force: tabulated bond type of atom %d unknown\n", p1->p.identity);
	return;
      }
      break;
#endif
#ifdef OVERLAPPED
    case BONDED_IA_OVERLAPPED:
      switch(iaparams->p.overlap.type) {
      case OVERLAP_BOND_LENGTH:
        bond_broken = calc_overlap_bond_force(p1, p2, iaparams, dx, force);
        break;
      case OVERLAP_BOND_ANGLE:
        bond_broken = calc_overlap_angle_force(p1, p2, p3, iaparams, force, force2);
        break;
      case OVERLAP_BOND_DIHEDRAL:
        bond_broken = calc_overlap_dihedral_force(p1, p2, p3, p4, iaparams, force, force2, force3);
        break;
      default:
        errtxt = runtime_error(128 + ES_INTEGER_SPACE);
        ERROR_SPRINTF(errtxt,"{081 add_bonded_force: overlapped bond type of atom %d unknown\n", p1->p.identity);
        return;
      }
      break;
#endif
#ifdef BOND_VIRTUAL
    case BONDED_IA_VIRTUAL_BOND:
      bond_broken = 0;
      force[0]=force[1]=force[2]=0.0;
      break;
#endif
    default :
      errtxt = runtime_error(128 + ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{082 add_bonded_force: bond type of atom %d unknown\n", p1->p.identity);
      return;
    }

    switch (n_partners) {
    case 1:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 2*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{083 bond broken between particles %d and %d} ",
		p1->p.identity, p2->p.identity); 
	continue;
      }
      
#ifdef ADRESS
      force_weight = adress_bonded_force_weight(p1,p2);
#endif

      for (j = 0; j < 3; j++) {
#ifdef ADRESS
        tmp=force_weight*force[j];
	p1->f.f[j] += tmp;
	p2->f.f[j] -= tmp;
#else // ADRESS

	switch (type) {
#ifdef BOND_ENDANGLEDIST
	case BONDED_IA_ENDANGLEDIST:
          p1->f.f[j] += force[j];
          p2->f.f[j] += force2[j];
	  break;
#endif // BOND_ENDANGLEDIST
	default:
	  p1->f.f[j] += force[j];
	  p2->f.f[j] -= force[j];
#ifdef ROTATION
	  p1->f.torque[j] += torque1[j];
	  p2->f.torque[j] += torque2[j];
#endif
	}
#endif // NOT ADRESS

#ifdef NPT
	if(integ_switch == INTEG_METHOD_NPT_ISO)
	  nptiso.p_vir[j] += force[j] * dx[j];
#endif
      }
      break;
    case 2:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 3*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{084 bond broken between particles %d, %d and %d} ",
		p1->p.identity, p2->p.identity, p3->p.identity); 
	continue;
      }
      
#ifdef ADRESS
      force_weight=adress_angle_force_weight(p1,p2,p3);
#endif
      for (j = 0; j < 3; j++) {
#ifdef ADRESS
	p1->f.f[j] += force_weight*force[j];
	p2->f.f[j] += force_weight*force2[j];
	p3->f.f[j] -= force_weight*(force[j] + force2[j]);
#else
	p1->f.f[j] += force[j];
	p2->f.f[j] += force2[j];
	p3->f.f[j] -= (force[j] + force2[j]);
#endif
      }
      break;
    case 3:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 4*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{085 bond broken between particles %d, %d, %d and %d} ",
		p1->p.identity, p2->p.identity, p3->p.identity, p4->p.identity); 
	continue;
      }
#ifdef ADRESS
      force_weight=adress_dihedral_force_weight(p1,p2,p3,p4);
#endif 
      for (j = 0; j < 3; j++) {
#ifdef ADRESS
	p1->f.f[j] += force_weight*force[j];
	p2->f.f[j] += force_weight*force2[j];
	p3->f.f[j] += force_weight*force3[j];
	p4->f.f[j] -= force_weight*(force[j] + force2[j] + force3[j]);
#else
	p1->f.f[j] += force[j];
	p2->f.f[j] += force2[j];
	p3->f.f[j] += force3[j];
	p4->f.f[j] -= (force[j] + force2[j] + force3[j]);
#endif
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

MDINLINE void check_particle_force(Particle *part)
{
  
  int i;
  for (i=0; i< 3; i++) {
    if isnan(part->f.f[i]) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{999 force on particle was NAN.} ");
    }
  }

#ifdef ROTATION
  for (i=0; i< 3; i++) {
    if isnan(part->f.torque[i]) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{999 force on particle was NAN.} ");
    }
  }
#endif
}

/*@}*/

#endif
