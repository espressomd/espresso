/*
  Copyright (C) 2010,2012,2013,2014 The ESPResSo project
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
#ifndef _FORCES_INLINE_HPP
#define _FORCES_INLINE_HPP

/*#include "forces.hpp"

#include <vector>
#include "utils.hpp"
#include "thermostat.hpp"
#ifdef MOLFORCES
#include "topology.hpp"
#endif
#include "npt.hpp"
#include "virtual_sites.hpp"

#include "p3m.hpp"
#include "p3m-dipolar.hpp"
#include "lj.hpp"
#include "ljgen.hpp"
#include "steppot.hpp"
#include "hertzian.hpp"
#include "gaussian.hpp"
#include "bmhtf-nacl.hpp"
#include "buckingham.hpp"
#include "soft_sphere.hpp"
#include "hat.hpp"
#include "maggs.hpp"
#include "tab.hpp"
#include "overlap.hpp"
#include "ljcos.hpp"
#include "ljcos2.hpp"
#include "ljangle.hpp"
#include "gb.hpp"
#include "fene.hpp"
#include "object-in-fluid/stretching_force.hpp"
#include "object-in-fluid/stretchlin_force.hpp"
#include "object-in-fluid/area_force_local.hpp"
#include "object-in-fluid/area_force_global.hpp"
#include "object-in-fluid/bending_force.hpp"
#include "object-in-fluid/volume_force.hpp"
#include "harmonic.hpp"
#include "subt_lj.hpp"
#include "angle.hpp"
#include "angle_harmonic.hpp"
#include "angle_cosine.hpp"
#include "angle_cossquare.hpp"
#include "angledist.hpp"
#include "dihedral.hpp"
#include "debye_hueckel.hpp"
#include "endangledist.hpp"
#include "reaction_field.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "comforce.hpp"
#include "comfixed.hpp"
#include "molforces.hpp"
#include "morse.hpp"
#include "elc.hpp"
#include "iccp3m.hpp"
#include "collision.hpp" 
#include "external_potential.hpp"*/

#ifdef MOLFORCES
#include "topology.hpp"
#endif

#include "magnetic_non_p3m_methods.hpp"
#include "mdlc_correction.hpp"
#include "constraint.hpp"
#include "EspressoSystemInterface.hpp"
#include "forces.hpp"

#include "npt.hpp"
#include "p3m-dipolar.hpp"
#include "lj.hpp"
#include "ljgen.hpp"
#include "steppot.hpp"
#include "hertzian.hpp"
#include "gaussian.hpp"
#include "bmhtf-nacl.hpp"
#include "buckingham.hpp"
#include "soft_sphere.hpp"
#include "hat.hpp"
#include "tab.hpp"
#include "overlap.hpp"
#include "ljcos.hpp"
#include "ljcos2.hpp"
#include "ljangle.hpp"
#include "gb.hpp"
#include "fene.hpp"
#include "object-in-fluid/stretching_force.hpp"
#include "object-in-fluid/stretchlin_force.hpp"
#include "object-in-fluid/area_force_local.hpp"
#include "object-in-fluid/area_force_global.hpp"
#include "object-in-fluid/bending_force.hpp"
#include "object-in-fluid/volume_force.hpp"
#include "harmonic.hpp"
#include "subt_lj.hpp"
#include "angle_harmonic.hpp"
#include "angle_cosine.hpp"
#include "angle_cossquare.hpp"
#include "angledist.hpp"
#include "debye_hueckel.hpp"
#include "endangledist.hpp"
#include "reaction_field.hpp"
#include "comforce.hpp"
#include "comfixed.hpp"
#include "molforces.hpp"
#include "morse.hpp"
#include "elc.hpp"
#include "collision.hpp"
#include "metadynamics.hpp"
#include "angle.hpp"

/** initialize the forces for a ghost particle */
inline void init_ghost_force(Particle *part)
{
  part->f.f[0] = 0;
  part->f.f[1] = 0;
  part->f.f[2] = 0;

#ifdef ROTATION
  {
    double scale;
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;

    /* and rescale quaternion, so it is exactly of unit length */
    scale = sqrt( SQR(part->r.quat[0]) + SQR(part->r.quat[1]) +
          SQR(part->r.quat[2]) + SQR(part->r.quat[3]));
    part->r.quat[0]/= scale;
    part->r.quat[1]/= scale;
    part->r.quat[2]/= scale;
    part->r.quat[3]/= scale;
  }
#endif
}

/** initialize the forces for a real particle */
inline void init_local_particle_force(Particle *part) {
  if ( thermo_switch & THERMO_LANGEVIN )
    friction_thermo_langevin(part);
  else {
    part->f.f[0] = 0;
    part->f.f[1] = 0;
    part->f.f[2] = 0;
  }

#ifdef EXTERNAL_FORCES
  if(part->p.ext_flag & PARTICLE_EXT_FORCE) {
    part->f.f[0] += part->p.ext_force[0];
    part->f.f[1] += part->p.ext_force[1];
    part->f.f[2] += part->p.ext_force[2];
  }
#endif

#ifdef ROTATION
  {
    double scale;
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;

    #ifdef EXTERNAL_FORCES
      if(part->p.ext_flag & PARTICLE_EXT_TORQUE) {
        part->f.torque[0] += part->p.ext_torque[0];
        part->f.torque[1] += part->p.ext_torque[1];
        part->f.torque[2] += part->p.ext_torque[2];
      }
    #endif

    /* and rescale quaternion, so it is exactly of unit length */
    scale = sqrt( SQR(part->r.quat[0]) + SQR(part->r.quat[1]) +
          SQR(part->r.quat[2]) + SQR(part->r.quat[3]));
    part->r.quat[0]/= scale;
    part->r.quat[1]/= scale;
    part->r.quat[2]/= scale;
    part->r.quat[3]/= scale;
  }
#endif
}

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
inline void force_calc()
{
  // Communication step: distribute ghost positions
  cells_update_ghosts();

  // VIRTUAL_SITES pos (and vel for DPD) update for security reason !!!
#ifdef VIRTUAL_SITES
  update_mol_vel_pos();
  ghost_communicator(&cell_structure.update_ghost_pos_comm);
#endif

#if defined(VIRTUAL_SITES_RELATIVE) && defined(LB)
  // This is on a workaround stage:
  // When using virtual sites relative and LB at the same time, it is necessary
  // to reassemble the cell lists after all position updates, also of virtual
  // particles.
  if ((lattice_switch & LATTICE_LB) && cell_structure.type == CELL_STRUCTURE_DOMDEC && (!dd.use_vList) )
    cells_update_ghosts();
#endif

#ifdef COLLISION_DETECTION
  prepare_collision_queue();
#endif

  espressoSystemInterface.update();

  // Compute the forces from the force objects
  for (ActorList::iterator actor = forceActors.begin();
          actor != forceActors.end(); ++actor)
      (*actor)->computeForces(espressoSystemInterface);

#ifdef LB_GPU
#ifdef SHANCHEN
  if (lattice_switch & LATTICE_LB_GPU && this_node == 0) lattice_boltzmann_calc_shanchen_gpu();
#endif // SHANCHEN

  // transfer_momentum_gpu check makes sure the LB fluid doesn't get updated on integrate 0
  // this_node==0 makes sure it is the master node where the gpu exists
  if (lattice_switch & LATTICE_LB_GPU && transfer_momentum_gpu && (this_node == 0) ) lb_calc_particle_lattice_ia_gpu();
#endif // LB_GPU

#ifdef ELECTROSTATICS
  if (iccp3m_initialized && iccp3m_cfg.set_flag)
    iccp3m_iteration();
#endif
  init_forces();

  calc_long_range_forces();

  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    layered_calculate_ia();
    break;
  case CELL_STRUCTURE_DOMDEC:
    if(dd.use_vList) {
      if (rebuild_verletlist)
    build_verlet_lists_and_calc_verlet_ia();
      else
    calculate_verlet_ia();
    }
    else
      calc_link_cell();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_calculate_ia();

  }

#ifdef VOLUME_FORCE
    double volume=0.;

    for (int i=0;i< MAX_OBJECTS_IN_FLUID;i++){
        calc_volume(&volume,i);
        if (volume<1e-100) break;
        add_volume_force(volume,i);
    }
#endif

#ifdef AREA_FORCE_GLOBAL
    double area=0.;

    for (int i=0;i< MAX_OBJECTS_IN_FLUID;i++){
        calc_area_global(&area,i);
        if (area<1e-100) break;
        add_area_global_force(area,i);
    }
#endif

#ifdef LB
  if (lattice_switch & LATTICE_LB) calc_particle_lattice_ia() ;
#endif

#ifdef COMFORCE
  calc_comforce();
#endif

#ifdef METADYNAMICS
  /* Metadynamics main function */
  meta_perform();
#endif

#ifdef CUDA
  copy_forces_from_GPU();
#endif

  // VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
  ghost_communicator(&cell_structure.collect_ghost_force_comm);
  init_forces_ghosts();
  distribute_mol_force();
#endif

  // Communication Step: ghost forces
  ghost_communicator(&cell_structure.collect_ghost_force_comm);

  // apply trap forces to trapped molecules
#ifdef MOLFORCES
  calc_and_apply_mol_constraints();
#endif

  // should be pretty late, since it needs to zero out the total force
#ifdef COMFIXED
  calc_comfixed();
#endif

  // mark that forces are now up-to-date
  recalc_forces = 0;

#ifdef COLLISION_DETECTION
  handle_collisions();
#endif
}


inline void 
calc_non_bonded_pair_force_parts(Particle *p1, Particle *p2, IA_parameters *ia_params,
                                 double d[3], double dist, double dist2, 
                                 double force[3], 
                                 double torque1[3] = NULL, double torque2[3] = NULL) {
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
  add_ljangle_pair_force(p1, p2, ia_params, d, dist);
#endif
  /* smooth step */
#ifdef SMOOTH_STEP
  add_SmSt_pair_force(p1, p2, ia_params, d, dist, dist2, force);
#endif
  /* Hertzian force */
#ifdef HERTZIAN
  add_hertzian_pair_force(p1, p2, ia_params, d, dist, dist2, force);
#endif
  /* Gaussian force */
#ifdef GAUSSIAN
  add_gaussian_pair_force(p1, p2, ia_params, d, dist, dist2, force);
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
  add_gb_pair_force(p1,p2,ia_params,d,dist,force,torque1,torque2);
#endif
#ifdef INTER_RF
  add_interrf_pair_force(p1,p2,ia_params,d,dist, force);
#endif
}

inline void
calc_non_bonded_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
                           double d[3], double dist, double dist2,
                           double force[3],
                           double torque1[3] = NULL, double torque2[3] = NULL) {
#ifdef MOL_CUT
   // You may want to put a correction factor and correction term for smoothing function else then theta
   if (checkIfParticlesInteractViaMolCut(p1,p2,ia_params)==1)
#endif
   {
      calc_non_bonded_pair_force_parts(p1, p2, ia_params, d, dist, dist2, 
                                       force, torque1, torque2);
   }
}

inline void 
calc_non_bonded_pair_force(Particle *p1, Particle *p2, 
                           double d[3], double dist, double dist2, 
                           double force[3]){
  IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
  calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2, force);
}


/** Calculate non bonded forces between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
inline void add_non_bonded_pair_force(Particle *p1, Particle *p2, 
					double d[3], double dist, double dist2)
{
  IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
  double force[3] = { 0., 0., 0. };
  double torque1[3] = { 0., 0., 0. };
  double torque2[3] = { 0., 0., 0. };
  int j;
  

#ifdef COLLISION_DETECTION
  if (collision_params.mode > 0)
    detect_collision(p1,p2);
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
  if ( thermo_switch & THERMO_INTER_DPD ) add_inter_dpd_pair_force(p1,p2,ia_params,d,dist,dist2);
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
  case COULOMB_P3M_GPU:
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
  case COULOMB_MMM1D:
	  if (q1q2) add_mmm1d_coulomb_pair_force(q1q2,d,dist2,dist,force);
	  break;
  case COULOMB_MMM2D:
	  if (q1q2) add_mmm2d_coulomb_pair_force(q1q2,d,dist2,dist,force);
	  break;
  default:
	  break;
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
    p1->f.f[j] += force[j];
    p2->f.f[j] -= force[j];
#ifdef ROTATION
    p1->f.torque[j] += torque1[j];
    p2->f.torque[j] += torque2[j];
#endif
  }
}

/** Calculate bonded forces for one particle.
    @param p1 particle for which to calculate forces
*/
inline void add_bonded_force(Particle *p1)
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
  int i, j, type_num, n_partners, bond_broken;
  BondedInteraction type;

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
    case BONDED_IA_STRETCHING_FORCE:
      bond_broken = calc_stretching_force_pair_force(p1, p2, iaparams, dx, force);
      break;
    case BONDED_IA_STRETCHLIN_FORCE:
      bond_broken = calc_stretchlin_force_pair_force(p1, p2, iaparams, dx, force);
      break;
    case BONDED_IA_AREA_FORCE_LOCAL:
      bond_broken = calc_area_force_local(p1, p2, p3, iaparams, force, force2, force3);
      break;
#ifdef AREA_FORCE_GLOBAL
    case BONDED_IA_AREA_FORCE_GLOBAL:
      bond_broken = 0;
      break;
#endif
    case BONDED_IA_BENDING_FORCE:
      bond_broken = calc_bending_force(p1, p2, p3, p4, iaparams, force, force2);
      break;
#ifdef VOLUME_FORCE
    case BONDED_IA_VOLUME_FORCE:
      bond_broken = 0;
      break;
#endif
#ifdef LENNARD_JONES
    case BONDED_IA_SUBT_LJ:
      bond_broken = calc_subt_lj_pair_force(p1, p2, iaparams, dx, force);
      break;
#endif
#ifdef BOND_ANGLE_OLD
	/* the first case is not needed and should not be called */ 
    case BONDED_IA_ANGLE_OLD:
      bond_broken = calc_angle_force(p1, p2, p3, iaparams, force, force2);
      break;
#endif
#ifdef BOND_ANGLE
    case BONDED_IA_ANGLE_HARMONIC:
      bond_broken = calc_angle_harmonic_force(p1, p2, p3, iaparams, force, force2);
      break;
    case BONDED_IA_ANGLE_COSINE:
      bond_broken = calc_angle_cosine_force(p1, p2, p3, iaparams, force, force2);
      break;
    case BONDED_IA_ANGLE_COSSQUARE:
      bond_broken = calc_angle_cossquare_force(p1, p2, p3, iaparams, force, force2);
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
      
      for (j = 0; j < 3; j++) {
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
      
      for (j = 0; j < 3; j++) {
switch (type) {
	case BONDED_IA_AREA_FORCE_LOCAL:
		p1->f.f[j] += force[j];
		p2->f.f[j] += force2[j];
		p3->f.f[j] += force3[j];
		break;
#ifdef AREA_FORCE_GLOBAL
	case BONDED_IA_AREA_FORCE_GLOBAL:
		break;
#endif
#ifdef VOLUME_FORCE
	case BONDED_IA_VOLUME_FORCE:
		break;
#endif
	default:
		p1->f.f[j] += force[j];
		p2->f.f[j] += force2[j];
		p3->f.f[j] -= (force[j] + force2[j]);
	}
      }
      break;
    case 3:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 4*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{085 bond broken between particles %d, %d, %d and %d} ",
		p1->p.identity, p2->p.identity, p3->p.identity, p4->p.identity); 
	continue;
      }
      for (j = 0; j < 3; j++) {
	switch (type) {
	case BONDED_IA_BENDING_FORCE:
		p1->f.f[j] -= (force[j]*0.5+force2[j]*0.5);
		p2->f.f[j] += force[j];
		p3->f.f[j] -= (force[j]*0.5+force2[j]*0.5);
		p4->f.f[j] += force2[j];
		break;
	default:
		p1->f.f[j] += force[j];
		p2->f.f[j] += force2[j];
		p3->f.f[j] += force3[j];
		p4->f.f[j] -= (force[j] + force2[j] + force3[j]);
	}
      }
      break;
    }
  }
}  

/** add force to another. This is used when collecting ghost forces. */
inline void add_force(ParticleForce *F_to, ParticleForce *F_add)
{
  for (int i = 0; i < 3; i++)
    F_to->f[i] += F_add->f[i];
#ifdef ROTATION
  for (int i = 0; i < 3; i++)
    F_to->torque[i] += F_add->torque[i];
#endif
}

inline void check_particle_force(Particle *part) {
  for (int i=0; i< 3; i++) {
    if (isnan(part->f.f[i])) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{999 force on particle %d was NAN.} ",
		    part->p.identity);
    }
  }

#ifdef ROTATION
  for (int i=0; i< 3; i++) {
    if (isnan(part->f.torque[i])) {
      char *errtext = runtime_error(128);
      ERROR_SPRINTF(errtext,"{999 force on particle %d was NAN.} ",
		    part->p.identity);
    }
  }
#endif
}

#endif
