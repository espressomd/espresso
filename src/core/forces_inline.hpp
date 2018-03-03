/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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

#include "config.hpp"

#ifdef MOLFORCES
#include "topology.hpp"
#endif

#include "thermostat.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "p3m.hpp"
#include "external_potential.hpp"
#include "angle_cosine.hpp"
#include "angle_cossquare.hpp"
#include "angle_harmonic.hpp"
#include "angledist.hpp"
#include "bmhtf-nacl.hpp"
#include "buckingham.hpp"
#include "collision.hpp"
#include "constraints.hpp"
#include "dihedral.hpp"
#include "thermalized_bond.hpp"
#include "elc.hpp"
#include "endangledist.hpp"
#include "fene.hpp"
#include "forces.hpp"
#include "gaussian.hpp"
#include "gb.hpp"
#include "harmonic.hpp"
#include "harmonic_dumbbell.hpp"
#include "hat.hpp"
#include "hertzian.hpp"
#include "hydrogen_bond.hpp"
#include "lj.hpp"
#include "ljangle.hpp"
#include "ljcos.hpp"
#include "ljcos2.hpp"
#include "ljgen.hpp"
#include "magnetic_non_p3m_methods.hpp"
#include "mdlc_correction.hpp"
#include "metadynamics.hpp"
#include "molforces.hpp"
#include "morse.hpp"
#include "npt.hpp"
#include "object-in-fluid/affinity.hpp"
#include "object-in-fluid/membrane_collision.hpp"
#include "object-in-fluid/oif_global_forces.hpp"
#include "object-in-fluid/oif_local_forces.hpp"
#include "object-in-fluid/out_direction.hpp"
#include "overlap.hpp"
#include "p3m-dipolar.hpp"
#include "quartic.hpp"
#include "soft_sphere.hpp"
#include "steppot.hpp"
#include "subt_lj.hpp"
#include "tab.hpp"
#include "thole.hpp"
#include "twist_stack.hpp"
#include "umbrella.hpp"
#ifdef ELECTROSTATICS
#include "bonded_coulomb.hpp"
#include "debye_hueckel.hpp"
#include "reaction_field.hpp"
#include "scafacos.hpp"
#endif
#ifdef P3M
#include "bonded_coulomb_p3m_sr.hpp"
#endif
#ifdef IMMERSED_BOUNDARY
#include "immersed_boundary/ibm_main.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_triel.hpp"
#include "immersed_boundary/ibm_volume_conservation.hpp"
#endif
#ifdef DPD
#include "dpd.hpp"
#endif

/** initialize the forces for a ghost particle */
inline void init_ghost_force(Particle *part) {
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
    scale = sqrt(Utils::sqr(part->r.quat[0]) + Utils::sqr(part->r.quat[1]) +
                 Utils::sqr(part->r.quat[2]) + Utils::sqr(part->r.quat[3]));
    part->r.quat[0] /= scale;
    part->r.quat[1] /= scale;
    part->r.quat[2] /= scale;
    part->r.quat[3] /= scale;
  }
#endif
}

/** initialize the forces for a real particle */
inline void init_local_particle_force(Particle *part) {
  if (thermo_switch & THERMO_LANGEVIN)
    friction_thermo_langevin(part);
  else {
    part->f.f[0] = 0;
    part->f.f[1] = 0;
    part->f.f[2] = 0;
  }

#ifdef EXTERNAL_FORCES
  if (part->p.ext_flag & PARTICLE_EXT_FORCE) {
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
    if (part->p.ext_flag & PARTICLE_EXT_TORQUE) {
      part->f.torque[0] += part->p.ext_torque[0];
      part->f.torque[1] += part->p.ext_torque[1];
      part->f.torque[2] += part->p.ext_torque[2];
    }
#endif

#ifdef ENGINE
    // apply a swimming force in the direction of
    // the particle's orientation axis
    if (part->swim.swimming) {
      part->f.f[0] += part->swim.f_swim * part->r.quatu[0];
      part->f.f[1] += part->swim.f_swim * part->r.quatu[1];
      part->f.f[2] += part->swim.f_swim * part->r.quatu[2];
    }
#endif

    /* and rescale quaternion, so it is exactly of unit length */
    scale = sqrt(Utils::sqr(part->r.quat[0]) + Utils::sqr(part->r.quat[1]) +
                 Utils::sqr(part->r.quat[2]) + Utils::sqr(part->r.quat[3]));
    part->r.quat[0] /= scale;
    part->r.quat[1] /= scale;
    part->r.quat[2] /= scale;
    part->r.quat[3] /= scale;
  }
#endif
}

inline void calc_non_bonded_pair_force_parts(
    const Particle *const p1, const Particle *const p2,
    IA_parameters *ia_params, double d[3], double dist, double dist2,
    double force[3], double torque1[3] = nullptr, double torque2[3] = nullptr) {
#ifdef NO_INTRA_NB
  if (p1->p.mol_id == p2->p.mol_id)
    return;
#endif
/* lennard jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* lennard jones generic */
#ifdef LENNARD_JONES_GENERIC
  add_ljgen_pair_force(p1, p2, ia_params, d, dist, force);
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
  add_BMHTF_pair_force(p1, p2, ia_params, d, dist, dist2, force);
#endif
/* buckingham*/
#ifdef BUCKINGHAM
  add_buck_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* morse*/
#ifdef MORSE
  add_morse_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/*soft-sphere potential*/
#ifdef SOFT_SPHERE
  add_soft_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/*repulsive membrane potential*/
#ifdef MEMBRANE_COLLISION
  add_membrane_collision_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/*hat potential*/
#ifdef HAT
  add_hat_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* lennard jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* lennard jones cosine */
#ifdef LJCOS2
  add_ljcos2_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* thole damping */
#ifdef THOLE
  add_thole_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* tabulated */
#ifdef TABULATED
  add_tabulated_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* Gay-Berne */
#ifdef GAY_BERNE
  add_gb_pair_force(p1, p2, ia_params, d, dist, force, torque1, torque2);
#endif
#ifdef INTER_RF
  add_interrf_pair_force(p1, p2, ia_params, d, dist, force);
#endif
}

inline void calc_non_bonded_pair_force(Particle *p1, Particle *p2,
                                       IA_parameters *ia_params, double d[3],
                                       double dist, double dist2,
                                       double force[3],
                                       double torque1[3] = nullptr,
                                       double torque2[3] = nullptr) {
    calc_non_bonded_pair_force_parts(p1, p2, ia_params, d, dist, dist2, force,
                                     torque1, torque2);
}

inline void calc_non_bonded_pair_force(Particle *p1, Particle *p2, double d[3],
                                       double dist, double dist2,
                                       double force[3]) {
  IA_parameters *ia_params = get_ia_param(p1->p.type, p2->p.type);
  calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2, force);
}

/** Calculate non bonded forces between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2.
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
inline void add_non_bonded_pair_force(Particle *p1, Particle *p2, double d[3],
                                      double dist, double dist2) {

  IA_parameters *ia_params = get_ia_param(p1->p.type, p2->p.type);
  double force[3] = {0., 0., 0.};
  double torque1[3] = {0., 0., 0.};
  double torque2[3] = {0., 0., 0.};
  int j;

/***********************************************/
/* bond creation and breaking                  */
/***********************************************/

#ifdef COLLISION_DETECTION
  if (collision_params.mode != COLLISION_MODE_OFF)
    detect_collision(p1, p2,dist);
#endif

/*affinity potential*/
#ifdef AFFINITY
  add_affinity_pair_force(p1, p2, ia_params, d, dist, force);
#endif

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node,
                      p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* non bonded pair potentials                  */
  /***********************************************/

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
    calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2, force,
                               torque1, torque2);

/***********************************************/
/* short range electrostatics                  */
/***********************************************/

#ifdef ELECTROSTATICS
  if (coulomb.method == COULOMB_DH)
    add_dh_coulomb_pair_force(p1, p2, d, dist, force);

  if (coulomb.method == COULOMB_RF)
    add_rf_coulomb_pair_force(p1, p2, d, dist, force);
#endif

/*********************************************************************/
/* everything before this contributes to the virial pressure in NpT, */
/* but nothing afterwards                                            */
/*********************************************************************/
#ifdef NPT
  for (j = 0; j < 3; j++)
    if (integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[j] += force[j] * d[j];
#endif

/***********************************************/
/* thermostat                                  */
/***********************************************/

/** The inter dpd force should not be part of the virial */
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    add_dpd_pair_force(p1, p2, ia_params, d, dist, dist2);
  }
#endif

/***********************************************/
/* semi-bonded multi-body potentials            */
/***********************************************/

/* Directional LJ */
#ifdef LJ_ANGLE
  /* This is a multi-body forces that changes the forces of 6 particles */
  add_ljangle_force(p1, p2, ia_params, d, dist);
#endif

/***********************************************/
/* long range electrostatics                   */
/***********************************************/

#ifdef ELECTROSTATICS

  /* real space coulomb */
  const double q1q2 = p1->p.q * p2->p.q;

  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M: {
    if (q1q2) {
      p3m_add_pair_force(q1q2, d, dist2, dist, force);

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
      double eng = p3m_add_pair_force(q1q2, d, dist2, dist, force);
      if (integ_switch == INTEG_METHOD_NPT_ISO)
        nptiso.p_vir[0] += eng;
    }
#else
    if (q1q2)
      p3m_add_pair_force(q1q2, d, dist2, dist, force);
#endif
    break;
  }
#endif
  case COULOMB_MMM1D:
    if (q1q2)
      add_mmm1d_coulomb_pair_force(q1q2, d, dist2, dist, force);
    break;
  case COULOMB_MMM2D:
    if (q1q2)
      add_mmm2d_coulomb_pair_force(q1q2, d, dist2, dist, force);
    break;
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    if (q1q2) {
      Scafacos::add_pair_force(p1, p2, d, dist, force);
    }
    break;
#endif
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
  case DIPOLAR_MDLC_P3M:
  // fall trough
  case DIPOLAR_P3M: {
#ifdef NPT
    double eng = dp3m_add_pair_force(p1, p2, d, dist2, dist, force);
    if (integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += eng;
#else
    dp3m_add_pair_force(p1, p2, d, dist2, dist, force);
#endif
    break;
  }
#endif /*ifdef DP3M */
  default:
    break;
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


/**New add_bonded_force function for bond classes**/
inline void add_bonded_force(Particle *p1)
{

  int bond_broken = bond_container.force_loop(p1);
  if(bond_broken == 2)
    return;

}

inline void check_particle_force(Particle *part) {
  for (int i = 0; i < 3; i++) {
    if (std::isnan(part->f.f[i])) {
      runtimeErrorMsg() << "force on particle " << part->p.identity
                        << " was NAN.";
    }
  }

#ifdef ROTATION
  for (int i = 0; i < 3; i++) {
    if (std::isnan(part->f.torque[i])) {
      runtimeErrorMsg() << "torque on particle " << part->p.identity
                        << " was NAN.";
    }
  }
#endif
}

inline void add_single_particle_force(Particle *p) {
  add_bonded_force(p);
#ifdef CONSTRAINTS
  add_constraints_forces(p);
#endif
#ifdef EXTERNAL_FORCES
  add_external_potential_forces(p);
#endif
}

#endif
