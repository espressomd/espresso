/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#include "bonded_interactions/angle_cosine.hpp"
#include "bonded_interactions/angle_cossquare.hpp"
#include "bonded_interactions/angle_harmonic.hpp"
#include "bonded_interactions/bonded_tab.hpp"
#include "bonded_interactions/dihedral.hpp"
#include "bonded_interactions/fene.hpp"
#include "bonded_interactions/harmonic.hpp"
#include "bonded_interactions/harmonic_dumbbell.hpp"
#include "bonded_interactions/quartic.hpp"
#include "bonded_interactions/subt_lj.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "bonded_interactions/umbrella.hpp"
#include "collision.hpp"
#include "forces.hpp"
#include "metadynamics.hpp"
#include "nonbonded_interactions/bmhtf-nacl.hpp"
#include "nonbonded_interactions/buckingham.hpp"
#include "nonbonded_interactions/gaussian.hpp"
#include "nonbonded_interactions/gb.hpp"
#include "nonbonded_interactions/hat.hpp"
#include "nonbonded_interactions/hertzian.hpp"
#include "nonbonded_interactions/lj.hpp"
#include "nonbonded_interactions/ljcos.hpp"
#include "nonbonded_interactions/ljcos2.hpp"
#include "nonbonded_interactions/ljgen.hpp"
#include "nonbonded_interactions/morse.hpp"
#include "nonbonded_interactions/nonbonded_tab.hpp"
#include "nonbonded_interactions/soft_sphere.hpp"
#include "nonbonded_interactions/steppot.hpp"
#include "nonbonded_interactions/thole.hpp"
#include "nonbonded_interactions/wca.hpp"
#include "npt.hpp"
#include "object-in-fluid/affinity.hpp"
#include "object-in-fluid/membrane_collision.hpp"
#include "object-in-fluid/oif_global_forces.hpp"
#include "object-in-fluid/oif_local_forces.hpp"
#include "object-in-fluid/out_direction.hpp"
#include "thermostat.hpp"

#ifdef DIPOLES
#include "electrostatics_magnetostatics/dipole_inline.hpp"
#endif

#ifdef ELECTROSTATICS
#include "bonded_interactions/bonded_coulomb.hpp"
#include "bonded_interactions/bonded_coulomb_sr.hpp"
#include "electrostatics_magnetostatics/coulomb_inline.hpp"
#endif
#ifdef IMMERSED_BOUNDARY
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_triel.hpp"
#endif
#ifdef DPD
#include "dpd.hpp"
#endif

/** Initialize the forces for a ghost particle */
inline void init_ghost_force(Particle *part) {
  part->f.f = Utils::Vector3d{};

#ifdef ROTATION
  part->f.torque[0] = 0;
  part->f.torque[1] = 0;
  part->f.torque[2] = 0;
#endif
}

/** Initialize the forces for a real particle */
inline void init_local_particle_force(Particle *part) {
  if (thermo_switch & THERMO_LANGEVIN)
    friction_thermo_langevin(part);
  else {
    part->f.f = Utils::Vector3d{};
  }

#ifdef EXTERNAL_FORCES
  // If individual coordinates are fixed, set force to 0.
  for (int j = 0; j < 3; j++)
    if (part->p.ext_flag & COORD_FIXED(j))
      part->f.f[j] = 0;
  // Add external force
  if (part->p.ext_flag & PARTICLE_EXT_FORCE)
    part->f.f += part->p.ext_force;
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
      part->f.f += part->swim.f_swim * part->r.calc_director();
    }
#endif

    /* and rescale quaternion, so it is exactly of unit length */
    scale = sqrt(Utils::sqr(part->r.quat[0]) + Utils::sqr(part->r.quat[1]) +
                 Utils::sqr(part->r.quat[2]) + Utils::sqr(part->r.quat[3]));
    if (scale == 0) {
      part->r.quat[0] = 1;
    } else {
      part->r.quat[0] /= scale;
      part->r.quat[1] /= scale;
      part->r.quat[2] /= scale;
      part->r.quat[3] /= scale;
    }
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
/* Lennard-Jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(ia_params, d, dist, force);
#endif
/* WCA */
#ifdef WCA
  add_wca_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* Lennard-Jones generic */
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
/* Buckingham*/
#ifdef BUCKINGHAM
  add_buck_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* Morse*/
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
/* Lennard-Jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* Lennard-Jones cosine */
#ifdef LJCOS2
  add_ljcos2_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* Thole damping */
#ifdef THOLE
  add_thole_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* tabulated */
#ifdef TABULATED
  add_tabulated_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* Gay-Berne */
#ifdef GAY_BERNE
  // The gb force function isn't inlined, probably due to its size
  if (dist < ia_params->GB_cut) {
    add_gb_pair_force(p1, p2, ia_params, d, dist, force, torque1, torque2);
  }
#endif
}

inline void calc_non_bonded_pair_force(const Particle *p1, const Particle *p2,
                                       IA_parameters *ia_params, double *d,
                                       double dist, double dist2, double *force,
                                       double *torque1 = nullptr,
                                       double *torque2 = nullptr) {
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
 *  @param p1        pointer to particle 1.
 *  @param p2        pointer to particle 2.
 *  @param d         vector between p1 and p2.
 *  @param dist      distance between p1 and p2.
 *  @param dist2     distance squared between p1 and p2.
 */
inline void add_non_bonded_pair_force(Particle *p1, Particle *p2, double d[3],
                                      double dist, double dist2) {

  IA_parameters *ia_params = get_ia_param(p1->p.type, p2->p.type);
  Utils::Vector3d force{};
  double torque1[3] = {0., 0., 0.};
  double torque2[3] = {0., 0., 0.};
  int j;

  // Early exit if there is no interactoin to calculate
  // The exception for MMM2d is there, because the method assumes that
  // pairs within a cell system layer but outside the cutoff are considered

#ifdef ELECTROSTATICS
  if (!(coulomb.method == COULOMB_MMM2D))
#endif
  {
    if (dist > max_cut)
      return;
  }

  /***********************************************/
  /* bond creation and breaking                  */
  /***********************************************/

#ifdef COLLISION_DETECTION
  if (collision_params.mode != COLLISION_MODE_OFF)
    detect_collision(p1, p2, dist);
#endif

/*affinity potential*/
#ifdef AFFINITY
  // Prevent jump to non-inlined function
  if (dist < ia_params->affinity_cut) {
    add_affinity_pair_force(p1, p2, ia_params, d, dist, force.data());
  }
#endif

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node,
                      p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* non bonded pair potentials                  */
  /***********************************************/

  if (dist < ia_params->max_cut) {
#ifdef EXCLUSIONS
    if (do_nonbonded(p1, p2))
#endif
      calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2,
                                 force.data(), torque1, torque2);
  }

  /***********************************************/
  /* short range electrostatics                  */
  /***********************************************/

#ifdef ELECTROSTATICS
  Coulomb::calc_pair_force(p1, p2, d, dist, force);
#endif

/*********************************************************************/
/* everything before this contributes to the virial pressure in NpT, */
/* but nothing afterwards                                            */
/*********************************************************************/
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    for (j = 0; j < 3; j++) {
      nptiso.p_vir[j] += force[j] * d[j];
    }
  }
#endif

/***********************************************/
/* thermostat                                  */
/***********************************************/

/** The inter dpd force should not be part of the virial */
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    force += dpd_pair_force(p1, p2, ia_params, d, dist, dist2);
  }
#endif

  /***********************************************/
  /* long range magnetostatics                   */
  /***********************************************/

#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  Dipole::calc_pair_force(p1, p2, d, dist, dist2, force);
#endif /* ifdef DIPOLES */

  /***********************************************/
  /* add total nonbonded forces to particle      */
  /***********************************************/

  p1->f.f += force;
  p2->f.f -= force;

  for (j = 0; j < 3; j++) {
#ifdef ROTATION
    p1->f.torque[j] += torque1[j];
    p2->f.torque[j] += torque2[j];
#endif
  }
}

inline int calc_bond_pair_force(Particle *p1, Particle *p2,
                                Bonded_ia_parameters *iaparams, double *dx_,
                                double *force) {
  int bond_broken = 0;

  auto const dx = Utils::Vector3d{dx_, dx_ + 3};

  switch (iaparams->type) {
  case BONDED_IA_FENE:
    bond_broken = calc_fene_pair_force(iaparams, dx, force);
    break;
#ifdef ROTATION
  case BONDED_IA_HARMONIC_DUMBBELL:
    bond_broken = calc_harmonic_dumbbell_pair_force(p1, iaparams, dx, force);
    break;
#endif
  case BONDED_IA_HARMONIC:
    bond_broken = calc_harmonic_pair_force(iaparams, dx, force);
    break;
  case BONDED_IA_QUARTIC:
    bond_broken = calc_quartic_pair_force(iaparams, dx, force);
    break;
#ifdef ELECTROSTATICS
  case BONDED_IA_BONDED_COULOMB:
    bond_broken = calc_bonded_coulomb_pair_force(p1, p2, iaparams, dx, force);
    break;
  case BONDED_IA_BONDED_COULOMB_SR:
    bond_broken = calc_bonded_coulomb_sr_pair_force(iaparams, dx, force);
    break;
#endif
#ifdef LENNARD_JONES
  case BONDED_IA_SUBT_LJ:
    bond_broken = calc_subt_lj_pair_force(p1, p2, iaparams, dx, force);
    break;
#endif
#ifdef TABULATED
  case BONDED_IA_TABULATED:
    if (iaparams->num == 1)
      bond_broken = calc_tab_bond_force(iaparams, dx, force);
    break;
#endif
#ifdef UMBRELLA
  case BONDED_IA_UMBRELLA:
    bond_broken = calc_umbrella_pair_force(p1, p2, iaparams, dx, force);
    break;
#endif
  default:
    bond_broken = 0;
    break;

  } // switch type

  return bond_broken;
}

/** Calculate bonded forces for one particle.
 *  @param p1   particle for which to calculate forces
 */
inline void add_bonded_force(Particle *p1) {
  Particle *p3 = nullptr, *p4 = nullptr;
  int bond_broken = 1;

  int i = 0;
  while (i < p1->bl.n) {
    double dx[3] = {0., 0., 0.};
    double force[3] = {0., 0., 0.};
    double force2[3] = {0., 0., 0.};
    double force3[3] = {0., 0., 0.};
#if defined(OIF_LOCAL_FORCES)
    double force4[3] = {0., 0., 0.};
#endif
    int type_num = p1->bl.e[i++];
    auto iaparams = &bonded_ia_params[type_num];
    int type = iaparams->type;
    int n_partners = iaparams->num;

    Particle *p2 = nullptr;
    if (n_partners) {
      p2 = local_particles[p1->bl.e[i++]];
      if (!p2) {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity;
        return;
      }

      /* fetch particle 3 eventually */
      if (n_partners >= 2) {
        p3 = local_particles[p1->bl.e[i++]];
        if (!p3) {
          runtimeErrorMsg()
              << "bond broken between particles " << p1->p.identity << ", "
              << p1->bl.e[i - 2] << " and " << p1->bl.e[i - 1]
              << " (particles are not stored on the same node)";
          return;
        }
      }

      /* fetch particle 4 eventually */
      if (n_partners >= 3) {
        p4 = local_particles[p1->bl.e[i++]];
        if (!p4) {
          runtimeErrorMsg()
              << "bond broken between particles " << p1->p.identity << ", "
              << p1->bl.e[i - 3] << ", " << p1->bl.e[i - 2] << " and "
              << p1->bl.e[i - 1] << " (particles not stored on the same node)";
          return;
        }
      }
    }

    if (n_partners == 1) {
      /* because of the NPT pressure calculation for pair forces, we need the
         1->2 distance vector here. For many body interactions this vector is
         not needed,
         and the pressure calculation not yet clear. */
      get_mi_vector(dx, p1->r.p, p2->r.p);
      bond_broken = calc_bond_pair_force(p1, p2, iaparams, dx, force);

#ifdef NPT
      if (integ_switch == INTEG_METHOD_NPT_ISO)
        for (int j = 0; j < 3; j++)
          nptiso.p_vir[j] += force[j] * dx[j];
#endif

      switch (type) {
      case BONDED_IA_THERMALIZED_DIST:
        bond_broken =
            calc_thermalized_bond_forces(p1, p2, iaparams, dx, force, force2);
        break;

      default:
        break;
      }
    } // 1 partner
    else if (n_partners == 2) {
      switch (type) {
      case BONDED_IA_ANGLE_HARMONIC:
        bond_broken = calc_angle_harmonic_force(p1, p2, p3, iaparams, force,
                                                force2, force3);
        break;
      case BONDED_IA_ANGLE_COSINE:
        bond_broken = calc_angle_cosine_force(p1, p2, p3, iaparams, force,
                                              force2, force3);
        break;
      case BONDED_IA_ANGLE_COSSQUARE:
        bond_broken = calc_angle_cossquare_force(p1, p2, p3, iaparams, force,
                                                 force2, force3);
        break;
#ifdef OIF_GLOBAL_FORCES
      case BONDED_IA_OIF_GLOBAL_FORCES:
        bond_broken = 0;
        break;
#endif
#ifdef TABULATED
      case BONDED_IA_TABULATED:
        if (iaparams->num == 2)
          bond_broken =
              calc_tab_angle_force(p1, p2, p3, iaparams, force, force2, force3);
        break;
#endif
#ifdef IMMERSED_BOUNDARY
      case BONDED_IA_IBM_TRIEL:
        bond_broken = IBM_Triel_CalcForce(p1, p2, p3, iaparams);
        break;
#endif
      default:
        runtimeErrorMsg() << "add_bonded_force: bond type of atom "
                          << p1->p.identity << " unknown " << type << ","
                          << n_partners << "\n";
        return;
      }
    } // 2 partners (angle bonds...)
    else if (n_partners == 3) {
      switch (type) {
#ifdef MEMBRANE_COLLISION
      case BONDED_IA_OIF_OUT_DIRECTION:
        bond_broken = calc_out_direction(p1, p2, p3, p4, iaparams);
        break;
#endif
#ifdef OIF_LOCAL_FORCES
      case BONDED_IA_OIF_LOCAL_FORCES:
        bond_broken = calc_oif_local(p1, p2, p3, p4, iaparams, force, force2,
                                     force3, force4);
        break;
#endif
// IMMERSED_BOUNDARY
#ifdef IMMERSED_BOUNDARY
      case BONDED_IA_IBM_TRIBEND: {
        IBM_Tribend_CalcForce(p1, p2, p3, p4, *iaparams);
        bond_broken = 0;

        break;
      }
#endif
      case BONDED_IA_DIHEDRAL:
        bond_broken = calc_dihedral_force(p1, p2, p3, p4, iaparams, force,
                                          force2, force3);
        break;
#ifdef TABULATED
      case BONDED_IA_TABULATED:
        if (iaparams->num == 3)
          bond_broken = calc_tab_dihedral_force(p1, p2, p3, p4, iaparams, force,
                                                force2, force3);
        break;
#endif
      default:
        runtimeErrorMsg() << "add_bonded_force: bond type of atom "
                          << p1->p.identity << " unknown " << type << ","
                          << n_partners << "\n";
        return;
      }
    } // 3 bond partners

    switch (n_partners) {
    case 1:
      if (bond_broken) {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << " and " << p2->p.identity
                          << ". Distance vector: " << dx[0] << " " << dx[1]
                          << " " << dx[2];
        continue;
      }

      for (int j = 0; j < 3; j++) {
        switch (type) {
        case BONDED_IA_THERMALIZED_DIST:
          p1->f.f[j] += force[j];
          p2->f.f[j] += force2[j];
          break;
        default:
          p1->f.f[j] += force[j];
          p2->f.f[j] -= force[j];
        }
      }
      break;
    case 2:
      if (bond_broken) {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p2->p.identity << " and "
                          << p3->p.identity;
        continue;
      }

      for (int j = 0; j < 3; j++) {
        switch (type) {
        default:
          p1->f.f[j] += force[j];
          p2->f.f[j] += force2[j];
          p3->f.f[j] += force3[j];
        }
      }
      break;
    case 3:
      if (bond_broken) {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p2->p.identity << ", " << p3->p.identity
                          << " and " << p4->p.identity;
        continue;
      }

      switch (type) {
      case BONDED_IA_DIHEDRAL:
        for (int j = 0; j < 3; j++) {
          p1->f.f[j] += force[j];
          p2->f.f[j] += force2[j];
          p3->f.f[j] += force3[j];
          p4->f.f[j] -= force[j] + force2[j] + force3[j];
        }
        break;

#ifdef OIF_LOCAL_FORCES
      case BONDED_IA_OIF_LOCAL_FORCES:
        for (int j = 0; j < 3; j++) {
          p1->f.f[j] += force2[j];
          p2->f.f[j] += force[j];
          p3->f.f[j] += force3[j];
          p4->f.f[j] += force4[j];
        }
        break;
#endif
      } // Switch type of 4-particle bond
      break;
    } // switch number of partners (add forces to particles)
  }   // loop over the particle's bond list
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
  if (p->bl.n) {
    add_bonded_force(p);
  }
}

#endif
