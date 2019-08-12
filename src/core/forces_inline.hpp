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
#include "forces.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_triel.hpp"
#include "metadynamics.hpp"
#include "nonbonded_interactions/bmhtf-nacl.hpp"
#include "nonbonded_interactions/buckingham.hpp"
#include "nonbonded_interactions/gaussian.hpp"
#include "nonbonded_interactions/gay_berne.hpp"
#include "nonbonded_interactions/hat.hpp"
#include "nonbonded_interactions/hertzian.hpp"
#include "nonbonded_interactions/lj.hpp"
#include "nonbonded_interactions/ljcos.hpp"
#include "nonbonded_interactions/ljcos2.hpp"
#include "nonbonded_interactions/ljgen.hpp"
#include "nonbonded_interactions/morse.hpp"
#include "nonbonded_interactions/nonbonded_tab.hpp"
#include "nonbonded_interactions/smooth_step.hpp"
#include "nonbonded_interactions/soft_sphere.hpp"
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
#ifdef DPD
#include "dpd.hpp"
#endif

/** Initialize the forces for a ghost particle */
inline ParticleForce init_ghost_force(const Particle *) { return {}; }

/** Initialize the forces for a real particle */
inline ParticleForce init_local_particle_force(const Particle *part) {
  auto f = (thermo_switch & THERMO_LANGEVIN) ? friction_thermo_langevin(part)
                                             : ParticleForce{};

#ifdef EXTERNAL_FORCES
  // If individual coordinates are fixed, set force to 0.
  for (int j = 0; j < 3; j++)
    if (part->p.ext_flag & COORD_FIXED(j))
      f.f[j] = 0;
  // Add external force
  if (part->p.ext_flag & PARTICLE_EXT_FORCE)
    f.f += part->p.ext_force;
#endif

#ifdef ROTATION
  {

#ifdef EXTERNAL_FORCES
    if (part->p.ext_flag & PARTICLE_EXT_TORQUE) {
      f.torque += part->p.ext_torque;
    }
#endif

#ifdef ENGINE
    // apply a swimming force in the direction of
    // the particle's orientation axis
    if (part->swim.swimming) {
      f.f += part->swim.f_swim * part->r.calc_director();
    }
#endif
  }
#endif

  return f;
}

inline Utils::Vector3d calc_non_bonded_pair_force_parts(
    Particle const *const p1, Particle const *const p2,
    IA_parameters const *const ia_params, Utils::Vector3d const &d,
    double const dist, double const dist2, Utils::Vector3d *torque1 = nullptr,
    Utils::Vector3d *torque2 = nullptr) {
#ifdef NO_INTRA_NB
  if (p1->p.mol_id == p2->p.mol_id)
    return {};
#endif
  Utils::Vector3d force{};
/* Lennard-Jones */
#ifdef LENNARD_JONES
  force += add_lj_pair_force(ia_params, d, dist);
#endif
/* WCA */
#ifdef WCA
  force += add_wca_pair_force(p1, p2, ia_params, d, dist);
#endif
/* Lennard-Jones generic */
#ifdef LENNARD_JONES_GENERIC
  force += add_ljgen_pair_force(p1, p2, ia_params, d, dist);
#endif
/* smooth step */
#ifdef SMOOTH_STEP
  force += add_SmSt_pair_force(p1, p2, ia_params, d, dist, dist2);
#endif
/* Hertzian force */
#ifdef HERTZIAN
  force += add_hertzian_pair_force(p1, p2, ia_params, d, dist, dist2);
#endif
/* Gaussian force */
#ifdef GAUSSIAN
  force += add_gaussian_pair_force(p1, p2, ia_params, d, dist, dist2);
#endif
/* BMHTF NaCl */
#ifdef BMHTF_NACL
  force += add_BMHTF_pair_force(p1, p2, ia_params, d, dist, dist2);
#endif
/* Buckingham*/
#ifdef BUCKINGHAM
  force += add_buck_pair_force(p1, p2, ia_params, d, dist);
#endif
/* Morse*/
#ifdef MORSE
  force += add_morse_pair_force(p1, p2, ia_params, d, dist);
#endif
/*soft-sphere potential*/
#ifdef SOFT_SPHERE
  force += add_soft_pair_force(p1, p2, ia_params, d, dist);
#endif
/*repulsive membrane potential*/
#ifdef MEMBRANE_COLLISION
  force += add_membrane_collision_pair_force(p1, p2, ia_params, d, dist);
#endif
/*hat potential*/
#ifdef HAT
  force += add_hat_pair_force(p1, p2, ia_params, d, dist);
#endif
/* Lennard-Jones cosine */
#ifdef LJCOS
  force += add_ljcos_pair_force(p1, p2, ia_params, d, dist);
#endif
/* Lennard-Jones cosine */
#ifdef LJCOS2
  force += add_ljcos2_pair_force(p1, p2, ia_params, d, dist);
#endif
/* Thole damping */
#ifdef THOLE
  force += add_thole_pair_force(p1, p2, ia_params, d, dist);
#endif
/* tabulated */
#ifdef TABULATED
  force += add_tabulated_pair_force(p1, p2, ia_params, d, dist);
#endif
/* Gay-Berne */
#ifdef GAY_BERNE
  // The gb force function isn't inlined, probably due to its size
  if (dist < ia_params->gay_berne.cut) {
    auto const forces =
        add_gb_pair_force(p1, p2, ia_params, d, dist, torque1, torque2);
    force += std::get<0>(forces);
    if (torque1) {
      *torque1 += std::get<1>(forces);
    }
    if (torque2) {
      *torque2 += std::get<2>(forces);
    }
  }
#endif
  return force;
}

inline Utils::Vector3d
calc_non_bonded_pair_force(Particle const *const p1, Particle const *const p2,
                           IA_parameters const *const ia_params,
                           Utils::Vector3d const &d, double dist, double dist2,
                           Utils::Vector3d *torque1 = nullptr,
                           Utils::Vector3d *torque2 = nullptr) {
  return calc_non_bonded_pair_force_parts(p1, p2, ia_params, d, dist, dist2,
                                          torque1, torque2);
}

inline Utils::Vector3d calc_non_bonded_pair_force(Particle const *const p1,
                                                  Particle const *const p2,
                                                  Utils::Vector3d const &d,
                                                  double dist, double dist2) {
  IA_parameters const *const ia_params = get_ia_param(p1->p.type, p2->p.type);
  return calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2);
}

/** Calculate non-bonded forces between a pair of particles and update their
 *  forces and torques.
 *  @param[in,out] p1   particle 1.
 *  @param[in,out] p2   particle 2.
 *  @param[in] d        vector between @p p1 and @p p2.
 *  @param dist         distance between @p p1 and @p p2.
 *  @param dist2        distance squared between @p p1 and @p p2.
 */
inline void add_non_bonded_pair_force(Particle *const p1, Particle *const p2,
                                      Utils::Vector3d const &d, double dist,
                                      double dist2) {
  IA_parameters const *const ia_params = get_ia_param(p1->p.type, p2->p.type);
  Utils::Vector3d force{};
  Utils::Vector3d *torque1 = nullptr;
  Utils::Vector3d *torque2 = nullptr;
#if defined(ROTATION) || defined(DIPOLES)
  Utils::Vector3d _torque1{};
  Utils::Vector3d _torque2{};
  torque1 = &_torque1;
  torque2 = &_torque2;
#endif

  /***********************************************/
  /* bond creation and breaking                  */
  /***********************************************/

#ifdef AFFINITY
  /* affinity potential */
  // Prevent jump to non-inlined function
  if (dist < ia_params->affinity.cut) {
    force += add_affinity_pair_force(p1, p2, ia_params, d, dist);
  }
#endif

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node,
                      p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* non-bonded pair potentials                  */
  /***********************************************/

  if (dist < ia_params->max_cut) {
#ifdef EXCLUSIONS
    if (do_nonbonded(p1, p2))
#endif
      force += calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2,
                                          torque1, torque2);
  }

  /***********************************************/
  /* short-range electrostatics                  */
  /***********************************************/

#ifdef ELECTROSTATICS
  {
    auto const forces = Coulomb::calc_pair_force(p1, p2, d, dist);
    force += std::get<0>(forces);
#ifdef P3M
    // forces from the virtual charges
    p1->f.f += std::get<1>(forces);
    p2->f.f += std::get<2>(forces);
#endif
  }
#endif

  /*********************************************************************/
  /* everything before this contributes to the virial pressure in NpT, */
  /* but nothing afterwards                                            */
  /*********************************************************************/
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    for (int j = 0; j < 3; j++) {
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
  /* long-range magnetostatics                   */
  /***********************************************/

#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  {
    auto const forces = Dipole::calc_pair_force(p1, p2, d, dist, dist2);
    force += std::get<0>(forces);
    *torque1 += std::get<1>(forces);
    *torque2 += std::get<2>(forces);
  }
#endif /* ifdef DIPOLES */

  /***********************************************/
  /* add total non-bonded forces to particles    */
  /***********************************************/

  p1->f.f += force;
  p2->f.f -= force;
#if defined(ROTATION) || defined(DIPOLES)
  p1->f.torque += *torque1;
  p2->f.torque += *torque2;
#endif
}

/** Compute the bonded interaction force between particle pairs.
 *
 *  @param[in] p1          First particle.
 *  @param[in] p2          Second particle.
 *  @param[in] iaparams    Bonded parameters for the interaction.
 *  @param[in] dx          Vector between @p p1 and @p p2.
 *  @param[out] torque     Torque on @p p1.
 */
inline boost::optional<Utils::Vector3d>
calc_bond_pair_force(Particle const *const p1, Particle const *const p2,
                     Bonded_ia_parameters const *const iaparams,
                     Utils::Vector3d const &dx, Utils::Vector3d &torque) {

  boost::optional<Utils::Vector3d> result;

  switch (iaparams->type) {
  case BONDED_IA_FENE:
    result = calc_fene_pair_force(iaparams, dx);
    break;
#ifdef ROTATION
  case BONDED_IA_HARMONIC_DUMBBELL: {
    auto values = calc_harmonic_dumbbell_pair_force(p1, iaparams, dx);
    if (values) {
      result = boost::optional<Utils::Vector3d>(std::get<0>(values.get()));
      torque = std::get<1>(values.get());
    } else {
      result = boost::optional<Utils::Vector3d>();
    }
    break;
  }
#endif
  case BONDED_IA_HARMONIC:
    result = calc_harmonic_pair_force(iaparams, dx);
    break;
  case BONDED_IA_QUARTIC:
    result = calc_quartic_pair_force(iaparams, dx);
    break;
#ifdef ELECTROSTATICS
  case BONDED_IA_BONDED_COULOMB:
    result = calc_bonded_coulomb_pair_force(p1, p2, iaparams, dx);
    break;
  case BONDED_IA_BONDED_COULOMB_SR:
    result = calc_bonded_coulomb_sr_pair_force(iaparams, dx);
    break;
#endif
#ifdef LENNARD_JONES
  case BONDED_IA_SUBT_LJ:
    result = calc_subt_lj_pair_force(p1, p2, iaparams, dx);
    break;
#endif
  case BONDED_IA_TABULATED_DISTANCE:
    result = calc_tab_bond_force(iaparams, dx);
    break;
#ifdef UMBRELLA
  case BONDED_IA_UMBRELLA:
    result = calc_umbrella_pair_force(p1, p2, iaparams, dx);
    break;
#endif
  default:
    result = boost::optional<Utils::Vector3d>(Utils::Vector3d{});
    break;

  } // switch type

  return result;
}

/** Calculate bonded forces for one particle.
 *  @param p1   particle for which to calculate forces
 */
inline void add_bonded_force(Particle *const p1) {
  int i = 0;
  while (i < p1->bl.n) {
    Utils::Vector3d force1{};
    Utils::Vector3d force2{};
    Utils::Vector3d force3{};
    Utils::Vector3d force4{};
    Utils::Vector3d torque1{};
    Particle *p2 = nullptr;
    Particle *p3 = nullptr;
    Particle *p4 = nullptr;
    int type_num = p1->bl.e[i++];
    Bonded_ia_parameters const *const iaparams = &bonded_ia_params[type_num];
    int type = iaparams->type;
    int n_partners = iaparams->num;
    bool bond_broken = true;

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
      auto const dx = get_mi_vector(p1->r.p, p2->r.p, box_geo);
      auto result = calc_bond_pair_force(p1, p2, iaparams, dx, torque1);
      if (result) {
        force1 = result.get();
        bond_broken = false;
      }

#ifdef NPT
      if (integ_switch == INTEG_METHOD_NPT_ISO) {
        for (int j = 0; j < 3; j++)
          nptiso.p_vir[j] += force1[j] * dx[j];
      }
#endif

      switch (type) {
      case BONDED_IA_THERMALIZED_DIST: {
        auto result = calc_thermalized_bond_forces(p1, p2, iaparams, dx);
        if (result) {
          std::tie(force1, force2) = result.get();
          bond_broken = false;
        }
        break;
      }
      default:
        break;
      }
    } // 1 partner
    else if (n_partners == 2) {
      switch (type) {
      case BONDED_IA_ANGLE_HARMONIC:
        std::tie(force1, force2, force3) =
            calc_angle_harmonic_force(p1, p2, p3, iaparams);
        bond_broken = false;
        break;
      case BONDED_IA_ANGLE_COSINE:
        std::tie(force1, force2, force3) =
            calc_angle_cosine_force(p1, p2, p3, iaparams);
        bond_broken = false;
        break;
      case BONDED_IA_ANGLE_COSSQUARE:
        std::tie(force1, force2, force3) =
            calc_angle_cossquare_force(p1, p2, p3, iaparams);
        bond_broken = false;
        break;
#ifdef OIF_GLOBAL_FORCES
      case BONDED_IA_OIF_GLOBAL_FORCES:
        bond_broken = false;
        break;
#endif
      case BONDED_IA_TABULATED_ANGLE:
        std::tie(force1, force2, force3) =
            calc_tab_angle_force(p1, p2, p3, iaparams);
        bond_broken = false;
        break;
      case BONDED_IA_IBM_TRIEL: {
        auto result = IBM_Triel_CalcForce(p1, p2, p3, iaparams);
        if (result) {
          std::tie(force1, force2, force3) = result.get();
          bond_broken = false;
        }
        break;
      }
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
      case BONDED_IA_OIF_OUT_DIRECTION: {
        p1->p.out_direction = calc_out_direction(p2, p3, p4, iaparams);
        bond_broken = false;
        break;
      }
#endif
#ifdef OIF_LOCAL_FORCES
      case BONDED_IA_OIF_LOCAL_FORCES:
        // in OIF nomenclature, particles p2 and p3 are common to both triangles
        std::tie(force1, force2, force3, force4) =
            calc_oif_local(p1, p2, p3, p4, iaparams);
        bond_broken = false;
        break;
#endif
      // IMMERSED_BOUNDARY
      case BONDED_IA_IBM_TRIBEND:
        std::tie(force1, force2, force3, force4) =
            IBM_Tribend_CalcForce(p1, p2, p3, p4, iaparams);
        bond_broken = false;
        break;
      case BONDED_IA_DIHEDRAL: {
        auto result = calc_dihedral_force(p1, p2, p3, p4, iaparams);
        if (result) {
          std::tie(force1, force2, force3) = result.get();
          force4 = -(force1 + force2 + force3);
          bond_broken = false;
        }
        break;
      }
      case BONDED_IA_TABULATED_DIHEDRAL: {
        auto result = calc_tab_dihedral_force(p1, p2, p3, p4, iaparams);
        if (result) {
          std::tie(force1, force2, force3) = result.get();
          force4 = -(force1 + force2 + force3);
          bond_broken = false;
        }
        break;
      }
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
                          << " and " << p2->p.identity;
        continue;
      }

      switch (type) {
      case BONDED_IA_THERMALIZED_DIST:
        p1->f.f += force1;
        p2->f.f += force2;
        break;
      default:
        p1->f.f += force1;
        p2->f.f -= force1;
      }
#ifdef ROTATION
      if (type == BONDED_IA_HARMONIC_DUMBBELL) {
        p1->f.torque += torque1;
      }
#endif
      break;
    case 2:
      if (bond_broken) {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p2->p.identity << " and "
                          << p3->p.identity;
        continue;
      }

      p1->f.f += force1;
      p2->f.f += force2;
      p3->f.f += force3;
      break;
    case 3:
      if (bond_broken) {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p2->p.identity << ", " << p3->p.identity
                          << " and " << p4->p.identity;
        continue;
      }

      switch (type) {
      case BONDED_IA_TABULATED_DIHEDRAL:
      case BONDED_IA_DIHEDRAL:
      case BONDED_IA_IBM_TRIBEND:
#ifdef OIF_LOCAL_FORCES
      case BONDED_IA_OIF_LOCAL_FORCES:
#endif
        p1->f.f += force1;
        p2->f.f += force2;
        p3->f.f += force3;
        p4->f.f += force4;
        break;
      } // Switch type of 4-particle bond
      break;
    } // switch number of partners (add forces to particles)
  }   // loop over the particle's bond list
}

inline void check_particle_force(Particle const *const part) {
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

inline void add_single_particle_force(Particle *const p) {
  if (p->bl.n) {
    add_bonded_force(p);
  }
}

#endif
