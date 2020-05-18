/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
#include "bonded_interactions/thermalized_bond.hpp"
#include "bonded_interactions/umbrella.hpp"
#include "errorhandling.hpp"
#include "exclusions.hpp"
#include "forces.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_triel.hpp"
#include "integrators/langevin_inline.hpp"
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
#include "object-in-fluid/oif_global_forces.hpp"
#include "object-in-fluid/oif_local_forces.hpp"
#include "rotation.hpp"
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
inline ParticleForce init_ghost_force(Particle const &) { return {}; }

/** External particle forces */
inline ParticleForce external_force(Particle const &p) {
  ParticleForce f = {};

#ifdef EXTERNAL_FORCES
  f.f += p.p.ext_force;
#ifdef ROTATION
  f.torque += p.p.ext_torque;
#endif
#endif

#ifdef ENGINE
  // apply a swimming force in the direction of
  // the particle's orientation axis
  if (p.p.swim.swimming) {
    f.f += p.p.swim.f_swim * p.r.calc_director();
  }
#endif

  return f;
}

inline ParticleForce thermostat_force(Particle const &p) {
  extern LangevinThermostat langevin;
  if (!(thermo_switch & THERMO_LANGEVIN)) {
    return {};
  }

#ifdef ROTATION
  return {friction_thermo_langevin(langevin, p),
          p.p.rotation ? convert_vector_body_to_space(
                             p, friction_thermo_langevin_rotation(langevin, p))
                       : Utils::Vector3d{}};
#else
  return friction_thermo_langevin(langevin, p);
#endif
}

/** Initialize the forces for a real particle */
inline ParticleForce init_local_particle_force(Particle const &part) {
  return thermostat_force(part) + external_force(part);
}

inline Utils::Vector3d calc_non_bonded_pair_force_parts(
    Particle const &p1, Particle const &p2, IA_parameters const &ia_params,
    Utils::Vector3d const &d, double const dist,
    Utils::Vector3d *torque1 = nullptr, Utils::Vector3d *torque2 = nullptr) {
#ifdef NO_INTRA_NB
  if (p1.p.mol_id == p2.p.mol_id)
    return {};
#endif
  Utils::Vector3d force{};
  double force_factor = 0;
/* Lennard-Jones */
#ifdef LENNARD_JONES
  force_factor += lj_pair_force_factor(ia_params, dist);
#endif
/* WCA */
#ifdef WCA
  force_factor += wca_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones generic */
#ifdef LENNARD_JONES_GENERIC
  force_factor += ljgen_pair_force_factor(ia_params, dist);
#endif
/* smooth step */
#ifdef SMOOTH_STEP
  force_factor += SmSt_pair_force_factor(ia_params, dist);
#endif
/* Hertzian force */
#ifdef HERTZIAN
  force_factor += hertzian_pair_force_factor(ia_params, dist);
#endif
/* Gaussian force */
#ifdef GAUSSIAN
  force_factor += gaussian_pair_force_factor(ia_params, dist);
#endif
/* BMHTF NaCl */
#ifdef BMHTF_NACL
  force_factor += BMHTF_pair_force_factor(ia_params, dist);
#endif
/* Buckingham*/
#ifdef BUCKINGHAM
  force_factor += buck_pair_force_factor(ia_params, dist);
#endif
/* Morse*/
#ifdef MORSE
  force_factor += morse_pair_force_factor(ia_params, dist);
#endif
/*soft-sphere potential*/
#ifdef SOFT_SPHERE
  force_factor += soft_pair_force_factor(ia_params, dist);
#endif
/*hat potential*/
#ifdef HAT
  force_factor += hat_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones cosine */
#ifdef LJCOS
  force_factor += ljcos_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones cosine */
#ifdef LJCOS2
  force_factor += ljcos2_pair_force_factor(ia_params, dist);
#endif
/* Thole damping */
#ifdef THOLE
  force += thole_pair_force(p1, p2, ia_params, d, dist);
#endif
/* tabulated */
#ifdef TABULATED
  force_factor += tabulated_pair_force_factor(ia_params, dist);
#endif
/* Gay-Berne */
#ifdef GAY_BERNE
  // The gb force function isn't inlined, probably due to its size
  if (dist < ia_params.gay_berne.cut) {
    auto const forces =
        gb_pair_force(p1.r.calc_director(), p2.r.calc_director(), ia_params, d,
                      dist, torque1, torque2);
    force += std::get<0>(forces);
    if (torque1) {
      *torque1 += std::get<1>(forces);
    }
    if (torque2) {
      *torque2 += std::get<2>(forces);
    }
  }
#endif
  force += force_factor * d;
  return force;
}

inline Utils::Vector3d calc_non_bonded_pair_force(
    Particle const &p1, Particle const &p2, IA_parameters const &ia_params,
    Utils::Vector3d const &d, double dist, Utils::Vector3d *torque1 = nullptr,
    Utils::Vector3d *torque2 = nullptr) {
  return calc_non_bonded_pair_force_parts(p1, p2, ia_params, d, dist, torque1,
                                          torque2);
}

inline Utils::Vector3d calc_non_bonded_pair_force(Particle const &p1,
                                                  Particle const &p2,
                                                  Utils::Vector3d const &d,
                                                  double dist) {
  IA_parameters const &ia_params = *get_ia_param(p1.p.type, p2.p.type);
  return calc_non_bonded_pair_force(p1, p2, ia_params, d, dist);
}

/** Calculate non-bonded forces between a pair of particles and update their
 *  forces and torques.
 *  @param[in,out] p1   particle 1.
 *  @param[in,out] p2   particle 2.
 *  @param[in] d        vector between @p p1 and @p p2.
 *  @param dist         distance between @p p1 and @p p2.
 *  @param dist2        distance squared between @p p1 and @p p2.
 */
inline void add_non_bonded_pair_force(Particle &p1, Particle &p2,
                                      Utils::Vector3d const &d, double dist,
                                      double dist2) {
  IA_parameters const &ia_params = *get_ia_param(p1.p.type, p2.p.type);
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
  /* non-bonded pair potentials                  */
  /***********************************************/

  if (dist < ia_params.max_cut) {
#ifdef EXCLUSIONS
    if (do_nonbonded(p1, p2))
#endif
      force += calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, torque1,
                                          torque2);
  }

  /***********************************************/
  /* short-range electrostatics                  */
  /***********************************************/

#ifdef ELECTROSTATICS
  {
    auto const forces = Coulomb::pair_force(p1, p2, d, dist);
    force += std::get<0>(forces);
#ifdef P3M
    // forces from the virtual charges
    p1.f.f += std::get<1>(forces);
    p2.f.f += std::get<2>(forces);
#endif
  }
#endif

  /*********************************************************************/
  /* everything before this contributes to the virial pressure in NpT, */
  /* but nothing afterwards                                            */
  /*********************************************************************/
#ifdef NPT
  npt_add_virial_contribution(force, d);
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
    auto const forces = Dipole::pair_force(p1, p2, d, dist, dist2);
    force += std::get<0>(forces);
    *torque1 += std::get<1>(forces);
    *torque2 += std::get<2>(forces);
  }
#endif /* ifdef DIPOLES */

  /***********************************************/
  /* add total non-bonded forces to particles    */
  /***********************************************/

  p1.f.f += force;
  p2.f.f -= force;
#if defined(ROTATION) || defined(DIPOLES)
  p1.f.torque += *torque1;
  p2.f.torque += *torque2;
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
calc_bond_pair_force(Particle const &p1, Particle const &p2,
                     Bonded_ia_parameters const &iaparams,
                     Utils::Vector3d const &dx, Utils::Vector3d &torque) {
  switch (iaparams.type) {
  case BONDED_IA_FENE:
    return fene_pair_force(iaparams, dx);
#ifdef ROTATION
  case BONDED_IA_HARMONIC_DUMBBELL: {
    auto values =
        harmonic_dumbbell_pair_force(p1.r.calc_director(), iaparams, dx);
    if (values) {
      torque = std::get<1>(values.get());
      return boost::optional<Utils::Vector3d>(std::get<0>(values.get()));
    }

    return {};
  }
#endif
  case BONDED_IA_HARMONIC:
    return harmonic_pair_force(iaparams, dx);
  case BONDED_IA_QUARTIC:
    return quartic_pair_force(iaparams, dx);
#ifdef ELECTROSTATICS
  case BONDED_IA_BONDED_COULOMB:
    return bonded_coulomb_pair_force(p1.p.q * p2.p.q, iaparams, dx);
  case BONDED_IA_BONDED_COULOMB_SR:
    return bonded_coulomb_sr_pair_force(iaparams, dx);
#endif
  case BONDED_IA_TABULATED_DISTANCE:
    return tab_bond_force(iaparams, dx);
#ifdef UMBRELLA
  case BONDED_IA_UMBRELLA:
    return umbrella_pair_force(iaparams, dx);
#endif
  case BONDED_IA_VIRTUAL_BOND:
  case BONDED_IA_RIGID_BOND:
    return Utils::Vector3d{};
  default:
    throw BondUnknownTypeError(iaparams.type);
  } // switch type
}

inline bool add_bonded_two_body_force(Bonded_ia_parameters const &iaparams,
                                      Particle &p1, Particle &p2) {
  auto const dx = get_mi_vector(p1.r.p, p2.r.p, box_geo);

  switch (iaparams.type) {
  case BONDED_IA_THERMALIZED_DIST: {
    auto result = thermalized_bond_forces(p1, p2, iaparams, dx);
    if (result) {
      using std::get;
      p1.f.f += get<0>(result.get());
      p2.f.f += get<1>(result.get());

      return false;
    }
  }
  default: {
    Utils::Vector3d torque1{};
    auto result = calc_bond_pair_force(p1, p2, iaparams, dx, torque1);
    if (result) {
      p1.f.f += result.get();
      p2.f.f -= result.get();

#ifdef ROTATION
      p1.f.torque += torque1;
#endif

#ifdef NPT
      npt_add_virial_contribution(result.get(), dx);
#endif
      return false;
    }
  }
  }

  return true;
}

inline boost::optional<
    std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
calc_bonded_three_body_force(Bonded_ia_parameters const &iaparams,
                             Particle const &p1, Particle const &p2,
                             Particle const &p3) {
  switch (iaparams.type) {
  case BONDED_IA_ANGLE_HARMONIC:
    return angle_harmonic_force(p1.r.p, p2.r.p, p3.r.p, iaparams);
  case BONDED_IA_ANGLE_COSINE:
    return angle_cosine_force(p1.r.p, p2.r.p, p3.r.p, iaparams);
  case BONDED_IA_ANGLE_COSSQUARE:
    return angle_cossquare_force(p1.r.p, p2.r.p, p3.r.p, iaparams);
  case BONDED_IA_TABULATED_ANGLE:
    return tab_angle_force(p1.r.p, p2.r.p, p3.r.p, iaparams);
  case BONDED_IA_IBM_TRIEL: {
    return IBM_Triel_CalcForce(p1, p2, p3, iaparams);
  }
  default:
    throw BondUnknownTypeError(iaparams.type);
  }
}

inline bool add_bonded_three_body_force(Bonded_ia_parameters const &iaparams,
                                        Particle &p1, Particle &p2,
                                        Particle &p3) {
  switch (iaparams.type) {
  case BONDED_IA_OIF_GLOBAL_FORCES:
    return false;
  default: {
    auto const result = calc_bonded_three_body_force(iaparams, p1, p2, p3);
    if (result) {
      using std::get;
      auto const &forces = result.get();

      p1.f.f += get<0>(forces);
      p2.f.f += get<1>(forces);
      p3.f.f += get<2>(forces);

      return false;
    }
    break;
  }
  }

  return true;
}

inline boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                                  Utils::Vector3d, Utils::Vector3d>>
calc_bonded_four_body_force(Bonded_ia_parameters const &iaparams,
                            Particle const &p1, Particle const &p2,
                            Particle const &p3, Particle const &p4) {
  switch (iaparams.type) {
  case BONDED_IA_OIF_LOCAL_FORCES:
    return calc_oif_local(p1, p2, p3, p4, iaparams);
  case BONDED_IA_IBM_TRIBEND:
    return IBM_Tribend_CalcForce(p1, p2, p3, p4, iaparams);
  case BONDED_IA_DIHEDRAL:
    return dihedral_force(p2.r.p, p1.r.p, p3.r.p, p4.r.p, iaparams);
  case BONDED_IA_TABULATED_DIHEDRAL:
    return tab_dihedral_force(p2.r.p, p1.r.p, p3.r.p, p4.r.p, iaparams);
  default:
    throw BondUnknownTypeError(iaparams.type);
  }
}

inline bool add_bonded_four_body_force(Bonded_ia_parameters const &iaparams,
                                       Particle &p1, Particle &p2, Particle &p3,
                                       Particle &p4) {
  auto const result = calc_bonded_four_body_force(iaparams, p1, p2, p3, p4);
  if (result) {
    using std::get;
    auto const &forces = result.get();

    p1.f.f += get<0>(forces);
    p2.f.f += get<1>(forces);
    p3.f.f += get<2>(forces);
    p4.f.f += get<3>(forces);

    return false;
  }

  return true;
}

inline bool add_bonded_force(Particle &p1, int bond_id,
                             Utils::Span<Particle *> partners) {
  auto const &iaparams = bonded_ia_params[bond_id];

  switch (iaparams.num) {
  case 0:
    return false;
  case 1:
    return add_bonded_two_body_force(iaparams, p1, *partners[0]);
  case 2:
    return add_bonded_three_body_force(iaparams, p1, *partners[0],
                                       *partners[1]);
  case 3:
    return add_bonded_four_body_force(iaparams, p1, *partners[0], *partners[1],
                                      *partners[2]);
  default:
    throw BondInvalidSizeError{iaparams.num};
  }
}

/** Calculate bonded forces for one particle.
 *  @param p particle for which to calculate forces
 */
inline void add_single_particle_force(Particle &p) {
  cell_structure.execute_bond_handler(p, add_bonded_force);
}
#endif
