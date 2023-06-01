/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef CORE_FORCES_INLINE_HPP
#define CORE_FORCES_INLINE_HPP
/** \file
 *  Force calculation.
 */

#include "config/config.hpp"

#include "forces.hpp"

#include "actor/visitors.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/thermalized_bond_kernel.hpp"
#include "electrostatics/coulomb_inline.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_triel.hpp"
#include "magnetostatics/dipoles_inline.hpp"
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
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "nonbonded_interactions/nonbonded_tab.hpp"
#include "nonbonded_interactions/smooth_step.hpp"
#include "nonbonded_interactions/soft_sphere.hpp"
#include "nonbonded_interactions/thole.hpp"
#include "nonbonded_interactions/wca.hpp"
#include "object-in-fluid/oif_global_forces.hpp"
#include "object-in-fluid/oif_local_forces.hpp"

#ifdef DPD
#include "dpd.hpp"
#endif

#include "Particle.hpp"
#include "bond_error.hpp"
#include "errorhandling.hpp"
#include "exclusions.hpp"
#include "thermostat.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <tuple>

inline ParticleForce calc_non_bonded_pair_force(
    Particle const &p1, Particle const &p2, IA_parameters const &ia_params,
    Utils::Vector3d const &d, double const dist,
    Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel) {

  ParticleForce pf{};
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
  pf.f += thole_pair_force(p1, p2, ia_params, d, dist, coulomb_kernel);
#endif
/* tabulated */
#ifdef TABULATED
  force_factor += tabulated_pair_force_factor(ia_params, dist);
#endif
/* Gay-Berne */
#ifdef GAY_BERNE
  pf += gb_pair_force(p1.quat(), p2.quat(), ia_params, d, dist);
#endif
  pf.f += force_factor * d;
  return pf;
}

inline ParticleForce calc_opposing_force(ParticleForce const &pf,
                                         Utils::Vector3d const &d) {
  ParticleForce out{-pf.f};
#ifdef ROTATION
  // if torque is a null vector, the opposing torque is a null vector too
  // (this check guards from returning a small yet non-null opposing
  // torque due to numerical imprecision)
  if (pf.torque[0] != 0. || pf.torque[1] != 0. || pf.torque[2] != 0.) {
    out.torque = -(pf.torque + vector_product(d, pf.f));
  }
#endif
  return out;
}

/** Calculate non-bonded forces between a pair of particles and update their
 *  forces and torques.
 *  @param[in,out] p1      particle 1.
 *  @param[in,out] p2      particle 2.
 *  @param[in] d           vector between @p p1 and @p p2.
 *  @param[in] dist        distance between @p p1 and @p p2.
 *  @param[in] dist2       distance squared between @p p1 and @p p2.
 *  @param[in] coulomb_kernel  %Coulomb force kernel.
 *  @param[in] dipoles_kernel  Dipolar force kernel.
 *  @param[in] elc_kernel      ELC force correction kernel.
 */
inline void add_non_bonded_pair_force(
    Particle &p1, Particle &p2, Utils::Vector3d const &d, double dist,
    double dist2,
    Coulomb::ShortRangeForceKernel::kernel_type const *coulomb_kernel,
    Dipoles::ShortRangeForceKernel::kernel_type const *dipoles_kernel,
    Coulomb::ShortRangeForceCorrectionsKernel::kernel_type const *elc_kernel) {
  auto const &ia_params = get_ia_param(p1.type(), p2.type());
  ParticleForce pf{};

  /***********************************************/
  /* non-bonded pair potentials                  */
  /***********************************************/

  if (dist < ia_params.max_cut) {
#ifdef EXCLUSIONS
    if (do_nonbonded(p1, p2))
#endif
      pf += calc_non_bonded_pair_force(p1, p2, ia_params, d, dist,
                                       coulomb_kernel);
  }

  /***********************************************/
  /* short-range electrostatics                  */
  /***********************************************/

#ifdef ELECTROSTATICS
  // real-space electrostatic charge-charge interaction
  auto const q1q2 = p1.q() * p2.q();
  if (q1q2 != 0. and coulomb_kernel != nullptr) {
    pf.f += (*coulomb_kernel)(q1q2, d, dist);
#ifdef P3M
    if (elc_kernel)
      (*elc_kernel)(p1, p2, q1q2);
#endif // P3M
  }
#endif // ELECTROSTATICS

  /*********************************************************************/
  /* everything before this contributes to the virial pressure in NpT, */
  /* but nothing afterwards                                            */
  /*********************************************************************/
#ifdef NPT
  npt_add_virial_force_contribution(pf.f, d);
#endif

  /***********************************************/
  /* thermostat                                  */
  /***********************************************/

  /* The inter dpd force should not be part of the virial */
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    auto const force = dpd_pair_force(p1, p2, ia_params, d, dist, dist2);
    p1.force() += force;
    p2.force() -= force;
  }
#endif

  /***********************************************/
  /* short-range magnetostatics                  */
  /***********************************************/

#ifdef DIPOLES
  // real-space magnetic dipole-dipole
  if (dipoles_kernel) {
    pf += (*dipoles_kernel)(p1, p2, d, dist, dist2);
  }
#endif

  /***********************************************/
  /* add total non-bonded forces to particles    */
  /***********************************************/

  p1.force_and_torque() += pf;
  p2.force_and_torque() += calc_opposing_force(pf, d);
}

/** Compute the bonded interaction force between particle pairs.
 *
 *  @param[in] p1          First particle.
 *  @param[in] p2          Second particle.
 *  @param[in] iaparams    Bonded parameters for the interaction.
 *  @param[in] dx          Vector between @p p1 and @p p2.
 *  @param[in] kernel      %Coulomb force kernel.
 */
inline boost::optional<Utils::Vector3d> calc_bond_pair_force(
    Particle const &p1, Particle const &p2,
    Bonded_IA_Parameters const &iaparams, Utils::Vector3d const &dx,
    Coulomb::ShortRangeForceKernel::kernel_type const *kernel) {
  if (auto const *iap = boost::get<FeneBond>(&iaparams)) {
    return iap->force(dx);
  }
  if (auto const *iap = boost::get<HarmonicBond>(&iaparams)) {
    return iap->force(dx);
  }
  if (auto const *iap = boost::get<QuarticBond>(&iaparams)) {
    return iap->force(dx);
  }
#ifdef ELECTROSTATICS
  if (auto const *iap = boost::get<BondedCoulomb>(&iaparams)) {
    return iap->force(p1.q() * p2.q(), dx);
  }
  if (auto const *iap = boost::get<BondedCoulombSR>(&iaparams)) {
    return iap->force(dx, *kernel);
  }
#endif
#ifdef BOND_CONSTRAINT
  if (boost::get<RigidBond>(&iaparams)) {
    return Utils::Vector3d{};
  }
#endif
#ifdef TABULATED
  if (auto const *iap = boost::get<TabulatedDistanceBond>(&iaparams)) {
    return iap->force(dx);
  }
#endif
  if (boost::get<VirtualBond>(&iaparams)) {
    return Utils::Vector3d{};
  }
  throw BondUnknownTypeError();
}

inline bool add_bonded_two_body_force(
    Bonded_IA_Parameters const &iaparams, Particle &p1, Particle &p2,
    Coulomb::ShortRangeForceKernel::kernel_type const *kernel) {
  auto const dx = box_geo.get_mi_vector(p1.pos(), p2.pos());

  if (auto const *iap = boost::get<ThermalizedBond>(&iaparams)) {
    auto result = iap->forces(p1, p2, dx);
    if (result) {
      auto const &forces = result.get();

      p1.force() += std::get<0>(forces);
      p2.force() += std::get<1>(forces);

      return false;
    }
  } else {
    auto result = calc_bond_pair_force(p1, p2, iaparams, dx, kernel);
    if (result) {
      p1.force() += result.get();
      p2.force() -= result.get();

#ifdef NPT
      npt_add_virial_force_contribution(result.get(), dx);
#endif
      return false;
    }
  }
  return true;
}

inline boost::optional<
    std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
calc_bonded_three_body_force(Bonded_IA_Parameters const &iaparams,
                             Particle const &p1, Particle const &p2,
                             Particle const &p3) {
  if (auto const *iap = boost::get<AngleHarmonicBond>(&iaparams)) {
    return iap->forces(p1.pos(), p2.pos(), p3.pos());
  }
  if (auto const *iap = boost::get<AngleCosineBond>(&iaparams)) {
    return iap->forces(p1.pos(), p2.pos(), p3.pos());
  }
  if (auto const *iap = boost::get<AngleCossquareBond>(&iaparams)) {
    return iap->forces(p1.pos(), p2.pos(), p3.pos());
  }
#ifdef TABULATED
  if (auto const *iap = boost::get<TabulatedAngleBond>(&iaparams)) {
    return iap->forces(p1.pos(), p2.pos(), p3.pos());
  }
#endif
  if (auto const *iap = boost::get<IBMTriel>(&iaparams)) {
    return iap->calc_forces(p1, p2, p3);
  }
  throw BondUnknownTypeError();
}

inline bool add_bonded_three_body_force(Bonded_IA_Parameters const &iaparams,
                                        Particle &p1, Particle &p2,
                                        Particle &p3) {
  if (boost::get<OifGlobalForcesBond>(&iaparams)) {
    return false;
  }
  auto const result = calc_bonded_three_body_force(iaparams, p1, p2, p3);
  if (result) {
    auto const &forces = result.get();

    p1.force() += std::get<0>(forces);
    p2.force() += std::get<1>(forces);
    p3.force() += std::get<2>(forces);

    return false;
  }
  return true;
}

inline boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                                  Utils::Vector3d, Utils::Vector3d>>
calc_bonded_four_body_force(Bonded_IA_Parameters const &iaparams,
                            Particle const &p1, Particle const &p2,
                            Particle const &p3, Particle const &p4) {
  if (auto const *iap = boost::get<OifLocalForcesBond>(&iaparams)) {
    return iap->calc_forces(p1, p2, p3, p4);
  }
  if (auto const *iap = boost::get<IBMTribend>(&iaparams)) {
    return iap->calc_forces(p1, p2, p3, p4);
  }
  if (auto const *iap = boost::get<DihedralBond>(&iaparams)) {
    return iap->forces(p2.pos(), p1.pos(), p3.pos(), p4.pos());
  }
#ifdef TABULATED
  if (auto const *iap = boost::get<TabulatedDihedralBond>(&iaparams)) {
    return iap->forces(p2.pos(), p1.pos(), p3.pos(), p4.pos());
  }
#endif
  throw BondUnknownTypeError();
}

inline bool add_bonded_four_body_force(Bonded_IA_Parameters const &iaparams,
                                       Particle &p1, Particle &p2, Particle &p3,
                                       Particle &p4) {
  auto const result = calc_bonded_four_body_force(iaparams, p1, p2, p3, p4);
  if (result) {
    auto const &forces = result.get();

    p1.force() += std::get<0>(forces);
    p2.force() += std::get<1>(forces);
    p3.force() += std::get<2>(forces);
    p4.force() += std::get<3>(forces);

    return false;
  }

  return true;
}

inline bool
add_bonded_force(Particle &p1, int bond_id, Utils::Span<Particle *> partners,
                 Coulomb::ShortRangeForceKernel::kernel_type const *kernel) {

  // Consider for bond breakage
  if (partners.size() == 1) { // pair bonds
    auto d = box_geo.get_mi_vector(p1.pos(), partners[0]->pos()).norm();
    if (BondBreakage::check_and_handle_breakage(
            p1.id(), {{partners[0]->id(), boost::none}}, bond_id, d)) {
      return false;
    }
  }
  if (partners.size() == 2) { // angle bond
    auto d =
        box_geo.get_mi_vector(partners[0]->pos(), partners[1]->pos()).norm();
    if (BondBreakage::check_and_handle_breakage(
            p1.id(), {{partners[0]->id(), partners[1]->id()}}, bond_id, d)) {
      return false;
    }
  }

  auto const &iaparams = *bonded_ia_params.at(bond_id);

  switch (number_of_partners(iaparams)) {
  case 0:
    return false;
  case 1:
    return add_bonded_two_body_force(iaparams, p1, *partners[0], kernel);
  case 2:
    return add_bonded_three_body_force(iaparams, p1, *partners[0],
                                       *partners[1]);
  case 3:
    return add_bonded_four_body_force(iaparams, p1, *partners[0], *partners[1],
                                      *partners[2]);
  default:
    throw BondInvalidSizeError{number_of_partners(iaparams)};
  }
}

#endif // CORE_FORCES_INLINE_HPP
