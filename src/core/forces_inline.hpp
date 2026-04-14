/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

/** \file
 *  Force calculation.
 */

#include <config/config.hpp>

#include "BoxGeometry.hpp"
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

#ifdef ESPRESSO_DPD
#include "dpd.hpp"
#endif

#include "Particle.hpp"
#include "bond_error.hpp"
#include "errorhandling.hpp"
#include "exclusions.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <optional>
#include <span>
#include <tuple>
#include <variant>

ESPRESSO_ATTR_ALWAYS_INLINE
inline Utils::Vector3d calc_central_radial_force(IA_parameters const &ia_params,
                                                 Utils::Vector3d const &d,
                                                 double const dist) {

  auto force_factor = 0.;
/* Lennard-Jones */
#ifdef ESPRESSO_LENNARD_JONES
  force_factor += lj_pair_force_factor(ia_params, dist);
#endif
/* WCA */
#ifdef ESPRESSO_WCA
  force_factor += wca_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones generic */
#ifdef ESPRESSO_LENNARD_JONES_GENERIC
  force_factor += ljgen_pair_force_factor(ia_params, dist);
#endif
/* smooth step */
#ifdef ESPRESSO_SMOOTH_STEP
  force_factor += SmSt_pair_force_factor(ia_params, dist);
#endif
/* Hertzian force */
#ifdef ESPRESSO_HERTZIAN
  force_factor += hertzian_pair_force_factor(ia_params, dist);
#endif
/* Gaussian force */
#ifdef ESPRESSO_GAUSSIAN
  force_factor += gaussian_pair_force_factor(ia_params, dist);
#endif
/* BMHTF NaCl */
#ifdef ESPRESSO_BMHTF_NACL
  force_factor += BMHTF_pair_force_factor(ia_params, dist);
#endif
/* Buckingham*/
#ifdef ESPRESSO_BUCKINGHAM
  force_factor += buck_pair_force_factor(ia_params, dist);
#endif
/* Morse*/
#ifdef ESPRESSO_MORSE
  force_factor += morse_pair_force_factor(ia_params, dist);
#endif
/*soft-sphere potential*/
#ifdef ESPRESSO_SOFT_SPHERE
  force_factor += soft_pair_force_factor(ia_params, dist);
#endif
/*hat potential*/
#ifdef ESPRESSO_HAT
  force_factor += hat_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones cosine */
#ifdef ESPRESSO_LJCOS
  force_factor += ljcos_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones cosine */
#ifdef ESPRESSO_LJCOS2
  force_factor += ljcos2_pair_force_factor(ia_params, dist);
#endif
/* tabulated */
#ifdef ESPRESSO_TABULATED
  force_factor += tabulated_pair_force_factor(ia_params, dist);
#endif
  return force_factor * d;
}

inline ParticleForce calc_non_central_force(Particle const &p1,
                                            Particle const &p2,
                                            IA_parameters const &ia_params,
                                            Utils::Vector3d const &d,
                                            double const dist) {

  ParticleForce pf{};
/* Gay-Berne */
#ifdef ESPRESSO_GAY_BERNE
  pf += gb_pair_force(p1.quat(), p2.quat(), ia_params, d, dist);
#endif
  return pf;
}

inline ParticleForce calc_opposing_force(ParticleForce const &pf,
                                         Utils::Vector3d const &d) {
  ParticleForce out{-pf.f};
#ifdef ESPRESSO_ROTATION
  // if torque is a null vector, the opposing torque is a null vector too
  // (this check guards from returning a small yet non-null opposing
  // torque due to numerical imprecision)
  if (pf.torque[0] != 0. || pf.torque[1] != 0. || pf.torque[2] != 0.) {
    out.torque = -(pf.torque + vector_product(d, pf.f));
  }
#endif
  return out;
}

/** Compute the bonded interaction force between particle pairs.
 *
 *  @param[in] iaparams    Bonded parameters for the interaction.
 *  @param[in] q1q2        Product of the particle charges.
 *  @param[in] dx          Vector between @p p1 and @p p2.
 *  @param[in] kernel      Coulomb force kernel.
 */
ESPRESSO_ATTR_ALWAYS_INLINE
inline std::optional<Utils::Vector3d> calc_bond_pair_force(
    Bonded_IA_Parameters const &iaparams, Utils::Vector3d const &dx,
    double const q1q2,
    Coulomb::ShortRangeForceKernel::kernel_type const *kernel) {
  if (auto const *iap = std::get_if<FeneBond>(&iaparams)) {
    return iap->force(dx);
  }
  if (auto const *iap = std::get_if<HarmonicBond>(&iaparams)) {
    return iap->force(dx);
  }
  if (auto const *iap = std::get_if<QuarticBond>(&iaparams)) {
    return iap->force(dx);
  }
#ifdef ESPRESSO_ELECTROSTATICS
  if (auto const *iap = std::get_if<BondedCoulomb>(&iaparams)) {
    return iap->force(q1q2, dx);
  }
  if (auto const *iap = std::get_if<BondedCoulombSR>(&iaparams)) {
    return iap->force(dx, *kernel);
  }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (std::get_if<RigidBond>(&iaparams)) {
    return Utils::Vector3d{};
  }
#endif
#ifdef ESPRESSO_TABULATED
  if (auto const *iap = std::get_if<TabulatedDistanceBond>(&iaparams)) {
    return iap->force(dx);
  }
#endif
  if (std::get_if<VirtualBond>(&iaparams)) {
    return Utils::Vector3d{};
  }
  throw BondUnknownTypeError();
}

ESPRESSO_ATTR_ALWAYS_INLINE
inline std::optional<
    std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
calc_bonded_three_body_force(Bonded_IA_Parameters const &iaparams,
                             Utils::Vector3d const &vec1,
                             Utils::Vector3d const &vec2) {
  if (auto const *iap = std::get_if<AngleHarmonicBond>(&iaparams)) {
    return iap->forces(vec1, vec2);
  }
  if (auto const *iap = std::get_if<AngleCosineBond>(&iaparams)) {
    return iap->forces(vec1, vec2);
  }
  if (auto const *iap = std::get_if<AngleCossquareBond>(&iaparams)) {
    return iap->forces(vec1, vec2);
  }
#ifdef ESPRESSO_TABULATED
  if (auto const *iap = std::get_if<TabulatedAngleBond>(&iaparams)) {
    return iap->forces(vec1, vec2);
  }
#endif
  if (auto const *iap = std::get_if<IBMTriel>(&iaparams)) {
    return iap->calc_forces(vec1, vec2);
  }
  throw BondUnknownTypeError();
}

ESPRESSO_ATTR_ALWAYS_INLINE
inline std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                                Utils::Vector3d, Utils::Vector3d>>
calc_bonded_four_body_force(
    Bonded_IA_Parameters const &iaparams, BoxGeometry const &box_geo,
    Utils::Vector3d const &pos1, Utils::Vector3d const &pos2,
    Utils::Vector3d const &pos3, Utils::Vector3d const &pos4,
    Utils::Vector3d const &vel1, Utils::Vector3d const &vel3,
    Utils::Vector3i const &image1) {
  if (auto const *iap = std::get_if<OifLocalForcesBond>(&iaparams)) {
    // note: particles in a dihedral bond are ordered as p2-p1-p3-p4
    auto const fp2 = box_geo.unfolded_position(pos1, image1);
    auto const fp1 = fp2 + box_geo.get_mi_vector(pos2, fp2);
    auto const fp3 = fp2 + box_geo.get_mi_vector(pos3, fp2);
    auto const fp4 = fp2 + box_geo.get_mi_vector(pos4, fp2);
    return iap->calc_forces(fp2, fp1, fp3, fp4, vel1, vel3);
  }
  if (auto const *iap = std::get_if<IBMTribend>(&iaparams)) {
    return iap->calc_forces(box_geo, pos1, pos2, pos3, pos4);
  }
  // note: particles in a dihedral bond are ordered as p2-p1-p3-p4
  auto const v12 = box_geo.get_mi_vector(pos1, pos2);
  auto const v23 = box_geo.get_mi_vector(pos3, pos1);
  auto const v34 = box_geo.get_mi_vector(pos4, pos3);
  if (auto const *iap = std::get_if<DihedralBond>(&iaparams)) {
    return iap->forces(v12, v23, v34);
  }
#ifdef ESPRESSO_TABULATED
  if (auto const *iap = std::get_if<TabulatedDihedralBond>(&iaparams)) {
    return iap->forces(v12, v23, v34);
  }
#endif
  throw BondUnknownTypeError();
}
