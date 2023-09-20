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
#ifndef CORE_ENERGY_INLINE_HPP
#define CORE_ENERGY_INLINE_HPP
/** \file
 *  Energy calculation.
 */

#include "config/config.hpp"

#include "energy.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "electrostatics/coulomb_inline.hpp"
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

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "bond_error.hpp"
#include "errorhandling.hpp"
#include "exclusions.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/variant.hpp>

/** Calculate non-bonded energies between a pair of particles.
 *  @param p1         particle 1.
 *  @param p2         particle 2.
 *  @param ia_params  the interaction parameters between the two particles
 *  @param d          vector between p1 and p2.
 *  @param dist       distance between p1 and p2.
 *  @param coulomb_kernel   %Coulomb energy kernel.
 *  @return the short-range interaction energy between the two particles
 */
inline double calc_non_bonded_pair_energy(
    Particle const &p1, Particle const &p2, IA_parameters const &ia_params,
    Utils::Vector3d const &d, double const dist,
    Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_kernel) {

  double ret = 0;

#ifdef LENNARD_JONES
  /* Lennard-Jones */
  ret += lj_pair_energy(ia_params, dist);
#endif
#ifdef WCA
  /* WCA */
  ret += wca_pair_energy(ia_params, dist);
#endif

#ifdef LENNARD_JONES_GENERIC
  /* Generic Lennard-Jones */
  ret += ljgen_pair_energy(ia_params, dist);
#endif

#ifdef SMOOTH_STEP
  /* smooth step */
  ret += SmSt_pair_energy(ia_params, dist);
#endif

#ifdef HERTZIAN
  /* Hertzian potential */
  ret += hertzian_pair_energy(ia_params, dist);
#endif

#ifdef GAUSSIAN
  /* Gaussian potential */
  ret += gaussian_pair_energy(ia_params, dist);
#endif

#ifdef BMHTF_NACL
  /* BMHTF NaCl */
  ret += BMHTF_pair_energy(ia_params, dist);
#endif

#ifdef MORSE
  /* Morse */
  ret += morse_pair_energy(ia_params, dist);
#endif

#ifdef BUCKINGHAM
  /* Buckingham */
  ret += buck_pair_energy(ia_params, dist);
#endif

#ifdef SOFT_SPHERE
  /* soft-sphere */
  ret += soft_pair_energy(ia_params, dist);
#endif

#ifdef HAT
  /* hat */
  ret += hat_pair_energy(ia_params, dist);
#endif

#ifdef LJCOS2
  /* Lennard-Jones */
  ret += ljcos2_pair_energy(ia_params, dist);
#endif

#ifdef THOLE
  /* Thole damping */
  ret += thole_pair_energy(p1, p2, ia_params, d, dist, coulomb_kernel);
#endif

#ifdef TABULATED
  /* tabulated */
  ret += tabulated_pair_energy(ia_params, dist);
#endif

#ifdef LJCOS
  /* Lennard-Jones cosine */
  ret += ljcos_pair_energy(ia_params, dist);
#endif

#ifdef GAY_BERNE
  /* Gay-Berne */
  ret += gb_pair_energy(p1.quat(), p2.quat(), ia_params, d, dist);
#endif

  return ret;
}

/** Add non-bonded and short-range Coulomb energies between a pair of particles
 *  to the energy observable.
 *  @param p1        particle 1.
 *  @param p2        particle 2.
 *  @param d         vector between p1 and p2.
 *  @param dist      distance between p1 and p2.
 *  @param dist2     distance squared between p1 and p2.
 *  @param[in] coulomb_kernel   %Coulomb energy kernel.
 *  @param[in] dipoles_kernel   Dipolar energy kernel.
 *  @param[in,out] obs_energy   energy observable.
 */
inline void add_non_bonded_pair_energy(
    Particle const &p1, Particle const &p2, Utils::Vector3d const &d,
    double const dist, double const dist2,
    Coulomb::ShortRangeEnergyKernel::kernel_type const *coulomb_kernel,
    Dipoles::ShortRangeEnergyKernel::kernel_type const *dipoles_kernel,
    Observable_stat &obs_energy) {
  auto const &ia_params = get_ia_param(p1.type(), p2.type());

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
    obs_energy.add_non_bonded_contribution(
        p1.type(), p2.type(),
        calc_non_bonded_pair_energy(p1, p2, ia_params, d, dist,
                                    coulomb_kernel));

#ifdef ELECTROSTATICS
  if (!obs_energy.coulomb.empty() and coulomb_kernel != nullptr) {
    auto const q1q2 = p1.q() * p2.q();
    obs_energy.coulomb[0] += (*coulomb_kernel)(p1, p2, q1q2, d, dist);
  }
#endif

#ifdef DIPOLES
  if (!obs_energy.dipolar.empty() and dipoles_kernel != nullptr)
    obs_energy.dipolar[0] += (*dipoles_kernel)(p1, p2, d, dist, dist2);
#endif
}

inline boost::optional<double>
calc_bonded_energy(Bonded_IA_Parameters const &iaparams, Particle const &p1,
                   Utils::Span<Particle *> partners, BoxGeometry const &box_geo,
                   Coulomb::ShortRangeEnergyKernel::kernel_type const *kernel) {
  auto const n_partners = static_cast<int>(partners.size());

  auto p2 = (n_partners > 0) ? partners[0] : nullptr;
  auto p3 = (n_partners > 1) ? partners[1] : nullptr;
  auto p4 = (n_partners > 2) ? partners[2] : nullptr;

  if (n_partners == 1) {
    auto const dx = box_geo.get_mi_vector(p1.pos(), p2->pos());
    if (auto const *iap = boost::get<FeneBond>(&iaparams)) {
      return iap->energy(dx);
    }
    if (auto const *iap = boost::get<HarmonicBond>(&iaparams)) {
      return iap->energy(dx);
    }
    if (auto const *iap = boost::get<QuarticBond>(&iaparams)) {
      return iap->energy(dx);
    }
#ifdef ELECTROSTATICS
    if (auto const *iap = boost::get<BondedCoulomb>(&iaparams)) {
      return iap->energy(p1.q() * p2->q(), dx);
    }
    if (auto const *iap = boost::get<BondedCoulombSR>(&iaparams)) {
      return iap->energy(p1, *p2, dx, *kernel);
    }
#endif
#ifdef BOND_CONSTRAINT
    if (boost::get<RigidBond>(&iaparams)) {
      return {0.};
    }
#endif
#ifdef TABULATED
    if (auto const *iap = boost::get<TabulatedDistanceBond>(&iaparams)) {
      return iap->energy(dx);
    }
#endif
    if (boost::get<VirtualBond>(&iaparams)) {
      return {0.};
    }
    throw BondUnknownTypeError();
  } // 1 partner
  if (n_partners == 2) {
    auto const vec1 = box_geo.get_mi_vector(p2->pos(), p1.pos());
    auto const vec2 = box_geo.get_mi_vector(p3->pos(), p1.pos());
    if (auto const *iap = boost::get<AngleHarmonicBond>(&iaparams)) {
      return iap->energy(vec1, vec2);
    }
    if (auto const *iap = boost::get<AngleCosineBond>(&iaparams)) {
      return iap->energy(vec1, vec2);
    }
    if (auto const *iap = boost::get<AngleCossquareBond>(&iaparams)) {
      return iap->energy(vec1, vec2);
    }
    if (auto const *iap = boost::get<TabulatedAngleBond>(&iaparams)) {
      return iap->energy(vec1, vec2);
    }
    if (boost::get<IBMTriel>(&iaparams)) {
      runtimeWarningMsg() << "Unsupported bond type " +
                                 std::to_string(iaparams.which()) +
                                 " in energy calculation.";
      return 0.;
    }
    throw BondUnknownTypeError();
  } // 2 partners
  if (n_partners == 3) {
    // note: particles in a dihedral bond are ordered as p2-p1-p3-p4
    auto const v12 = box_geo.get_mi_vector(p1.pos(), p2->pos());
    auto const v23 = box_geo.get_mi_vector(p3->pos(), p1.pos());
    auto const v34 = box_geo.get_mi_vector(p4->pos(), p3->pos());
    if (auto const *iap = boost::get<DihedralBond>(&iaparams)) {
      return iap->energy(v12, v23, v34);
    }
    if (auto const *iap = boost::get<TabulatedDihedralBond>(&iaparams)) {
      return iap->energy(v12, v23, v34);
    }
    if (boost::get<IBMTribend>(&iaparams)) {
      runtimeWarningMsg() << "Unsupported bond type " +
                                 std::to_string(iaparams.which()) +
                                 " in energy calculation.";
      return 0.;
    }
    throw BondUnknownTypeError();
  } // 3 partners
  if (n_partners == 0) {
    return 0.;
  }

  throw BondInvalidSizeError(n_partners);
}

/** Calculate kinetic energies from translation for one particle.
 *  @param p   particle for which to calculate energies
 */
inline double translational_kinetic_energy(Particle const &p) {
  return p.is_virtual() ? 0. : 0.5 * p.mass() * p.v().norm2();
}

/** Calculate kinetic energies from rotation for one particle.
 *  @param p   particle for which to calculate energies
 */
inline double rotational_kinetic_energy(Particle const &p) {
#ifdef ROTATION
  return p.can_rotate()
             ? 0.5 * (hadamard_product(p.omega(), p.omega()) * p.rinertia())
             : 0.0;
#else
  return 0.0;
#endif
}

/** Calculate kinetic energies for one particle.
 *  @param p   particle for which to calculate energies
 */
inline double calc_kinetic_energy(Particle const &p) {
  return translational_kinetic_energy(p) + rotational_kinetic_energy(p);
}

#endif // CORE_ENERGY_INLINE_HPP
