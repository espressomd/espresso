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

#pragma once

#include "config/config.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "errorhandling.hpp"
#include "exclusions.hpp"
#include "forces_inline.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <string>
#include <tuple>

/** Calculate non-bonded energies between a pair of particles.
 *  @param p1        pointer to particle 1.
 *  @param p2        pointer to particle 2.
 *  @param d         vector between p1 and p2.
 *  @param dist      distance between p1 and p2.
 *  @param ia_params              non-bonded interaction kernels.
 *  @param kernel_forces          Coulomb force kernel.
 *  @param kernel_pressure        Coulomb pressure kernel.
 *  @param[in,out] obs_pressure   pressure observable.
 */
inline void add_non_bonded_pair_virials(
    Particle const &p1, Particle const &p2, Utils::Vector3d const &d,
    double dist, IA_parameters const &ia_params,
    Coulomb::ShortRangeForceKernel::kernel_type const *kernel_forces,
    Coulomb::ShortRangePressureKernel::kernel_type const *kernel_pressure,
    Observable_stat &obs_pressure) {
#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
  {
    auto const force = calc_central_radial_force(p1, p2, ia_params, d, dist).f +
                       calc_central_radial_charge_force(p1, p2, ia_params, d,
                                                        dist, kernel_forces)
                           .f +
                       calc_non_central_force(p1, p2, ia_params, d, dist).f;
    auto const stress = Utils::tensor_product(d, force);
    obs_pressure.add_non_bonded_contribution(p1.type(), p2.type(), p1.mol_id(),
                                             p2.mol_id(), flatten(stress));
  }

#ifdef ELECTROSTATICS
  if (!obs_pressure.coulomb.empty() and kernel_pressure != nullptr) {
    /* real space Coulomb */
    auto const p_coulomb = (*kernel_pressure)(p1.q() * p2.q(), d, dist);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        obs_pressure.coulomb[i * 3 + j] += p_coulomb(i, j);
      }
    }
  }
#endif // ELECTROSTATICS

#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  if (Dipoles::get_dipoles().impl->solver) {
    fprintf(stderr, "calculating pressure for magnetostatics which doesn't "
                    "have it implemented\n");
  }
#endif // DIPOLES
}

inline boost::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_virial_pressure_tensor(
    Bonded_IA_Parameters const &iaparams, Particle const &p1,
    Particle const &p2, BoxGeometry const &box_geo,
    Coulomb::ShortRangeForceKernel::kernel_type const *kernel) {
  auto const dx = box_geo.get_mi_vector(p1.pos(), p2.pos());
  auto const result = calc_bond_pair_force(p1, p2, iaparams, dx, kernel);
  if (result) {
    auto const &force = result.get();

    return Utils::tensor_product(force, dx);
  }

  return {};
}

inline boost::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_three_body_pressure_tensor(Bonded_IA_Parameters const &iaparams,
                                       Particle const &p1, Particle const &p2,
                                       Particle const &p3,
                                       BoxGeometry const &box_geo) {
  if ((boost::get<AngleHarmonicBond>(&iaparams) != nullptr) ||
      (boost::get<AngleCosineBond>(&iaparams) != nullptr) ||
#ifdef TABULATED
      (boost::get<TabulatedAngleBond>(&iaparams) != nullptr) ||
#endif
      (boost::get<AngleCossquareBond>(&iaparams) != nullptr)) {
    auto const dx21 = -box_geo.get_mi_vector(p1.pos(), p2.pos());
    auto const dx31 = box_geo.get_mi_vector(p3.pos(), p1.pos());

    auto const result =
        calc_bonded_three_body_force(iaparams, box_geo, p1, p2, p3);
    if (result) {
      Utils::Vector3d force2, force3;
      std::tie(std::ignore, force2, force3) = result.get();

      return Utils::tensor_product(force2, dx21) +
             Utils::tensor_product(force3, dx31);
    }
  } else {
    runtimeWarningMsg() << "Unsupported bond type " +
                               std::to_string(iaparams.which()) +
                               " in pressure calculation.";
    return Utils::Matrix<double, 3, 3>{};
  }

  return {};
}

inline boost::optional<Utils::Matrix<double, 3, 3>> calc_bonded_pressure_tensor(
    Bonded_IA_Parameters const &iaparams, Particle const &p1,
    Utils::Span<Particle *> partners, BoxGeometry const &box_geo,
    Coulomb::ShortRangeForceKernel::kernel_type const *kernel) {
  switch (number_of_partners(iaparams)) {
  case 1:
    return calc_bonded_virial_pressure_tensor(iaparams, p1, *partners[0],
                                              box_geo, kernel);
  case 2:
    return calc_bonded_three_body_pressure_tensor(iaparams, p1, *partners[0],
                                                  *partners[1], box_geo);
  default:
    runtimeWarningMsg() << "Unsupported bond type " +
                               std::to_string(iaparams.which()) +
                               " in pressure calculation.";
    return Utils::Matrix<double, 3, 3>{};
  }
}

/** Calculate kinetic pressure (aka energy) for one particle.
 *  @param[in]   p1            particle for which to calculate pressure
 *  @param[out]  obs_pressure  pressure observable
 */
inline void add_kinetic_virials(Particle const &p1,
                                Observable_stat &obs_pressure) {
  if (p1.is_virtual())
    return;

  /* kinetic pressure */
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      obs_pressure.kinetic[k * 3 + l] += p1.v()[k] * p1.v()[l] * p1.mass();
}
