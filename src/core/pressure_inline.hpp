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

#include <config/config.hpp>

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "errorhandling.hpp"
#include "exclusions.hpp"
#include "forces_inline.hpp"

#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>

#include <cstdio>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <variant>

// Overload for Cabana kernels: takes positions and charge product directly
inline std::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_virial_pressure_tensor(
    Bonded_IA_Parameters const &iaparams, Utils::Vector3d const &pos1,
    Utils::Vector3d const &pos2, BoxGeometry const &box_geo,
    Coulomb::ShortRangeForceKernel::kernel_type const *kernel, double q1q2) {
  auto const dx = box_geo.get_mi_vector(pos1, pos2);
  auto const pair_force = calc_bond_pair_force(iaparams, dx,
#ifdef ESPRESSO_ELECTROSTATICS
                                               q1q2, kernel
#else
                                               0.0, nullptr
#endif
  );
  std::optional<Utils::Matrix<double, 3, 3>> pressure{std::nullopt};
  if (pair_force) {
    pressure = Utils::tensor_product(dx, *pair_force);
  }
  return pressure;
}

inline std::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_virial_pressure_tensor(
    Bonded_IA_Parameters const &iaparams, Particle const &p1,
    Particle const &p2, BoxGeometry const &box_geo,
    Coulomb::ShortRangeForceKernel::kernel_type const *kernel) {
  return calc_bonded_virial_pressure_tensor(iaparams, p1.pos(), p2.pos(),
                                            box_geo, kernel,
#ifdef ESPRESSO_ELECTROSTATICS
                                            p1.q() * p2.q()
#else
                                            0.0
#endif
  );
}

// Overload for Cabana kernels: takes positions directly
inline std::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_three_body_pressure_tensor(Bonded_IA_Parameters const &iaparams,
                                       Utils::Vector3d const &pos1,
                                       Utils::Vector3d const &pos2,
                                       Utils::Vector3d const &pos3,
                                       BoxGeometry const &box_geo) {
  if (std::holds_alternative<AngleHarmonicBond>(iaparams) or
      std::holds_alternative<AngleCosineBond>(iaparams) or
#ifdef ESPRESSO_TABULATED
      std::holds_alternative<TabulatedAngleBond>(iaparams) or
#endif
      std::holds_alternative<AngleCossquareBond>(iaparams)) {
    auto const dx21 = -box_geo.get_mi_vector(pos1, pos2);
    auto const dx31 = box_geo.get_mi_vector(pos3, pos1);

    auto const result = calc_bonded_three_body_force(iaparams, dx21, dx31);
    if (result) {
      Utils::Vector3d force2, force3;
      std::tie(std::ignore, force2, force3) = result.value();
      return Utils::tensor_product(dx21, force2) +
             Utils::tensor_product(dx31, force3);
    }
  } else {
    runtimeWarningMsg() << "Unsupported bond type " +
                               std::to_string(iaparams.index()) +
                               " in pressure calculation.";
    return Utils::Matrix<double, 3, 3>{};
  }
  return {};
}

inline std::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_four_body_pressure_tensor(Bonded_IA_Parameters const &iaparams,
                                      Utils::Vector3d const &pos1,
                                      Utils::Vector3d const &pos2,
                                      Utils::Vector3d const &pos3,
                                      Utils::Vector3d const &pos4,
                                      BoxGeometry const &box_geo) {
  if (std::holds_alternative<DihedralBond>(iaparams)
#ifdef ESPRESSO_TABULATED
      or std::holds_alternative<TabulatedDihedralBond>(iaparams)
#endif
  ) {
    auto const v12 = box_geo.get_mi_vector(pos1, pos2);
    auto const v23 = box_geo.get_mi_vector(pos3, pos1);
    auto const v34 = box_geo.get_mi_vector(pos4, pos3);

    auto const result = calc_bonded_dihedral_force(iaparams, v12, v23, v34);

    if (result) {
      Utils::Vector3d force2, force3, force4;
      std::tie(std::ignore, force2, force3, force4) = result.value();

      return -Utils::tensor_product(v12, force2) +
             Utils::tensor_product(v23, force3) +
             Utils::tensor_product(v23 + v34, force4);
    }
  } else {
    runtimeWarningMsg() << "Unsupported bond type " +
                               std::to_string(iaparams.index()) +
                               " in pressure calculation.";
    return Utils::Matrix<double, 3, 3>{};
  }

  return {};
}
