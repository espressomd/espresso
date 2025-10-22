/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "BoxGeometry.hpp"
#include "cell_system/CellStructure.hpp"

#include <utils/Vector.hpp>

#include <optional>
#include <tuple>

enum class tElasticLaw { NeoHookean, Skalak };

/** Parameters for IBM elastic triangle (triel) */
struct IBMTriel {
  // These values encode the reference state
  double l0;
  double lp0;
  double sinPhi0;
  double cosPhi0;
  double area0;

  // These values are cache values to speed up computation
  double a1;
  double a2;
  double b1;
  double b2;

  // These are interaction parameters
  // k1 is used for Neo-Hookean
  // k1 and k2 are used Skalak
  double maxDist;
  tElasticLaw elasticLaw;
  double k1;
  double k2;

  /** Particle ids */
  std::tuple<int, int, int> p_ids;
  bool is_initialized;

  double cutoff() const { return maxDist; }

  static constexpr int num = 2;

  /** Set the IBM Triel parameters.
   *  Also calculate and store the reference state.
   */
  void initialize(BoxGeometry const &box_geo,
                  CellStructure const &cell_structure);

  IBMTriel(int ind1, int ind2, int ind3, double maxDist, tElasticLaw elasticLaw,
           double k1, double k2)
      : l0{0.}, lp0{0.}, sinPhi0{0.}, cosPhi0{0.}, area0{0.}, a1{0.}, a2{0.},
        b1{0.}, b2{0.}, maxDist{maxDist}, elasticLaw{elasticLaw}, k1{k1},
        k2{k2}, p_ids{ind1, ind2, ind3}, is_initialized{false} {}

  /** Calculate the forces.
   *  The equations can be found in Appendix C of @cite kruger12a.
   *  @return the forces on @p p1, @p p2, @p p3
   */
  std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
  calc_forces(Utils::Vector3d const &vec1, Utils::Vector3d const &vec2) const;
};
