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
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"

#include <utils/Vector.hpp>

#include <tuple>

/** Parameters for IBM tribend */
struct IBMTribend {
  /**
   * @brief Bare bending modulus.
   *
   * If triangle pairs appear only once, the total bending force should get a
   * factor 2. For the numerical model, a factor @f$ \sqrt(3) @f$ should be
   * added, see @cite gompper96a and @cite kruger12a. This is an approximation,
   * it holds strictly only for a sphere
   */
  double kb;

  /** Reference angle */
  double theta0;

  /** Particle ids */
  std::tuple<int, int, int, int> p_ids;

  bool flat;
  bool is_initialized;

  double cutoff() const { return 0.; }

  static constexpr int num = 3;

  /** Set the IBM Tribend parameters.
   *  Also calculate and store the reference state.
   *  See details in @cite gompper96a and @cite kruger12a.
   */
  void initialize(BoxGeometry const &box_geo,
                  CellStructure const &cell_structure);

  IBMTribend(int ind1, int ind2, int ind3, int ind4, double kb, bool flat)
      : kb{kb}, theta0{0.}, p_ids{ind1, ind2, ind3, ind4}, flat{flat},
        is_initialized{false} {}

  /** Calculate the forces
   *  The equations can be found in Appendix C of @cite kruger12a.
   *  @return forces on @p p1, @p p2, @p p3, @p p4
   */
  std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
  calc_forces(BoxGeometry const &box_geo, Particle const &p1,
              Particle const &p2, Particle const &p3, Particle const &p4) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar & kb;
    ar & theta0;
  }
};
