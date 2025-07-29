/*
 * Copyright (C) 2025 The ESPResSo project
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

#include "BoxGeometry.hpp"
#include "lb/Solver.hpp"

#include <utils/Vector.hpp>

#include <cstddef>

namespace Observables {

/** @brief Bookkeeping of the box gometry and LB object memory address. */
class SanityChecksLB {
  Utils::Vector3d m_box_l = Utils::Vector3d::broadcast(0.);
  double m_agrid = 0.;
  std::size_t m_lb_obj_tag = 0ul;

  std::size_t get_lb_obj_tag(LB::Solver const &solver) const {
    std::size_t lb_obj_tag{0u};
    solver.connect([&lb_obj_tag](auto const &lb_obj) {
      lb_obj_tag = reinterpret_cast<std::size_t>(&lb_obj);
    });
    return lb_obj_tag;
  }

public:
  SanityChecksLB() = default;
  SanityChecksLB(BoxGeometry const &box_geo, LB::Solver const &solver) {
    m_box_l = box_geo.length();
    m_agrid = solver.get_agrid();
    m_lb_obj_tag = get_lb_obj_tag(solver);
  }

  bool mismatch(BoxGeometry const &box_geo, LB::Solver const &solver) const {
    if (m_lb_obj_tag == 0u) {
      return true;
    }
    auto const same_lb_obj = get_lb_obj_tag(solver) == m_lb_obj_tag;
    auto const same_box_l = m_box_l == box_geo.length();
    auto const same_agrid = m_agrid == solver.get_agrid();
    return not(same_lb_obj and same_box_l and same_agrid);
  }
};

} // namespace Observables
