/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include <boost/mpi/communicator.hpp>

#include <unordered_map>

class ExclusionRadius {
  boost::mpi::communicator const &m_comm;
  double m_max_exclusion_range = 0.;
  void recalc_derived_parameters();

public:
  using map_type = std::unordered_map<int, double>;
  explicit ExclusionRadius(boost::mpi::communicator const &comm)
      : m_comm{comm} {}
  void set_exclusion_range(double range);
  void set_exclusion_radius_per_type(map_type const &map);
  bool check_exclusion_range(int p_id, int p_type);
  bool check_exclusion_range(int pid);

  bool neighbor_search_order_n = true;
  double exclusion_range = 0.;
  map_type exclusion_radius_per_type{};
};
