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

#include "config/config.hpp"

#if defined(SCAFACOS) or defined(SCAFACOS_DIPOLES)

#include "scafacos/ScafacosContext.hpp"

#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives.hpp>

#include <cstddef>
#include <functional>

namespace detail {
std::tuple<Utils::Vector3d const &, Utils::Vector3i, std::size_t>
get_system_params() {
  auto &cell_structure = *System::get_system().cell_structure;
  auto periodicity = Utils::Vector3i{static_cast<int>(box_geo.periodic(0)),
                                     static_cast<int>(box_geo.periodic(1)),
                                     static_cast<int>(box_geo.periodic(2))};
  auto const n_part = boost::mpi::all_reduce(
      comm_cart, cell_structure.local_particles().size(), std::plus<>());
  return {box_geo.length(), periodicity, n_part};
}
} // namespace detail

#endif // SCAFACOS or SCAFACOS_DIPOLES
