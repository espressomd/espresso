/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include "ibm_common.hpp"

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/serialization/optional.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <optional>
#include <stdexcept>

Utils::Vector3d get_ibm_particle_position(CellStructure const &cell_structure,
                                          int pid) {
  auto const p = cell_structure.get_local_particle(pid);
  // A Particle is a non-owning view (its store pointer is meaningless on
  // another rank), so reduce the position value across ranks. The owning rank
  // contributes its position; every other rank contributes nullopt.
  std::optional<Utils::Vector3d> opt_pos{std::nullopt};

  if (p and not p->is_ghost()) {
    opt_pos = p->pos();
  }
  opt_pos =
      boost::mpi::all_reduce(comm_cart, opt_pos,
                             [](std::optional<Utils::Vector3d> const &acc,
                                std::optional<Utils::Vector3d> const &item) {
                               if (acc) {
                                 return acc;
                               }
                               return item;
                             });
  if (opt_pos)
    return opt_pos.value();
  throw std::runtime_error("Immersed Boundary: Particle not found");
}
