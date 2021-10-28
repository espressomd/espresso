/*
 * Copyright (C) 2010-2021 The ESPResSo project
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

#include "cells.hpp"
#include "communication.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/optional.hpp>

#include <stdexcept>

Utils::Vector3d get_ibm_particle_position(int pid) {
  auto *p = cell_structure.get_local_particle(pid);
  boost::optional<Particle> opt_part{boost::none};

  if (p and not p->l.ghost) {
    opt_part = *p;
  }
  opt_part = boost::mpi::all_reduce(comm_cart, opt_part,
                                    [](boost::optional<Particle> const &acc,
                                       boost::optional<Particle> const &item) {
                                      if (acc) {
                                        return acc;
                                      }
                                      return item;
                                    });
  if (opt_part)
    return opt_part.get().r.p;
  throw std::runtime_error("Immersed Boundary: Particle not found");
}