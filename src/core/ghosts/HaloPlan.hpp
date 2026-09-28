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

#include "cell_system/Cell.hpp"
#include <boost/mpi/communicator.hpp>
#include <optional>
#include <utils/Vector.hpp>
#include <vector>

namespace GhostComm {
enum class Direction { Push, Reduce }; // real->ghost ; ghost->real
enum class Combine { Overwrite, Add };
struct ExchangeOp {
  Direction direction;
  Combine combine;
};

struct SendRegion {
  Cell *cell;
  Utils::Vector3d shift;
};
struct NeighborComm {
  int peer;
  std::vector<SendRegion> send;
  std::vector<Cell *> recv; // recv[k] <-> peer.send[k]
};
struct LocalComm {
  Cell *src;
  Cell *dst;
  Utils::Vector3d shift;
};

enum class CollectivePattern { None, Broadcast, ReduceSum };
struct CollectiveSection {
  CollectivePattern pattern;
  std::vector<Cell *> cells;
};

struct HaloPlan {
  boost::mpi::communicator comm;
  std::vector<NeighborComm> neighbors;
  std::vector<LocalComm> local;
  std::optional<CollectiveSection> collective;
  /**
   * @brief Lazily-set flag for the once-per-plan symmetry check.
   *
   * Set to true on the first call to halo_exchange_start for this plan
   * (see HaloExchange.cpp, ESPRESSO_ADDITIONAL_CHECKS guard).  The plan
   * object is rebuilt on every decomposition change, so the flag naturally
   * rearms: a new plan always starts with symmetry_validated = false.
   *
   * mutable so halo_exchange_start can set it through a const HaloPlan&
   * without requiring a non-const overload.
   */
  mutable bool symmetry_validated = false;
};
} // namespace GhostComm
