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

#include "communication.hpp"

#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"

#ifdef WALBERLA
#include <walberla_bridge/walberla_init.hpp>
#endif

#include <utils/mpi/cart_comm.hpp>

#include <boost/mpi.hpp>

#include <mpi.h>

#include <cassert>
#include <memory>
#include <utility>

boost::mpi::communicator comm_cart;

namespace Communication {
static auto const &mpi_datatype_cache =
    boost::mpi::detail::mpi_datatype_cache();
static std::shared_ptr<boost::mpi::environment> mpi_env;
static std::unique_ptr<MpiCallbacks> m_callbacks;

/* We use a singleton callback class for now. */
MpiCallbacks &mpiCallbacks() {
  assert(m_callbacks && "Mpi not initialized!");

  return *m_callbacks;
}
} // namespace Communication

using Communication::mpiCallbacks;

int this_node = -1;
int n_nodes = -1;

namespace Communication {
void init(std::shared_ptr<boost::mpi::environment> mpi_env) {
  Communication::mpi_env = std::move(mpi_env);

  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  node_grid = Utils::Mpi::dims_create<3>(n_nodes);

  comm_cart =
      Utils::Mpi::cart_create(comm_cart, node_grid, /* reorder */ false);

  this_node = comm_cart.rank();

  Communication::m_callbacks =
      std::make_unique<Communication::MpiCallbacks>(comm_cart);

  ErrorHandling::init_error_handling(mpiCallbacks());

#ifdef WALBERLA
  walberla::mpi_init();
#endif

  on_program_start();
}
} // namespace Communication

std::shared_ptr<boost::mpi::environment> mpi_init(int argc, char **argv) {
  return std::make_shared<boost::mpi::environment>(argc, argv);
}

void mpi_loop() {
  if (this_node != 0)
    mpiCallbacks().loop();
}
