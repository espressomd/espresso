/*
 * Copyright (C) 2023 The ESPResSo project
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

#ifdef WALBERLA

#include <walberla_bridge/electrokinetics/EKWalberlaNodeState.hpp>

#include "EKSpecies.hpp"
#include "WalberlaCheckpoint.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"

#include "core/grid_based_algorithms/LBCheckpointFile.hpp"

#include "script_interface/communication.hpp"

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/optional.hpp>

namespace ScriptInterface::walberla {

boost::optional<EKWalberlaNodeState>
EKSpecies::get_node_checkpoint(Utils::Vector3i const &ind) const {
  auto const density = get_ekinstance()->get_node_density(ind);
  auto const is_b_d = get_ekinstance()->get_node_is_density_boundary(ind);
  auto const dens_b = get_ekinstance()->get_node_density_at_boundary(ind);
  auto const is_b_f = get_ekinstance()->get_node_is_flux_boundary(ind);
  auto const flux_b = get_ekinstance()->get_node_flux_at_boundary(ind);

  if (density and is_b_d and is_b_f and
      ((*is_b_d) ? dens_b.has_value() : true) and
      ((*is_b_f) ? flux_b.has_value() : true)) {
    EKWalberlaNodeState cpnode;
    cpnode.density = *density;
    cpnode.is_boundary_density = *is_b_d;
    if (*is_b_d) {
      cpnode.density_boundary = *dens_b;
    }
    cpnode.is_boundary_flux = *is_b_f;
    if (*is_b_f) {
      cpnode.flux_boundary = *flux_b;
    }
    return cpnode;
  }
  return {boost::none};
}

void EKSpecies::load_checkpoint(std::string const &filename, int mode) {
  auto const err_msg = std::string("Error while reading EK checkpoint: ");
  auto const binary = mode == static_cast<int>(CptMode::binary);
  auto const &comm = context()->get_comm();
  auto const is_head_node = context()->is_head_node();

  // open file and set exceptions
  LBCheckpointFile cpfile(filename, std::ios_base::in, binary);
  if (!cpfile.stream) {
    if (is_head_node) {
      throw std::runtime_error(err_msg + "could not open file " + filename);
    }
    return;
  }
  cpfile.stream.exceptions(std::ios_base::failbit | std::ios_base::badbit);

  // check the grid size in the checkpoint header matches the current grid size
  auto const check_header = [&](Utils::Vector3i const &expected_grid_size) {
    Utils::Vector3i grid_size;
    cpfile.read(grid_size);
    if (grid_size != expected_grid_size) {
      std::stringstream message;
      message << " grid dimensions mismatch,"
              << " read [" << grid_size << "],"
              << " expected [" << expected_grid_size << "].";
      throw std::runtime_error(err_msg + message.str());
    }
  };

  try {
    auto const grid_size = get_lattice()->lattice()->get_grid_dimensions();
    check_header(grid_size);

    EKWalberlaNodeState cpnode;
    for (int i = 0; i < grid_size[0]; i++) {
      for (int j = 0; j < grid_size[1]; j++) {
        for (int k = 0; k < grid_size[2]; k++) {
          auto const ind = Utils::Vector3i{{i, j, k}};
          cpfile.read(cpnode.density);
          cpfile.read(cpnode.is_boundary_density);
          if (cpnode.is_boundary_density) {
            cpfile.read(cpnode.density_boundary);
          }
          cpfile.read(cpnode.is_boundary_flux);
          if (cpnode.is_boundary_flux) {
            cpfile.read(cpnode.flux_boundary);
          }

          get_ekinstance()->set_node_density(ind, cpnode.density);
          if (cpnode.is_boundary_density) {
            get_ekinstance()->set_node_density_boundary(
                ind, cpnode.density_boundary);
          }
          if (cpnode.is_boundary_flux) {
            get_ekinstance()->set_node_flux_boundary(ind, cpnode.flux_boundary);
          }
        }
      }
    }
    comm.barrier();
    get_ekinstance()->ghost_communication();
    // check EOF
    if (!binary) {
      if (cpfile.stream.peek() == '\n') {
        static_cast<void>(cpfile.stream.get());
      }
    }
    if (cpfile.stream.peek() != EOF) {
      throw std::runtime_error(err_msg + "extra data found, expected EOF.");
    }
  } catch (std::ios_base::failure const &) {
    auto const eof_error = cpfile.stream.eof();
    cpfile.stream.close();
    if (eof_error) {
      if (is_head_node) {
        throw std::runtime_error(err_msg + "EOF found.");
      }
      return;
    }
    if (is_head_node) {
      throw std::runtime_error(err_msg + "incorrectly formatted data.");
    }
    return;
  } catch (std::runtime_error const &) {
    cpfile.stream.close();
    if (is_head_node) {
      throw;
    }
    return;
  }
}

void EKSpecies::save_checkpoint(std::string const &filename, int mode) {
  auto const err_msg = std::string("Error while writing EK checkpoint: ");
  auto const binary = mode == static_cast<int>(CptMode::binary);
  auto const unit_test_mode = (mode != static_cast<int>(CptMode::ascii)) and
                              (mode != static_cast<int>(CptMode::binary));
  auto const &comm = context()->get_comm();
  auto const is_head_node = context()->is_head_node();
  bool failure = false;

  // open file and set exceptions
  std::shared_ptr<LBCheckpointFile> cpfile;
  if (is_head_node) {
    cpfile = std::make_shared<LBCheckpointFile>(filename, std::ios_base::out,
                                                binary);
    failure = !cpfile->stream;
    boost::mpi::broadcast(comm, failure, 0);
    if (failure) {
      throw std::runtime_error(err_msg + "could not open file " + filename);
    }
    cpfile->stream.exceptions(std::ios_base::failbit | std::ios_base::badbit);
    if (!binary) {
      cpfile->stream.precision(16);
      cpfile->stream << std::fixed;
    }
  } else {
    boost::mpi::broadcast(comm, failure, 0);
    if (failure) {
      return;
    }
  }

  try {
    auto const grid_size = get_lattice()->lattice()->get_grid_dimensions();
    if (is_head_node) {
      cpfile->write(grid_size);
      unit_test_handle(mode);
    }

    EKWalberlaNodeState cpnode;
    for (int i = 0; i < grid_size[0]; i++) {
      for (int j = 0; j < grid_size[1]; j++) {
        for (int k = 0; k < grid_size[2]; k++) {
          auto const ind = Utils::Vector3i{{i, j, k}};
          auto const result = get_node_checkpoint(ind);
          if (!unit_test_mode) {
            assert(1 == boost::mpi::all_reduce(comm, static_cast<int>(!!result),
                                               std::plus<>()) &&
                   "Incorrect number of return values");
          }
          if (is_head_node) {
            if (result) {
              cpnode = *result;
            } else {
              comm.recv(boost::mpi::any_source, 42, cpnode);
            }
            cpfile->write(cpnode.density);
            cpfile->write(cpnode.is_boundary_density);
            if (cpnode.is_boundary_density) {
              cpfile->write(cpnode.density_boundary);
            }
            cpfile->write(cpnode.is_boundary_flux);
            if (cpnode.is_boundary_flux) {
              cpfile->write(cpnode.flux_boundary);
            }
            boost::mpi::broadcast(comm, failure, 0);
          } else {
            if (result) {
              comm.send(0, 42, *result);
            }
            boost::mpi::broadcast(comm, failure, 0);
            if (failure) {
              return;
            }
          }
        }
      }
    }
  } catch (std::exception const &error) {
    if (is_head_node) {
      failure = true;
      boost::mpi::broadcast(comm, failure, 0);
      cpfile->stream.close();
      if (dynamic_cast<std::ios_base::failure const *>(&error)) {
        throw std::runtime_error(err_msg + "could not write to " + filename);
      }
      throw;
    }
  }
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
