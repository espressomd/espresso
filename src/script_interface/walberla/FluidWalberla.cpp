/*
 * Copyright (C) 2021 The ESPResSo project
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
#include "config.hpp"

#ifdef LB_WALBERLA

#include <walberla_bridge/LBWalberlaNodeState.hpp>

#include "LatticeWalberla.hpp"

#include "core/communication.hpp"
#include "core/grid_based_algorithms/lb_walberla_instance.hpp"
#include "core/grid_based_algorithms/LBCheckpointFile.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

void do_reallocate_ubb_field() { lb_walberla()->reallocate_ubb_field(); }

REGISTER_CALLBACK(do_reallocate_ubb_field)

void do_ghost_communication() { lb_walberla()->ghost_communication(); }

REGISTER_CALLBACK(do_ghost_communication)

void set_node_from_checkpoint(Utils::Vector3i ind, LBWalberlaNodeState cpt) {
  lb_walberla()->set_node_pop(ind, cpt.populations);
  lb_walberla()->set_node_last_applied_force(ind, cpt.last_applied_force);
  if (cpt.is_boundary) {
    lb_walberla()->set_node_velocity_at_boundary(ind, cpt.slip_velocity, false);
  }
}

REGISTER_CALLBACK(set_node_from_checkpoint)

boost::optional<LBWalberlaNodeState> get_node_checkpoint(Utils::Vector3i ind) {
  auto const pop = lb_walberla()->get_node_pop(ind);
  auto const laf = lb_walberla()->get_node_last_applied_force(ind);
  auto const lbb = lb_walberla()->get_node_is_boundary(ind);
  auto const vbb = lb_walberla()->get_node_velocity_at_boundary(ind);
  if (pop and laf and lbb and ((*lbb) ? vbb.has_value() : true)) {
    LBWalberlaNodeState cpnode;
    cpnode.populations = *pop;
    cpnode.last_applied_force = *laf;
    cpnode.is_boundary = *lbb;
    if (*lbb) {
      cpnode.slip_velocity = *vbb;
    }
    return cpnode;
  }
  return {boost::none};
}

REGISTER_CALLBACK_ONE_RANK(get_node_checkpoint)

namespace ScriptInterface::walberla {

void load_checkpoint(const std::string &filename, bool binary) {
  auto const err_msg = std::string("Error while reading LB checkpoint: ");

  // open file and set exceptions
  LBCheckpointFile cpfile(filename, std::ios_base::in, binary);
  if (!cpfile.stream) {
    throw std::runtime_error(err_msg + "could not open file " + filename);
  }
  cpfile.stream.exceptions(std::ios_base::failbit | std::ios_base::badbit);

  // check the grid size in the checkpoint header matches the current grid size
  auto const check_header = [&](Utils::Vector3i const &expected_grid_size,
                                std::size_t expected_pop_size) {
    Utils::Vector3i grid_size;
    std::size_t pop_size;
    cpfile.read(grid_size);
    cpfile.read(pop_size);
    if (grid_size != expected_grid_size) {
      std::stringstream message;
      message << " grid dimensions mismatch,"
              << " read [" << grid_size << "],"
              << " expected [" << expected_grid_size << "].";
      throw std::runtime_error(err_msg + message.str());
    }
    if (pop_size != expected_pop_size) {
      throw std::runtime_error(err_msg + "population size mismatch, read " +
                               std::to_string(pop_size) + ", expected " +
                               std::to_string(expected_pop_size) + ".");
    }
  };

  try {
      auto const grid_size = lb_walberla()->lattice().get_grid_dimensions();
      auto const pop_size = lb_walberla()->stencil_size();
      check_header(grid_size, pop_size);

      LBWalberlaNodeState cpnode;
      cpnode.populations.resize(pop_size);
      for (int i = 0; i < grid_size[0]; i++) {
        for (int j = 0; j < grid_size[1]; j++) {
          for (int k = 0; k < grid_size[2]; k++) {
            auto const ind = Utils::Vector3i{{i, j, k}};
            cpfile.read(cpnode.populations);
            cpfile.read(cpnode.last_applied_force);
            cpfile.read(cpnode.is_boundary);
            if (cpnode.is_boundary) {
              cpfile.read(cpnode.slip_velocity);
            }
            ::Communication::mpiCallbacks().call_all(
                set_node_from_checkpoint, ind, cpnode);
          }
        }
      }
      ::Communication::mpiCallbacks().call_all(
          do_reallocate_ubb_field);
      ::Communication::mpiCallbacks().call_all(
          do_ghost_communication);
    // check EOF
    if (!binary) {
      if (cpfile.stream.peek() == '\n') {
        static_cast<void>(cpfile.stream.get());
      }
    }
    if (cpfile.stream.peek() != EOF) {
      throw std::runtime_error(err_msg + "extra data found, expected EOF.");
    }
  } catch (std::ios_base::failure const &fail) {
    auto const eof_error = cpfile.stream.eof();
    cpfile.stream.close();
    if (eof_error) {
      throw std::runtime_error(err_msg + "EOF found.");
    }
    throw std::runtime_error(err_msg + "incorrectly formatted data.");
  } catch (std::runtime_error const &fail) {
    cpfile.stream.close();
    throw;
  }
}

void save_checkpoint(const std::string &filename, bool binary) {
  auto const err_msg = std::string("Error while writing LB checkpoint: ");

  // open file and set exceptions
  LBCheckpointFile cpfile(filename, std::ios_base::out, binary);
  if (!cpfile.stream) {
    throw std::runtime_error(err_msg + "could not open file " + filename);
  }
  cpfile.stream.exceptions(std::ios_base::failbit | std::ios_base::badbit);

  try {
      if (!binary) {
        cpfile.stream.precision(16);
        cpfile.stream << std::fixed;
      }

      auto const grid_size = lb_walberla()->lattice().get_grid_dimensions();
      auto const pop_size = lb_walberla()->stencil_size();
      cpfile.write(grid_size);
      cpfile.write(pop_size);

      for (int i = 0; i < grid_size[0]; i++) {
        for (int j = 0; j < grid_size[1]; j++) {
          for (int k = 0; k < grid_size[2]; k++) {
            auto const ind = Utils::Vector3i{{i, j, k}};
            auto const cpnode = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                get_node_checkpoint, ind);
            cpfile.write(cpnode.populations);
            cpfile.write(cpnode.last_applied_force);
            cpfile.write(cpnode.is_boundary);
            if (cpnode.is_boundary) {
              cpfile.write(cpnode.slip_velocity);
            }
          }
        }
      }
  } catch (std::ios_base::failure const &fail) {
    cpfile.stream.close();
    throw std::runtime_error(err_msg + "could not write data to " + filename);
  } catch (std::runtime_error const &fail) {
    cpfile.stream.close();
    throw;
  }
}

} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
