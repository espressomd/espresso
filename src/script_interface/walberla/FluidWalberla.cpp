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
#include "config/config.hpp"

#ifdef WALBERLA

#include <walberla_bridge/lattice_boltzmann/LBWalberlaNodeState.hpp>

#include "FluidWalberla.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"

#include "core/BoxGeometry.hpp"
#include "core/grid.hpp"
#include "core/grid_based_algorithms/LBCheckpointFile.hpp"

#include "script_interface/communication.hpp"

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/optional.hpp>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

boost::optional<LBWalberlaNodeState>
FluidWalberla::get_node_checkpoint(Utils::Vector3i const &ind) const {
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

/** Inject code for unit tests. */
inline void unit_test_handle(int mode) {
  switch (mode) {
  case static_cast<int>(CptMode::ascii):
  case static_cast<int>(CptMode::binary):
    return;
  case static_cast<int>(CptMode::unit_test_runtime_error):
    throw std::runtime_error("unit test error");
  case static_cast<int>(CptMode::unit_test_ios_failure):
    throw std::ios_base::failure("unit test error");
  default:
    throw std::domain_error("Unknown mode " + std::to_string(mode));
  }
}

void FluidWalberla::load_checkpoint(std::string const &filename, int mode) {
  auto const err_msg = std::string("Error while reading LB checkpoint: ");
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
    auto const grid_size = m_lb_fluid->get_lattice().get_grid_dimensions();
    auto const pop_size = m_lb_fluid->stencil_size();
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
          m_lb_fluid->set_node_pop(ind, cpnode.populations);
          m_lb_fluid->set_node_last_applied_force(ind,
                                                  cpnode.last_applied_force);
          if (cpnode.is_boundary) {
            m_lb_fluid->set_node_velocity_at_boundary(ind,
                                                      cpnode.slip_velocity);
          }
        }
      }
    }
    comm.barrier();
    m_lb_fluid->reallocate_ubb_field();
    m_lb_fluid->ghost_communication();
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

void FluidWalberla::save_checkpoint(std::string const &filename, int mode) {
  auto const err_msg = std::string("Error while writing LB checkpoint: ");
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
    auto const grid_size = m_lb_fluid->get_lattice().get_grid_dimensions();
    auto const pop_size = m_lb_fluid->stencil_size();
    if (is_head_node) {
      cpfile->write(grid_size);
      cpfile->write(pop_size);
      unit_test_handle(mode);
    }

    LBWalberlaNodeState cpnode;
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
            cpfile->write(cpnode.populations);
            cpfile->write(cpnode.last_applied_force);
            cpfile->write(cpnode.is_boundary);
            if (cpnode.is_boundary) {
              cpfile->write(cpnode.slip_velocity);
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

std::vector<Variant> FluidWalberla::get_average_pressure_tensor() const {
  auto const local = m_lb_fluid->get_pressure_tensor() / m_conv_press;
  auto const tensor_flat = mpi_reduce_sum(context()->get_comm(), local);
  auto tensor = Utils::Matrix<double, 3, 3>{};
  std::copy(tensor_flat.begin(), tensor_flat.end(), tensor.m_data.begin());
  return std::vector<Variant>{tensor.row<0>().as_vector(),
                              tensor.row<1>().as_vector(),
                              tensor.row<2>().as_vector()};
}

Variant
FluidWalberla::get_interpolated_velocity(Utils::Vector3d const &pos) const {
  auto const lb_pos = folded_position(pos, box_geo) * m_conv_dist;
  auto const result = m_lb_fluid->get_velocity_at_pos(lb_pos);
  return mpi_reduce_optional(context()->get_comm(), result) / m_conv_speed;
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
