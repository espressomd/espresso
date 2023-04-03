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
#include "WalberlaCheckpoint.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"

#include "core/BoxGeometry.hpp"
#include "core/grid.hpp"

#include "script_interface/communication.hpp"

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/optional.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {

void FluidWalberla::load_checkpoint(std::string const &filename, int mode) {
  auto &lb_obj = *m_lb_fluid;

  auto const read_metadata = [&lb_obj](CheckpointFile &cpfile) {
    auto const expected_grid_size = lb_obj.get_lattice().get_grid_dimensions();
    auto const expected_pop_size = lb_obj.stencil_size();
    Utils::Vector3i read_grid_size;
    std::size_t read_pop_size;
    cpfile.read(read_grid_size);
    cpfile.read(read_pop_size);
    if (read_grid_size != expected_grid_size) {
      std::stringstream message;
      message << "grid dimensions mismatch, "
              << "read [" << read_grid_size << "], "
              << "expected [" << expected_grid_size << "].";
      throw std::runtime_error(message.str());
    }
    if (read_pop_size != expected_pop_size) {
      throw std::runtime_error("population size mismatch, read " +
                               std::to_string(read_pop_size) + ", expected " +
                               std::to_string(expected_pop_size) + ".");
    }
  };

  auto const read_data = [&lb_obj](CheckpointFile &cpfile) {
    auto const grid_size = lb_obj.get_lattice().get_grid_dimensions();
    auto const i_max = grid_size[0];
    auto const j_max = grid_size[1];
    auto const k_max = grid_size[2];
    LBWalberlaNodeState cpnode;
    cpnode.populations.resize(lb_obj.stencil_size());
    for (int i = 0; i < i_max; i++) {
      for (int j = 0; j < j_max; j++) {
        for (int k = 0; k < k_max; k++) {
          auto const ind = Utils::Vector3i{{i, j, k}};
          cpfile.read(cpnode.populations);
          cpfile.read(cpnode.last_applied_force);
          cpfile.read(cpnode.is_boundary);
          if (cpnode.is_boundary) {
            cpfile.read(cpnode.slip_velocity);
          }
          lb_obj.set_node_population(ind, cpnode.populations);
          lb_obj.set_node_last_applied_force(ind, cpnode.last_applied_force);
          if (cpnode.is_boundary) {
            lb_obj.set_node_velocity_at_boundary(ind, cpnode.slip_velocity);
          }
        }
      }
    }
  };

  auto const on_success = [&lb_obj]() {
    lb_obj.reallocate_ubb_field();
    lb_obj.ghost_communication();
  };

  load_checkpoint_common(*context(), "LB", filename, mode, read_metadata,
                         read_data, on_success);
}

void FluidWalberla::save_checkpoint(std::string const &filename, int mode) {
  auto &lb_obj = *m_lb_fluid;

  auto const write_metadata = [&lb_obj,
                               mode](std::shared_ptr<CheckpointFile> cpfile_ptr,
                                     Context const &context) {
    auto const grid_size = lb_obj.get_lattice().get_grid_dimensions();
    auto const pop_size = lb_obj.stencil_size();
    if (context.is_head_node()) {
      cpfile_ptr->write(grid_size);
      cpfile_ptr->write(pop_size);
      unit_test_handle(mode);
    }
  };

  auto const on_failure = [](std::shared_ptr<CheckpointFile> const &,
                             Context const &context) {
    if (context.is_head_node()) {
      auto failure = true;
      boost::mpi::broadcast(context.get_comm(), failure, 0);
    }
  };

  auto const write_data = [&lb_obj,
                           mode](std::shared_ptr<CheckpointFile> cpfile_ptr,
                                 Context const &context) {
    auto const get_node_checkpoint = [&](Utils::Vector3i const &ind)
        -> boost::optional<LBWalberlaNodeState> {
      auto const pop = lb_obj.get_node_population(ind);
      auto const laf = lb_obj.get_node_last_applied_force(ind);
      auto const lbb = lb_obj.get_node_is_boundary(ind);
      auto const vbb = lb_obj.get_node_velocity_at_boundary(ind);
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
    };

    auto failure = false;
    auto const &comm = context.get_comm();
    auto const is_head_node = context.is_head_node();
    auto const unit_test_mode = (mode != static_cast<int>(CptMode::ascii)) and
                                (mode != static_cast<int>(CptMode::binary));
    auto const grid_size = lb_obj.get_lattice().get_grid_dimensions();
    auto const i_max = grid_size[0];
    auto const j_max = grid_size[1];
    auto const k_max = grid_size[2];
    LBWalberlaNodeState cpnode;
    for (int i = 0; i < i_max; i++) {
      for (int j = 0; j < j_max; j++) {
        for (int k = 0; k < k_max; k++) {
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
            auto &cpfile = *cpfile_ptr;
            cpfile.write(cpnode.populations);
            cpfile.write(cpnode.last_applied_force);
            cpfile.write(cpnode.is_boundary);
            if (cpnode.is_boundary) {
              cpfile.write(cpnode.slip_velocity);
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
  };

  save_checkpoint_common(*context(), "LB", filename, mode, write_metadata,
                         write_data, on_failure);
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
