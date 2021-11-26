/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include "grid_based_algorithms/lb_walberla_interface.hpp"

#include "BoxGeometry.hpp"
#include "MpiCallbacks.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>

#include <boost/serialization/vector.hpp>

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

ActiveLB lattice_switch = ActiveLB::NONE;

struct NoLBActive : public std::exception {
  const char *what() const noexcept override { return "LB not activated"; }
};

void lb_lbfluid_init() {}

void lb_lbfluid_propagate() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    lb_walberla()->integrate();
#endif
  }
}

void lb_lbfluid_sanity_checks(double time_step) {
  if (lattice_switch == ActiveLB::NONE)
    return;

  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    lb_sanity_checks(*lb_walberla(), *lb_walberla_params(), time_step);
#endif
  }
}

double lb_lbfluid_get_agrid() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    return lb_walberla_params()->get_agrid();
#endif
  }
  throw NoLBActive();
}

void check_tau_time_step_consistency(double tau, double time_step) {
  auto const eps = std::numeric_limits<float>::epsilon();
  if ((tau - time_step) / (tau + time_step) < -eps)
    throw std::invalid_argument("LB tau (" + std::to_string(tau) +
                                ") must be >= MD time_step (" +
                                std::to_string(time_step) + ")");
  auto const factor = tau / time_step;
  if (fabs(round(factor) - factor) / factor > eps)
    throw std::invalid_argument("LB tau (" + std::to_string(tau) +
                                ") must be an integer multiple of the "
                                "MD time_step (" +
                                std::to_string(time_step) + "). Factor is " +
                                std::to_string(factor));
}

double lb_lbfluid_get_tau() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return lb_walberla_params()->get_tau();
  }
#endif
  throw NoLBActive();
}

double lb_lbfluid_get_kT() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    return lb_walberla()->get_kT();
#endif
  }
  throw NoLBActive();
}

double lb_lbfluid_get_lattice_speed() {
  return lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
}

void lb_lbfluid_save_checkpoint(const std::string &filename, bool binary) {
  auto const err_msg = std::string("Error while writing LB checkpoint: ");

  // open file and set exceptions
  auto flags = std::ios_base::out;
  if (binary)
    flags |= std::ios_base::binary;
  std::fstream cpfile;
  cpfile.open(filename, flags);
  if (!cpfile) {
    throw std::runtime_error(err_msg + "could not open file " + filename);
  }
  cpfile.exceptions(std::ios_base::failbit | std::ios_base::badbit);

#ifdef LB_WALBERLA
  // write the grid size in the checkpoint header
  auto const write_header = [&](Utils::Vector3i const &grid_size,
                                std::size_t pop_size) {
    if (!binary) {
      cpfile << Utils::Vector3i::formatter(" ") << grid_size << "\n";
      cpfile << pop_size << "\n";
    } else {
      cpfile.write(reinterpret_cast<const char *>(grid_size.data()),
                   3 * sizeof(int));
      cpfile.write(reinterpret_cast<const char *>(&pop_size), sizeof(pop_size));
    }
  };
#endif // WALBERLA

  try {
    if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
      if (!binary) {
        cpfile.precision(16);
        cpfile << std::fixed;
      }

      auto const gridsize = lb_walberla()->lattice().get_grid_dimensions();
      auto const pop_size = lb_walberla()->stencil_size();
      write_header(gridsize, pop_size);

      for (int i = 0; i < gridsize[0]; i++) {
        for (int j = 0; j < gridsize[1]; j++) {
          for (int k = 0; k < gridsize[2]; k++) {
            Utils::Vector3i const ind{{i, j, k}};
            auto const pop = lb_lbnode_get_pop(ind);
            auto const laf = lb_lbnode_get_last_applied_force(ind);
            if (!binary) {
              for (auto const p : pop) {
                cpfile << p << "\n";
              }
              cpfile << Utils::Vector3d::formatter(" ") << laf << "\n";
            } else {
              cpfile.write(reinterpret_cast<const char *>(pop.data()),
                           pop_size * sizeof(double));
              cpfile.write(reinterpret_cast<const char *>(laf.data()),
                           3 * sizeof(double));
            }
          }
        }
      }
#endif // WALBERLA
    }
  } catch (std::ios_base::failure const &fail) {
    cpfile.close();
    throw std::runtime_error(err_msg + "could not write data to " + filename);
  } catch (std::runtime_error const &fail) {
    cpfile.close();
    throw;
  }
  cpfile.close();
}

void lb_lbfluid_load_checkpoint(const std::string &filename, bool binary) {
  auto const err_msg = std::string("Error while reading LB checkpoint: ");

  // open file and set exceptions
  auto flags = std::ios_base::in;
  if (binary)
    flags |= std::ios_base::binary;
  std::fstream cpfile;
  cpfile.open(filename, flags);
  if (!cpfile) {
    throw std::runtime_error(err_msg + "could not open file " + filename);
  }
  cpfile.exceptions(std::ios_base::failbit | std::ios_base::badbit);

#ifdef LB_WALBERLA
  // check the grid size in the checkpoint header matches the current grid size
  auto const check_header = [&](Utils::Vector3i const &expected_grid_size,
                                std::size_t expected_pop_size) {
    Utils::Vector3i grid_size;
    std::size_t pop_size;
    if (!binary) {
      for (auto &n : grid_size) {
        cpfile >> n;
      }
      cpfile >> pop_size;
    } else {
      cpfile.read(reinterpret_cast<char *>(grid_size.data()), 3 * sizeof(int));
      cpfile.read(reinterpret_cast<char *>(&pop_size), sizeof(pop_size));
    }
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
#endif // WALBERLA

  try {
    if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
      auto const gridsize = lb_walberla()->lattice().get_grid_dimensions();
      auto const pop_size = lb_walberla()->stencil_size();
      check_header(gridsize, pop_size);

      Utils::Vector3d laf;
      std::vector<double> pop(pop_size);
      for (int i = 0; i < gridsize[0]; i++) {
        for (int j = 0; j < gridsize[1]; j++) {
          for (int k = 0; k < gridsize[2]; k++) {
            Utils::Vector3i const ind{{i, j, k}};
            if (!binary) {
              for (auto &p : pop) {
                cpfile >> p;
              }
              for (auto &l : laf) {
                cpfile >> l;
              }
            } else {
              cpfile.read(reinterpret_cast<char *>(pop.data()),
                          pop_size * sizeof(double));
              cpfile.read(reinterpret_cast<char *>(laf.data()),
                          3 * sizeof(double));
            }
            ::Communication::mpiCallbacks().call_all(
                Walberla::set_node_from_checkpoint, ind, pop, laf);
          }
        }
      }
      ::Communication::mpiCallbacks().call_all(
          Walberla::do_ghost_communication);
#endif // WALBERLA
    } else {
      throw std::runtime_error(
          "To load an LB checkpoint one needs to have already "
          "initialized the LB fluid with the same grid size.");
    }
    // check EOF
    if (!binary) {
      if (cpfile.peek() == '\n') {
        static_cast<void>(cpfile.get());
      }
    }
    if (cpfile.peek() != EOF) {
      throw std::runtime_error(err_msg + "extra data found, expected EOF.");
    }
  } catch (std::ios_base::failure const &fail) {
    auto const eof_error = cpfile.eof();
    cpfile.close();
    if (eof_error) {
      throw std::runtime_error(err_msg + "EOF found.");
    }
    throw std::runtime_error(err_msg + "incorrectly formatted data.");
  } catch (std::runtime_error const &fail) {
    cpfile.close();
    throw;
  }
  cpfile.close();
}

Utils::Vector3i lb_lbfluid_get_shape() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return lb_walberla()->lattice().get_grid_dimensions();
  }
#endif
  throw NoLBActive();
}

bool lb_lbnode_is_index_valid(Utils::Vector3i const &ind) {
  auto const limit = lb_lbfluid_get_shape();
  return ind < limit && ind >= Utils::Vector3i{};
}

double lb_lbnode_get_density(const Utils::Vector3i &ind) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, Walberla::get_node_density, ind);
  }
#endif
  throw NoLBActive();
}

const Utils::Vector3d lb_lbnode_get_velocity(const Utils::Vector3i &ind) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, Walberla::get_node_velocity, ind);
  }
#endif
  throw NoLBActive();
}

const Utils::Vector3d
lb_lbnode_get_velocity_at_boundary(const Utils::Vector3i &ind) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank,
        Walberla::get_node_velocity_at_boundary, ind);
  }
#endif
  throw NoLBActive();
}

const Utils::Vector3d
lb_lbnode_get_last_applied_force(const Utils::Vector3i &ind) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank,
        Walberla::get_node_last_applied_force, ind);
  }
#endif

  throw NoLBActive();
}

#ifdef LB_WALBERLA
/** Revert the correction done by waLBerla on off-diagonal terms. */
inline void walberla_off_diagonal_correction(Utils::Vector6d &tensor) {
  auto const visc = lb_walberla()->get_viscosity();
  auto const revert_factor = visc / (visc + 1.0 / 6.0);
  tensor[1] *= revert_factor;
  tensor[3] *= revert_factor;
  tensor[4] *= revert_factor;
}
#endif

const Utils::Vector6d
lb_lbnode_get_pressure_tensor(const Utils::Vector3i &ind) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    Utils::Vector6d tensor = ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, Walberla::get_node_pressure_tensor,
        ind);

    walberla_off_diagonal_correction(tensor);
    return tensor;
  }
#endif
  throw NoLBActive();
}

Utils::Vector6d lb_lbfluid_get_pressure_tensor_local() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    auto const gridsize = lb_walberla()->lattice().get_grid_dimensions();
    Utils::Vector6d tensor{};
    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          const Utils::Vector3i node{{i, j, k}};
          auto const node_tensor = Walberla::get_node_pressure_tensor(node);
          if (node_tensor) {
            tensor += *node_tensor;
          }
        }
      }
    }
    return tensor;
  }
#endif
  return {};
}

REGISTER_CALLBACK_REDUCTION(lb_lbfluid_get_pressure_tensor_local,
                            std::plus<Utils::Vector6d>())

const Utils::Vector6d lb_lbfluid_get_pressure_tensor() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    auto const gridsize = lb_walberla()->lattice().get_grid_dimensions();
    auto const number_of_nodes = gridsize[0] * gridsize[1] * gridsize[2];
    auto tensor = ::Communication::mpiCallbacks().call(
        ::Communication::Result::reduction, std::plus<Utils::Vector6d>(),
        lb_lbfluid_get_pressure_tensor_local);
    tensor /= static_cast<double>(number_of_nodes);

    walberla_off_diagonal_correction(tensor);
    return tensor;
  }
#endif
  throw NoLBActive();
}

bool lb_lbnode_is_boundary(const Utils::Vector3i &ind) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, Walberla::get_node_is_boundary, ind);
  }
#endif
  throw NoLBActive();
}

const Utils::Vector3d lb_lbnode_get_boundary_force(const Utils::Vector3i &ind) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, Walberla::get_node_boundary_force,
        ind);
  }
#endif
  throw NoLBActive();
}

void lb_lbnode_remove_from_boundary(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(
        Walberla::remove_node_from_boundary, ind);
#endif
  } else {
    throw NoLBActive();
  }
}

void lb_lbfluid_clear_boundaries() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(Walberla::clear_boundaries);
#endif
  } else {
    throw NoLBActive();
  }
}

void lb_lbfluid_update_boundary_from_shape(
    std::vector<int> const &raster_flat,
    std::vector<double> const &slip_velocity_flat) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(
        Walberla::update_boundary_from_shape, raster_flat, slip_velocity_flat);
#endif
  } else {
    throw NoLBActive();
  }
}

void lb_lbfluid_update_boundary_from_list(std::vector<int> const &nodes_flat,
                                          std::vector<double> const &vel_flat) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(
        Walberla::update_boundary_from_list, nodes_flat, vel_flat);
#endif
  } else {
    throw NoLBActive();
  }
}

const std::vector<double> lb_lbnode_get_pop(const Utils::Vector3i &ind) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, Walberla::get_node_pop, ind);
  }
#endif
  throw NoLBActive();
}

void lb_lbnode_set_density(const Utils::Vector3i &ind, double p_density) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    ::Communication::mpiCallbacks().call_all(Walberla::set_node_density, ind,
                                             p_density);
  } else
#endif
  {
    throw NoLBActive();
  }
}

void lb_lbnode_set_velocity(const Utils::Vector3i &ind,
                            const Utils::Vector3d &u) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(Walberla::set_node_velocity, ind,
                                             u);
#endif
  } else {
    throw NoLBActive();
  }
}

void lb_lbnode_set_velocity_at_boundary(const Utils::Vector3i &ind,
                                        const Utils::Vector3d &u) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(
        Walberla::set_node_velocity_at_boundary, ind, u);
#endif
  } else {
    throw NoLBActive();
  }
}

void lb_lbnode_set_last_applied_force(const Utils::Vector3i &ind,
                                      const Utils::Vector3d &f) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(
        Walberla::set_node_last_applied_force, ind, f);
#endif
  } else {
    throw NoLBActive();
  }
}

void lb_lbnode_set_pop(const Utils::Vector3i &ind,
                       const std::vector<double> &p_pop) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(Walberla::set_node_pop, ind,
                                             p_pop);
#endif
  } else {
    throw NoLBActive();
  }
}

ActiveLB lb_lbfluid_get_lattice_switch() { return lattice_switch; }

Utils::Vector3d lb_lbfluid_calc_fluid_momentum() {
  Utils::Vector3d fluid_momentum{};
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    fluid_momentum = ::Communication::mpiCallbacks().call(
        ::Communication::Result::Reduction(), std::plus<>(),
        Walberla::get_momentum);
#endif
  } else
    throw NoLBActive();

  return fluid_momentum;
}

const Utils::Vector3d
lb_lbfluid_get_interpolated_velocity(const Utils::Vector3d &pos) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    auto const folded_pos = folded_position(pos, box_geo);
    return mpi_call(::Communication::Result::one_rank,
                    Walberla::get_velocity_at_pos,
                    folded_pos / lb_lbfluid_get_agrid());
#endif
  }
  throw NoLBActive();
}

const Utils::Vector3d
lb_lbfluid_get_force_to_be_applied(const Utils::Vector3d &pos) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    auto const agrid = lb_lbfluid_get_agrid();
    auto const ind = Utils::Vector3i{static_cast<int>(pos[0] / agrid),
                                     static_cast<int>(pos[1] / agrid),
                                     static_cast<int>(pos[2] / agrid)};
    auto const res = lb_walberla()->get_node_force_to_be_applied(ind);
    if (!res) {
      std::cout << this_node << ": position: [" << pos << "]\n";
      throw std::runtime_error(
          "Force to be applied could not be obtained from Walberla");
    }
    return *res;
#endif
  }
  throw NoLBActive();
}

void lb_lbfluid_add_force_at_pos(const Utils::Vector3d &pos,
                                 const Utils::Vector3d &f) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    auto const folded_pos = folded_position(pos, box_geo);
    ::Communication::mpiCallbacks().call_all(
        Walberla::add_force_at_pos, folded_pos / lb_lbfluid_get_agrid(), f);
#endif
  } else {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_interpolated_density(const Utils::Vector3d &pos) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    auto const folded_pos = folded_position(pos, box_geo);
    return mpi_call(::Communication::Result::one_rank,
                    Walberla::get_interpolated_density_at_pos,
                    folded_pos / lb_lbfluid_get_agrid());
#endif
  }
  throw NoLBActive();
}
