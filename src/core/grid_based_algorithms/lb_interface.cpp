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
#include "grid_based_algorithms/lb_boundaries.hpp"
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

#include <algorithm>
#include <cmath>
#include <cstdio>
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

void lb_lbfluid_integrate() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    lb_walberla()->integrate();
#endif
  } else
    throw NoLBActive();
}

void lb_lbfluid_propagate() {
  if (lattice_switch != ActiveLB::NONE) {
    lb_lbfluid_integrate();
  }
}

/**
 * @brief Check the boundary velocities.
 * Sanity check if the velocity defined at LB boundaries is within the Mach
 * number limits of the scheme i.e. u < 0.3.
 */
void lb_boundary_mach_check() {
  // Boundary velocities are stored in MD units, therefore we need to scale them
  // in order to get lattice units.
  auto const conv_fac = lb_lbfluid_get_tau() / lb_lbfluid_get_agrid();
  double constexpr mach_limit = 0.3;
  using LBBoundaries::lbboundaries;
  if (std::any_of(lbboundaries.begin(), lbboundaries.end(),
                  [conv_fac, mach_limit](auto const &b) {
                    return (b->velocity() * conv_fac).norm() >= mach_limit;
                  })) {
    runtimeErrorMsg() << "Lattice velocity exceeds the Mach number limit";
  }
}

void lb_lbfluid_sanity_checks(double time_step) {
  if (lattice_switch == ActiveLB::NONE)
    return;

  // LB GPU interface functions only work on the head node.
  lb_boundary_mach_check();
  if (time_step > 0.)
    check_tau_time_step_consistency(lb_lbfluid_get_tau(), time_step);

  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    // Make sure, Walberla and Espresso agree on domain decomposition
    auto walberla_domain = lb_walberla()->get_local_domain();
    // Unit conversion
    auto const agrid = lb_lbfluid_get_agrid();
    walberla_domain.first *= agrid;
    walberla_domain.second *= agrid;

    auto const my_left = local_geo.my_left();
    auto const my_right = local_geo.my_right();
    auto const tol = lb_lbfluid_get_agrid() / 1E6;
    if ((walberla_domain.first - my_left).norm2() > tol or
        (walberla_domain.second - my_right).norm2() > tol) {
      std::cout << this_node << ": left ESPResSo: [" << my_left << "], "
                << "left waLBerla: [" << walberla_domain.first << "]\n";
      std::cout << this_node << ": right ESPResSo: [" << my_right << "], "
                << "right waLBerla: [" << walberla_domain.second << "]\n";
      throw std::runtime_error(
          "waLBerla and ESPResSo disagree about domain decomposition.");
    }
#endif
  }
}

uint64_t lb_lbfluid_get_rng_state() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return Walberla::get_rng_state();
  }
#endif
  throw NoLBActive();
}

void lb_lbfluid_set_rng_state(uint64_t counter) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    Walberla::set_rng_state(counter);
  } else
#endif
  {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_viscosity() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return lb_walberla()->get_viscosity();
  }
#endif
  throw NoLBActive();
}

double lb_lbfluid_get_bulk_viscosity() {
  if (lattice_switch == ActiveLB::WALBERLA)
    throw std::runtime_error("Getting bulk viscosity not implemented.");
  throw NoLBActive();
}

double lb_lbfluid_get_agrid() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    return lb_walberla_params()->get_agrid();
#endif
  }
  throw NoLBActive();
}

void lb_lbfluid_set_ext_force_density(const Utils::Vector3d &force_density) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    ::Communication::mpiCallbacks().call_all(Walberla::set_ext_force_density,
                                             force_density);
#endif
  } else {
    throw NoLBActive();
  }
}

const Utils::Vector3d lb_lbfluid_get_ext_force_density() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    return lb_walberla()->get_external_force();
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
                                ") must be integer multiple of "
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

void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch) {
  switch (local_lattice_switch) {
  case ActiveLB::NONE:
  case ActiveLB::WALBERLA:
    break;
  default:
    throw std::invalid_argument("Invalid lattice switch.");
  }
  mpi_set_lattice_switch(local_lattice_switch);
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
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    std::fstream cpfile;
    if (binary) {
      cpfile.open(filename, std::ios::out | std::ios::binary);
    } else {
      cpfile.open(filename, std::ios::out);
      cpfile.precision(16);
      cpfile << std::fixed;
    }

    auto const gridsize = lb_walberla()->get_grid_dimensions();
    auto const pop_size = lb_walberla()->stencil_size();
    std::vector<double> pop(pop_size);

    if (!binary) {
      cpfile << gridsize[0] << " " << gridsize[1] << " " << gridsize[2] << "\n";
      cpfile << pop_size << "\n";
    } else {
      cpfile.write(reinterpret_cast<const char *>(gridsize.data()),
                   3 * sizeof(gridsize[0]));
      cpfile.write(reinterpret_cast<const char *>(&pop_size), sizeof(pop_size));
    }

    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          Utils::Vector3i const ind{{i, j, k}};
          auto const pop = lb_lbnode_get_pop(ind);
          auto const laf = lb_lbnode_get_last_applied_force(ind);
          if (!binary) {
            for (std::size_t n = 0; n < pop_size; n++) {
              cpfile << pop[n] << " ";
            }
            cpfile << "\n";
            for (std::size_t n = 0; n < 3; n++) {
              cpfile << laf[n] << " ";
            }
            cpfile << "\n";
          } else {
            cpfile.write(reinterpret_cast<const char *>(pop.data()),
                         pop_size * sizeof(double));
            cpfile.write(reinterpret_cast<const char *>(laf.data()),
                         3 * sizeof(double));
          }
        }
      }
    }
    cpfile.close();
  }
#endif
}

void lb_lbfluid_load_checkpoint(const std::string &filename, bool binary) {
  int res;
  std::string err_msg = "Error while reading LB checkpoint: ";
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    FILE *cpfile;
    cpfile = fopen(filename.c_str(), "r");
    if (!cpfile) {
      throw std::runtime_error(err_msg + "could not open file for reading.");
    }

    auto const pop_size = lb_walberla()->stencil_size();
    std::size_t saved_pop_size;
    Utils::Vector3d laf;
    auto const gridsize = lb_walberla()->get_grid_dimensions();
    Utils::Vector3i saved_gridsize;
    if (!binary) {
      res = fscanf(cpfile, "%i %i %i\n%zu\n", &saved_gridsize[0],
                   &saved_gridsize[1], &saved_gridsize[2], &saved_pop_size);
      if (res == EOF) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "EOF found.");
      }
      if (res != 4) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
    } else {
      if (fread(&saved_gridsize[0], sizeof(int), 3, cpfile) != 3) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
      if (fread(&saved_pop_size, sizeof(std::size_t), 1, cpfile) != 1) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
    }
    if (saved_gridsize != gridsize) {
      fclose(cpfile);
      std::stringstream message;
      message << " grid dimensions mismatch,"
              << " read [" << saved_gridsize << "],"
              << " expected [" << gridsize << "].";
      throw std::runtime_error(err_msg + message.str());
    }
    if (saved_pop_size != pop_size) {
      fclose(cpfile);
      throw std::runtime_error(err_msg + "population size mismatch, read " +
                               std::to_string(saved_pop_size) + ", expected " +
                               std::to_string(pop_size) + ".");
    }

    std::vector<double> pop(saved_pop_size);
    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          Utils::Vector3i const ind{{i, j, k}};
          if (!binary) {
            for (std::size_t f = 0; f < saved_pop_size; ++f) {
              res = fscanf(cpfile, "%lf ", &pop[f]);
              if (res == EOF) {
                fclose(cpfile);
                throw std::runtime_error(err_msg + "EOF found.");
              }
              if (res != 1) {
                fclose(cpfile);
                throw std::runtime_error(err_msg +
                                         "incorrectly formatted data.");
              }
            }
            res = fscanf(cpfile, "\n");
            if (res == EOF) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "EOF found.");
            }
            res = fscanf(cpfile, "%lf %lf %lf \n", &laf[0], &laf[1], &laf[2]);
            if (res == EOF) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "EOF found.");
            }
            if (res != 3) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "incorrectly formatted data.");
            }
          } else {
            if (fread(pop.data(), sizeof(double), saved_pop_size, cpfile) !=
                saved_pop_size) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "incorrectly formatted data.");
            }
            if (fread(laf.data(), sizeof(double), 3, cpfile) != 3) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "incorrectly formatted data.");
            }
          }
          ::Communication::mpiCallbacks().call_all(
              Walberla::set_node_from_checkpoint, ind, pop, laf);
        }
      }
    }
    ::Communication::mpiCallbacks().call_all(Walberla::do_ghost_communication);
    if (!binary) {
      // skip spaces
      for (int n = 0; n < 2; ++n) {
        res = fgetc(cpfile);
        if (res != (int)' ' && res != (int)'\n')
          break;
      }
    } else {
      res = fgetc(cpfile);
    }
    if (res != EOF) {
      fclose(cpfile);
      throw std::runtime_error(err_msg + "extra data found, expected EOF.");
    }
    fclose(cpfile);
  } else
#endif
  {
    throw std::runtime_error(
        "To load an LB checkpoint one needs to have already "
        "initialized the LB fluid with the same grid size.");
  }
}

Utils::Vector3i lb_lbfluid_get_shape() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return lb_walberla()->get_grid_dimensions();
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

void lb_lbfluid_create_vtk(unsigned delta_N, unsigned initial_count,
                           unsigned flag_observables,
                           std::string const &identifier,
                           std::string const &base_folder,
                           std::string const &prefix) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    ::Communication::mpiCallbacks().call_all(Walberla::create_vtk, delta_N,
                                             initial_count, flag_observables,
                                             identifier, base_folder, prefix);
    return;
  }
#endif
  throw NoLBActive();
}

void lb_lbfluid_write_vtk(std::string const &vtk_uid) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    ::Communication::mpiCallbacks().call_all(Walberla::write_vtk, vtk_uid);
    return;
  }
#endif
  throw NoLBActive();
}

void lb_lbfluid_switch_vtk(std::string const &vtk_uid, int status) {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    ::Communication::mpiCallbacks().call_all(Walberla::switch_vtk, vtk_uid,
                                             status);
    return;
  }
#endif
  throw NoLBActive();
}

/** Revert the correction done by waLBerla on off-diagonal terms. */
inline void walberla_off_diagonal_correction(Utils::Vector6d &tensor) {
  auto const visc = lb_lbfluid_get_viscosity();
  auto const revert_factor = visc / (visc + 1.0 / 6.0);
  tensor[1] *= revert_factor;
  tensor[3] *= revert_factor;
  tensor[4] *= revert_factor;
}

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
    auto const gridsize = lb_walberla()->get_grid_dimensions();
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
    auto const gridsize = lb_walberla()->get_grid_dimensions();
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

void mpi_set_lattice_switch_local(ActiveLB lattice_switch) {
  ::lattice_switch = lattice_switch;
}

REGISTER_CALLBACK(mpi_set_lattice_switch_local)

void mpi_set_lattice_switch(ActiveLB lattice_switch) {
  mpi_call_all(mpi_set_lattice_switch_local, lattice_switch);
}
