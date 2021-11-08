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
#include "lb_interface.hpp"
#include "BoxGeometry.hpp"
#include "MpiCallbacks.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "electrokinetics.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "halo.hpp"
#include "lb-d3q19.hpp"
#include "lb.hpp"
#include "lb_boundaries.hpp"
#include "lb_collective_interface.hpp"
#include "lb_constants.hpp"
#include "lb_interpolation.hpp"
#include "lbgpu.hpp"

#include <utils/Vector.hpp>

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

ActiveLB lattice_switch = ActiveLB::NONE;

struct NoLBActive : public std::exception {
  const char *what() const noexcept override { return "LB not activated"; }
};

void lb_lbfluid_integrate() {
  if (lattice_switch == ActiveLB::CPU) {
    lb_integrate();
  } else if (lattice_switch == ActiveLB::GPU and this_node == 0) {
#ifdef CUDA
#ifdef ELECTROKINETICS
    if (ek_initialized) {
      ek_integrate();
    } else {
#endif
      lb_integrate_GPU();
#ifdef ELECTROKINETICS
    }
#endif
#endif
  }
}

void lb_lbfluid_propagate() {
  if (lattice_switch != ActiveLB::NONE) {
    lb_lbfluid_integrate();
    if (lb_lbfluid_get_kT() > 0.0) {
      if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
        rng_counter_fluid_gpu->increment();
#endif
      } else if (lattice_switch == ActiveLB::CPU) {
        rng_counter_fluid->increment();
      }
    }
  }
}

/**
 * @brief Check the boundary velocities.
 */
inline void lb_boundary_mach_check() {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
  if (LBBoundaries::sanity_check_mach_limit()) {
    runtimeErrorMsg() << "Lattice velocity exceeds the Mach number limit";
  }
#endif
}

void lb_lbfluid_sanity_checks(double time_step) {
  if (lattice_switch == ActiveLB::GPU && this_node == 0) {
#ifdef CUDA
    lb_GPU_sanity_checks();
    lb_boundary_mach_check();
    if (time_step > 0.)
      check_tau_time_step_consistency(lb_lbfluid_get_tau(), time_step);
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    lb_sanity_checks(lbpar);
    lb_boundary_mach_check();
    if (time_step > 0.)
      check_tau_time_step_consistency(lb_lbfluid_get_tau(), time_step);
  }
}

void lb_lbfluid_on_integration_start() {
  if (lattice_switch == ActiveLB::CPU) {
    halo_communication(update_halo_comm,
                       reinterpret_cast<char *>(lbfluid[0].data()));
  }
}

/** (Re-)initialize the fluid. */
void lb_lbfluid_reinit_parameters() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    if (this_node == 0)
      lb_reinit_parameters_gpu();
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lb_reinit_parameters(lbpar);
  }
}

/** Perform a full initialization of the lattice Boltzmann system.
 *  All derived parameters and the fluid are reset to their default values.
 */
void lb_lbfluid_init() {
  if (lattice_switch == ActiveLB::GPU && this_node == 0) {
#ifdef CUDA
    lb_init_gpu();
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lb_init(lbpar);
  }
}

uint64_t lb_lbfluid_get_rng_state() {
  if (lattice_switch == ActiveLB::CPU) {
    return lb_fluid_get_rng_state();
  }
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return lb_fluid_get_rng_state_gpu();
#endif
  }
  throw NoLBActive();
}

void lb_lbfluid_set_rng_state(uint64_t counter) {
  if (lattice_switch == ActiveLB::CPU) {
    lb_fluid_set_rng_state(counter);
  } else if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lb_fluid_set_rng_state_gpu(counter);
#endif
  } else {
    throw NoLBActive();
  }
}

void lb_lbfluid_set_density(double density) {
  if (density <= 0)
    throw std::invalid_argument("Density has to be > 0. but got " +
                                std::to_string(density));
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lbpar_gpu.rho = static_cast<float>(density);
    lb_reinit_fluid_gpu();
    lb_lbfluid_reinit_parameters();
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.density = density;
    mpi_bcast_lb_params(LBParam::DENSITY);
  } else {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_density() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return static_cast<double>(lbpar_gpu.rho);
#else
    return {};
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.density;
  }
  throw NoLBActive();
}

void lb_lbfluid_set_viscosity(double viscosity) {
  if (viscosity <= 0)
    throw std::invalid_argument("Viscosity has to be >0.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lbpar_gpu.viscosity = static_cast<float>(viscosity);
    lb_lbfluid_reinit_parameters();
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.viscosity = viscosity;
    mpi_bcast_lb_params(LBParam::VISCOSITY);
  } else {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_viscosity() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return static_cast<double>(lbpar_gpu.viscosity);
#else
    return {};
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.viscosity;
  }
  throw NoLBActive();
}

void lb_lbfluid_set_bulk_viscosity(double bulk_viscosity) {
  if (bulk_viscosity <= 0)
    throw std::invalid_argument("Bulk viscosity has to be >0. but got " +
                                std::to_string(bulk_viscosity));
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lbpar_gpu.bulk_viscosity = static_cast<float>(bulk_viscosity);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_reinit_parameters();
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.bulk_viscosity = bulk_viscosity;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(LBParam::BULKVISC);
  } else {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_bulk_viscosity() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return lbpar_gpu.bulk_viscosity;
#else
    return {};
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.bulk_viscosity;
  }
  throw NoLBActive();
}

void lb_lbfluid_set_gamma_odd(double gamma_odd) {
  if (fabs(gamma_odd) > 1)
    throw std::invalid_argument("Gamma odd has to be <= 1.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lbpar_gpu.gamma_odd = static_cast<float>(gamma_odd);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_reinit_parameters();
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.gamma_odd = gamma_odd;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(LBParam::GAMMA_ODD);
  } else {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_gamma_odd() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return lbpar_gpu.gamma_odd;
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.gamma_odd;
  }
  throw NoLBActive();
}

void lb_lbfluid_set_gamma_even(double gamma_even) {
  if (fabs(gamma_even) > 1)
    throw std::invalid_argument("gamma_even has to be <= 1.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lbpar_gpu.gamma_even = static_cast<float>(gamma_even);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_reinit_parameters();
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.gamma_even = gamma_even;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(LBParam::DENSITY);
  } else {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_gamma_even() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return lbpar_gpu.gamma_even;
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.gamma_even;
  }
  throw NoLBActive();
}

void lb_lbfluid_set_agrid(double agrid) {
  if (agrid <= 0)
    throw std::invalid_argument("agrid has to be > 0.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lb_set_agrid_gpu(agrid);
    lb_init_gpu();
#if defined(LB_BOUNDARIES_GPU)
    LBBoundaries::lb_init_boundaries();
#endif
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.agrid = agrid;
    mpi_bcast_lb_params(LBParam::AGRID);
  } else {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_agrid() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return lbpar_gpu.agrid;
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.agrid;
  }
  throw NoLBActive();
}

void lb_lbfluid_set_ext_force_density(const Utils::Vector3d &force_density) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lbpar_gpu.ext_force_density[0] = static_cast<float>(force_density[0]);
    lbpar_gpu.ext_force_density[1] = static_cast<float>(force_density[1]);
    lbpar_gpu.ext_force_density[2] = static_cast<float>(force_density[2]);
    if (force_density[0] != 0 || force_density[1] != 0 ||
        force_density[2] != 0) {
      lbpar_gpu.external_force_density = 1;
    } else {
      lbpar_gpu.external_force_density = 0;
    }
    lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);

#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.ext_force_density = force_density;
    mpi_bcast_lb_params(LBParam::EXT_FORCE_DENSITY);
  } else {
    throw NoLBActive();
  }
}

const Utils::Vector3d lb_lbfluid_get_ext_force_density() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return {static_cast<double>(lbpar_gpu.ext_force_density[0]),
            static_cast<double>(lbpar_gpu.ext_force_density[1]),
            static_cast<double>(lbpar_gpu.ext_force_density[2])};
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.ext_force_density;
  }
  throw NoLBActive();
}

void lb_lbfluid_set_tau(double tau) {
  if (tau <= 0.)
    throw std::invalid_argument("LB tau has to be positive.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lbpar_gpu.tau = static_cast<float>(tau);
    lb_lbfluid_reinit_parameters();
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.tau = tau;
    mpi_bcast_lb_params(LBParam::TAU);
  } else {
    throw NoLBActive();
  }
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
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return lbpar_gpu.tau;
#else
    return {};
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.tau;
  }
  throw NoLBActive();
}

void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch) {
  switch (local_lattice_switch) {
  case ActiveLB::NONE:
  case ActiveLB::CPU:
  case ActiveLB::GPU:
    break;
  default:
    throw std::invalid_argument("Invalid lattice switch.");
  }
  mpi_set_lattice_switch(local_lattice_switch);
}

void lb_lbfluid_set_kT(double kT) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lbpar_gpu.kT = static_cast<float>(kT);
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.kT = kT;
    mpi_bcast_lb_params(LBParam::KT);
  } else {
    throw NoLBActive();
  }
}

double lb_lbfluid_get_kT() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return static_cast<double>(lbpar_gpu.kT);
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.kT;
  }
  throw NoLBActive();
}

double lb_lbfluid_get_lattice_speed() {
  return lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
}

void lb_lbfluid_print_vtk_boundary(const std::string &filename) {
  std::fstream cpfile;
  cpfile.open(filename, std::ios::out);

  if (!cpfile) {
    throw std::runtime_error("Could not open '" + filename + "' for writing.");
  }

  auto const vtk_writer = [&](std::string const &label,
                              auto const &write_boundaries) {
    using Utils::Vector3d;
    cpfile.precision(6);
    cpfile << std::fixed;
    auto constexpr vtk_format = Vector3d::formatter(" ");
    auto const agrid = lb_lbfluid_get_agrid();
    auto const grid_size = lb_lbfluid_get_shape();
    auto const origin = Vector3d::broadcast(0.5) * agrid;
    cpfile << "# vtk DataFile Version 2.0\n"
           << label << "\n"
           << "ASCII\n"
           << "DATASET STRUCTURED_POINTS\n"
           << "DIMENSIONS " << vtk_format << grid_size << "\n"
           << "ORIGIN " << vtk_format << origin << "\n"
           << "SPACING " << vtk_format << Vector3d::broadcast(agrid) << "\n"
           << "POINT_DATA " << Utils::product(grid_size) << "\n"
           << "SCALARS boundary float 1\n"
           << "LOOKUP_TABLE default\n";
    write_boundaries();
  };

  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    std::vector<unsigned int> bound_array(lbpar_gpu.number_of_nodes);
    lb_get_boundary_flags_GPU(bound_array.data());
    vtk_writer("lbboundaries", [&]() {
      for (unsigned int j = 0; j < lbpar_gpu.number_of_nodes; ++j) {
        cpfile << bound_array[j] << "\n";
      }
    });
#endif //  CUDA
  } else {
    vtk_writer("lbboundaries", [&]() {
      auto const grid_size = lb_lbfluid_get_shape();
      Utils::Vector3i pos;
      for (pos[2] = 0; pos[2] < grid_size[2]; pos[2]++)
        for (pos[1] = 0; pos[1] < grid_size[1]; pos[1]++)
          for (pos[0] = 0; pos[0] < grid_size[0]; pos[0]++)
            cpfile << lb_lbnode_get_boundary(pos) << "\n";
    });
  }
  cpfile.close();
}

void lb_lbfluid_print_vtk_velocity(const std::string &filename,
                                   std::vector<int> bb1, std::vector<int> bb2) {
  std::fstream cpfile;
  cpfile.open(filename, std::ios::out);

  if (!cpfile) {
    throw std::runtime_error("Could not open '" + filename + "' for writing.");
  }

  auto bb_low = Utils::Vector3i{};
  auto bb_high = lb_lbfluid_get_shape();

  auto const vtk_writer = [&](std::string const &label, auto const &get_vel) {
    using Utils::Vector3d;
    cpfile.precision(6);
    cpfile << std::fixed;
    auto constexpr vtk_format = Vector3d::formatter(" ");
    auto const agrid = lb_lbfluid_get_agrid();
    auto const bb_dim = bb_high - bb_low;
    auto const origin = (bb_low + Vector3d::broadcast(0.5)) * agrid;
    auto const lattice_speed = lb_lbfluid_get_lattice_speed();
    cpfile << "# vtk DataFile Version 2.0\n"
           << label << "\n"
           << "ASCII\n"
           << "DATASET STRUCTURED_POINTS\n"
           << "DIMENSIONS " << vtk_format << bb_dim << "\n"
           << "ORIGIN " << vtk_format << origin << "\n"
           << "SPACING " << vtk_format << Vector3d::broadcast(agrid) << "\n"
           << "POINT_DATA " << Utils::product(bb_dim) << "\n"
           << "SCALARS velocity float 3\n"
           << "LOOKUP_TABLE default\n";

    Utils::Vector3i pos;
    for (pos[2] = bb_low[2]; pos[2] < bb_high[2]; pos[2]++)
      for (pos[1] = bb_low[1]; pos[1] < bb_high[1]; pos[1]++)
        for (pos[0] = bb_low[0]; pos[0] < bb_high[0]; pos[0]++)
          cpfile << vtk_format << get_vel(pos) * lattice_speed << "\n";
  };

  int it = 0;
  for (auto val1 = bb1.begin(), val2 = bb2.begin();
       val1 != bb1.end() && val2 != bb2.end(); ++val1, ++val2) {
    if (*val1 == -1 || *val2 == -1) {
      break;
    }
    auto const lower = std::min(*val1, *val2);
    auto const upper = std::max(*val1, *val2);
    if (lower < 0 or upper >= bb_high[it]) {
      throw std::runtime_error(
          "Tried to access index " + std::to_string(lower) + " and index " +
          std::to_string(upper) + " on dimension " + std::to_string(it) +
          " that has size " + std::to_string(bb_high[it]));
    }
    bb_low[it] = lower;
    bb_high[it] = upper;
    it++;
  }

  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    host_values.resize(lbpar_gpu.number_of_nodes);
    lb_get_values_GPU(host_values.data());
    auto const box_l_x = lb_lbfluid_get_shape()[0];
    vtk_writer("lbfluid_gpu", [box_l_x](Utils::Vector3i const &pos) {
      auto const j = box_l_x * box_l_x * pos[2] + box_l_x * pos[1] + pos[0];
      return Utils::Vector3d{host_values[j].v};
    });
#endif //  CUDA
  } else {
    vtk_writer("lbfluid_cpu", lb_lbnode_get_velocity);
  }
  cpfile.close();
}

void lb_lbfluid_print_boundary(const std::string &filename) {
  std::fstream cpfile;
  cpfile.open(filename, std::ios::out);

  if (!cpfile) {
    throw std::runtime_error("Could not open '" + filename + "' for writing.");
  }

  using Utils::Vector3d;
  auto constexpr vtk_format = Vector3d::formatter(" ");
  cpfile.precision(6);
  cpfile << std::fixed;

  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    std::vector<unsigned int> bound_array(lbpar_gpu.number_of_nodes);
    lb_get_boundary_flags_GPU(bound_array.data());
    auto const agrid = lb_lbfluid_get_agrid();
    Utils::Vector3d pos;
    for (unsigned int j = 0; j < lbpar_gpu.number_of_nodes; ++j) {
      auto const k = j / lbpar_gpu.dim[0];
      auto const l = k / lbpar_gpu.dim[1];
      pos[0] = (static_cast<double>(j % lbpar_gpu.dim[0]) + 0.5) * agrid;
      pos[1] = (static_cast<double>(k % lbpar_gpu.dim[1]) + 0.5) * agrid;
      pos[2] = (static_cast<double>(l) + 0.5) * agrid;
      cpfile << vtk_format << pos << " " << bound_array[j] << "\n";
    }
#endif //  CUDA
  } else {
    auto constexpr shift = Vector3d{0.5, 0.5, 0.5};
    auto const agrid = lb_lbfluid_get_agrid();
    auto const grid_size = lb_lbfluid_get_shape();
    Utils::Vector3i pos;
    for (pos[2] = 0; pos[2] < grid_size[2]; pos[2]++)
      for (pos[1] = 0; pos[1] < grid_size[1]; pos[1]++)
        for (pos[0] = 0; pos[0] < grid_size[0]; pos[0]++) {
          auto const flag = (lb_lbnode_get_boundary(pos) != 0) ? 1 : 0;
          cpfile << vtk_format << (pos + shift) * agrid << " " << flag << "\n";
        }
  }
  cpfile.close();
}

void lb_lbfluid_print_velocity(const std::string &filename) {
  std::fstream cpfile;
  cpfile.open(filename, std::ios::out);

  if (!cpfile) {
    throw std::runtime_error("Could not open '" + filename + "' for writing.");
  }

  using Utils::Vector3d;
  auto constexpr vtk_format = Vector3d::formatter(" ");
  cpfile.precision(6);
  cpfile << std::fixed;

  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    std::vector<LB_rho_v_pi_gpu> host_values(lbpar_gpu.number_of_nodes);
    lb_get_values_GPU(host_values.data());
    auto const agrid = lb_lbfluid_get_agrid();
    auto const lattice_speed = lb_lbfluid_get_lattice_speed();
    Utils::Vector3d pos;
    for (unsigned int j = 0; j < lbpar_gpu.number_of_nodes; ++j) {
      auto const k = j / lbpar_gpu.dim[0];
      auto const l = k / lbpar_gpu.dim[1];
      pos[0] = (static_cast<double>(j % lbpar_gpu.dim[0]) + 0.5) * agrid;
      pos[1] = (static_cast<double>(k % lbpar_gpu.dim[1]) + 0.5) * agrid;
      pos[2] = (static_cast<double>(l) + 0.5) * agrid;
      auto const velocity = Utils::Vector3f(host_values[j].v) * lattice_speed;
      cpfile << vtk_format << pos << " " << vtk_format << velocity << "\n";
    }
#endif //  CUDA
  } else {
    auto constexpr shift = Vector3d{0.5, 0.5, 0.5};
    auto const agrid = lb_lbfluid_get_agrid();
    auto const grid_size = lb_lbfluid_get_shape();
    auto const lattice_speed = lb_lbfluid_get_lattice_speed();
    Utils::Vector3i pos;
    for (pos[2] = 0; pos[2] < grid_size[2]; pos[2]++)
      for (pos[1] = 0; pos[1] < grid_size[1]; pos[1]++)
        for (pos[0] = 0; pos[0] < grid_size[0]; pos[0]++)
          cpfile << vtk_format << (pos + shift) * agrid << " " << vtk_format
                 << lb_lbnode_get_velocity(pos) * lattice_speed << "\n";
  }

  cpfile.close();
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

  // write the grid size in the checkpoint header
  auto const write_header = [&](Utils::Vector3i const &grid_size) {
    if (!binary) {
      cpfile << Utils::Vector3i::formatter(" ") << grid_size << "\n";
    } else {
      cpfile.write(reinterpret_cast<const char *>(grid_size.data()),
                   3 * sizeof(int));
    }
  };

  try {
    if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
      if (!binary) {
        cpfile.precision(8);
        cpfile << std::fixed;
      }

      auto const gridsize = lb_lbfluid_get_shape();
      auto const data_length = lbpar_gpu.number_of_nodes * D3Q19::n_vel;
      write_header(gridsize);

      std::vector<float> host_checkpoint_vd(data_length);
      lb_save_checkpoint_GPU(host_checkpoint_vd.data());
      if (!binary) {
        for (auto const p : host_checkpoint_vd) {
          cpfile << p << "\n";
        }
      } else {
        cpfile.write(reinterpret_cast<char *>(host_checkpoint_vd.data()),
                     data_length * sizeof(float));
      }
#endif //  CUDA
    } else if (lattice_switch == ActiveLB::CPU) {
      if (!binary) {
        cpfile.precision(16);
        cpfile << std::fixed;
      }

      auto const gridsize = lb_lbfluid_get_shape();
      write_header(gridsize);

      for (int i = 0; i < gridsize[0]; i++) {
        for (int j = 0; j < gridsize[1]; j++) {
          for (int k = 0; k < gridsize[2]; k++) {
            Utils::Vector3i const ind{{i, j, k}};
            auto const pop = mpi_call(::Communication::Result::one_rank,
                                      mpi_lb_get_populations, ind);
            if (!binary) {
              for (auto const p : pop) {
                cpfile << p << "\n";
              }
            } else {
              cpfile.write(reinterpret_cast<const char *>(pop.data()),
                           D3Q19::n_vel * sizeof(double));
            }
          }
        }
      }
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

  // check the grid size in the checkpoint header matches the current grid size
  auto const check_header = [&](Utils::Vector3i const &expected_grid_size) {
    Utils::Vector3i grid_size;
    if (!binary) {
      for (auto &n : grid_size) {
        cpfile >> n;
      }
    } else {
      cpfile.read(reinterpret_cast<char *>(grid_size.data()), 3 * sizeof(int));
    }
    if (grid_size != expected_grid_size) {
      std::stringstream message;
      message << " grid dimensions mismatch,"
              << " read [" << grid_size << "],"
              << " expected [" << expected_grid_size << "].";
      throw std::runtime_error(err_msg + message.str());
    }
  };

  try {
    if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
      auto const gridsize = lb_lbfluid_get_shape();
      auto const data_length = lbpar_gpu.number_of_nodes * D3Q19::n_vel;
      std::vector<float> host_checkpoint_vd(data_length);
      check_header(gridsize);

      if (!binary) {
        for (auto &p : host_checkpoint_vd) {
          cpfile >> p;
        }
      } else {
        cpfile.read(reinterpret_cast<char *>(host_checkpoint_vd.data()),
                    data_length * sizeof(float));
      }
      lb_load_checkpoint_GPU(host_checkpoint_vd.data());
#endif //  CUDA
    } else if (lattice_switch == ActiveLB::CPU) {
      auto const gridsize = lb_lbfluid_get_shape();
      mpi_bcast_lb_params(LBParam::DENSITY);
      check_header(gridsize);

      Utils::Vector19d pop;
      for (int i = 0; i < gridsize[0]; i++) {
        for (int j = 0; j < gridsize[1]; j++) {
          for (int k = 0; k < gridsize[2]; k++) {
            Utils::Vector3i const ind{{i, j, k}};
            if (!binary) {
              for (auto &p : pop) {
                cpfile >> p;
              }
            } else {
              cpfile.read(reinterpret_cast<char *>(pop.data()),
                          D3Q19::n_vel * sizeof(double));
            }
            lb_lbnode_set_pop(ind, pop);
          }
        }
      }
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
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return {static_cast<int>(lbpar_gpu.dim[0]),
            static_cast<int>(lbpar_gpu.dim[1]),
            static_cast<int>(lbpar_gpu.dim[2])};
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lblattice.global_grid;
  }
  throw NoLBActive();
}

bool lb_lbnode_is_index_valid(Utils::Vector3i const &ind) {
  auto const limit = lb_lbfluid_get_shape();
  return ind < limit && ind >= Utils::Vector3i{};
}

double lb_lbnode_get_density(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    auto const single_nodeindex = calculate_node_index(lbpar_gpu, ind);
    static LB_rho_v_pi_gpu host_print_values;
    lb_print_node_GPU(single_nodeindex, &host_print_values);
    return host_print_values.rho;
#else
    return {};
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, mpi_lb_get_density, ind);
  }
  throw NoLBActive();
}

const Utils::Vector3d lb_lbnode_get_velocity(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    static LB_rho_v_pi_gpu host_print_values;
    auto const single_nodeindex = calculate_node_index(lbpar_gpu, ind);
    lb_print_node_GPU(single_nodeindex, &host_print_values);
    return {static_cast<double>(host_print_values.v[0]),
            static_cast<double>(host_print_values.v[1]),
            static_cast<double>(host_print_values.v[2])};
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    auto const density = ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, mpi_lb_get_density, ind);
    auto const momentum_density = ::Communication::mpiCallbacks().call(
        ::Communication::Result::one_rank, mpi_lb_get_momentum_density, ind);
    return momentum_density / density;
  }
  throw NoLBActive();
}

const Utils::Vector6d
lb_lbnode_get_pressure_tensor(const Utils::Vector3i &ind) {
  // Add equilibrium pressure to the diagonal (in LB units)
  auto const p0 = lb_lbfluid_get_density() * D3Q19::c_sound_sq<double>;

  auto tensor = lb_lbnode_get_pressure_tensor_neq(ind);
  tensor[0] += p0;
  tensor[2] += p0;
  tensor[5] += p0;

  return tensor;
}

const Utils::Vector6d
lb_lbnode_get_pressure_tensor_neq(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    Utils::Vector6d tensor{};
    static LB_rho_v_pi_gpu host_print_values;
    auto const single_nodeindex = calculate_node_index(lbpar_gpu, ind);
    lb_print_node_GPU(single_nodeindex, &host_print_values);
    for (int i = 0; i < 6; i++) {
      tensor[i] = static_cast<double>(host_print_values.pi[i]);
    }
    return tensor;
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return mpi_call(::Communication::Result::one_rank,
                    mpi_lb_get_pressure_tensor, ind);
  }
  throw NoLBActive();
}

const Utils::Vector6d lb_lbfluid_get_pressure_tensor() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    auto const stress_tmp = stress_tensor_GPU();
    Utils::Vector6d tensor(stress_tmp.begin(), stress_tmp.end());

    // Normalize
    tensor /= static_cast<double>(lbpar_gpu.number_of_nodes);

    // Add equilibrium pressure to the diagonal (in LB units)
    double const p0 = lb_lbfluid_get_density() * D3Q19::c_sound_sq<double>;

    tensor[0] += p0;
    tensor[2] += p0;
    tensor[5] += p0;
    return tensor;
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    auto const grid_size = lb_lbfluid_get_shape();
    Utils::Vector6d tensor{};
    for (int i = 0; i < grid_size[0]; i++) {
      for (int j = 0; j < grid_size[1]; j++) {
        for (int k = 0; k < grid_size[2]; k++) {
          const Utils::Vector3i node{{i, j, k}};
          tensor += lb_lbnode_get_pressure_tensor(node);
        }
      }
    }

    tensor /= static_cast<double>(Utils::product(grid_size));
    return tensor;
  }
  throw NoLBActive();
}

int lb_lbnode_get_boundary(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    unsigned int host_flag;
    auto const single_nodeindex = calculate_node_index(lbpar_gpu, ind);
    lb_get_boundary_flag_GPU(single_nodeindex, &host_flag);
    return static_cast<int>(host_flag);
#else
    return {};
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return mpi_call(::Communication::Result::one_rank, mpi_lb_get_boundary_flag,
                    ind);
  }
  throw NoLBActive();
}

const Utils::Vector19d lb_lbnode_get_pop(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    float population[D3Q19::n_vel];

    lb_lbfluid_get_population(ind, population);
    Utils::Vector19d p_pop;
    for (std::size_t i = 0; i < D3Q19::n_vel; ++i)
      p_pop[i] = static_cast<double>(population[i]);
    return p_pop;
#else
    return {};
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return mpi_call(::Communication::Result::one_rank, mpi_lb_get_populations,
                    ind);
  }
  throw NoLBActive();
}

void lb_lbnode_set_density(const Utils::Vector3i &ind, double p_density) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    auto const single_nodeindex = calculate_node_index(lbpar_gpu, ind);
    auto const host_density = static_cast<float>(p_density);
    lb_set_node_rho_GPU(single_nodeindex, host_density);
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    auto const tensor = lb_lbnode_get_pressure_tensor(ind);
    auto const momentum_density =
        lb_lbnode_get_velocity(ind) * lb_lbnode_get_density(ind);
    auto const population =
        lb_get_population_from_density_momentum_density_stress(
            p_density, momentum_density, tensor);
    mpi_call_all(mpi_lb_set_population, ind, population);
  } else {
    throw NoLBActive();
  }
}

void lb_lbnode_set_velocity(const Utils::Vector3i &ind,
                            const Utils::Vector3d &u) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    float host_velocity[3];
    host_velocity[0] = static_cast<float>(u[0]);
    host_velocity[1] = static_cast<float>(u[1]);
    host_velocity[2] = static_cast<float>(u[2]);
    auto const single_nodeindex = calculate_node_index(lbpar_gpu, ind);
    lb_set_node_velocity_GPU(single_nodeindex, host_velocity);
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    auto const density = lb_lbnode_get_density(ind);
    auto const momentum_density = u * density;
    auto const tensor = lb_lbnode_get_pressure_tensor(ind);
    auto const population =
        lb_get_population_from_density_momentum_density_stress(
            density, momentum_density, tensor);
    mpi_call_all(mpi_lb_set_population, ind, population);
    mpi_call_all(mpi_lb_set_force_density, ind, Utils::Vector3d{});
  } else {
    throw NoLBActive();
  }
}

void lb_lbnode_set_pop(const Utils::Vector3i &ind,
                       const Utils::Vector19d &p_pop) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    float population[D3Q19::n_vel];

    for (std::size_t i = 0; i < D3Q19::n_vel; ++i)
      population[i] = static_cast<float>(p_pop[i]);

    lb_lbfluid_set_population(ind, population);
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    mpi_call_all(mpi_lb_set_population, ind, p_pop);
  } else {
    throw NoLBActive();
  }
}

const Lattice &lb_lbfluid_get_lattice() { return lblattice; }

ActiveLB lb_lbfluid_get_lattice_switch() { return lattice_switch; }

static void mpi_lb_lbfluid_calc_fluid_momentum_local() {
  lb_calc_fluid_momentum(nullptr, lbpar, lbfields, lblattice);
}

REGISTER_CALLBACK(mpi_lb_lbfluid_calc_fluid_momentum_local)

Utils::Vector3d lb_lbfluid_calc_fluid_momentum() {
  Utils::Vector3d fluid_momentum{};
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lb_calc_fluid_momentum_GPU(fluid_momentum.data());
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    mpi_call(mpi_lb_lbfluid_calc_fluid_momentum_local);
    lb_calc_fluid_momentum(fluid_momentum.data(), lbpar, lbfields, lblattice);
  }
  return fluid_momentum;
}

const Utils::Vector3d
lb_lbfluid_get_interpolated_velocity(const Utils::Vector3d &pos) {
  auto const folded_pos = folded_position(pos, box_geo);
  auto const interpolation_order = lb_lbinterpolation_get_interpolation_order();
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    Utils::Vector3d interpolated_u{};
    switch (interpolation_order) {
    case (InterpolationOrder::linear):
      lb_get_interpolated_velocity_gpu<8>(folded_pos.data(),
                                          interpolated_u.data(), 1);
      break;
    case (InterpolationOrder::quadratic):
      lb_get_interpolated_velocity_gpu<27>(folded_pos.data(),
                                           interpolated_u.data(), 1);
      break;
    }
    return interpolated_u;
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    switch (interpolation_order) {
    case (InterpolationOrder::quadratic):
      throw std::runtime_error("The non-linear interpolation scheme is not "
                               "implemented for the CPU LB.");
    case (InterpolationOrder::linear):
      return mpi_call(::Communication::Result::one_rank,
                      mpi_lb_get_interpolated_velocity, folded_pos);
    }
  }
  throw NoLBActive();
}

double lb_lbfluid_get_interpolated_density(const Utils::Vector3d &pos) {
  auto const folded_pos = folded_position(pos, box_geo);
  auto const interpolation_order = lb_lbinterpolation_get_interpolation_order();
  if (lattice_switch == ActiveLB::GPU) {
    throw std::runtime_error(
        "Density interpolation is not implemented for the GPU LB.");
  }
  if (lattice_switch == ActiveLB::CPU) {
    switch (interpolation_order) {
    case (InterpolationOrder::quadratic):
      throw std::runtime_error("The non-linear interpolation scheme is not "
                               "implemented for the CPU LB.");
    case (InterpolationOrder::linear):
      return mpi_call(::Communication::Result::one_rank,
                      mpi_lb_get_interpolated_density, folded_pos);
    }
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
