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
#include "MpiCallbacks.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "electrokinetics.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "lb-d3q19.hpp"
#include "lb.hpp"
#include "lb_collective_interface.hpp"
#include "lb_interpolation.hpp"
#include "lbgpu.hpp"

#include <utils/index.hpp>
using Utils::get_linear_index;

#include <cstdio>
#include <fstream>

ActiveLB lattice_switch = ActiveLB::NONE;

struct NoLBActive : public std::exception {
  const char *what() const noexcept override { return "LB not activated"; }
};

void lb_lbfluid_update() {
  if (lattice_switch == ActiveLB::CPU) {
    lattice_boltzmann_update();
  } else if (lattice_switch == ActiveLB::GPU and this_node == 0) {
#ifdef CUDA
#ifdef ELECTROKINETICS
    if (ek_initialized) {
      ek_integrate();
    } else {
#endif
      lattice_boltzmann_update_gpu();
#ifdef ELECTROKINETICS
    }
#endif
#endif
  }
}

void lb_lbfluid_propagate() {
  if (lattice_switch != ActiveLB::NONE) {
    lb_lbfluid_update();
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

void lb_lbfluid_sanity_checks() {
  extern double time_step;
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
  lb_lbfluid_sanity_checks();
  if (lattice_switch == ActiveLB::CPU) {
    halo_communication(&update_halo_comm,
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

/** (Re-)initialize the fluid according to the value of density. */

void lb_lbfluid_reinit_fluid() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lb_reinit_fluid_gpu();
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lb_reinit_fluid(lbfields, lblattice, lbpar);
  } else {
    throw NoLBActive();
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
    lb_lbfluid_reinit_fluid();
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
#else
    return {};
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
    return {{lbpar_gpu.ext_force_density[0], lbpar_gpu.ext_force_density[1],
             lbpar_gpu.ext_force_density[2]}};
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

void check_tau_time_step_consistency(double tau, double time_s) {
  auto const eps = std::numeric_limits<float>::epsilon();
  if ((tau - time_s) / (tau + time_s) < -eps)
    throw std::invalid_argument("LB tau (" + std::to_string(tau) +
                                ") must be >= MD time_step (" +
                                std::to_string(time_s) + ")");
  auto const factor = tau / time_s;
  if (fabs(round(factor) - factor) / factor > eps)
    throw std::invalid_argument("LB tau (" + std::to_string(tau) +
                                ") must be integer multiple of "
                                "MD time_step (" +
                                std::to_string(time_s) + "). Factor is " +
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
  lattice_switch = local_lattice_switch;
  mpi_bcast_parameter(FIELD_LATTICE_SWITCH);
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
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
  }

  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    std::vector<unsigned int> bound_array(lbpar_gpu.number_of_nodes);
    lb_get_boundary_flags_GPU(bound_array.data());

    fprintf(fp,
            "# vtk DataFile Version 2.0\nlbboundaries\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %u %u %u\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %u\n"
            "SCALARS boundary float 1\nLOOKUP_TABLE default\n",
            lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z,
            lbpar_gpu.agrid * 0.5, lbpar_gpu.agrid * 0.5, lbpar_gpu.agrid * 0.5,
            lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid,
            lbpar_gpu.number_of_nodes);
    for (int j = 0; j < int(lbpar_gpu.number_of_nodes); ++j) {
      fprintf(fp, "%d \n", bound_array[j]);
    }
#endif //  CUDA
  } else {
    Utils::Vector3i pos;
    auto const grid_size = lblattice.global_grid;

    fprintf(fp,
            "# vtk DataFile Version 2.0\nlbboundaries\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\n"
            "SCALARS boundary float 1\nLOOKUP_TABLE default\n",
            grid_size[0], grid_size[1], grid_size[2], lblattice.agrid * 0.5,
            lblattice.agrid * 0.5, lblattice.agrid * 0.5, lblattice.agrid,
            lblattice.agrid, lblattice.agrid,
            grid_size[0] * grid_size[1] * grid_size[2]);

    for (pos[2] = 0; pos[2] < grid_size[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < grid_size[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < grid_size[0]; pos[0]++) {
          auto boundary = lb_lbnode_get_boundary(pos);
          fprintf(fp, "%d \n", boundary);
        }
      }
    }
  }
  fclose(fp);
}

void lb_lbfluid_print_vtk_velocity(const std::string &filename,
                                   std::vector<int> bb1, std::vector<int> bb2) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
  }

  std::vector<int> bb_low;
  std::vector<int> bb_high;

  for (auto val1 = bb1.begin(), val2 = bb2.begin();
       val1 != bb1.end() && val2 != bb2.end(); ++val1, ++val2) {
    if (*val1 == -1 || *val2 == -1) {
      bb_low = {0, 0, 0};
      if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
        bb_high = {static_cast<int>(lbpar_gpu.dim_x) - 1,
                   static_cast<int>(lbpar_gpu.dim_y) - 1,
                   static_cast<int>(lbpar_gpu.dim_z) - 1};
#endif //  CUDA
      } else {
        bb_high = {lblattice.global_grid[0] - 1, lblattice.global_grid[1] - 1,
                   lblattice.global_grid[2] - 1};
      }
      break;
    }

    bb_low.push_back(std::min(*val1, *val2));
    bb_high.push_back(std::max(*val1, *val2));
  }

  Utils::Vector3i pos;
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    host_values.resize(lbpar_gpu.number_of_nodes);
    lb_get_values_GPU(host_values.data());
    auto const lattice_speed = lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
    fprintf(fp,
            "# vtk DataFile Version 2.0\nlbfluid_gpu\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\n"
            "SCALARS velocity float 3\nLOOKUP_TABLE default\n",
            bb_high[0] - bb_low[0] + 1, bb_high[1] - bb_low[1] + 1,
            bb_high[2] - bb_low[2] + 1, (bb_low[0] + 0.5) * lbpar_gpu.agrid,
            (bb_low[1] + 0.5) * lbpar_gpu.agrid,
            (bb_low[2] + 0.5) * lbpar_gpu.agrid, lbpar_gpu.agrid,
            lbpar_gpu.agrid, lbpar_gpu.agrid,
            (bb_high[0] - bb_low[0] + 1) * (bb_high[1] - bb_low[1] + 1) *
                (bb_high[2] - bb_low[2] + 1));
    for (pos[2] = bb_low[2]; pos[2] <= bb_high[2]; pos[2]++)
      for (pos[1] = bb_low[1]; pos[1] <= bb_high[1]; pos[1]++)
        for (pos[0] = bb_low[0]; pos[0] <= bb_high[0]; pos[0]++) {
          auto const j =
              static_cast<int>(lbpar_gpu.dim_y * lbpar_gpu.dim_x * pos[2] +
                               lbpar_gpu.dim_x * pos[1] + pos[0]);
          fprintf(fp, "%f %f %f\n", host_values[j].v[0] * lattice_speed,
                  host_values[j].v[1] * lattice_speed,
                  host_values[j].v[2] * lattice_speed);
        }
#endif //  CUDA
  } else {
    fprintf(fp,
            "# vtk DataFile Version 2.0\nlbfluid_cpu\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\n"
            "SCALARS velocity float 3\nLOOKUP_TABLE default\n",
            bb_high[0] - bb_low[0] + 1, bb_high[1] - bb_low[1] + 1,
            bb_high[2] - bb_low[2] + 1, (bb_low[0] + 0.5) * lblattice.agrid,
            (bb_low[1] + 0.5) * lblattice.agrid,
            (bb_low[2] + 0.5) * lblattice.agrid, lblattice.agrid,
            lblattice.agrid, lblattice.agrid,
            (bb_high[0] - bb_low[0] + 1) * (bb_high[1] - bb_low[1] + 1) *
                (bb_high[2] - bb_low[2] + 1));

    for (pos[2] = bb_low[2]; pos[2] <= bb_high[2]; pos[2]++)
      for (pos[1] = bb_low[1]; pos[1] <= bb_high[1]; pos[1]++)
        for (pos[0] = bb_low[0]; pos[0] <= bb_high[0]; pos[0]++) {
          auto u = lb_lbnode_get_velocity(pos) * lb_lbfluid_get_lattice_speed();
          fprintf(fp, "%f %f %f\n", u[0], u[1], u[2]);
        }
  }
  fclose(fp);
}

void lb_lbfluid_print_boundary(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
  }

  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    std::vector<unsigned int> bound_array(lbpar_gpu.number_of_nodes);
    lb_get_boundary_flags_GPU(bound_array.data());

    Utils::Vector3i xyz;
    for (int j = 0; j < static_cast<int>(lbpar_gpu.number_of_nodes); ++j) {
      xyz[0] = j % lbpar_gpu.dim_x;
      auto k = j / lbpar_gpu.dim_x;
      xyz[1] = k % lbpar_gpu.dim_y;
      k /= lbpar_gpu.dim_y;
      xyz[2] = k;
      fprintf(fp, "%f %f %f %u\n", (xyz[0] + 0.5) * lbpar_gpu.agrid,
              (xyz[1] + 0.5) * lbpar_gpu.agrid,
              (xyz[2] + 0.5) * lbpar_gpu.agrid, bound_array[j]);
    }
#endif //  CUDA
  } else {
    Utils::Vector3i pos;

    for (pos[2] = 0; pos[2] < lblattice.global_grid[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < lblattice.global_grid[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < lblattice.global_grid[0]; pos[0]++) {
          auto boundary = lb_lbnode_get_boundary(pos);
          boundary = (boundary != 0 ? 1 : 0);
          fprintf(fp, "%f %f %f %d\n", (pos[0] + 0.5) * lblattice.agrid,
                  (pos[1] + 0.5) * lblattice.agrid,
                  (pos[2] + 0.5) * lblattice.agrid, boundary);
        }
      }
    }
  }
  fclose(fp);
}

void lb_lbfluid_print_velocity(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
  }

  auto const lattice_speed = lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
  auto const agrid = lb_lbfluid_get_agrid();
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    std::vector<LB_rho_v_pi_gpu> host_values(lbpar_gpu.number_of_nodes);
    lb_get_values_GPU(host_values.data());
    Utils::Vector3i xyz;
    int j;
    for (j = 0; j < int(lbpar_gpu.number_of_nodes); ++j) {
      xyz[0] = j % lbpar_gpu.dim_x;
      auto k = j / lbpar_gpu.dim_x;
      xyz[1] = k % lbpar_gpu.dim_y;
      k /= lbpar_gpu.dim_y;
      xyz[2] = k;
      fprintf(fp, "%f %f %f %f %f %f\n", (xyz[0] + 0.5) * agrid,
              (xyz[1] + 0.5) * agrid, (xyz[2] + 0.5) * agrid,
              host_values[j].v[0] * lattice_speed,
              host_values[j].v[1] * lattice_speed,
              host_values[j].v[2] * lattice_speed);
    }
#endif //  CUDA
  } else {
    Utils::Vector3i pos;

    for (pos[2] = 0; pos[2] < lblattice.global_grid[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < lblattice.global_grid[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < lblattice.global_grid[0]; pos[0]++) {
          auto const u = lb_lbnode_get_velocity(pos) * lattice_speed;
          fprintf(fp, "%f %f %f %f %f %f\n", (pos[0] + 0.5) * agrid,
                  (pos[1] + 0.5) * agrid, (pos[2] + 0.5) * agrid, u[0], u[1],
                  u[2]);
        }
      }
    }
  }

  fclose(fp);
}

void lb_lbfluid_save_checkpoint(const std::string &filename, bool binary) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    std::vector<float> host_checkpoint_vd(19 * lbpar_gpu.number_of_nodes);
    lb_save_checkpoint_GPU(host_checkpoint_vd.data());
    if (!binary) {
      std::fstream cpfile(filename, std::ios::out);
      cpfile << std::fixed;
      cpfile.precision(8);
      cpfile << lbpar_gpu.dim_x << " ";
      cpfile << lbpar_gpu.dim_y << " ";
      cpfile << lbpar_gpu.dim_z << "\n";
      for (int n = 0; n < (19 * int(lbpar_gpu.number_of_nodes)); n++) {
        cpfile << host_checkpoint_vd[n] << "\n";
      }
      cpfile.close();
    } else {
      std::fstream cpfile(filename, std::ios::out | std::ios::binary);
      cpfile.write(reinterpret_cast<char *>(&lbpar_gpu.dim_x),
                   sizeof(lbpar_gpu.dim_x));
      cpfile.write(reinterpret_cast<char *>(&lbpar_gpu.dim_y),
                   sizeof(lbpar_gpu.dim_y));
      cpfile.write(reinterpret_cast<char *>(&lbpar_gpu.dim_z),
                   sizeof(lbpar_gpu.dim_z));
      cpfile.write(reinterpret_cast<char *>(host_checkpoint_vd.data()),
                   19 * sizeof(float) * lbpar_gpu.number_of_nodes);
      cpfile.close();
    }
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    std::fstream cpfile;
    if (binary) {
      cpfile.open(filename, std::ios::out | std::ios::binary);
    } else {
      cpfile.open(filename, std::ios::out);
      cpfile.precision(16);
      cpfile << std::fixed;
    }

    auto const gridsize = lblattice.global_grid;

    if (!binary) {
      cpfile << gridsize[0] << " " << gridsize[1] << " " << gridsize[2] << "\n";
    } else {
      cpfile.write(reinterpret_cast<const char *>(gridsize.data()),
                   3 * sizeof(gridsize[0]));
    }

    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          Utils::Vector3i ind{{i, j, k}};
          auto const pop = mpi_call(::Communication::Result::one_rank,
                                    mpi_lb_get_populations, ind);
          if (!binary) {
            for (auto const &p : pop) {
              cpfile << p << "\n";
            }
          } else {
            cpfile.write(reinterpret_cast<const char *>(&pop[0]),
                         pop.size() * sizeof(double));
          }
        }
      }
    }
    cpfile.close();
  }
}

void lb_lbfluid_load_checkpoint(const std::string &filename, bool binary) {
  int res;
  std::string err_msg = "Error while reading LB checkpoint: ";
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    FILE *cpfile;
    cpfile = fopen(filename.c_str(), "r");
    if (!cpfile) {
      throw std::runtime_error(err_msg + "could not open file for reading.");
    }
    std::vector<float> host_checkpoint_vd(lbpar_gpu.number_of_nodes * 19);
    if (!binary) {
      int saved_gridsize[3];
      for (int &n : saved_gridsize) {
        res = fscanf(cpfile, "%i", &n);
        if (res == EOF) {
          fclose(cpfile);
          throw std::runtime_error(err_msg + "EOF found.");
        }
        if (res != 1) {
          fclose(cpfile);
          throw std::runtime_error(err_msg + "incorrectly formatted data.");
        }
      }
      if (saved_gridsize[0] != lbpar_gpu.dim_x ||
          saved_gridsize[1] != lbpar_gpu.dim_y ||
          saved_gridsize[2] != lbpar_gpu.dim_z) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "grid dimensions mismatch, read [" +
                                 std::to_string(saved_gridsize[0]) + ' ' +
                                 std::to_string(saved_gridsize[1]) + ' ' +
                                 std::to_string(saved_gridsize[2]) +
                                 "], expected [" +
                                 std::to_string(lbpar_gpu.dim_x) + ' ' +
                                 std::to_string(lbpar_gpu.dim_y) + ' ' +
                                 std::to_string(lbpar_gpu.dim_z) + "].");
      }
      for (int n = 0; n < (19 * int(lbpar_gpu.number_of_nodes)); n++) {
        res = fscanf(cpfile, "%f", &host_checkpoint_vd[n]);
        if (res == EOF) {
          fclose(cpfile);
          throw std::runtime_error(err_msg + "EOF found.");
        }
        if (res != 1) {
          fclose(cpfile);
          throw std::runtime_error(err_msg + "incorrectly formatted data.");
        }
      }
    } else {
      Utils::Vector3i saved_gridsize;
      if (fread(saved_gridsize.data(), sizeof(int), 3, cpfile) != 3) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
      if (saved_gridsize[0] != lbpar_gpu.dim_x ||
          saved_gridsize[1] != lbpar_gpu.dim_y ||
          saved_gridsize[2] != lbpar_gpu.dim_z) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "grid dimensions mismatch, read [" +
                                 std::to_string(saved_gridsize[0]) + ' ' +
                                 std::to_string(saved_gridsize[1]) + ' ' +
                                 std::to_string(saved_gridsize[2]) +
                                 "], expected [" +
                                 std::to_string(lbpar_gpu.dim_x) + ' ' +
                                 std::to_string(lbpar_gpu.dim_y) + ' ' +
                                 std::to_string(lbpar_gpu.dim_z) + "].");
      }
      if (fread(host_checkpoint_vd.data(), sizeof(float),
                19 * int(lbpar_gpu.number_of_nodes),
                cpfile) != (unsigned int)(19 * lbpar_gpu.number_of_nodes)) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
    }
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
    lb_load_checkpoint_GPU(host_checkpoint_vd.data());
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    FILE *cpfile;
    cpfile = fopen(filename.c_str(), "r");
    if (!cpfile) {
      throw std::runtime_error(err_msg + "could not open file for reading.");
    }

    auto const gridsize = lblattice.global_grid;
    int saved_gridsize[3];
    mpi_bcast_lb_params(LBParam::DENSITY);

    if (!binary) {
      res = fscanf(cpfile, "%i %i %i\n", &saved_gridsize[0], &saved_gridsize[1],
                   &saved_gridsize[2]);
      if (res == EOF) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "EOF found.");
      }
      if (res != 3) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
    } else {
      if (fread(&saved_gridsize[0], sizeof(int), 3, cpfile) != 3) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
    }
    if (saved_gridsize[0] != gridsize[0] || saved_gridsize[1] != gridsize[1] ||
        saved_gridsize[2] != gridsize[2]) {
      fclose(cpfile);
      throw std::runtime_error(err_msg + "grid dimensions mismatch, read [" +
                               std::to_string(saved_gridsize[0]) + ' ' +
                               std::to_string(saved_gridsize[1]) + ' ' +
                               std::to_string(saved_gridsize[2]) +
                               "], expected [" + std::to_string(gridsize[0]) +
                               ' ' + std::to_string(gridsize[1]) + ' ' +
                               std::to_string(gridsize[2]) + "].");
    }

    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          Utils::Vector3i ind{{i, j, k}};
          Utils::Vector19d pop;
          if (!binary) {
            res = fscanf(cpfile,
                         "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
                         "%lf %lf %lf %lf %lf %lf \n",
                         &pop[0], &pop[1], &pop[2], &pop[3], &pop[4], &pop[5],
                         &pop[6], &pop[7], &pop[8], &pop[9], &pop[10], &pop[11],
                         &pop[12], &pop[13], &pop[14], &pop[15], &pop[16],
                         &pop[17], &pop[18]);
            if (res == EOF) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "EOF found.");
            }
            if (res != 19) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "incorrectly formatted data.");
            }
          } else {
            if (fread(pop.data(), sizeof(double), 19, cpfile) != 19) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "incorrectly formatted data.");
            }
          }
          lb_lbnode_set_pop(ind, pop);
        }
      }
    }
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
  } else {
    throw std::runtime_error(
        "To load an LB checkpoint one needs to have already "
        "initialized the LB fluid with the same grid size.");
  }
}

Utils::Vector3i lb_lbfluid_get_shape() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return {static_cast<int>(lbpar_gpu.dim_x),
            static_cast<int>(lbpar_gpu.dim_y),
            static_cast<int>(lbpar_gpu.dim_z)};
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
    auto const single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                                  ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
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
    auto const single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                                  ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_print_node_GPU(single_nodeindex, &host_print_values);
    return {{host_print_values.v[0], host_print_values.v[1],
             host_print_values.v[2]}};
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

const Utils::Vector6d lb_lbnode_get_stress(const Utils::Vector3i &ind) {
  // Add equilibrium stress to the diagonal (in LB units)
  auto const p0 = lb_lbfluid_get_density() * D3Q19::c_sound_sq<double>;

  auto stress = lb_lbnode_get_stress_neq(ind);
  stress[0] += p0;
  stress[2] += p0;
  stress[5] += p0;

  return stress;
}

const Utils::Vector6d lb_lbnode_get_stress_neq(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    Utils::Vector6d stress{};
    static LB_rho_v_pi_gpu host_print_values;
    auto const single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                                  ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_print_node_GPU(single_nodeindex, &host_print_values);
    for (int i = 0; i < 6; i++) {
      stress[i] = host_print_values.pi[i];
    }
    return stress;
#endif //  CUDA
  }
  if (lattice_switch == ActiveLB::CPU) {
    return mpi_call(::Communication::Result::one_rank, mpi_lb_get_stress, ind);
  }
  throw NoLBActive();
}

const Utils::Vector6d lb_lbfluid_get_stress() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    Utils::Vector6d stress{};
    // Copy observable data from gpu
    host_values.resize(lbpar_gpu.number_of_nodes);
    lb_get_values_GPU(host_values.data());
    std::for_each(host_values.begin(), host_values.end(),
                  [&stress](LB_rho_v_pi_gpu &v) {
                    for (int i = 0; i < 6; i++)
                      stress[i] += v.pi[i];
                  });

    // Normalize
    stress /= static_cast<double>(lbpar_gpu.number_of_nodes);

    // Add equilibrium stress to the diagonal (in LB units)
    double const p0 = lb_lbfluid_get_density() * D3Q19::c_sound_sq<double>;

    stress[0] += p0;
    stress[2] += p0;
    stress[5] += p0;
    return stress;
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    Utils::Vector6d stress{};
    for (int i = 0; i < lblattice.global_grid[0]; i++) {
      for (int j = 0; j < lblattice.global_grid[1]; j++) {
        for (int k = 0; k < lblattice.global_grid[2]; k++) {
          const Utils::Vector3i node{{i, j, k}};
          stress += lb_lbnode_get_stress(node);
        }
      }
    }

    auto const number_of_nodes = lblattice.global_grid[0] *
                                 lblattice.global_grid[1] *
                                 lblattice.global_grid[2];

    stress /= static_cast<double>(number_of_nodes);
    return stress;
  }
  throw NoLBActive();
}

int lb_lbnode_get_boundary(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    unsigned int host_flag;
    auto const single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                                  ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
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
    float population[19];

    lb_lbfluid_get_population(ind, population);
    Utils::Vector19d p_pop;
    for (int i = 0; i < LBQ; ++i)
      p_pop[i] = population[i];
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
    auto const single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                                  ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    auto const host_density = static_cast<float>(p_density);
    lb_set_node_rho_GPU(single_nodeindex, host_density);
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    auto const stress = lb_lbnode_get_stress(ind);
    auto const momentum_density =
        lb_lbnode_get_velocity(ind) * lb_lbnode_get_density(ind);
    auto const population =
        lb_get_population_from_density_momentum_density_stress(
            p_density, momentum_density, stress);
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
    auto const single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                                  ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_set_node_velocity_GPU(single_nodeindex, host_velocity);
#endif //  CUDA
  } else if (lattice_switch == ActiveLB::CPU) {
    auto const density = lb_lbnode_get_density(ind);
    auto const momentum_density = u * density;
    auto const stress = lb_lbnode_get_stress(ind);
    auto const population =
        lb_get_population_from_density_momentum_density_stress(
            density, momentum_density, stress);
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
    float population[19];

    for (int i = 0; i < LBQ; ++i)
      population[i] = p_pop[i];

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

Utils::Vector3d lb_lbfluid_calc_fluid_momentum() {
  Utils::Vector3d fluid_momentum{};
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lb_calc_fluid_momentum_GPU(fluid_momentum.data());
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    mpi_gather_stats(6, fluid_momentum.data(), nullptr, nullptr, nullptr);
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
