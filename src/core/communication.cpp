/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <mpi.h>
#ifdef OPEN_MPI
#include <dlfcn.h>
#endif
#include <cassert>

#include "communication.hpp"

#include "errorhandling.hpp"

#include "CellStructure.hpp"
#include "EspressoSystemInterface.hpp"
#include "cells.hpp"
#include "cuda_interface.hpp"
#include "energy.hpp"
#include "event.hpp"
#include "galilei.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/icc.hpp"

#include <boost/mpi.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <utils/mpi/cart_comm.hpp>

using namespace std;

namespace Communication {
auto const &mpi_datatype_cache = boost::mpi::detail::mpi_datatype_cache();
std::shared_ptr<boost::mpi::environment> mpi_env;
} // namespace Communication

boost::mpi::communicator comm_cart;

namespace Communication {
std::unique_ptr<MpiCallbacks> m_callbacks;

/* We use a singleton callback class for now. */
MpiCallbacks &mpiCallbacks() {
  assert(m_callbacks && "Mpi not initialized!");

  return *m_callbacks;
}
} // namespace Communication

using Communication::mpiCallbacks;

int this_node = -1;
int n_nodes = -1;

// if you want to add a callback, add it here, and here only
#define CALLBACK_LIST                                                          \
  CB(mpi_who_has_slave)                                                        \
  CB(mpi_gather_stats_slave)                                                   \
  CB(mpi_bcast_coulomb_params_slave)                                           \
  CB(mpi_bcast_cuda_global_part_vars_slave)                                    \
  CB(mpi_resort_particles_slave)                                               \
  CB(mpi_get_pairs_slave)                                                      \
  CB(mpi_get_particles_slave)                                                  \
  CB(mpi_rotate_system_slave)                                                  \
  CB(mpi_update_particle_slave)

// create the forward declarations
#define CB(name) void name(int node, int param);
#ifndef DOXYGEN
/* this conditional on DOXYGEN prevents an interaction in Doxygen between
 * CALLBACK_LIST and whatever follows next, e.g. a function "int foo();"
 * would otherwise become "CALLBACK_LIST int foo();" */
CALLBACK_LIST
#endif

#undef CB

/**********************************************
 * procedures
 **********************************************/

#if defined(OPEN_MPI)
namespace {
/** Workaround for "Read -1, expected XXXXXXX, errno = 14" that sometimes
 *  appears when CUDA is used. This is a bug in OpenMPI 2.0-2.1.2 and 3.0.0
 *  according to
 *  https://www.mail-archive.com/users@lists.open-mpi.org/msg32357.html,
 *  so we set btl_vader_single_copy_mechanism = none.
 */
void openmpi_fix_vader() {
  if ((OMPI_MAJOR_VERSION == 2 && OMPI_MINOR_VERSION == 1 &&
       OMPI_RELEASE_VERSION < 3) or
      (OMPI_MAJOR_VERSION == 3 && OMPI_MINOR_VERSION == 0 &&
       OMPI_RELEASE_VERSION == 0)) {
    setenv("OMPI_MCA_btl_vader_single_copy_mechanism", "none", 0);
  }
}

/**
 * @brief Assure that openmpi is loaded to the global namespace.
 *
 * This was originally inspired by mpi4py
 * (https://github.com/mpi4py/mpi4py/blob/4e3f47b6691c8f5a038e73f84b8d43b03f16627f/src/lib-mpi/compat/openmpi.h).
 * It's needed because OpenMPI dlopens its submodules. These are unable to find
 * the top-level OpenMPI library if that was dlopened itself, i.e. when the
 * Python interpreter dlopens a module that is linked against OpenMPI. It's
 * about some weird two-level symbol namespace thing.
 */
void openmpi_global_namespace() {
  if (OMPI_MAJOR_VERSION >= 3)
    return;
#ifdef RTLD_NOLOAD
  const int mode = RTLD_NOW | RTLD_GLOBAL | RTLD_NOLOAD;
#else
  const int mode = RTLD_NOW | RTLD_GLOBAL;
#endif

  const void *_openmpi_symbol = dlsym(RTLD_DEFAULT, "MPI_Init");
  if (!_openmpi_symbol) {
    fprintf(stderr, "Aborting because unable to find OpenMPI symbol.\n");
    errexit();
  }

  Dl_info _openmpi_info;
  dladdr(_openmpi_symbol, &_openmpi_info);

  const void *handle = dlopen(_openmpi_info.dli_fname, mode);

  if (!handle) {
    fprintf(stderr, "Aborting because unable to load libmpi into the "
                    "global symbol space.\n");
    errexit();
  }
}
} // namespace
#endif

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

#define CB(name) Communication::m_callbacks->add(&(name));
  CALLBACK_LIST
#undef CB

  ErrorHandling::init_error_handling(mpiCallbacks());

  on_program_start();
}
} // namespace Communication

std::shared_ptr<boost::mpi::environment> mpi_init() {
#ifdef OPEN_MPI
  openmpi_fix_vader();
  openmpi_global_namespace();
#endif

  return std::make_shared<boost::mpi::environment>();
}

/********************* STEEPEST DESCENT ********/
static int mpi_steepest_descent_slave(int steps, int) {
  return integrate(steps, -1);
}
REGISTER_CALLBACK_MASTER_RANK(mpi_steepest_descent_slave)

int mpi_steepest_descent(int steps) {
  return mpi_call(Communication::Result::master_rank,
                  mpi_steepest_descent_slave, steps, 0);
}

/********************* INTEGRATE ********/
static int mpi_integrate_slave(int n_steps, int reuse_forces) {
  integrate(n_steps, reuse_forces);

  return check_runtime_errors_local();
}
REGISTER_CALLBACK_REDUCTION(mpi_integrate_slave, std::plus<int>())

int mpi_integrate(int n_steps, int reuse_forces) {
  return mpi_call(Communication::Result::reduction, std::plus<int>(),
                  mpi_integrate_slave, n_steps, reuse_forces);
}

/*************** GATHER ************/
void mpi_gather_stats(GatherStats job, double *result) {
  auto job_slave = static_cast<int>(job);
  switch (job) {
  case GatherStats::energy:
    mpi_call(mpi_gather_stats_slave, -1, job_slave);
    energy_calc(sim_time);
    break;
  case GatherStats::pressure:
    mpi_call(mpi_gather_stats_slave, -1, job_slave);
    pressure_calc();
    break;
  case GatherStats::lb_fluid_momentum:
    mpi_call(mpi_gather_stats_slave, -1, job_slave);
    lb_calc_fluid_momentum(result, lbpar, lbfields, lblattice);
    break;
#ifdef LB_BOUNDARIES
  case GatherStats::lb_boundary_forces:
    mpi_call(mpi_gather_stats_slave, -1, job_slave);
    lb_collect_boundary_forces(result);
    break;
#endif
  default:
    fprintf(
        stderr,
        "%d: INTERNAL ERROR: illegal request %d for mpi_gather_stats_slave\n",
        this_node, job_slave);
    errexit();
  }
}

void mpi_gather_stats_slave(int, int job_slave) {
  auto job = static_cast<GatherStats>(job_slave);
  switch (job) {
  case GatherStats::energy:
    energy_calc(sim_time);
    break;
  case GatherStats::pressure:
    pressure_calc();
    break;
  case GatherStats::lb_fluid_momentum:
    lb_calc_fluid_momentum(nullptr, lbpar, lbfields, lblattice);
    break;
#ifdef LB_BOUNDARIES
  case GatherStats::lb_boundary_forces:
    lb_collect_boundary_forces(nullptr);
    break;
#endif
  default:
    fprintf(
        stderr,
        "%d: INTERNAL ERROR: illegal request %d for mpi_gather_stats_slave\n",
        this_node, job_slave);
    errexit();
  }
}

/*************** TIME STEP ************/
void mpi_set_time_step_slave(double dt) {
  time_step = dt;

  on_parameter_change(FIELD_TIMESTEP);
}
REGISTER_CALLBACK(mpi_set_time_step_slave)

void mpi_set_time_step(double time_s) {
  if (time_s <= 0.)
    throw std::invalid_argument("time_step must be > 0.");
  if (lb_lbfluid_get_lattice_switch() != ActiveLB::NONE)
    check_tau_time_step_consistency(lb_lbfluid_get_tau(), time_s);
  mpi_call_all(mpi_set_time_step_slave, time_s);
}

/*************** BCAST COULOMB ************/
void mpi_bcast_coulomb_params() {
#if defined(ELECTROSTATICS) || defined(DIPOLES)
  mpi_call(mpi_bcast_coulomb_params_slave, 1, 0);
  mpi_bcast_coulomb_params_slave(-1, 0);
#endif
}

void mpi_bcast_coulomb_params_slave(int, int) {

#ifdef ELECTROSTATICS
  MPI_Bcast(&coulomb, sizeof(Coulomb_parameters), MPI_BYTE, 0, comm_cart);

  Coulomb::bcast_coulomb_params();
#endif

#ifdef DIPOLES
  MPI_Bcast(&dipole, sizeof(Dipole_parameters), MPI_BYTE, 0, comm_cart);

  Dipole::set_method_local(dipole.method);

  Dipole::bcast_params(comm_cart);
#endif

#if defined(ELECTROSTATICS) || defined(DIPOLES)
  on_coulomb_change();
#endif
}

/*************** BCAST CELL STRUCTURE *****************/

REGISTER_CALLBACK(cells_re_init)

void mpi_bcast_cell_structure(int cs) { mpi_call_all(cells_re_init, cs); }

REGISTER_CALLBACK(cells_set_use_verlet_lists)

void mpi_set_use_verlet_lists(bool use_verlet_lists) {
  mpi_call_all(cells_set_use_verlet_lists, use_verlet_lists);
}

/******************* BCAST CUDA GLOBAL PART VARS ********************/

void mpi_bcast_cuda_global_part_vars() {
#ifdef CUDA
  mpi_call(mpi_bcast_cuda_global_part_vars_slave, 1,
           0); // third parameter is meaningless
  mpi_bcast_cuda_global_part_vars_slave(-1, 0);
#endif
}

void mpi_bcast_cuda_global_part_vars_slave(int, int) {
#ifdef CUDA
  MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
            sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart);
  espressoSystemInterface.requestParticleStructGpu();
#endif
}

/********************* SET EXCLUSION ********/
#ifdef EXCLUSIONS
void mpi_send_exclusion_slave(int part1, int part2, int _delete) {
  local_change_exclusion(part1, part2, _delete);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_send_exclusion_slave)

void mpi_send_exclusion(int part1, int part2, int _delete) {
  mpi_call(mpi_send_exclusion_slave, part1, part2, _delete);
  mpi_send_exclusion_slave(part1, part2, _delete);
}
#endif

/********************* ICCP3M INIT ********/
#ifdef ELECTROSTATICS
void mpi_iccp3m_init_slave(const iccp3m_struct &iccp3m_cfg_) {
  iccp3m_cfg = iccp3m_cfg_;

  on_particle_charge_change();
  check_runtime_errors(comm_cart);
}

REGISTER_CALLBACK(mpi_iccp3m_init_slave)

int mpi_iccp3m_init() {
  mpi_call(mpi_iccp3m_init_slave, iccp3m_cfg);

  on_particle_charge_change();
  return check_runtime_errors(comm_cart);
}
#endif

/***** GALILEI TRANSFORM AND ASSOCIATED FUNCTIONS ****/
void mpi_kill_particle_motion_slave(int rotation) {
  local_kill_particle_motion(rotation, cell_structure.local_particles());
  on_particle_change();
}

REGISTER_CALLBACK(mpi_kill_particle_motion_slave)

void mpi_kill_particle_motion(int rotation) {
  mpi_call_all(mpi_kill_particle_motion_slave, rotation);
}

void mpi_kill_particle_forces_slave(int torque) {
  local_kill_particle_forces(torque, cell_structure.local_particles());
  on_particle_change();
}

REGISTER_CALLBACK(mpi_kill_particle_forces_slave)

void mpi_kill_particle_forces(int torque) {
  mpi_call_all(mpi_kill_particle_forces_slave, torque);
}

struct pair_sum {
  template <class T, class U>
  auto operator()(std::pair<T, U> l, std::pair<T, U> r) const {
    return std::pair<T, U>{l.first + r.first, l.second + r.second};
  }
};

Utils::Vector3d mpi_system_CMS() {
  auto const data =
      mpi_call(Communication::Result::reduction, pair_sum{}, local_system_CMS);
  return data.first / data.second;
}

REGISTER_CALLBACK_REDUCTION(local_system_CMS_velocity, pair_sum{})

Utils::Vector3d mpi_system_CMS_velocity() {
  auto const data = mpi_call(Communication::Result::reduction, pair_sum{},
                             local_system_CMS_velocity);
  return data.first / data.second;
}

REGISTER_CALLBACK_REDUCTION(local_system_CMS, pair_sum{})

void mpi_galilei_transform_slave(Utils::Vector3d const &cmsvel) {
  local_galilei_transform(cmsvel);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_galilei_transform_slave)

void mpi_galilei_transform() {
  auto const cmsvel = mpi_system_CMS_velocity();

  mpi_call_all(mpi_galilei_transform_slave, cmsvel);
}

/*********************** MAIN LOOP for slaves ****************/

void mpi_loop() {
  if (this_node != 0)
    mpiCallbacks().loop();
}

std::vector<int> mpi_resort_particles(int global_flag) {
  mpi_call(mpi_resort_particles_slave, global_flag, 0);
  cell_structure.resort_particles(global_flag);

  clear_particle_node();

  std::vector<int> n_parts;
  boost::mpi::gather(comm_cart,
                     static_cast<int>(cell_structure.local_particles().size()),
                     n_parts, 0);

  return n_parts;
}

void mpi_resort_particles_slave(int global_flag, int) {
  cell_structure.resort_particles(global_flag);

  boost::mpi::gather(
      comm_cart, static_cast<int>(cell_structure.local_particles().size()), 0);
}
