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
#include <cstring>
#include <memory>
#include <mpi.h>
#ifdef OPEN_MPI
#include <dlfcn.h>
#endif
#include <cassert>

#include "communication.hpp"

#include "errorhandling.hpp"

#include "EspressoSystemInterface.hpp"
#include "Particle.hpp"
#include "bonded_interactions/bonded_tab.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "cuda_interface.hpp"
#include "energy.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "galilei.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_interpolation.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "integrate.hpp"
#include "integrators/steepest_descent.hpp"
#include "io/mpiio/mpiio.hpp"
#include "nonbonded_interactions/nonbonded_tab.hpp"
#include "npt.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include "rotation.hpp"
#include "statistics.hpp"
#include "virtual_sites.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/icc.hpp"
#include "electrostatics_magnetostatics/mdlc_correction.hpp"

#include "serialization/IA_parameters.hpp"

#include <utils/Counter.hpp>
#include <utils/u32_to_u64.hpp>

#include <boost/mpi.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <utils/mpi/cart_comm.hpp>

using namespace std;

namespace Communication {
auto const &mpi_datatype_cache = boost::mpi::detail::mpi_datatype_cache();
std::unique_ptr<boost::mpi::environment> mpi_env;
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
  CB(mpi_place_particle_slave)                                                 \
  CB(mpi_bcast_ia_params_slave)                                                \
  CB(mpi_gather_stats_slave)                                                   \
  CB(mpi_bcast_coulomb_params_slave)                                           \
  CB(mpi_remove_particle_slave)                                                \
  CB(mpi_rescale_particles_slave)                                              \
  CB(mpi_bcast_cell_structure_slave)                                           \
  CB(mpi_bcast_nptiso_geom_slave)                                              \
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

int mpi_check_runtime_errors();

/**********************************************
 * procedures
 **********************************************/

#if defined(OPEN_MPI)
/** Workaround for "Read -1, expected XXXXXXX, errno = 14" that sometimes
 *  appears when CUDA is used. This is a bug in OpenMPI 2.0-2.1.2 and 3.0.0
 *  according to
 *  https://www.mail-archive.com/users@lists.open-mpi.org/msg32357.html,
 *  so we set btl_vader_single_copy_mechanism = none.
 */
static void openmpi_fix_vader() {
  if (OMPI_MAJOR_VERSION < 2 || OMPI_MAJOR_VERSION > 3)
    return;
  if (OMPI_MAJOR_VERSION == 2 && OMPI_MINOR_VERSION == 1 &&
      OMPI_RELEASE_VERSION >= 3)
    return;
  if (OMPI_MAJOR_VERSION == 3 &&
      (OMPI_MINOR_VERSION > 0 || OMPI_RELEASE_VERSION > 0))
    return;

  std::string varname = "btl_vader_single_copy_mechanism";
  std::string varval = "none";

  setenv((std::string("OMPI_MCA_") + varname).c_str(), varval.c_str(), 0);
}
#endif

void mpi_init() {
#ifdef OPEN_MPI
  openmpi_fix_vader();

  void *handle = nullptr;
  int mode = RTLD_NOW | RTLD_GLOBAL;
#ifdef RTLD_NOLOAD
  mode |= RTLD_NOLOAD;
#endif
  void *_openmpi_symbol = dlsym(RTLD_DEFAULT, "MPI_Init");
  if (!_openmpi_symbol) {
    fprintf(stderr, "%d: Aborting because unable to find OpenMPI symbol.\n",
            this_node);
    errexit();
  }
  Dl_info _openmpi_info;
  dladdr(_openmpi_symbol, &_openmpi_info);

  if (!handle)
    handle = dlopen(_openmpi_info.dli_fname, mode);

  if (!handle) {
    fprintf(stderr,
            "%d: Aborting because unable to load libmpi into the "
            "global symbol space.\n",
            this_node);
    errexit();
  }
#endif

#ifdef BOOST_MPI_HAS_NOARG_INITIALIZATION
  Communication::mpi_env = std::make_unique<boost::mpi::environment>();
#else
  int argc{};
  char **argv{};
  Communication::mpi_env =
      std::make_unique<boost::mpi::environment>(argc, argv);
#endif

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

/****************** PLACE/PLACE NEW PARTICLE ************/

void mpi_place_particle(int node, int id, const Utils::Vector3d &pos) {
  mpi_call(mpi_place_particle_slave, node, id);

  if (node == this_node)
    local_place_particle(id, pos, 0);
  else {
    comm_cart.send(node, SOME_TAG, pos);
  }

  cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();
}

void mpi_place_particle_slave(int pnode, int part) {
  if (pnode == this_node) {
    Utils::Vector3d pos;
    comm_cart.recv(0, SOME_TAG, pos);
    local_place_particle(part, pos, 0);
  }

  cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();
}

boost::optional<int> mpi_place_new_particle_slave(int part,
                                                  Utils::Vector3d const &pos) {
  auto p = local_place_particle(part, pos, 1);

  on_particle_change();

  if (p) {
    return comm_cart.rank();
  }

  return {};
}

REGISTER_CALLBACK_ONE_RANK(mpi_place_new_particle_slave)

int mpi_place_new_particle(int id, const Utils::Vector3d &pos) {
  return mpi_call(Communication::Result::one_rank, mpi_place_new_particle_slave,
                  id, pos);
}

/****************** REMOVE PARTICLE ************/
void mpi_remove_particle(int pnode, int part) {
  mpi_call_all(mpi_remove_particle_slave, pnode, part);
}

void mpi_remove_particle_slave(int pnode, int part) {
  if (part != -1) {
    cell_structure.remove_particle(part);
  } else {
    cell_structure.remove_all_particles();
  }

  on_particle_change();
}

/********************* STEEPEST DESCENT ********/
static int mpi_steepest_descent_slave(int steps, int) {
  return steepest_descent(steps);
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

/*************** BCAST IA ************/
static void mpi_bcast_all_ia_params_slave() {
  boost::mpi::broadcast(comm_cart, ia_params, 0);
}

REGISTER_CALLBACK(mpi_bcast_all_ia_params_slave)

void mpi_bcast_all_ia_params() { mpi_call_all(mpi_bcast_all_ia_params_slave); }

void mpi_bcast_ia_params(int i, int j) {
  mpi_call(mpi_bcast_ia_params_slave, i, j);

  if (j >= 0) {
    /* non-bonded interaction parameters */
    boost::mpi::broadcast(comm_cart, *get_ia_param(i, j), 0);
  } else {
    /* bonded interaction parameters */
    MPI_Bcast(&(bonded_ia_params[i]), sizeof(Bonded_ia_parameters), MPI_BYTE, 0,
              comm_cart);
    /* For tabulated potentials we have to send the tables extra */
    if (bonded_ia_params[i].type == BONDED_IA_TABULATED_DISTANCE or
        bonded_ia_params[i].type == BONDED_IA_TABULATED_ANGLE or
        bonded_ia_params[i].type == BONDED_IA_TABULATED_DIHEDRAL) {
      boost::mpi::broadcast(comm_cart, *bonded_ia_params[i].p.tab.pot, 0);
    }
  }

  on_short_range_ia_change();
}

void mpi_bcast_ia_params_slave(int i, int j) {
  if (j >= 0) { /* non-bonded interaction parameters */

    boost::mpi::broadcast(comm_cart, *get_ia_param(i, j), 0);
  } else {                   /* bonded interaction parameters */
    make_bond_type_exist(i); /* realloc bonded_ia_params on slave nodes! */
    MPI_Bcast(&(bonded_ia_params[i]), sizeof(Bonded_ia_parameters), MPI_BYTE, 0,
              comm_cart);
    /* For tabulated potentials we have to send the tables extra */
    if (bonded_ia_params[i].type == BONDED_IA_TABULATED_DISTANCE or
        bonded_ia_params[i].type == BONDED_IA_TABULATED_ANGLE or
        bonded_ia_params[i].type == BONDED_IA_TABULATED_DIHEDRAL) {
      auto *tab_pot = new TabulatedPotential();
      boost::mpi::broadcast(comm_cart, *tab_pot, 0);

      bonded_ia_params[i].p.tab.pot = tab_pot;
    }
  }

  on_short_range_ia_change();
}

/*************** BCAST IA SIZE ************/

REGISTER_CALLBACK(realloc_ia_params)
void mpi_bcast_max_seen_particle_type(int ns) {
  mpi_call_all(realloc_ia_params, ns);
}

/*************** GATHER ************/
void mpi_gather_stats(int job, void *result, void *result_t, void *result_nb,
                      void *result_t_nb) {
  switch (job) {
  case 1:
    mpi_call(mpi_gather_stats_slave, -1, 1);
    energy_calc((double *)result, sim_time);
    break;
  case 2:
    /* calculate and reduce (sum up) virials for 'analyze pressure' or
       'analyze stress_tensor' */
    mpi_call(mpi_gather_stats_slave, -1, 2);
    pressure_calc((double *)result, (double *)result_t, (double *)result_nb,
                  (double *)result_t_nb, 0);
    break;
  case 3:
    mpi_call(mpi_gather_stats_slave, -1, 3);
    pressure_calc((double *)result, (double *)result_t, (double *)result_nb,
                  (double *)result_t_nb, 1);
    break;
  case 6:
    mpi_call(mpi_gather_stats_slave, -1, 6);
    lb_calc_fluid_momentum((double *)result, lbpar, lbfields, lblattice);
    break;
  case 7:
    break;
#ifdef LB_BOUNDARIES
  case 8:
    mpi_call(mpi_gather_stats_slave, -1, 8);
    lb_collect_boundary_forces((double *)result);
    break;
#endif
  default:
    fprintf(
        stderr,
        "%d: INTERNAL ERROR: illegal request %d for mpi_gather_stats_slave\n",
        this_node, job);
    errexit();
  }
}

void mpi_gather_stats_slave(int, int job) {
  switch (job) {
  case 1:
    /* calculate and reduce (sum up) energies */
    energy_calc(nullptr, sim_time);
    break;
  case 2:
    /* calculate and reduce (sum up) virials for 'analyze pressure' or 'analyze
     * stress_tensor'*/
    pressure_calc(nullptr, nullptr, nullptr, nullptr, 0);
    break;
  case 3:
    /* calculate and reduce (sum up) virials, revert velocities half a timestep
     * for 'analyze p_inst' */
    pressure_calc(nullptr, nullptr, nullptr, nullptr, 1);
    break;
  case 6:
    lb_calc_fluid_momentum(nullptr, lbpar, lbfields, lblattice);
    break;
  case 7:
    break;
#ifdef LB_BOUNDARIES
  case 8:
    lb_collect_boundary_forces(nullptr);
    break;
#endif
  default:
    fprintf(
        stderr,
        "%d: INTERNAL ERROR: illegal request %d for mpi_gather_stats_slave\n",
        this_node, job);
    errexit();
  }
}

/*************** TIME STEP ************/
void mpi_set_time_step_slave(double dt) {
  time_step = dt;
  time_step_squared = time_step * time_step;
  time_step_squared_half = time_step_squared / 2.;
  time_step_half = time_step / 2.;

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

/****************** RESCALE PARTICLES ************/

void mpi_rescale_particles(int dir, double scale) {
  int pnode;

  mpi_call(mpi_rescale_particles_slave, -1, dir);
  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      local_rescale_particles(dir, scale);
    } else {
      MPI_Send(&scale, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
    }
  }
  on_particle_change();
}

void mpi_rescale_particles_slave(int, int dir) {
  double scale = 0.0;
  MPI_Recv(&scale, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  local_rescale_particles(dir, scale);
  on_particle_change();
}

/*************** BCAST CELL STRUCTURE *****************/

void mpi_bcast_cell_structure(int cs) {
  mpi_call(mpi_bcast_cell_structure_slave, -1, cs);
  cells_re_init(cs, cell_structure.min_range);
}

void mpi_bcast_cell_structure_slave(int, int cs) {
  cells_re_init(cs, cell_structure.min_range);
}

/*************** BCAST NPTISO GEOM *****************/

void mpi_bcast_nptiso_geom() {
  mpi_call(mpi_bcast_nptiso_geom_slave, -1, 0);
  mpi_bcast_nptiso_geom_slave(-1, 0);
}

void mpi_bcast_nptiso_geom_slave(int, int) {
  MPI_Bcast(&nptiso.geometry, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&nptiso.dimension, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&nptiso.cubic_box, 1, MPI_LOGICAL, 0, comm_cart);
  MPI_Bcast(&nptiso.non_const_dim, 1, MPI_INT, 0, comm_cart);
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

/********************* CALC MU MAX ********/

#ifdef DP3M
REGISTER_CALLBACK(calc_mu_max)

void mpi_bcast_max_mu() { mpi_call_all(calc_mu_max); }
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
  cells_resort_particles(global_flag);

  clear_particle_node();

  std::vector<int> n_parts;
  boost::mpi::gather(comm_cart, cells_get_n_particles(), n_parts, 0);

  return n_parts;
}

void mpi_resort_particles_slave(int global_flag, int) {
  cells_resort_particles(global_flag);

  boost::mpi::gather(comm_cart, cells_get_n_particles(), 0);
}
