/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#ifdef OPEN_MPI
#include <dlfcn.h>
#endif
#include <cassert>

#include "communication.hpp"

#include "errorhandling.hpp"

#include "EspressoSystemInterface.hpp"
#include "buckingham.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "cuda_interface.hpp"
#include "debye_hueckel.hpp"
#include "elc.hpp"
#include "energy.hpp"
#include "forces.hpp"
#include "galilei.hpp"
#include "gb.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "iccp3m.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "io/mpiio/mpiio.hpp"
#include "lb.hpp"
#include "lbboundaries.hpp"
#include "lbboundaries/LBBoundary.hpp"
#include "lj.hpp"
#include "ljcos.hpp"
#include "ljcos2.hpp"
#include "ljgen.hpp"
#include "maggs.hpp"
#include "mdlc_correction.hpp"
#include "minimize_energy.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "molforces.hpp"
#include "morse.hpp"
#include "npt.hpp"
#include "p3m-dipolar.hpp"
#include "p3m.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include "reaction_field.hpp"
#include "rotation.hpp"
#include "scafacos.hpp"
#include "statistics.hpp"
#include "statistics_chain.hpp"
#include "statistics_fluid.hpp"
#include "swimmer_reaction.hpp"
#include "tab.hpp"
#include "topology.hpp"
#include "virtual_sites.hpp"

#include "utils.hpp"
#include "utils/make_unique.hpp"
#include "utils/serialization/IA_parameters.hpp"
#include "utils/serialization/Particle.hpp"

#include <boost/mpi.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>

using namespace std;

namespace Communication {
auto const &mpi_datatype_cache = boost::mpi::detail::mpi_datatype_cache();
std::unique_ptr<boost::mpi::environment> mpi_env;
}

boost::mpi::communicator comm_cart;

namespace Communication {
std::unique_ptr<MpiCallbacks> m_callbacks;

/* We use a singelton callback class for now. */
MpiCallbacks &mpiCallbacks() {
  assert(m_callbacks && "Mpi not initialized!");

  return *m_callbacks;
}
}

using Communication::mpiCallbacks;

int this_node = -1;
int n_nodes = -1;

int graceful_exit = 0;
/* whether there is already a termination going on. */
static int terminated = 0;

// if you want to add a callback, add it here, and here only
#define CALLBACK_LIST                                                          \
  CB(mpi_bcast_parameter_slave)                                                \
  CB(mpi_who_has_slave)                                                        \
  CB(mpi_bcast_event_slave)                                                    \
  CB(mpi_place_particle_slave)                                                 \
  CB(mpi_send_v_slave)                                                         \
  CB(mpi_send_swimming_slave)                                                  \
  CB(mpi_send_f_slave)                                                         \
  CB(mpi_send_q_slave)                                                         \
  CB(mpi_send_type_slave)                                                      \
  CB(mpi_send_bond_slave)                                                      \
  CB(mpi_recv_part_slave)                                                      \
  CB(mpi_integrate_slave)                                                      \
  CB(mpi_bcast_ia_params_slave)                                                \
  CB(mpi_bcast_all_ia_params_slave)                                            \
  CB(mpi_bcast_max_seen_particle_type_slave)                                   \
  CB(mpi_gather_stats_slave)                                                   \
  CB(mpi_set_time_step_slave)                                                  \
  CB(mpi_bcast_coulomb_params_slave)                                           \
  CB(mpi_bcast_collision_params_slave)                                         \
  CB(mpi_send_ext_force_slave)                                                 \
  CB(mpi_send_ext_torque_slave)                                                \
  CB(mpi_place_new_particle_slave)                                             \
  CB(mpi_remove_particle_slave)                                                \
  CB(mpi_rescale_particles_slave)                                              \
  CB(mpi_bcast_cell_structure_slave)                                           \
  CB(mpi_send_quat_slave)                                                      \
  CB(mpi_send_omega_slave)                                                     \
  CB(mpi_send_torque_slave)                                                    \
  CB(mpi_send_mol_id_slave)                                                    \
  CB(mpi_bcast_nptiso_geom_slave)                                              \
  CB(mpi_update_mol_ids_slave)                                                 \
  CB(mpi_sync_topo_part_info_slave)                                            \
  CB(mpi_send_mass_slave)                                                      \
  CB(mpi_send_solvation_slave)                                                 \
  CB(mpi_send_exclusion_slave)                                                 \
  CB(mpi_bcast_lb_params_slave)                                                \
  CB(mpi_bcast_cuda_global_part_vars_slave)                                    \
  CB(mpi_send_dip_slave)                                                       \
  CB(mpi_send_dipm_slave)                                                      \
  CB(mpi_send_fluid_slave)                                                     \
  CB(mpi_recv_fluid_slave)                                                     \
  CB(mpi_local_stress_tensor_slave)                                            \
  CB(mpi_send_virtual_slave)                                                   \
  CB(mpi_iccp3m_iteration_slave)                                               \
  CB(mpi_iccp3m_init_slave)                                                    \
  CB(mpi_send_rotational_inertia_slave)                                        \
  CB(mpi_send_affinity_slave)                                                  \
  CB(mpi_rotate_particle_slave)                                                \
  CB(mpi_send_out_direction_slave)                                             \
  CB(mpi_send_mu_E_slave)                                                      \
  CB(mpi_bcast_max_mu_slave)                                                   \
  CB(mpi_send_vs_quat_slave)                                                   \
  CB(mpi_send_vs_relative_slave)                                               \
  CB(mpi_recv_fluid_populations_slave)                                         \
  CB(mpi_send_fluid_populations_slave)                                         \
  CB(mpi_recv_fluid_boundary_flag_slave)                                       \
  CB(mpi_set_particle_temperature_slave)                                       \
  CB(mpi_set_particle_gamma_slave)                                             \
  CB(mpi_set_particle_gamma_rot_slave)                                         \
  CB(mpi_kill_particle_motion_slave)                                           \
  CB(mpi_kill_particle_forces_slave)                                           \
  CB(mpi_system_CMS_slave)                                                     \
  CB(mpi_system_CMS_velocity_slave)                                            \
  CB(mpi_galilei_transform_slave)                                              \
  CB(mpi_setup_reaction_slave)                                                 \
  CB(mpi_send_rotation_slave)                                                  \
  CB(mpi_check_runtime_errors_slave)                                           \
  CB(mpi_minimize_energy_slave)                                                \
  CB(mpi_gather_cuda_devices_slave)                                            \
  CB(mpi_scafacos_set_parameters_slave)                                        \
  CB(mpi_scafacos_set_r_cut_and_tune_slave)                                    \
  CB(mpi_scafacos_free_slave)                                                  \
  CB(mpi_resort_particles_slave)                                               \
  CB(mpi_get_pairs_slave)                                                      \
  CB(mpi_get_particles_slave)                                                  \
  CB(mpi_rotate_system_slave)

// create the forward declarations
#define CB(name) void name(int node, int param);
CALLBACK_LIST

namespace {

// create the list of callbacks
#undef CB
#define CB(name) name,
std::vector<SlaveCallback *> slave_callbacks{CALLBACK_LIST};

/** The callback name list is only used for
    debugging.
*/
#ifdef COMM_DEBUG
// create the list of names
#undef CB
#define CB(name) #name,

std::vector<std::string> names{CALLBACK_LIST};
#endif
}

/** Forward declarations */

int mpi_check_runtime_errors(void);

/**********************************************
 * procedures
 **********************************************/

void mpi_init() {
#ifdef OPEN_MPI
  void *handle = 0;
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
    fprintf(stderr, "%d: Aborting because unable to load libmpi into the "
                    "global symbol space.\n",
            this_node);
    errexit();
  }
#endif

#ifdef BOOST_MPI_HAS_NOARG_INITIALIZATION
  Communication::mpi_env = Utils::make_unique<boost::mpi::environment>();
#else
  int argc{};
  char **argv{};
  Communication::mpi_env =
      Utils::make_unique<boost::mpi::environment>(argc, argv);
#endif

  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, node_grid);

  mpi_reshape_communicator({{node_grid[0], node_grid[1], node_grid[2]}},
                           /* periodicity */ {{1, 1, 1}});
  MPI_Cart_coords(comm_cart, this_node, 3, node_pos);

  Communication::m_callbacks =
      Utils::make_unique<Communication::MpiCallbacks>(comm_cart);

  for (int i = 0; i < slave_callbacks.size(); ++i) {
    mpiCallbacks().add(slave_callbacks[i]);
  }

  ErrorHandling::init_error_handling(mpiCallbacks());
  partCfg(Utils::make_unique<PartCfg>(mpiCallbacks(), GetLocalParts()));

  on_program_start();
}

void mpi_reshape_communicator(std::array<int, 3> const &node_grid,
                              std::array<int, 3> const &periodicity) {
  MPI_Comm temp_comm;
  MPI_Cart_create(MPI_COMM_WORLD, 3, const_cast<int *>(node_grid.data()),
                  const_cast<int *>(periodicity.data()), 0, &temp_comm);
  comm_cart =
      boost::mpi::communicator(temp_comm, boost::mpi::comm_take_ownership);

  this_node = comm_cart.rank();
}

void mpi_call(SlaveCallback cb, int node, int param) {
#ifdef COMM_DEBUG
  auto it = std::find(slave_callbacks.begin(), slave_callbacks.end(), cb);

  if (it != slave_callbacks.end()) {
    auto const id = it - slave_callbacks.begin();
    COMM_TRACE(fprintf(stderr, "%d: issuing %s %d %d\n", this_node,
                       names[id].c_str(), node, param));
  }
#endif /* COMM_DEBUG */
  mpiCallbacks().call(cb, node, param);

  COMM_TRACE(fprintf(stderr, "%d: finished sending.\n", this_node));
}

/**************** REQ_CHTOPL ***********/
void mpi_bcast_event(int event) {
  mpi_call(mpi_bcast_event_slave, -1, event);
  mpi_bcast_event_slave(-1, event);
}

void mpi_bcast_event_slave(int node, int event) {
  switch (event) {
#ifdef ELECTROSTATICS
#ifdef P3M
  case P3M_COUNT_CHARGES:
    p3m_count_charged_particles();
    break;
#endif
  case MAGGS_COUNT_CHARGES:
    maggs_count_charged_particles();
    break;
#endif
  case CHECK_PARTICLES:
    check_particles();
    break;

#ifdef DP3M
  case P3M_COUNT_DIPOLES:
    dp3m_count_magnetic_particles();
    break;
#endif

  default:;
  }
}

/****************** REQ_PLACE/REQ_PLACE_NEW ************/

void mpi_place_particle(int pnode, int part, double p[3]) {
  mpi_call(mpi_place_particle_slave, pnode, part);

  if (pnode == this_node)
    local_place_particle(part, p, 0);
  else
    MPI_Send(p, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);

  set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();
}

void mpi_place_particle_slave(int pnode, int part) {

  if (pnode == this_node) {
    double p[3];
    MPI_Recv(p, 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    local_place_particle(part, p, 0);
  }

  set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();
}

void mpi_place_new_particle(int pnode, int part, double p[3]) {
  mpi_call(mpi_place_new_particle_slave, pnode, part);
  added_particle(part);

  if (pnode == this_node)
    local_place_particle(part, p, 1);
  else
    MPI_Send(p, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_place_new_particle_slave(int pnode, int part) {

  added_particle(part);

  if (pnode == this_node) {
    double p[3];
    MPI_Recv(p, 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    local_place_particle(part, p, 1);
  }

  on_particle_change();
}

/****************** REQ_SET_V ************/
void mpi_send_v(int pnode, int part, double* v) {
  mpi_call(mpi_send_v_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];

    p->m.v = {v[0],v[1],v[2]};
  } else
    MPI_Send(v, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_send_v_slave(int pnode, int part) {
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->m.v.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
}

/****************** REQ_SET_SWIMMING ************/
void mpi_send_swimming(int pnode, int part, ParticleParametersSwimming swim) {
#ifdef ENGINE
  mpi_call(mpi_send_swimming_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->swim = swim;
  } else {
    MPI_Send(&swim, sizeof(ParticleParametersSwimming), MPI_BYTE, pnode,
             SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_swimming_slave(int pnode, int part) {
#ifdef ENGINE
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&p->swim, sizeof(ParticleParametersSwimming), MPI_BYTE, 0,
             SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/****************** REQ_SET_F ************/
void mpi_send_f(int pnode, int part, const Vector3d &F) {
  mpi_call(mpi_send_f_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->f.f = F;
  } else {
    comm_cart.send(pnode, SOME_TAG, F);
  }

  on_particle_change();
}

void mpi_send_f_slave(int pnode, int part) {
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    Vector3d F;
    comm_cart.recv(0, SOME_TAG, F);

    p->f.f = F;
  }

  on_particle_change();
}

/********************* REQ_SET_Q ********/
void mpi_send_q(int pnode, int part, double q) {
#ifdef ELECTROSTATICS
  mpi_call(mpi_send_q_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.q = q;
  } else {
    MPI_Send(&q, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_charge_change();
#endif
}

void mpi_send_q_slave(int pnode, int part) {
#ifdef ELECTROSTATICS
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&p->p.q, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_charge_change();
#endif
}

/********************* REQ_SET_MU_E ********/
void mpi_send_mu_E(int pnode, int part, double mu_E[3]) {
#ifdef LB_ELECTROHYDRODYNAMICS
  mpi_call(mpi_send_mu_E_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.mu_E[0] = mu_E[0];
    p->p.mu_E[1] = mu_E[1];
    p->p.mu_E[2] = mu_E[2];
  } else {
    MPI_Send(mu_E, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_mu_E_slave(int pnode, int part) {
#ifdef LB_ELECTROHYDRODYNAMICS
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->p.mu_E.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_SOLV ********/
void mpi_send_solvation(int pnode, int part, double *solvation) {
#ifdef SHANCHEN
  mpi_call(mpi_send_solvation_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    for (int ii = 0; ii < 2 * LB_COMPONENTS; ii++)
      p->p.solvation[ii] = solvation[ii];
  } else {
    MPI_Send(&solvation, LB_COMPONENTS, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_solvation_slave(int pnode, int part) {
#ifdef SHANCHEN
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&p->p.solvation, 2 * LB_COMPONENTS, MPI_DOUBLE, 0, SOME_TAG,
             comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_M ********/
void mpi_send_mass(int pnode, int part, double mass) {
#ifdef MASS
  mpi_call(mpi_send_mass_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.mass = mass;
  } else {
    MPI_Send(&mass, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_mass_slave(int pnode, int part) {
#ifdef MASS
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&p->p.mass, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_RINERTIA ********/

void mpi_send_rotational_inertia(int pnode, int part, double rinertia[3]) {
#ifdef ROTATIONAL_INERTIA
  mpi_call(mpi_send_rotational_inertia_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.rinertia[0] = rinertia[0];
    p->p.rinertia[1] = rinertia[1];
    p->p.rinertia[2] = rinertia[2];
  } else {
    MPI_Send(rinertia, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_rotational_inertia_slave(int pnode, int part) {
#ifdef ROTATIONAL_INERTIA
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->p.rinertia.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

void mpi_rotate_particle(int pnode, int part, double axis[3], double angle) {
#ifdef ROTATION
  mpi_call(mpi_rotate_particle_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    local_rotate_particle(p, axis, angle);
  } else {
    MPI_Send(axis, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
    MPI_Send(&angle, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_rotate_particle_slave(int pnode, int part) {
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    double axis[3], angle;
    MPI_Recv(axis, 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    MPI_Recv(&angle, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    local_rotate_particle(p, axis, angle);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_BOND_SITE ********/

void mpi_send_affinity_slave(int pnode, int part) {
#ifdef AFFINITY
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->p.bond_site.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

void mpi_send_affinity(int pnode, int part, double bond_site[3]) {
#ifdef AFFINITY
  mpi_call(mpi_send_affinity_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.bond_site[0] = bond_site[0];
    p->p.bond_site[1] = bond_site[1];
    p->p.bond_site[2] = bond_site[2];
  } else {
    MPI_Send(bond_site, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_OUT_DIRECTION ********/

void mpi_send_out_direction(int pnode, int part, double out_direction[3]) {
#ifdef MEMBRANE_COLLISION
  mpi_call(mpi_send_out_direction_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.out_direction[0] = out_direction[0];
    p->p.out_direction[1] = out_direction[1];
    p->p.out_direction[2] = out_direction[2];
  } else {
    MPI_Send(out_direction, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_out_direction_slave(int pnode, int part) {
#ifdef MEMBRANE_COLLISION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->p.out_direction.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_TYPE ********/
void mpi_send_type(int pnode, int part, int type) {
  mpi_call(mpi_send_type_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.type = type;
  } else
    MPI_Send(&type, 1, MPI_INT, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_send_type_slave(int pnode, int part) {
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&p->p.type, 1, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
}

/********************* REQ_SET_MOLID ********/
void mpi_send_mol_id(int pnode, int part, int mid) {
  mpi_call(mpi_send_mol_id_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.mol_id = mid;
  } else
    MPI_Send(&mid, 1, MPI_INT, pnode, SOME_TAG, comm_cart);

  on_particle_change();
}

void mpi_send_mol_id_slave(int pnode, int part) {
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&p->p.mol_id, 1, MPI_INT, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
}

/********************* REQ_SET_QUAT ********/

void mpi_send_quat(int pnode, int part, double quat[4]) {
#ifdef ROTATION
  mpi_call(mpi_send_quat_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->r.quat[0] = quat[0];
    p->r.quat[1] = quat[1];
    p->r.quat[2] = quat[2];
    p->r.quat[3] = quat[3];
    convert_quat_to_quatu(p->r.quat, p->r.quatu);
#ifdef DIPOLES
    convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
  } else {
    MPI_Send(quat, 4, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_quat_slave(int pnode, int part) {
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->r.quat.data(), 4, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
    convert_quat_to_quatu(p->r.quat, p->r.quatu);
#ifdef DIPOLES
    convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_OMEGA ********/

void mpi_send_omega(int pnode, int part, double omega[3]) {
#ifdef ROTATION
  mpi_call(mpi_send_omega_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /*  memmove(p->omega, omega, 3*sizeof(double));*/
    p->m.omega[0] = omega[0];
    p->m.omega[1] = omega[1];
    p->m.omega[2] = omega[2];
  } else {
    MPI_Send(omega, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_omega_slave(int pnode, int part) {
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->m.omega.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_TORQUE ********/

void mpi_send_torque(int pnode, int part, double torque[3]) {
#ifdef ROTATION
  mpi_call(mpi_send_torque_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->f.torque[0] = torque[0];
    p->f.torque[1] = torque[1];
    p->f.torque[2] = torque[2];
  } else {
    MPI_Send(torque, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_torque_slave(int pnode, int part) {
#ifdef ROTATION
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->f.torque.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_DIP ********/

void mpi_send_dip(int pnode, int part, double dip[3]) {
#ifdef DIPOLES
  mpi_call(mpi_send_dip_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->r.dip[0] = dip[0];
    p->r.dip[1] = dip[1];
    p->r.dip[2] = dip[2];
#ifdef ROTATION
    convert_dip_to_quat(p->r.dip, p->r.quat, &p->p.dipm);
    convert_quat_to_quatu(p->r.quat, p->r.quatu);
#else
    p->p.dipm = sqrt(p->r.dip[0] * p->r.dip[0] + p->r.dip[1] * p->r.dip[1] +
                     p->r.dip[2] * p->r.dip[2]);
#endif
  } else {
    MPI_Send(dip, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_dip_slave(int pnode, int part) {
#ifdef DIPOLES
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->r.dip.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
#ifdef ROTATION
    convert_dip_to_quat(p->r.dip, p->r.quat, &p->p.dipm);
    convert_quat_to_quatu(p->r.quat, p->r.quatu);
#else
    p->p.dipm = sqrt(p->r.dip[0] * p->r.dip[0] + p->r.dip[1] * p->r.dip[1] +
                     p->r.dip[2] * p->r.dip[2]);
#endif
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_DIPM ********/

void mpi_send_dipm(int pnode, int part, double dipm) {
#ifdef DIPOLES
  mpi_call(mpi_send_dipm_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.dipm = dipm;
#ifdef ROTATION
    convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
  } else {
    MPI_Send(&dipm, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_dipm_slave(int pnode, int part) {
#ifdef DIPOLES
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&p->p.dipm, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
#ifdef ROTATION
    convert_quatu_to_dip(p->r.quatu, p->p.dipm, p->r.dip);
#endif
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_ISVI ********/

void mpi_send_virtual(int pnode, int part, int isVirtual) {
#ifdef VIRTUAL_SITES
  mpi_call(mpi_send_virtual_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.isVirtual = isVirtual;
  } else {
    MPI_Send(&isVirtual, 1, MPI_INT, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_virtual_slave(int pnode, int part) {
#ifdef VIRTUAL_SITES
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&(p->p.isVirtual), 1, MPI_INT, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/********************* REQ_SET_BOND ********/
void mpi_send_vs_quat(int pnode, int part, double *vs_quat) {
#ifdef VIRTUAL_SITES_RELATIVE
  mpi_call(mpi_send_vs_quat_slave, pnode, part);
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    for (int i = 0; i < 4; ++i) {
      p->p.vs_quat[i] = vs_quat[i];
    }
  } else {
    MPI_Send(vs_quat, 4, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}
void mpi_send_vs_quat_slave(int pnode, int part) {
#ifdef VIRTUAL_SITES_RELATIVE
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(p->p.vs_quat, 4, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

void mpi_send_vs_relative(int pnode, int part, int vs_relative_to,
                          double vs_distance, double *rel_ori) {
#ifdef VIRTUAL_SITES_RELATIVE
  mpi_call(mpi_send_vs_relative_slave, pnode, part);

  // If the particle is on the node on which this function was called
  // set the values locally
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.vs_relative_to_particle_id = vs_relative_to;
    p->p.vs_relative_distance = vs_distance;
    for (int i = 0; i < 4; i++) {
      p->p.vs_relative_rel_orientation[i] = rel_ori[i];
    }
  } else {
    MPI_Send(&vs_relative_to, 1, MPI_INT, pnode, SOME_TAG, comm_cart);
    MPI_Send(&vs_distance, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
    MPI_Send(rel_ori, 4, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_vs_relative_slave(int pnode, int part) {
#ifdef VIRTUAL_SITES_RELATIVE
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Recv(&p->p.vs_relative_to_particle_id, 1, MPI_INT, 0, SOME_TAG,
             comm_cart, MPI_STATUS_IGNORE);
    MPI_Recv(&p->p.vs_relative_distance, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
    MPI_Recv(p->p.vs_relative_rel_orientation, 4, MPI_DOUBLE, 0, SOME_TAG,
             comm_cart, MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

// ********************************

void mpi_send_rotation(int pnode, int part, short int rot) {
  mpi_call(mpi_send_rotation_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    p->p.rotation = rot;
  } else {
    MPI_Send(&rot, 1, MPI_SHORT, pnode, SOME_TAG, MPI_COMM_WORLD);
  }

  on_particle_change();
}

void mpi_send_rotation_slave(int pnode, int part) {
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&p->p.rotation, 1, MPI_SHORT, 0, SOME_TAG, MPI_COMM_WORLD,
             &status);
  }

  on_particle_change();
}

/********************* REQ_SET_BOND ********/

int mpi_send_bond(int pnode, int part, int *bond, int _delete) {
  int bond_size, stat = 0;

  mpi_call(mpi_send_bond_slave, pnode, part);

  bond_size = (bond) ? bonded_ia_params[bond[0]].num + 1 : 0;

  if (pnode == this_node) {
    stat = local_change_bond(part, bond, _delete);
    on_particle_change();
    return stat;
  }
  /* else */
  MPI_Send(&bond_size, 1, MPI_INT, pnode, SOME_TAG, comm_cart);
  if (bond_size)
    MPI_Send(bond, bond_size, MPI_INT, pnode, SOME_TAG, comm_cart);
  MPI_Send(&_delete, 1, MPI_INT, pnode, SOME_TAG, comm_cart);
  MPI_Recv(&stat, 1, MPI_INT, pnode, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  on_particle_change();
  return stat;
}

void mpi_send_bond_slave(int pnode, int part) {
  int bond_size = 0, _delete = 0, stat;

  if (pnode == this_node) {
    MPI_Recv(&bond_size, 1, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    int *bond;
    if (bond_size) {
      bond = (int *)Utils::malloc(bond_size * sizeof(int));
      MPI_Recv(bond, bond_size, MPI_INT, 0, SOME_TAG, comm_cart,
               MPI_STATUS_IGNORE);
    } else
      bond = nullptr;
    MPI_Recv(&_delete, 1, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    stat = local_change_bond(part, bond, _delete);
    if (bond)
      free(bond);
    MPI_Send(&stat, 1, MPI_INT, 0, SOME_TAG, comm_cart);
  }

  on_particle_change();
}

/****************** REQ_GET_PART ************/
Particle mpi_recv_part(int pnode, int part) {
  Particle ret;

  mpi_call(mpi_recv_part_slave, pnode, part);
  comm_cart.recv(pnode, SOME_TAG, ret);

  return ret;
}

void mpi_recv_part_slave(int pnode, int part) {
  if (pnode != this_node)
    return;

  assert(local_particles[part]);
  comm_cart.send(0, SOME_TAG, *local_particles[part]);
}

/****************** REQ_REM_PART ************/
void mpi_remove_particle(int pnode, int part) {
  mpi_call(mpi_remove_particle_slave, pnode, part);
  mpi_remove_particle_slave(pnode, part);
}

void mpi_remove_particle_slave(int pnode, int part) {
  if (part != -1) {
    n_part--;

    if (pnode == this_node)
      local_remove_particle(part);

    remove_all_bonds_to(part);
  } else
    local_remove_all_particles();

  on_particle_change();
}

/********************* REQ_MIN_ENERGY ********/

int mpi_minimize_energy(void) {
  mpi_call(mpi_minimize_energy_slave, 0, 0);
  return minimize_energy();
}

void mpi_minimize_energy_slave(int a, int b) { minimize_energy(); }

/********************* REQ_INTEGRATE ********/
int mpi_integrate(int n_steps, int reuse_forces) {
  mpi_call(mpi_integrate_slave, n_steps, reuse_forces);
  integrate_vv(n_steps, reuse_forces);
  COMM_TRACE(
      fprintf(stderr, "%d: integration task %d done.\n", this_node, n_steps));
  return mpi_check_runtime_errors();
}

void mpi_integrate_slave(int n_steps, int reuse_forces) {
  integrate_vv(n_steps, reuse_forces);
  COMM_TRACE(fprintf(
      stderr, "%d: integration for %d n_steps with %d reuse_forces done.\n",
      this_node, n_steps, reuse_forces));
}

/*************** REQ_BCAST_IA ************/
void mpi_bcast_all_ia_params() {
  mpi_call(mpi_bcast_all_ia_params_slave, -1, -1);
  boost::mpi::broadcast(comm_cart, ia_params, 0);
}

void mpi_bcast_all_ia_params_slave(int a, int b) {
  boost::mpi::broadcast(comm_cart, ia_params, 0);
}

void mpi_bcast_ia_params(int i, int j) {
  mpi_call(mpi_bcast_ia_params_slave, i, j);

  if (j >= 0) {
    /* non-bonded interaction parameters */
    boost::mpi::broadcast(comm_cart, *get_ia_param(i, j), 0);

    *get_ia_param(j, i) = *get_ia_param(i, j);
  } else {
    /* bonded interaction parameters */
    MPI_Bcast(&(bonded_ia_params[i]), sizeof(Bonded_ia_parameters), MPI_BYTE, 0,
              comm_cart);
#ifdef TABULATED
    /* For tabulated potentials we have to send the tables extra */
    if (bonded_ia_params[i].type == BONDED_IA_TABULATED) {
      boost::mpi::broadcast(comm_cart, *bonded_ia_params[i].p.tab.pot, 0);
    }
#endif
  }

  on_short_range_ia_change();
}

void mpi_bcast_ia_params_slave(int i, int j) {
  if (j >= 0) { /* non-bonded interaction parameters */

    boost::mpi::broadcast(comm_cart, *get_ia_param(i, j), 0);

    *get_ia_param(j, i) = *get_ia_param(i, j);

  } else {                   /* bonded interaction parameters */
    make_bond_type_exist(i); /* realloc bonded_ia_params on slave nodes! */
    MPI_Bcast(&(bonded_ia_params[i]), sizeof(Bonded_ia_parameters), MPI_BYTE, 0,
              comm_cart);
#ifdef TABULATED
    /* For tabulated potentials we have to send the tables extra */
    if (bonded_ia_params[i].type == BONDED_IA_TABULATED) {
      auto *tab_pot = new TabulatedPotential();
      boost::mpi::broadcast(comm_cart, *tab_pot, 0);

      bonded_ia_params[i].p.tab.pot = tab_pot;
    }
#endif
  }

  on_short_range_ia_change();
}

/*************** REQ_BCAST_IA_SIZE ************/

void mpi_bcast_max_seen_particle_type(int ns) {
  mpi_call(mpi_bcast_max_seen_particle_type_slave, -1, ns);
  mpi_bcast_max_seen_particle_type_slave(-1, ns);
}

void mpi_bcast_max_seen_particle_type_slave(int pnode, int ns) {
  realloc_ia_params(ns);
}

/*************** REQ_GATHER ************/
void mpi_gather_stats(int job, void *result, void *result_t, void *result_nb,
                      void *result_t_nb) {
  switch (job) {
  case 1:
    mpi_call(mpi_gather_stats_slave, -1, 1);
    energy_calc((double *)result);
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
  case 4:
    mpi_call(mpi_gather_stats_slave, -1, 4);
    predict_momentum_particles((double *)result);
    break;
#ifdef LB
  case 5:
    mpi_call(mpi_gather_stats_slave, -1, 5);
    lb_calc_fluid_mass((double *)result);
    break;
  case 6:
    mpi_call(mpi_gather_stats_slave, -1, 6);
    lb_calc_fluid_momentum((double *)result);
    break;
  case 7:
    mpi_call(mpi_gather_stats_slave, -1, 7);
    lb_calc_fluid_temp((double *)result);
    break;
#ifdef LB_BOUNDARIES
  case 8:
    mpi_call(mpi_gather_stats_slave, -1, 8);
    lb_collect_boundary_forces((double *)result);
    break;
#endif
#endif
  default:
    fprintf(
        stderr,
        "%d: INTERNAL ERROR: illegal request %d for mpi_gather_stats_slave\n",
        this_node, job);
    errexit();
  }
}

void mpi_gather_stats_slave(int ana_num, int job) {
  switch (job) {
  case 1:
    /* calculate and reduce (sum up) energies */
    energy_calc(nullptr);
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
  case 4:
    predict_momentum_particles(nullptr);
    break;
#ifdef LB
  case 5:
    lb_calc_fluid_mass(nullptr);
    break;
  case 6:
    lb_calc_fluid_momentum(nullptr);
    break;
  case 7:
    lb_calc_fluid_temp(nullptr);
    break;
#ifdef LB_BOUNDARIES
  case 8:
    lb_collect_boundary_forces(nullptr);
    break;
#endif
#endif
  default:
    fprintf(
        stderr,
        "%d: INTERNAL ERROR: illegal request %d for mpi_gather_stats_slave\n",
        this_node, job);
    errexit();
  }
}

/*************** REQ_GET_LOCAL_STRESS_TENSOR ************/
void mpi_local_stress_tensor(DoubleList *TensorInBin, int bins[3],
                             int periodic[3], double range_start[3],
                             double range[3]) {

  PTENSOR_TRACE(fprintf(stderr, "%d: mpi_local_stress_tensor: Broadcasting "
                                "local_stress_tensor parameters\n",
                        this_node));

  mpi_call(mpi_local_stress_tensor_slave, -1, 0);

  MPI_Bcast(bins, 3, MPI_INT, 0, comm_cart);
  MPI_Bcast(periodic, 3, MPI_INT, 0, comm_cart);
  MPI_Bcast(range_start, 3, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(range, 3, MPI_DOUBLE, 0, comm_cart);

  PTENSOR_TRACE(fprintf(
      stderr, "%d: mpi_local_stress_tensor: Call local_stress_tensor_calc\n",
      this_node));
  local_stress_tensor_calc(TensorInBin, bins, periodic, range_start, range);

  PTENSOR_TRACE(fprintf(stderr, "%d: mpi_local_stress_tensor: Reduce local "
                                "stress tensors with MPI_Reduce\n",
                        this_node));
  for (int i = 0; i < bins[0] * bins[1] * bins[2]; i++) {
    MPI_Reduce(MPI_IN_PLACE, TensorInBin[i].e, 9, MPI_DOUBLE, MPI_SUM, 0,
               comm_cart);
  }
}

void mpi_local_stress_tensor_slave(int ana_num, int job) {
  int bins[3] = {0, 0, 0};
  int periodic[3] = {0, 0, 0};
  double range_start[3] = {0, 0, 0};
  double range[3] = {0, 0, 0};
  int i, j;

  MPI_Bcast(bins, 3, MPI_INT, 0, comm_cart);
  MPI_Bcast(periodic, 3, MPI_INT, 0, comm_cart);
  MPI_Bcast(range_start, 3, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(range, 3, MPI_DOUBLE, 0, comm_cart);

  auto TensorInBin =
      std::vector<DoubleList>(bins[0] * bins[1] * bins[2], DoubleList(9, 0.0));

  local_stress_tensor_calc(TensorInBin.data(), bins, periodic, range_start,
                           range);

  for (i = 0; i < bins[0] * bins[1] * bins[2]; i++) {
    MPI_Reduce(TensorInBin[i].e, nullptr, 9, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
    PTENSOR_TRACE(fprintf(
        stderr, "%d: mpi_local_stress_tensor: Tensor sent in bin %d is {",
        this_node, i));
    for (j = 0; j < 9; j++) {
      PTENSOR_TRACE(fprintf(stderr, "%f ", TensorInBin[i].e[j]));
    }
    PTENSOR_TRACE(fprintf(stderr, "}\n"));
  }
}

/*************** REQ_SET_TIME_STEP ************/
void mpi_set_time_step(double time_s) {
  double old_ts = time_step;

  mpi_call(mpi_set_time_step_slave, -1, 0);

  time_step = time_s;

  time_step_squared = time_step * time_step;
  time_step_squared_half = time_step_squared / 2.;
  time_step_half = time_step / 2.;

  MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, comm_cart);

  rescale_velocities(time_step / old_ts);
  on_parameter_change(FIELD_TIMESTEP);
}

void mpi_set_time_step_slave(int node, int i) {
  double old_ts = time_step;

  MPI_Bcast(&time_step, 1, MPI_DOUBLE, 0, comm_cart);
  rescale_velocities(time_step / old_ts);
  on_parameter_change(FIELD_TIMESTEP);
  time_step_squared = time_step * time_step;
  time_step_squared_half = time_step_squared / 2.;
  time_step_half = time_step / 2.;
}

int mpi_check_runtime_errors(void) {
  mpi_call(mpi_check_runtime_errors_slave, 0, 0);
  return check_runtime_errors();
}

void mpi_check_runtime_errors_slave(int a, int b) { check_runtime_errors(); }

/*************** REQ_BCAST_COULOMB ************/
void mpi_bcast_coulomb_params() {
#if defined(ELECTROSTATICS) || defined(DIPOLES)
  mpi_call(mpi_bcast_coulomb_params_slave, 1, 0);
  mpi_bcast_coulomb_params_slave(-1, 0);
#endif
}

void mpi_bcast_coulomb_params_slave(int node, int parm) {

#if defined(ELECTROSTATICS) || defined(DIPOLES)
  MPI_Bcast(&coulomb, sizeof(Coulomb_parameters), MPI_BYTE, 0, comm_cart);

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE:
  // fall through, scafacos has internal parameter propagation
  case COULOMB_SCAFACOS:
    break;
#ifdef P3M
  case COULOMB_ELC_P3M:
    MPI_Bcast(&elc_params, sizeof(ELC_struct), MPI_BYTE, 0, comm_cart);
  // fall through
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    MPI_Bcast(&p3m.params, sizeof(p3m_parameter_struct), MPI_BYTE, 0,
              comm_cart);
    break;
#endif
  case COULOMB_DH:
    MPI_Bcast(&dh_params, sizeof(Debye_hueckel_params), MPI_BYTE, 0, comm_cart);
    break;
  case COULOMB_MMM1D:
  case COULOMB_MMM1D_GPU:
    MPI_Bcast(&mmm1d_params, sizeof(MMM1D_struct), MPI_BYTE, 0, comm_cart);
    break;
  case COULOMB_MMM2D:
    MPI_Bcast(&mmm2d_params, sizeof(MMM2D_struct), MPI_BYTE, 0, comm_cart);
    break;
  case COULOMB_MAGGS:
    MPI_Bcast(&maggs, sizeof(MAGGS_struct), MPI_BYTE, 0, comm_cart);
    break;
  case COULOMB_RF:
  case COULOMB_INTER_RF:
    MPI_Bcast(&rf_params, sizeof(Reaction_field_params), MPI_BYTE, 0,
              comm_cart);
    break;
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: cannot bcast coulomb params for "
                    "unknown method %d\n",
            this_node, coulomb.method);
    errexit();
  }
#endif

#ifdef DIPOLES
  set_dipolar_method_local(coulomb.Dmethod);

  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE:
    break;
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    MPI_Bcast(&dlc_params, sizeof(DLC_struct), MPI_BYTE, 0, comm_cart);
  // fall through
  case DIPOLAR_P3M:
    MPI_Bcast(&dp3m.params, sizeof(p3m_parameter_struct), MPI_BYTE, 0,
              comm_cart);
    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    break;
  case DIPOLAR_MDLC_DS:
  // fall trough
  case DIPOLAR_DS:
    break;
  case DIPOLAR_DS_GPU:
    break;
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
    break;
#endif
  case DIPOLAR_SCAFACOS:
    break;
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: cannot bcast dipolar params for "
                    "unknown method %d\n",
            this_node, coulomb.Dmethod);
    errexit();
  }

#endif

  on_coulomb_change();
#endif
}

/*************** REQ_BCAST_COULOMB ************/
void mpi_bcast_collision_params() {
#ifdef COLLISION_DETECTION
  mpi_call(mpi_bcast_collision_params_slave, 1, 0);
  mpi_bcast_collision_params_slave(-1, 0);
#endif
}

void mpi_bcast_collision_params_slave(int node, int parm) {
#ifdef COLLISION_DETECTION
  MPI_Bcast(&collision_params, sizeof(Collision_parameters), MPI_BYTE, 0,
            comm_cart);

  recalc_forces = 1;
#endif
}

/****************** REQ_SET_PERM ************/

void mpi_send_permittivity_slave(int node, int index) {
#ifdef ELECTROSTATICS
  if (node == this_node) {
    double data[3];
    int indices[3];
    MPI_Recv(data, 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    MPI_Recv(indices, 3, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    for (int d = 0; d < 3; d++) {
      maggs_set_permittivity(indices[0], indices[1], indices[2], d, data[d]);
    }
  }
#endif
}

void mpi_send_permittivity(int node, int index, int *indices,
                           double *permittivity) {
#ifdef ELECTROSTATICS
  if (node == this_node) {
    for (int d = 0; d < 3; d++) {
      maggs_set_permittivity(indices[0], indices[1], indices[2], d,
                             permittivity[d]);
    }
  } else {
    mpi_call(mpi_send_permittivity_slave, node, index);
    MPI_Send(permittivity, 3, MPI_DOUBLE, node, SOME_TAG, comm_cart);
    MPI_Send(indices, 3, MPI_INT, node, SOME_TAG, comm_cart);
  }
#endif
}

/****************** REQ_SET_EXT ************/

void mpi_send_ext_torque(int pnode, int part, int flag, int mask,
                         double torque[3]) {
#ifdef EXTERNAL_FORCES
#ifdef ROTATION
  mpi_call(mpi_send_ext_torque_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /* mask out old flags */
    p->p.ext_flag &= ~mask;
    /* set new values */
    p->p.ext_flag |= flag;

    if (mask & PARTICLE_EXT_TORQUE)
      p->p.ext_torque ={torque[0],torque[1],torque[2]};
  } else {
    int s_buf[2];
    s_buf[0] = flag;
    s_buf[1] = mask;
    MPI_Send(s_buf, 2, MPI_INT, pnode, SOME_TAG, comm_cart);
    if (mask & PARTICLE_EXT_TORQUE)
      MPI_Send(torque, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
#endif
}

void mpi_send_ext_torque_slave(int pnode, int part) {
#ifdef EXTERNAL_FORCES
#ifdef ROTATION
  if (pnode == this_node) {
    int s_buf[2] = {0, 0};
    Particle *p = local_particles[part];
    MPI_Recv(s_buf, 2, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    /* mask out old flags */
    p->p.ext_flag &= ~s_buf[1];
    /* set new values */
    p->p.ext_flag |= s_buf[0];

    if (s_buf[1] & PARTICLE_EXT_TORQUE)
      MPI_Recv(p->p.ext_torque.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
               MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
#endif
}

void mpi_send_ext_force(int pnode, int part, int flag, int mask,
                        double force[3]) {
#ifdef EXTERNAL_FORCES
  mpi_call(mpi_send_ext_force_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /* mask out old flags */
    p->p.ext_flag &= ~mask;
    /* set new values */
    p->p.ext_flag |= flag;
    if (mask & PARTICLE_EXT_FORCE)
      p->p.ext_force ={force[0],force[1],force[2]};
  } else {
    int s_buf[2];
    s_buf[0] = flag;
    s_buf[1] = mask;
    MPI_Send(s_buf, 2, MPI_INT, pnode, SOME_TAG, comm_cart);
    if (mask & PARTICLE_EXT_FORCE)
      MPI_Send(force, 3, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
#endif
}

void mpi_send_ext_force_slave(int pnode, int part) {
#ifdef EXTERNAL_FORCES
  if (pnode == this_node) {
    int s_buf[2] = {0, 0};
    Particle *p = local_particles[part];
    MPI_Recv(s_buf, 2, MPI_INT, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    /* mask out old flags */
    p->p.ext_flag &= ~s_buf[1];
    /* set new values */
    p->p.ext_flag |= s_buf[0];

    if (s_buf[1] & PARTICLE_EXT_FORCE)
      MPI_Recv(p->p.ext_force.data(), 3, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
               MPI_STATUS_IGNORE);
  }

  on_particle_change();
#endif
}

/****************** REQ_RESCALE_PART ************/

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

void mpi_rescale_particles_slave(int pnode, int dir) {
  double scale = 0.0;
  MPI_Recv(&scale, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  local_rescale_particles(dir, scale);
  on_particle_change();
}

/*************** REQ_BCAST_CS *****************/

void mpi_bcast_cell_structure(int cs) {
  mpi_call(mpi_bcast_cell_structure_slave, -1, cs);
  cells_re_init(cs);
}

void mpi_bcast_cell_structure_slave(int pnode, int cs) { cells_re_init(cs); }

/*************** REQ_BCAST_NPTISO_GEOM *****************/

void mpi_bcast_nptiso_geom() {
  mpi_call(mpi_bcast_nptiso_geom_slave, -1, 0);
  mpi_bcast_nptiso_geom_slave(-1, 0);
}

void mpi_bcast_nptiso_geom_slave(int node, int parm) {
  MPI_Bcast(&nptiso.geometry, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&nptiso.dimension, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&nptiso.cubic_box, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&nptiso.non_const_dim, 1, MPI_INT, 0, comm_cart);
}

/***************REQ_UPDATE_MOL_IDS *********************/

void mpi_update_mol_ids() {
  mpi_call(mpi_update_mol_ids_slave, -1, 0);
  mpi_update_mol_ids_slave(-1, 0);
}

void mpi_update_mol_ids_slave(int node, int parm) {
  update_mol_ids_setchains();
}

/******************* REQ_SYNC_TOPO ********************/
int mpi_sync_topo_part_info() {
  int i;
  int molsize = 0;
  int moltype = 0;

  mpi_call(mpi_sync_topo_part_info_slave, -1, 0);
  int n_mols = topology.size();
  MPI_Bcast(&n_mols, 1, MPI_INT, 0, comm_cart);

  for (i = 0; i < n_mols; i++) {
    molsize = topology[i].part.n;
    moltype = topology[i].type;

#ifdef MOLFORCES
    MPI_Bcast(&(topology[i].trap_flag), 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(topology[i].trap_center, 3, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&(topology[i].trap_spring_constant), 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&(topology[i].drag_constant), 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&(topology[i].noforce_flag), 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(&(topology[i].isrelative), 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(&(topology[i].favcounter), 1, MPI_INT, 0, comm_cart);
    if (topology[i].favcounter == -1)
      MPI_Bcast(topology[i].fav, 3, MPI_DOUBLE, 0, comm_cart);
    /* check if any molecules are trapped */
    if ((topology[i].trap_flag != 32) && (topology[i].noforce_flag != 32)) {
      IsTrapped = 1;
    }
#endif

    MPI_Bcast(&molsize, 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(&moltype, 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(topology[i].part.e, topology[i].part.n, MPI_INT, 0, comm_cart);
    MPI_Bcast(&topology[i].type, 1, MPI_INT, 0, comm_cart);
  }

  sync_topo_part_info();

  return 1;
}

void mpi_sync_topo_part_info_slave(int node, int parm) {
  int i;
  int molsize = 0;
  int moltype = 0;
  int n_mols = 0;

  MPI_Bcast(&n_mols, 1, MPI_INT, 0, comm_cart);
  realloc_topology(n_mols);
  for (i = 0; i < n_mols; i++) {

#ifdef MOLFORCES
    MPI_Bcast(&(topology[i].trap_flag), 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(topology[i].trap_center, 3, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&(topology[i].trap_spring_constant), 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&(topology[i].drag_constant), 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&(topology[i].noforce_flag), 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(&(topology[i].isrelative), 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(&(topology[i].favcounter), 1, MPI_INT, 0, comm_cart);
    if (topology[i].favcounter == -1)
      MPI_Bcast(topology[i].fav, 3, MPI_DOUBLE, 0, comm_cart);
    /* check if any molecules are trapped */
    if ((topology[i].trap_flag != 32) && (topology[i].noforce_flag != 32)) {
      IsTrapped = 1;
    }
#endif

    MPI_Bcast(&molsize, 1, MPI_INT, 0, comm_cart);
    MPI_Bcast(&moltype, 1, MPI_INT, 0, comm_cart);
    topology[i].type = moltype;
    topology[i].part.resize(molsize);

    MPI_Bcast(topology[i].part.e, topology[i].part.n, MPI_INT, 0, comm_cart);
    MPI_Bcast(&topology[i].type, 1, MPI_INT, 0, comm_cart);
  }

  sync_topo_part_info();
}

/******************* REQ_BCAST_LBPAR ********************/

void mpi_bcast_lb_params(int field, int value) {
#ifdef LB
  mpi_call(mpi_bcast_lb_params_slave, field, value);
  mpi_bcast_lb_params_slave(field, value);
#endif
}

void mpi_bcast_lb_params_slave(int field, int value) {
#ifdef LB
  MPI_Bcast(&lbpar, sizeof(LB_Parameters), MPI_BYTE, 0, comm_cart);
  on_lb_params_change(field);
#endif
}

/******************* REQ_BCAST_CUDA_GLOBAL_PART_VARS ********************/

void mpi_bcast_cuda_global_part_vars() {
#ifdef CUDA
  mpi_call(mpi_bcast_cuda_global_part_vars_slave, 1,
           0); // third parameter is meaningless
  mpi_bcast_cuda_global_part_vars_slave(-1, 0);
#endif
}

void mpi_bcast_cuda_global_part_vars_slave(int node, int dummy) {
#ifdef CUDA
  MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
            sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart);
  espressoSystemInterface.requestParticleStructGpu();
#endif
}

/********************* REQ_SET_EXCL ********/
void mpi_send_exclusion(int part1, int part2, int _delete) {
#ifdef EXCLUSIONS
  mpi_call(mpi_send_exclusion_slave, part1, part2);

  MPI_Bcast(&_delete, 1, MPI_INT, 0, comm_cart);
  local_change_exclusion(part1, part2, _delete);
  on_particle_change();
#endif
}

void mpi_send_exclusion_slave(int part1, int part2) {
#ifdef EXCLUSIONS
  int _delete = 0;
  MPI_Bcast(&_delete, 1, MPI_INT, 0, comm_cart);
  local_change_exclusion(part1, part2, _delete);
  on_particle_change();
#endif
}

/************** REQ_SET_FLUID **************/
void mpi_send_fluid(int node, int index, double rho, double *j, double *pi) {
#ifdef LB
  if (node == this_node) {
    lb_calc_n_from_rho_j_pi(index, rho, j, pi);
  } else {
    double data[10] = {rho,   j[0],  j[1],  j[2],  pi[0],
                       pi[1], pi[2], pi[3], pi[4], pi[5]};
    mpi_call(mpi_send_fluid_slave, node, index);
    MPI_Send(data, 10, MPI_DOUBLE, node, SOME_TAG, comm_cart);
  }
#endif
}

void mpi_send_fluid_slave(int node, int index) {
#ifdef LB
  if (node == this_node) {
    double data[10];
    MPI_Recv(data, 10, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);

    lb_calc_n_from_rho_j_pi(index, data[0], &data[1], &data[4]);
  }
#endif
}

/************** REQ_GET_FLUID **************/
void mpi_recv_fluid(int node, int index, double *rho, double *j, double *pi) {
#ifdef LB
  if (node == this_node) {
    lb_calc_local_fields(index, rho, j, pi);
  } else {
    double data[10];
    mpi_call(mpi_recv_fluid_slave, node, index);
    MPI_Recv(data, 10, MPI_DOUBLE, node, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
    *rho = data[0];
    j[0] = data[1];
    j[1] = data[2];
    j[2] = data[3];
    pi[0] = data[4];
    pi[1] = data[5];
    pi[2] = data[6];
    pi[3] = data[7];
    pi[4] = data[8];
    pi[5] = data[9];
  }
#endif
}

void mpi_recv_fluid_slave(int node, int index) {
#ifdef LB
  if (node == this_node) {
    double data[10];
    lb_calc_local_fields(index, &data[0], &data[1], &data[4]);
    MPI_Send(data, 10, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
  }
#endif
}

/************** REQ_LB_GET_BOUNDARY_FLAG **************/
void mpi_recv_fluid_boundary_flag(int node, int index, int *boundary) {
#ifdef LB_BOUNDARIES
  if (node == this_node) {
    lb_local_fields_get_boundary_flag(index, boundary);
  } else {
    int data = 0;
    mpi_call(mpi_recv_fluid_boundary_flag_slave, node, index);
    MPI_Recv(&data, 1, MPI_INT, node, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    *boundary = data;
  }
#endif
}

void mpi_recv_fluid_boundary_flag_slave(int node, int index) {
#ifdef LB_BOUNDARIES
  if (node == this_node) {
    int data;
    lb_local_fields_get_boundary_flag(index, &data);
    MPI_Send(&data, 1, MPI_INT, 0, SOME_TAG, comm_cart);
  }
#endif
}

/********************* REQ_ICCP3M_ITERATION ********/
int mpi_iccp3m_iteration(int dummy) {
#ifdef ELECTROSTATICS
  mpi_call(mpi_iccp3m_iteration_slave, -1, 0);

  iccp3m_iteration();

  COMM_TRACE(fprintf(stderr, "%d: iccp3m iteration task %d done.\n", this_node,
                     dummy));

  return check_runtime_errors();
#else
  return 0;
#endif
}

void mpi_iccp3m_iteration_slave(int dummy, int dummy2) {
#ifdef ELECTROSTATICS
  iccp3m_iteration();
  COMM_TRACE(
      fprintf(stderr, "%d: iccp3m iteration task %d done.\n", dummy, dummy2));

  check_runtime_errors();
#endif
}

/********************* REQ_ICCP3M_INIT********/
int mpi_iccp3m_init(int n_induced_charges) {
#ifdef ELECTROSTATICS
  /* nothing has to be done on the master node, this
   * passes only the number of induced charges, in order for
   * slaves to allocate memory */

  mpi_call(mpi_iccp3m_init_slave, -1, n_induced_charges);

  bcast_iccp3m_cfg();

  COMM_TRACE(fprintf(stderr, "%d: iccp3m init task %d done.\n", this_node,
                     n_induced_charges));

  return check_runtime_errors();
#else
  return 0;
#endif
}

void mpi_iccp3m_init_slave(int node, int dummy) {
#ifdef ELECTROSTATICS
  COMM_TRACE(fprintf(stderr, "%d: iccp3m iteration task %d done.\n", this_node,
                     dummy));

  if (iccp3m_initialized == 0) {
    iccp3m_init();
    iccp3m_initialized = 1;
  }

  bcast_iccp3m_cfg();

  check_runtime_errors();
#endif
}

void mpi_recv_fluid_populations(int node, int index, double *pop) {
#ifdef LB
  if (node == this_node) {
    lb_get_populations(index, pop);
  } else {
    mpi_call(mpi_recv_fluid_populations_slave, node, index);
    MPI_Recv(pop, 19 * LB_COMPONENTS, MPI_DOUBLE, node, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
  }
  lbpar.resend_halo = 1;
#endif
}

void mpi_recv_fluid_populations_slave(int node, int index) {
#ifdef LB
  if (node == this_node) {
    double data[19 * LB_COMPONENTS];
    lb_get_populations(index, data);
    MPI_Send(data, 19 * LB_COMPONENTS, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
  }
  lbpar.resend_halo = 1;
#endif
}

void mpi_send_fluid_populations(int node, int index, double *pop) {
#ifdef LB
  if (node == this_node) {
    lb_set_populations(index, pop);
  } else {
    mpi_call(mpi_send_fluid_populations_slave, node, index);
    MPI_Send(pop, 19 * LB_COMPONENTS, MPI_DOUBLE, node, SOME_TAG, comm_cart);
  }
#endif
}

void mpi_send_fluid_populations_slave(int node, int index) {
#ifdef LB
  if (node == this_node) {
    double data[19 * LB_COMPONENTS];
    MPI_Recv(data, 19 * LB_COMPONENTS, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
    lb_set_populations(index, data);
  }
#endif
}

/****************************************************/

void mpi_bcast_max_mu() {
#ifdef DIPOLES
  mpi_call(mpi_bcast_max_mu_slave, -1, 0);

  calc_mu_max();

#endif
}

void mpi_bcast_max_mu_slave(int node, int dummy) {
#ifdef DIPOLES

  calc_mu_max();

#endif
}

#ifdef LANGEVIN_PER_PARTICLE

/******************** REQ_SEND_PARTICLE_T ********************/
void mpi_set_particle_temperature(int pnode, int part, double _T) {
  mpi_call(mpi_set_particle_temperature_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /* here the setting actually happens, if the particle belongs to the local
     * node */
    p->p.T = _T;
  } else {
    MPI_Send(&_T, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
  }

  on_particle_change();
}
#endif

void mpi_set_particle_temperature_slave(int pnode, int part) {
#ifdef LANGEVIN_PER_PARTICLE
  double s_buf = 0.;
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    MPI_Status status;
    MPI_Recv(&s_buf, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, &status);
    /* here the setting happens for nonlocal nodes */
    p->p.T = s_buf;
  }

  on_particle_change();
#endif
}

#ifdef LANGEVIN_PER_PARTICLE
#ifndef PARTICLE_ANISOTROPY
void mpi_set_particle_gamma(int pnode, int part, double gamma) {
#else
void mpi_set_particle_gamma(int pnode, int part, Vector3d gamma) {
#endif
  mpi_call(mpi_set_particle_gamma_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /* here the setting actually happens, if the particle belongs to the local
     * node */
    p->p.gamma = gamma;
  } else {
    comm_cart.send(pnode, SOME_TAG, gamma);
  }

  on_particle_change();
}
#endif

void mpi_set_particle_gamma_slave(int pnode, int part) {
#ifdef LANGEVIN_PER_PARTICLE
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    comm_cart.recv(0, SOME_TAG, p->p.gamma);
  }

  on_particle_change();
#endif
}

#if defined(LANGEVIN_PER_PARTICLE) && defined(ROTATION)
#ifndef PARTICLE_ANISOTROPY
void mpi_set_particle_gamma_rot(int pnode, int part, double gamma_rot)
#else
void mpi_set_particle_gamma_rot(int pnode, int part, Vector3d gamma_rot)
#endif
{
  mpi_call(mpi_set_particle_gamma_rot_slave, pnode, part);

  if (pnode == this_node) {
    Particle *p = local_particles[part];
    /* here the setting actually happens, if the particle belongs to the local
     * node */
    p->p.gamma_rot = gamma_rot;
  } else {
    comm_cart.send(pnode, SOME_TAG, gamma_rot);
  }

  on_particle_change();
}
#endif

void mpi_set_particle_gamma_rot_slave(int pnode, int part) {
#if defined(LANGEVIN_PER_PARTICLE) && defined(ROTATION)
  if (pnode == this_node) {
    Particle *p = local_particles[part];
    comm_cart.recv(0, SOME_TAG, p->p.gamma_rot);
  }

  on_particle_change();
#endif
}

/***** GALILEI TRANSFORM AND ASSOCIATED FUNCTIONS ****/

void mpi_kill_particle_motion(int rotation) {
  mpi_call(mpi_kill_particle_motion_slave, -1, rotation);
  local_kill_particle_motion(rotation);
  on_particle_change();
}

void mpi_kill_particle_motion_slave(int pnode, int rotation) {
  local_kill_particle_motion(rotation);
  on_particle_change();
}

void mpi_kill_particle_forces(int torque) {
  mpi_call(mpi_kill_particle_forces_slave, -1, torque);
  local_kill_particle_forces(torque);
  on_particle_change();
}

void mpi_kill_particle_forces_slave(int pnode, int torque) {
  local_kill_particle_forces(torque);
  on_particle_change();
}

void mpi_system_CMS() {
  int pnode;
  double data[4];
  double rdata[4];
  double *pdata = rdata;

  data[0] = 0.0;
  data[1] = 0.0;
  data[2] = 0.0;
  data[3] = 0.0;

  mpi_call(mpi_system_CMS_slave, -1, 0);

  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      local_system_CMS(pdata);
      data[0] += rdata[0];
      data[1] += rdata[1];
      data[2] += rdata[2];
      data[3] += rdata[3];
    } else {
      MPI_Recv(rdata, 4, MPI_DOUBLE, MPI_ANY_SOURCE, SOME_TAG, comm_cart,
               MPI_STATUS_IGNORE);
      data[0] += rdata[0];
      data[1] += rdata[1];
      data[2] += rdata[2];
      data[3] += rdata[3];
    }
  }

  gal.cms[0] = data[0] / data[3];
  gal.cms[1] = data[1] / data[3];
  gal.cms[2] = data[2] / data[3];
}

void mpi_system_CMS_slave(int node, int index) {
  double rdata[4];
  double *pdata = rdata;
  local_system_CMS(pdata);
  MPI_Send(rdata, 4, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
}

void mpi_system_CMS_velocity() {
  int pnode;
  double data[4];
  double rdata[4];
  double *pdata = rdata;

  data[0] = 0.0;
  data[1] = 0.0;
  data[2] = 0.0;
  data[3] = 0.0;

  mpi_call(mpi_system_CMS_velocity_slave, -1, 0);

  for (pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      local_system_CMS_velocity(pdata);
      data[0] += rdata[0];
      data[1] += rdata[1];
      data[2] += rdata[2];
      data[3] += rdata[3];
    } else {
      MPI_Recv(rdata, 4, MPI_DOUBLE, MPI_ANY_SOURCE, SOME_TAG, comm_cart,
               MPI_STATUS_IGNORE);
      data[0] += rdata[0];
      data[1] += rdata[1];
      data[2] += rdata[2];
      data[3] += rdata[3];
    }
  }

  gal.cms_vel[0] = data[0] / data[3];
  gal.cms_vel[1] = data[1] / data[3];
  gal.cms_vel[2] = data[2] / data[3];
}

void mpi_system_CMS_velocity_slave(int node, int index) {
  double rdata[4];
  double *pdata = rdata;
  local_system_CMS_velocity(pdata);
  MPI_Send(rdata, 4, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
}

void mpi_galilei_transform() {
  double cmsvel[3];

  mpi_system_CMS_velocity();
  memmove(cmsvel, gal.cms_vel, 3 * sizeof(double));

  mpi_call(mpi_galilei_transform_slave, -1, 0);
  MPI_Bcast(cmsvel, 3, MPI_DOUBLE, 0, comm_cart);

  local_galilei_transform(cmsvel);

  on_particle_change();
}

void mpi_galilei_transform_slave(int pnode, int i) {
  double cmsvel[3];
  MPI_Bcast(cmsvel, 3, MPI_DOUBLE, 0, comm_cart);

  local_galilei_transform(cmsvel);
  on_particle_change();
}

/******************** REQ_SWIMMER_REACTIONS ********************/

void mpi_setup_reaction() {
#ifdef SWIMMER_REACTIONS
  mpi_call(mpi_setup_reaction_slave, -1, 0);
  local_setup_reaction();
#endif
}

void mpi_setup_reaction_slave(int pnode, int i) {
#ifdef SWIMMER_REACTIONS
  local_setup_reaction();
#endif
}

/*********************** MAIN LOOP for slaves ****************/

void mpi_loop() {
  if (this_node != 0)
    mpiCallbacks().loop();
}

/*********************** error abort ****************/

void mpi_abort() {
  if (terminated)
    return;

  terminated = 1;
  MPI_Abort(comm_cart, -1);
}

/*********************** other stuff ****************/

#ifdef CUDA
std::vector<EspressoGpuDevice> mpi_gather_cuda_devices() {
  mpi_call(mpi_gather_cuda_devices_slave, 0, 0);
  return cuda_gather_gpus();
}
#endif

void mpi_gather_cuda_devices_slave(int dummy1, int dummy2) {
#ifdef CUDA
  cuda_gather_gpus();
#endif
}

std::vector<int> mpi_resort_particles(int global_flag) {
  mpi_call(mpi_resort_particles_slave, global_flag, 0);
  cells_resort_particles(global_flag);

  std::vector<int> n_parts;
  boost::mpi::gather(comm_cart, cells_get_n_particles(), n_parts, 0);

  return n_parts;
}

void mpi_resort_particles_slave(int global_flag, int) {
  cells_resort_particles(global_flag);

  boost::mpi::gather(comm_cart, cells_get_n_particles(), 0);
}
