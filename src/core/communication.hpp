/*
  Copyright (C) 2010-2018 The ESPResSo project
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
#ifndef _COMMUNICATION_HPP
#define _COMMUNICATION_HPP
/** \file
 *  This file contains the asynchronous MPI communication.
 *
 *  It is the header file for communication.cpp.
 *
 *  The asynchronous MPI communication is used during the script
 *  evaluation. Except for the master node that interprets the Tcl
 *  script, all other nodes wait in mpi_loop() for the master node to
 *  issue an action using mpi_call(). mpi_loop() immediately
 *  executes an MPI_Bcast and therefore waits for the master node to
 *  broadcast a command, which is done by mpi_call(). The request
 *  consists of three integers, the first one describing the action
 *  issued, the second and third an arbitrary parameter depending on
 *  the action issued. If applicable, the second parameter is the node
 *  number of the slave this request is dedicated to.
 *
 *  To add new actions (e.g. to implement new Tcl commands), do the
 *  following:
 *  - write the mpi_* function that is executed on the master
 *  - write the mpi_*_slave function
 *  - Add your slave function to \ref CALLBACK_LIST in communication.cpp
 *
 *  After this your procedure is free to do anything. However, it has
 *  to be in (MPI) sync with what your new mpi_*_slave does. This
 *  procedure is called immediately after the broadcast with the
 *  arbitrary integer as parameter. To this aim it has also to be added
 *  to \ref CALLBACK_LIST. A debug message will be created automatically
 *  in \ref anonymous_namespace{communication.cpp}::names "names".
 */

#include <boost/mpi/communicator.hpp>

#include <array>

#include "MpiCallbacks.hpp"

/** Included needed by callbacks. */
#include "cuda_init.hpp"
#include "grid_based_algorithms/lb_constants.hpp"
#include "particle_data.hpp"

#include "utils/serialization/array.hpp"

/**************************************************
 * exported variables
 **************************************************/

/** \name Exported Variables */
/*@{*/
/** The number of this node. */
extern int this_node;
/** The total number of nodes. */
extern int n_nodes;
// extern MPI_Comm comm_cart;
extern boost::mpi::communicator comm_cart;
/*@}*/

/**
 * Default MPI tag used by callbacks.
 */
#ifndef SOME_TAG
#define SOME_TAG 42
#endif

namespace Communication {
/**
 * @brief Returns a reference to the global callback class instance.
 */
MpiCallbacks &mpiCallbacks();
} // namespace Communication

/**************************************************
 * for every procedure requesting a MPI negotiation
 * a slave exists which processes this request on
 * the slave nodes. It is denoted by *_slave.
 **************************************************/

/**********************************************
 * slave callbacks.
 **********************************************/
typedef void(SlaveCallback)(int node, int param);

/** \name Exported Functions */
/*@{*/
/** Initialize MPI and determine \ref n_nodes and \ref this_node. */
void mpi_init();

/* Call a slave function. */
template <class... Args, class... ArgRef>
void mpi_call(void (*fp)(Args...), ArgRef &&... args) {
  Communication::mpiCallbacks().call(fp, std::forward<ArgRef>(args)...);
}

/** Process requests from master node. Slave nodes main loop. */
void mpi_loop();

/**
 * @brief Replace the MPI communicator by a new one with the given periodicity
 * and node grid.
 */
void mpi_reshape_communicator(std::array<int, 3> const &node_grid,
                              std::array<int, 3> const &periodicity = {
                                  {1, 1, 1}});

/** Issue REQ_PLACE: move particle to a position on a node.
 *  Also calls \ref on_particle_change.
 *  \param id    the particle to move.
 *  \param node  the node to attach it to.
 *  \param pos   the particles position.
 */
void mpi_place_particle(int node, int id, double pos[3]);

/** Issue REQ_PLACE: create particle at a position on a node.
 *  Also calls \ref on_particle_change.
 *  \param id    the particle to create.
 *  \param node  the node to attach it to.
 *  \param pos   the particles position.
 */
void mpi_place_new_particle(int node, int id, double pos[3]);

/** Issue REQ_SET_EXCLUSION: send exclusions.
 *  Also calls \ref on_particle_change.
 *  \param part     identity of first particle of the exclusion.
 *  \param part2    identity of second particle of the exclusion.
 *  \param _delete  if true, do not add the exclusion, rather delete it if
 *                  found
 */
void mpi_send_exclusion(int part, int part2, int _delete);

/** Issue REQ_REM_PART: remove a particle.
 *  Also calls \ref on_particle_change.
 *  \param id    the particle to remove.
 *  \param node  the node it is attached to.
 */
void mpi_remove_particle(int node, int id);

/** Issue REQ_GET_PART: recv particle data. The data has to be freed later
 *  using \ref free_particle, otherwise the dynamically allocated parts, bonds
 *  and exclusions are left over.
 *  \param part  the particle.
 *  \param node  the node it is attached to.
 *  \note Gets a copy of the particle data not a pointer to the actual particle
 *  used in integration
 */
Particle mpi_recv_part(int node, int part);

/** Issue REQ_INTEGRATE: start integrator.
 *  @param n_steps       how many steps to do.
 *  @param reuse_forces  whether to trust the old forces for the first half step
 *  @return nonzero on error
 */
int mpi_integrate(int n_steps, int reuse_forces);

/** Issue REQ_MIN_ENERGY: start energy minimization.
 *  @return nonzero on error
 */
int mpi_minimize_energy();

void mpi_bcast_all_ia_params();

/** Issue REQ_BCAST_IA: send new ia params.
 *  Also calls \ref on_short_range_ia_change.
 *
 *  mpi_bcast_ia_params is used for both, bonded and non-bonded
 *  interaction parameters. Therefore i and j are used depending on
 *  their value:
 *
 *  \param i   particle type for non bonded interaction parameters /
 *             bonded interaction type number.
 *  \param j   if not negative: particle type for non bonded interaction
 *             parameters / if negative: flag for bonded interaction
 */
void mpi_bcast_ia_params(int i, int j);

/** Issue REQ_BCAST_IA_SIZE: send new size of \ref ia_params.
 *  \param s   the new size for \ref ia_params.
 */
void mpi_bcast_max_seen_particle_type(int s);

/** Issue REQ_GATHER: gather data for analysis in analyze.
 *  \todo update parameter descriptions
 *  \param job what to do:
 *      \arg \c 1 calculate and reduce (sum up) energies,
 *           using \ref energy_calc.
 *      \arg \c 2 calculate and reduce (sum up) pressure, stress tensor,
 *           using \ref pressure_calc.
 *      \arg \c 3 calculate and reduce (sum up) instantaneous pressure,
 *           using \ref pressure_calc.
 *      \arg \c 4 use \ref predict_momentum_particles
 *      \arg \c 6 use \ref lb_calc_fluid_momentum
 *      \arg \c 8 use \ref lb_collect_boundary_forces
 *  \param result where to store the gathered value(s):
 *      \arg for \c job=1 unused (the results are stored in a global
 *           energy array of type \ref Observable_stat)
 *      \arg for \c job=2 unused (the results are stored in a global
 *           virials array of type \ref Observable_stat)
 *      \arg for \c job=3 unused (the results are stored in a global
 *           virials array of type \ref Observable_stat)
 *  \param result_t where to store the gathered value(s):
 *      \arg for \c job=1 unused (the results are stored in a global
 *           energy array of type \ref Observable_stat)
 *      \arg for \c job=2 unused (the results are stored in a global
 *           p_tensor tensor of type \ref Observable_stat)
 *      \arg for \c job=3 unused (the results are stored in a global
 *           p_tensor tensor of type \ref Observable_stat)
 *  \param result_nb where to store the gathered value(s):
 *      \arg for \c job=1 unused (the results are stored in a global
 *           energy array of type \ref Observable_stat_non_bonded)
 *      \arg for \c job=2 unused (the results are stored in a global
 *           virials_non_bonded array of type \ref Observable_stat_non_bonded)
 *      \arg for \c job=3 unused (the results are stored in a global
 *           virials_non_bonded array of type \ref Observable_stat_non_bonded)
 *  \param result_t_nb where to store the gathered value(s):
 *      \arg for \c job=1 unused (the results are stored in a global
 *           energy array of type \ref Observable_stat_non_bonded)
 *      \arg for \c job=2 unused (the results are stored in a global
 *           p_tensor_non_bonded tensor of type \ref Observable_stat_non_bonded)
 *      \arg for \c job=3 unused (the results are stored in a global
 *           p_tensor_non_bonded tensor of type \ref Observable_stat_non_bonded)
 */
void mpi_gather_stats(int job, void *result, void *result_t, void *result_nb,
                      void *result_t_nb);

/** Issue REQ_SET_TIME_STEP: send new \ref time_step and rescale the
 *  velocities accordingly.
 */
void mpi_set_time_step(double time_step);

/** Issue REQ_BCAST_COULOMB: send new Coulomb parameters. */
void mpi_bcast_coulomb_params();

/** Issue REQ_RESCALE_PART: rescales all particle positions in direction 'dir'
 *  by a factor 'scale'.
 */
void mpi_rescale_particles(int dir, double scale);

/** Issue REQ_BCAST_CS: change the cell structure on all nodes. */
void mpi_bcast_cell_structure(int cs);

/** Issue REQ_BCAST_NPTISO_GEOM: broadcast nptiso geometry parameter to all
 *  nodes.
 */
void mpi_bcast_nptiso_geom();

/** Issue REQ_UPDATE_MOL_IDS: Update the molecule ids so that they are
 *  in sync with the topology. Note that this only makes sense if you
 *  have a simple topology such that each particle can only belong to
 *  a single molecule */
void mpi_update_mol_ids();

void mpi_bcast_lb_particle_coupling();

Utils::Vector3d mpi_recv_lb_interpolated_velocity(int node,
                                                  Utils::Vector3d const &pos);

/** Issue REQ_BCAST_cuda_global_part_vars: Broadcast a parameter for CUDA */
void mpi_bcast_cuda_global_part_vars();

/** Issue REQ_ICCP3M_ITERATION: performs iccp3m iteration.
 *  @return nonzero on error
 */
int mpi_iccp3m_iteration();

/** Issue REQ_ICCP3M_INIT: performs iccp3m initialization
 *  @return nonzero on error
 */
int mpi_iccp3m_init();

/** Part of MDLC */
void mpi_bcast_max_mu();

/** Galilei and other: set all particle velocities and rotational inertias to
 *                     zero.
 *                     set all forces and torques on the particles to zero
 *                     calculate the centre of mass (CMS)
 *                     calculate the velocity of the CMS
 *                     remove the CMS velocity from the system
 */
void mpi_kill_particle_motion(int rotation);
void mpi_kill_particle_forces(int torque);
void mpi_system_CMS();
void mpi_system_CMS_velocity();
void mpi_galilei_transform();
void mpi_observable_lb_radial_velocity_profile();

/** Issue REQ_SWIMMER_REACTIONS: notify the system of changes to the reaction
 *  parameters
 */
void mpi_setup_reaction();

#ifdef CUDA
/** Gather CUDA devices from all nodes */
std::vector<EspressoGpuDevice> mpi_gather_cuda_devices();
#endif

/**
 * @brief Resort the particles.
 *
 * This function resorts the particles on the nodes.
 *
 * @param global_flag If true a global resort is done,
 *        if false particles are only exchanges between neighbors.
 * @return The number of particles on the nodes after the resort.
 */
std::vector<int> mpi_resort_particles(int global_flag);

/*@}*/

#endif
