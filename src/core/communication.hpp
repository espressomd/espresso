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
#ifndef _COMMUNICATION_HPP
#define _COMMUNICATION_HPP
/** \file
 *  This file contains the asynchronous MPI communication.
 *
 *  It is the header file for communication.cpp.
 *
 *  The asynchronous MPI communication is used during the script
 *  evaluation. Except for the master node that interprets the interface
 *  script, all other nodes wait in mpi_loop() for the master node to
 *  issue an action using mpi_call(). mpi_loop() immediately
 *  executes an MPI_Bcast and therefore waits for the master node to
 *  broadcast a command, which is done by mpi_call(). The request
 *  consists of a callback function and two arbitrary integers. If
 *  applicable, the first integer is the node number of the slave
 *  this request is dedicated to.
 *
 *  To add new actions (e.g. to implement new interface functionality), do the
 *  following:
 *  - write the @c mpi_* function that is executed on the master
 *  - write the @c mpi_*_slave function
 *  - in communication.cpp add your slave function to \ref CALLBACK_LIST or
 *    register it with one of the @c REGISTER_CALLBACK macros
 *
 *  After this, your procedure is free to do anything. However, it has
 *  to be in (MPI) sync with what your new @c mpi_*_slave does. This
 *  procedure is called immediately after the broadcast with the
 *  arbitrary integer as parameter. To this aim it has also to be added
 *  to \ref CALLBACK_LIST.
 */

#include "MpiCallbacks.hpp"

/* Includes needed by callbacks. */
#include "Particle.hpp"
#include "cuda_init.hpp"
#include "grid_based_algorithms/lb_constants.hpp"

#include <boost/mpi/communicator.hpp>
#include <vector>

/** The number of this node. */
extern int this_node;
/** The total number of nodes. */
extern int n_nodes;
/** The communicator */
extern boost::mpi::communicator comm_cart;
/** Statistics to calculate */
enum class GatherStats : int {
  energy,
  pressure,
  lb_fluid_momentum,
  lb_boundary_forces
};

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

/** Initialize MPI and determine \ref n_nodes and \ref this_node. */
void mpi_init();

/** @brief Call a slave function.
 *  @tparam Args   Slave function argument types
 *  @tparam ArgRef Slave function argument types
 *  @param fp      Slave function
 *  @param args    Slave function arguments
 */
template <class... Args, class... ArgRef>
void mpi_call(void (*fp)(Args...), ArgRef &&... args) {
  Communication::mpiCallbacks().call(fp, std::forward<ArgRef>(args)...);
}

/** @brief Call a slave function.
 *  @tparam Args   Slave function argument types
 *  @tparam ArgRef Slave function argument types
 *  @param fp      Slave function
 *  @param args    Slave function arguments
 */
template <class... Args, class... ArgRef>
void mpi_call_all(void (*fp)(Args...), ArgRef &&... args) {
  Communication::mpiCallbacks().call_all(fp, std::forward<ArgRef>(args)...);
}

/** @brief Call a slave function.
 *  @tparam Tag    Any tag type defined in @ref Communication::Result
 *  @tparam R      Return type of the slave function
 *  @tparam Args   Slave function argument types
 *  @tparam ArgRef Slave function argument types
 *  @param tag     Reduction strategy
 *  @param fp      Slave function
 *  @param args    Slave function arguments
 */
template <class Tag, class R, class... Args, class... ArgRef>
auto mpi_call(Tag tag, R (*fp)(Args...), ArgRef &&... args) {
  return Communication::mpiCallbacks().call(tag, fp,
                                            std::forward<ArgRef>(args)...);
}

/** @brief Call a slave function.
 *  @tparam Tag    Any tag type defined in @ref Communication::Result
 *  @tparam TagArg Types of arguments to @p Tag
 *  @tparam R      Return type of the slave function
 *  @tparam Args   Slave function argument types
 *  @tparam ArgRef Slave function argument types
 *  @param tag     Reduction strategy
 *  @param tag_arg Arguments to the reduction strategy
 *  @param fp      Slave function
 *  @param args    Slave function arguments
 */
template <class Tag, class TagArg, class R, class... Args, class... ArgRef>
auto mpi_call(Tag tag, TagArg &&tag_arg, R (*fp)(Args...), ArgRef &&... args) {
  return Communication::mpiCallbacks().call(tag, std::forward<TagArg>(tag_arg),
                                            fp, std::forward<ArgRef>(args)...);
}

/** Process requests from master node. Slave nodes main loop. */
void mpi_loop();

/** Move particle to a position on a node.
 *  Also calls \ref on_particle_change.
 *  \param id    the particle to move.
 *  \param node  the node to attach it to.
 *  \param pos   the particles position.
 */
void mpi_place_particle(int node, int id, const Utils::Vector3d &pos);

/** Create particle at a position on a node.
 *  Also calls \ref on_particle_change.
 *  \param id    the particle to create.
 *  \param pos   the particles position.
 */
int mpi_place_new_particle(int id, const Utils::Vector3d &pos);

/** Send exclusions.
 *  Also calls \ref on_particle_change.
 *  \param part     identity of first particle of the exclusion.
 *  \param part2    identity of second particle of the exclusion.
 *  \param _delete  if true, do not add the exclusion, rather delete it if
 *                  found
 */
void mpi_send_exclusion(int part, int part2, int _delete);

/** Remove a particle.
 *  Also calls \ref on_particle_change.
 *  \param id    the particle to remove.
 *  \param node  the node it is attached to.
 */
void mpi_remove_particle(int node, int id);

/** Start integrator.
 *  @param n_steps       how many steps to do.
 *  @param reuse_forces  whether to trust the old forces for the first half step
 *  @return nonzero on error
 */
int mpi_integrate(int n_steps, int reuse_forces);

/** Start steepest descent. */
int mpi_steepest_descent(int steps);

void mpi_bcast_all_ia_params();

/** Send new IA params.
 *  Also calls \ref on_short_range_ia_change.
 *
 *  Used for both bonded and non-bonded interaction parameters. Therefore
 *  @p i and @p j are used depending on their value:
 *
 *  \param i   particle type for non bonded-interaction parameters /
 *             bonded interaction type number.
 *  \param j   if not negative: particle type for non-bonded interaction
 *             parameters / if negative: flag for bonded interaction
 */
void mpi_bcast_ia_params(int i, int j);

/** Send new size of \ref ia_params.
 *  \param s   the new size for \ref ia_params.
 */
void mpi_bcast_max_seen_particle_type(int s);

/** Gather data for analysis.
 *  \param[in] job what to do:
 *      \arg for \ref GatherStats::energy, calculate and reduce (sum up)
 *           energies, using \ref energy_calc.
 *      \arg for \ref GatherStats::pressure, calculate and reduce (sum up)
 *           pressure, using \ref pressure_calc.
 *      \arg for \ref GatherStats::lb_fluid_momentum, use
 *           \ref lb_calc_fluid_momentum.
 *      \arg for \ref GatherStats::lb_boundary_forces, use
 *           \ref lb_collect_boundary_forces.
 *  \param[out] result where to store values gathered by
 *      \ref GatherStats::lb_fluid_momentum,
 *      \ref GatherStats::lb_boundary_forces
 */
void mpi_gather_stats(GatherStats job, double *result = nullptr);

/** Send new \ref time_step and rescale the velocities accordingly. */
void mpi_set_time_step(double time_step);

/** Send new Coulomb parameters. */
void mpi_bcast_coulomb_params();

/** Rescale all particle positions in direction @p dir by a factor @p scale. */
void mpi_rescale_particles(int dir, double scale);

/** Change the cell structure on all nodes. */
void mpi_bcast_cell_structure(int cs);

void mpi_set_use_verlet_lists(bool use_verlet_lists);

/** Broadcast nptiso geometry parameter to all nodes. */
void mpi_bcast_nptiso_geom();

/** Broadcast @ref CUDA_global_part_vars structure */
void mpi_bcast_cuda_global_part_vars();

/** Perform iccp3m initialization.
 *  @return nonzero on error
 */
int mpi_iccp3m_init();

/** Calculate the maximal dipole moment in the system (part of MDLC) */
void mpi_bcast_max_mu();

/** @name Galilei and other
 *  - set all particle velocities and rotational inertias to zero
 *  - set all forces and torques on the particles to zero
 *  - calculate the centre of mass (CMS)
 *  - calculate the velocity of the CMS
 *  - remove the CMS velocity from the system
 */
/*@{*/
void mpi_kill_particle_motion(int rotation);
void mpi_kill_particle_forces(int torque);
Utils::Vector3d mpi_system_CMS();
Utils::Vector3d mpi_system_CMS_velocity();
void mpi_galilei_transform();
/*@}*/

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

#endif
