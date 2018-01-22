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
#ifndef _COMMUNICATION_HPP
#define _COMMUNICATION_HPP
/** \file communication.hpp
    This file contains the asynchronous MPI communication.

    It is the header file for \ref communication.cpp "communication.c".

    The asynchronous MPI communication is used during the script
    evaluation. Except for the master node that interpretes the Tcl
    script, all other nodes wait in mpi_loop() for the master node to
    issue an action using mpi_call(). \ref mpi_loop immediately
    executes an MPI_Bcast and therefore waits for the master node to
    broadcast a command, which is done by mpi_call(). The request
    consists of three integers, the first one describing the action
    issued, the second and third an arbitrary parameter depending on
    the action issued. If applicable, the second parameter is the node
    number of the slave this request is dedicated to.

    To add new actions (e. g. to implement new Tcl commands), do the
    following:
    - write the mpi_* function that is executed on the master
    - write the mpi_*_slave function
    - Add your slave function to CALLBACK_LIST in communication.cpp

    After this your procedure is free to do anything. However, it has
    to be in (MPI) sync with what your new mpi_*_slave does.  This
    procedure is called immediately after the broadcast with the
    arbitrary integer as parameter.  To this aim it has also to be
    added to \ref CALLBACK_LIST "callbacks".  Last but not least for
    debugging purposes you can add a nice name to \ref #names in the
    same way.  */

#include <array>
#include <mpi.h>

#include <boost/mpi/communicator.hpp>

#include "MpiCallbacks.hpp"

/** Included needed by callbacks. */
#include "cuda_init.hpp"
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
 *
 */
MpiCallbacks &mpiCallbacks();
}

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
void mpi_call(SlaveCallback cb, int node, int param);

/** Process requests from master node. Slave nodes main loop. */
void mpi_loop();

/** Abort Espresso using MPI_Abort. */
void mpi_abort();

/** Finalize MPI. Called by all nodes upon exit */
void mpi_finalize();

/**
 * @brief Replace the MPI communicator by a new one with the given periodicity
 * and node grid.
 */
void mpi_reshape_communicator(std::array<int, 3> const &node_grid,
                              std::array<int, 3> const &periodicity = {{1, 1,
                                                                       1}});

/** Issue REQ_EVENT: tells all clients of some system change.
    The events are:
    <ul>
    <li> PARTICLE_CHANGED
    <li> INTERACTION_CHANGED
    </ul>
    Then all nodes execute the respective on_* procedure from \ref
   initialize.cpp
    Note that not all of these codes are used. Since some actions (like placing
   a
    particle) include communication anyways, this is handled by the way.
*/
void mpi_bcast_event(int event);

/** Issue REQ_PLACE: move particle to a position on a node.
    Also calls \ref on_particle_change.
    \param id   the particle to move.
    \param node the node to attach it to.
    \param pos  the particles position.
*/
void mpi_place_particle(int node, int id, double pos[3]);

/** Issue REQ_PLACE: create particle at a position on a node.
    Also calls \ref on_particle_change.
    \param id   the particle to create.
    \param node the node to attach it to.
    \param pos  the particles position.
*/
void mpi_place_new_particle(int node, int id, double pos[3]);

/** Issue REQ_SET_V: send particle velocity.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param v its new velocity.
*/
void mpi_send_v(int node, int part, double v[3]);

/** Issue REQ_SET_SWIMMING: send particle swimming properties.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param swim struct containing swimming parameters
*/
void mpi_send_swimming(int node, int part, ParticleParametersSwimming swim);

/** Issue REQ_SET_F: send particle force.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param F its new force.
*/
void mpi_send_f(int node, int part, double F[3]);

/** issue req_set_solv: send particle solvation free energy
    also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param solvation its new solvation free energy.
*/
void mpi_send_solvation(int node, int part, double *solvation);

/** Issue REQ_SET_M: send particle mass.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param mass its new mass.
*/
void mpi_send_mass(int node, int part, double mass);

/** Issue REQ_SET_Q: send particle charge.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param q its new charge.
*/
void mpi_send_q(int node, int part, double q);

/** Issue REQ_SET_MU_E: send particle electrophoretic mobility.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param mu_E its new mobility.
*/
void mpi_send_mu_E(int node, int part, double mu_E[3]);

#ifdef ROTATIONAL_INERTIA
/** Issue REQ_SET_ROTATIONAL_INERTIA: send particle rotational inertia.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param rinertia its new rotational inertia.
*/
void mpi_send_rotational_inertia(int node, int part, double rinertia[3]);
#endif

#ifdef AFFINITY
/** Issue REQ_SET_AFFINITY: send particle affinity.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param bond_site its new site of the affinity bond.
*/
void mpi_send_affinity(int node, int part, double bond_site[3]);
#endif

#ifdef MEMBRANE_COLLISION
/** Issue REQ_SET_MEMBRANE_COLLISION: send outward direction of the particle.
 Also calls \ref on_particle_change.
 \param part the particle.
 \param node the node it is attached to.
 \param out_direction its new outward direction.
 */
void mpi_send_out_direction(int node, int part, double out_direction[3]);
#endif

#ifdef ROTATION
/** Issue REQ_SET_QUAT: send particle orientation.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param quat its new quaternions.
*/
void mpi_send_quat(int node, int part, double quat[4]);

/** Issue REQ_SET_ROTATION: send particle rotation flag
    Also calls \ref on_particle_change.
    \param part the particle.
    \param pnode the node it is attached to.
    \param rot the rotation flag
*/
void mpi_send_rotation(int pnode, int part, short int rot);

/* Issue REQ_SET_LAMBDA: send particle angular velocity.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param omega its new angular velocity.
*/
void mpi_send_omega(int node, int part, double omega[3]);

/** Issue REQ_SET_TORQUE: send particle torque.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param torque its new torque.
*/
void mpi_send_torque(int node, int part, double torque[3]);
#endif

#ifdef DIPOLES
/** Issue REQ_SET_DIP: send particle dipole orientation.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param dip its new dipole orientation.
*/
void mpi_send_dip(int node, int part, double dip[3]);
/** Issue REQ_SET_DIPM: send particle dipole moment.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param dipm its new dipole moment (absolut value).
*/
void mpi_send_dipm(int node, int part, double dipm);
#endif

#ifdef VIRTUAL_SITES
/** Issue REQ_SET_DIPM: send particle dipole moment.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param isVirtual its new isVirtual.
*/
void mpi_send_virtual(int node, int part, int isVirtual);
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void mpi_send_vs_relative(int node, int part, int vs_relative_to,
                          double vs_distance, double *rel_ori);
#endif

#ifdef MULTI_TIMESTEP
/** Issue REQ_SET_SMALLER_TIMESTEP: send smaller time step value.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param smaller_timestep its new smaller_timestep.
*/
void mpi_send_smaller_timestep_flag(int node, int part,
                                    int smaller_timestep_flag);
#endif

/** Issue REQ_SET_TYPE: send particle type.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param type its new type.
*/
void mpi_send_type(int node, int part, int type);

/** Issue REQ_SET_MOL_ID: send molecule id.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param mid its new mol_id.
*/
void mpi_send_mol_id(int node, int part, int mid);

/** Issue REQ_SET_BOND: send bond.
    Also calls \ref on_particle_change.
    \param pnode    node it is attached to.
    \param part     identity of principal atom of the bond.
    \param bond     field containing the bond type number and the identity of
   all bond partners (secundary atoms of the bond).
    \param _delete   if true, do not add the bond, rather delete it if found
    \return 1 on success or 0 if not (e. g. bond to delete does not exist)
*/
int mpi_send_bond(int pnode, int part, int *bond, int _delete);

/** Issue REQ_SET_EXCLUSION: send exclusions.
    Also calls \ref on_particle_change.
    \param part     identity of first particle of the exclusion.
    \param part2    identity of secnd particle of the exclusion.
    \param _delete   if true, do not add the exclusion, rather delete it if
   found
*/
void mpi_send_exclusion(int part, int part2, int _delete);

/** Issue REQ_REM_PART: remove a particle.
    Also calls \ref on_particle_change.
    \param id   the particle to remove.
    \param node the node it is attached to.
*/
void mpi_remove_particle(int node, int id);

/** Issue REQ_GET_PART: recv particle data. The data has to be freed later
    using \ref free_particle, otherwise the dynamically allocated parts, bonds
    and exclusions are left over.
    \param part the particle.
    \param node the node it is attached to.
    \param part_data where to store the received data.
    \note Gets a copy of the particle data not a pointer to the actual particle
    used in integration
*/
Particle mpi_recv_part(int node, int part);

/** Issue REQ_INTEGRATE: start integrator.
    @param n_steps how many steps to do.
    @param reuse_forces whether to trust the old forces for the first half step
    @return nonzero on error
*/
int mpi_integrate(int n_steps, int reuse_forces);

/** Issue REQ_MIN_ENERGY: start energy minimization.
    @return nonzero on error
 */
int mpi_minimize_energy(void);

/** Issue REQ_BCAST_IA: send new ia params.
    Also calls \ref on_short_range_ia_change.

    mpi_bcast_ia_params is used for both, bonded and non-bonded
    interaction parameters. Therefor i and j are used depending on
    their value:

    \param i   particle type for non bonded interaction parameters /
               bonded interaction type number.
    \param j   if not negative: particle type for non bonded interaction
   parameters /
               if negative: flag for bonded interaction */
void mpi_bcast_ia_params(int i, int j);

/** Issue REQ_BCAST_IA_SIZE: send new size of \ref ia_params.
    \param s the new size for \ref ia_params.
*/
void mpi_bcast_n_particle_types(int s);

/** Issue REQ_GATHER: gather data for analysis in analyze.
    \param job what to do:
    <ul>
        <li> 1 calculate and reduce (sum up) energies, using \ref energy_calc.
        <li> 2 calculate and reduce (sum up) pressure, stress tensor, using \ref
   pressure_calc.
        <li> 3 calculate and reduce (sum up) instantaneous pressure, using \ref
   pressure_calc.
    </ul>
    \param result where to store the gathered value(s):
    <ul><li> job=1 unused (the results are stored in a global
             energy array of type \ref Observable_stat)
        <li> job=2 unused (the results are stored in a global
             virials array of type \ref Observable_stat)
        <li> job=3 unused (the results are stored in a global
             virials array of type \ref Observable_stat)
    \param result_t where to store the gathered value(s):
    <ul><li> job=1 unused (the results are stored in a global
             energy array of type \ref Observable_stat)
        <li> job=2 unused (the results are stored in a global
             p_tensor tensor of type \ref Observable_stat)
        <li> job=3 unused (the results are stored in a global
             p_tensor tensor of type \ref Observable_stat)
    \param result_nb where to store the gathered value(s):
    <ul><li> job=1 unused (the results are stored in a global
             energy array of type \ref Observable_stat_non_bonded)
        <li> job=2 unused (the results are stored in a global
             virials_non_bonded array of type \ref Observable_stat_non_bonded)
        <li> job=3 unused (the results are stored in a global
             virials_non_bonded array of type \ref Observable_stat_non_bonded)
    \param result_t_nb where to store the gathered value(s):
    <ul><li> job=1 unused (the results are stored in a global
             energy array of type \ref Observable_stat_non_bonded)
        <li> job=2 unused (the results are stored in a global
             p_tensor_non_bonded tensor of type \ref Observable_stat_non_bonded)
        <li> job=3 unused (the results are stored in a global
             p_tensor_non_bonded tensor of type \ref Observable_stat_non_bonded)
    </ul>
*/
void mpi_gather_stats(int job, void *result, void *result_t, void *result_nb,
                      void *result_t_nb);

/** Issue GET_LOCAL_STRESS_TENSOR: gather the contribution to the local stress
   tensors from
    each node.
 */

void mpi_local_stress_tensor(DoubleList *TensorInBin, int bins[3],
                             int periodic[3], double range_start[3],
                             double range[3]);

/** Issue REQ_SET_TIME_STEP: send new \ref time_step and rescale the
    velocities accordingly.
*/
void mpi_set_time_step(double time_step);

#ifdef MULTI_TIMESTEP
/** Issue REQ_SET_SMALLER_TIME_STEP: send new \ref smaller_time_step.
    Requires MULTI_TIMESTEP feature. */
void mpi_set_smaller_time_step(double smaller_time_step);
#endif

/** Issue REQ_BCAST_COULOMB: send new coulomb parameters. */
void mpi_bcast_coulomb_params();

/** send new collision parameters. */
void mpi_bcast_collision_params();

/** Issue REQ_SEND_EXT_FORCE: send nex external flag and external force. */
void mpi_send_ext_force(int pnode, int part, int flag, int mask,
                        double force[3]);

/** Issue REQ_SEND_EXT_TORQUE: send nex external flag and external torque. */
void mpi_send_ext_torque(int pnode, int part, int flag, int mask,
                         double torque[3]);

#ifdef LANGEVIN_PER_PARTICLE
/** Issue REQ_SEND_PARTICLE_T: send particle type specific temperature. */
void mpi_set_particle_temperature(int pnode, int part, double _T);

/** Issue REQ_SEND_PARTICLE_T: send particle type specific frictional
 * coefficient. */
#ifndef PARTICLE_ANISOTROPY
void mpi_set_particle_gamma(int pnode, int part, double gamma);
#else
void mpi_set_particle_gamma(int pnode, int part, Vector3d gamma);
#endif

#ifdef ROTATION
#ifndef PARTICLE_ANISOTROPY
void mpi_set_particle_gamma_rot(int pnode, int part, double gamma_rot);
#else
void mpi_set_particle_gamma_rot(int pnode, int part, Vector3d gamma_rot);
#endif // PARTICLE_ANISOTROPY
#endif
#endif // LANGEVIN_PER_PARTICLE

#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
/** Issue REQ_LB_BOUNDARY: set up walls for lb fluid */
void mpi_bcast_lbboundary(int del_num);
#endif

/** Issue REQ_RESCALE_PART: rescales all particle positions in direction 'dir'
 * by a factor 'scale'. */
void mpi_rescale_particles(int dir, double scale);

/** Issue REQ_BCAST_CS: change the cell structure on all nodes. */
void mpi_bcast_cell_structure(int cs);

/** Issue REQ_BCAST_NPTISO_GEOM: broadcast nptiso geometry parameter to all
 * nodes. */
void mpi_bcast_nptiso_geom(void);

/** Issue REQ_UPDATE_MOL_IDS: Update the molecule ids so that they are
    in sync with the topology.  Note that this only makes sense if you
    have a simple topology such that each particle can only belong to
    a single molecule */
void mpi_update_mol_ids(void);

/** Issue REQ_SYNC_TOPO: Update the molecules ids to that they correspond to the
 * topology */
int mpi_sync_topo_part_info(void);

/** Issue REQ_BCAST_LBPAR: Broadcast a parameter for Lattice Boltzmann.
 * @param field References the parameter field to be broadcasted. The references
 * are defined in \ref lb.hpp "lb.hpp"
 */
void mpi_bcast_lb_params(int field, int value = -1);

/** Issue REQ_BCAST_cuda_global_part_vars: Broadcast a parameter for CUDA
 */
void mpi_bcast_cuda_global_part_vars();

/** Issue REQ_SEND_FLUID: Send a single lattice site to a processor.
 * @param node  processor to send to
 * @param index index of the lattice site
 * @param rho   local fluid density
 * @param j     local fluid velocity
 * @param pi    local fluid pressure
 */
void mpi_send_fluid(int node, int index, double rho, double *j, double *pi);

/** Issue REQ_GET_FLUID: Receive a single lattice site from a processor.
 * @param node  processor to send to
 * @param index index of the lattice site
 * @param rho   local fluid density
 * @param j     local fluid velocity
 * @param pi    local fluid pressure
 */
void mpi_recv_fluid(int node, int index, double *rho, double *j, double *pi);

/** Issue REQ_LB_GET_BOUNDARY_FLAG: Receive a single lattice sites boundary flag
 * from a processor.
 * @param node     processor to send to
 * @param index    index of the lattice site
 * @param boundary local boundary flag
 */
void mpi_recv_fluid_boundary_flag(int node, int index, int *boundary);

/** Issue REQ_ICCP3M_ITERATION: performs iccp3m iteration.
    @return nonzero on error
*/
int mpi_iccp3m_iteration(int dummy);

/** Issue REQ_ICCP3M_INIT: performs iccp3m initialization
    @return nonzero on error
*/
int mpi_iccp3m_init(int dummy);

/** Issue REQ_RECV_FLUID_POPULATIONS: Send a single lattice site to a processor.
 * @param node  processor to send to
 * @param index index of the lattice site
 * @param pop   local fluid population
 */
void mpi_recv_fluid_populations(int node, int index, double *pop);

/** Issue REQ_SEND_FLUID_POPULATIONS: Send a single lattice site to a processor.
 * @param node  processor to send to
 * @param index index of the lattice site
 * @param pop   local fluid population
 */
void mpi_send_fluid_populations(int node, int index, double *pop);

/** Part of MDLC
 */
void mpi_bcast_max_mu();

/** Galilei and other: set all particle velocities and rotational inertias to
   zero.
                       set all forces and torques on the particles to zero
                       calculate the centre of mass (CMS)
                       calculate the velocity of the CMS
                       remove the CMS velocity from the system
 */
void mpi_kill_particle_motion(int rotation);
void mpi_kill_particle_forces(int torque);
void mpi_system_CMS();
void mpi_system_CMS_velocity();
void mpi_galilei_transform();
void mpi_observable_lb_radial_velocity_profile();

/** Issue REQ_CATALYTIC_REACTIONS: notify the system of changes to the reaction
 * parameters
 */
void mpi_setup_reaction();

void mpi_external_potential_broadcast(int number);
void mpi_external_potential_broadcast_slave(int node, int number);
void mpi_external_potential_tabulated_read_potential_file(int number);
void mpi_external_potential_sum_energies();
void mpi_external_potential_sum_energies_slave();

#ifdef CUDA
/** Gather CUDA devices from all nodes */
std::vector<EspressoGpuDevice> mpi_gather_cuda_devices();
#endif

/** MPI-IO output function.
 *  \param filename Filename prefix for the created files. Must be
 * null-terminated.
 *  \param fields Fields to dump.
 *  \param write 1 to write, 0 to read
 */
void mpi_mpiio(const char *filename, unsigned fields, int write);

/**
 * @brief Resort the particles.
 *
 * This function resorts the particles on the nodes.
 *
 * @param global_flag If true a global resort is done,
 *        if false particles are only exchanges between
 *        neighbors.
 * @return The number of particles on the nodes after
 *         the resort.
 */
std::vector<int> mpi_resort_particles(int global_flag);

/*@}*/

/** \name Event codes for \ref mpi_bcast_event
    These codes are used by \ref mpi_bcast_event to notify certain changes
    of doing something now.
*/
/*@{*/
#define P3M_COUNT_CHARGES 0
#define CHECK_PARTICLES 2
#define MAGGS_COUNT_CHARGES 3
#define P3M_COUNT_DIPOLES 5
/*@}*/

#endif
