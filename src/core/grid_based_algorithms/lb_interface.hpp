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
#ifndef CORE_LB_INTERFACE
#define CORE_LB_INTERFACE

#include "config.hpp"

#include <utils/Vector.hpp>

#include <cstdint>
#include <string>
#include <vector>

/** @brief LB implementation currently active. */
enum class ActiveLB : int { NONE, WALBERLA };

/** @brief Switch determining the type of lattice dynamics. */
extern ActiveLB lattice_switch;

/**
 * @brief Propagate the LB fluid.
 */
void lb_lbfluid_propagate();

/**
 * @brief Perform a full initialization of the lattice-Boltzmann system.
 * All derived parameters and the fluid are reset to their default values.
 */
void lb_lbfluid_init();

/**
 * @brief Get the current counter of the Philox RNG.
 */
uint64_t lb_lbfluid_get_rng_state();

/**
 * @brief Set the current counter of the Philox RNG.
 */
void lb_lbfluid_set_rng_state(uint64_t counter);

/**
 * @brief Get the global variable @ref lattice_switch.
 */
ActiveLB lb_lbfluid_get_lattice_switch();

/**
 * @brief Set the global variable @ref lattice_switch.
 */
void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch);

/**
 * @brief Check if tau is an integer multiple of time_step, throws if not
 */
void check_tau_time_step_consistency(double tau, double time_s);

/**
 * @brief Set the external force density acting on the LB fluid.
 */
void lb_lbfluid_set_ext_force_density(const Utils::Vector3d &force_density);

/**
 * @brief Perform LB parameter and boundary velocity checks.
 */
void lb_lbfluid_sanity_checks(double time_step);

/**
 * @brief Set the LB density for a single node.
 */
void lb_lbnode_set_density(const Utils::Vector3i &ind, double density);

/**
 * @brief Set the LB fluid velocity for a single node.
 */
void lb_lbnode_set_velocity(const Utils::Vector3i &ind,
                            const Utils::Vector3d &u);

/**
 * @brief Set the LB fluid velocity for a single boundary node.
 */
void lb_lbnode_set_velocity_at_boundary(const Utils::Vector3i &ind,
                                        const Utils::Vector3d &u);

/**
 * @brief Set the LB fluid populations for a single node.
 */
void lb_lbnode_set_pop(const Utils::Vector3i &ind,
                       const std::vector<double> &pop);

/**
 * @brief Set force applied on an lb node during the previous integration step
 */
void lb_lbnode_set_last_applied_force(const Utils::Vector3i &ind,
                                      const Utils::Vector3d &force);

/**
 * @brief Get the LB time step.
 */
double lb_lbfluid_get_tau();

/**
 * @brief Get the LB grid spacing.
 */
double lb_lbfluid_get_agrid();

/**
 * @brief Get the global LB bulk viscosity.
 */
double lb_lbfluid_get_bulk_viscosity();

/**
 * @brief Get the global LB viscosity.
 */
double lb_lbfluid_get_viscosity();

/**
 * @brief Get the external force density acting on the LB fluid.
 */
const Utils::Vector3d lb_lbfluid_get_ext_force_density();

/**
 * @brief Get the thermal energy parameter of the LB fluid.
 */
double lb_lbfluid_get_kT();

/**
 * @brief Get the lattice speed (agrid/tau).
 */
double lb_lbfluid_get_lattice_speed();

/**
 * @brief Create a VTK observable.
 */
void lb_lbfluid_create_vtk(unsigned delta_N, unsigned initial_count,
                           unsigned flag_observables,
                           std::string const &identifier,
                           std::string const &base_folder,
                           std::string const &prefix);

/**
 * @brief Write a VTK observable to disk.
 */
void lb_lbfluid_write_vtk(std::string const &vtk_uid);

/**
 * @brief Toggle a VTK observable on/off.
 */
void lb_lbfluid_switch_vtk(std::string const &vtk_uid, int status);

/**
 * @brief Get the LB fluid density for a single node.
 */
double lb_lbnode_get_density(const Utils::Vector3i &ind);

/**
 * @brief Get the LB fluid velocity for a single node.
 */
const Utils::Vector3d lb_lbnode_get_velocity(const Utils::Vector3i &ind);

/**
 * @brief Get the LB fluid velocity for a single node.
 */
const Utils::Vector3d
lb_lbnode_get_velocity_at_boundary(const Utils::Vector3i &ind);

/**
 * @brief Get the LB fluid pressure tensor for a single node.
 */
const Utils::Vector6d lb_lbnode_get_pressure_tensor(const Utils::Vector3i &ind);

/**
 * @brief Get force applied on an lb node during the previous integration step
 */
const Utils::Vector3d
lb_lbnode_get_last_applied_force(const Utils::Vector3i &ind);

/** @brief Calculate the average pressure tensor of all nodes by accumulating
 *  over all nodes and dividing by the number of nodes.
 *  Returns the lower triangle of the LB pressure tensor.
 */
const Utils::Vector6d lb_lbfluid_get_pressure_tensor();

/**
 * @brief Get the LB fluid boundary bool for a single node.
 */
bool lb_lbnode_is_boundary(const Utils::Vector3i &ind);

/**
 * @brief Clear boundaries
 */
void lb_lbfluid_clear_boundaries();

/**
 * @brief Add a boundary.
 */
void lb_lbfluid_update_boundary_from_shape(
    std::vector<int> const &raster_flat,
    std::vector<double> const &slip_velocity_flat);

/**
 * @brief Update a boundary slip velocity.
 */
void lb_lbfluid_update_boundary_from_list(std::vector<int> const &nodes_flat,
                                          std::vector<double> const &vel_flat);

/**
 * @brief Get the LB fluid velocity for a single node.
 */
const Utils::Vector3d lb_lbnode_get_boundary_force(const Utils::Vector3i &ind);

/**
 * @brief Remove single node from boundary
 */
void lb_lbnode_remove_from_boundary(const Utils::Vector3i &ind);

/**
 * @brief Get the LB fluid populations for a single node.
 */
const std::vector<double> lb_lbnode_get_pop(const Utils::Vector3i &ind);

void lb_lbfluid_save_checkpoint(const std::string &filename, bool binary);
void lb_lbfluid_load_checkpoint(const std::string &filename, bool binary);

/**
 * @brief Checks whether the given node index is within the LB lattice.
 */
bool lb_lbnode_is_index_valid(const Utils::Vector3i &ind);

/**
 * @brief returns the shape of the LB fluid lattice
 */
Utils::Vector3i lb_lbfluid_get_shape();

Utils::Vector3d lb_lbfluid_calc_fluid_momentum();

/**
 * @brief Calculates the interpolated fluid velocity on the master process.
 * @param pos Position at which the velocity is to be calculated.
 * @retval interpolated fluid velocity.
 */
const Utils::Vector3d
lb_lbfluid_get_interpolated_velocity(const Utils::Vector3d &pos);

/**
 * @brief Calculates the interpolated fluid density on the master process.
 * @param pos Position at which the density is to be calculated.
 * @retval interpolated fluid density.
 */
double lb_lbfluid_get_interpolated_density(const Utils::Vector3d &pos);

/**
 * @brief Calculates the interpolated force to be applied on the master process.
 * @param pos Position at which the force is to be calculated.
 * @retval interpolated force to be applied.
 */
const Utils::Vector3d
lb_lbfluid_get_force_to_be_applied(const Utils::Vector3d &pos);

/**
 * @brief Distributes a force at a position which will be applied during
 *        the next integration loop
 * @param pos Position at which the force is beeing applied.
 * @param f   The force vector that is beeing applied.
 */
void lb_lbfluid_add_force_at_pos(const Utils::Vector3d &pos,
                                 const Utils::Vector3d &f);

void mpi_set_lattice_switch(ActiveLB lattice_switch);

#endif
