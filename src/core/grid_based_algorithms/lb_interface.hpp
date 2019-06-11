#ifndef CORE_LB_INTERFACE
#define CORE_LB_INTERFACE

#include "config.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include "grid_based_algorithms/lb_constants.hpp"
#include <utils/Vector.hpp>

/** @brief LB implementation currently active. */
enum class ActiveLB { NONE, CPU, GPU };

/** @brief Switch determining the type of lattice dynamics. */
extern ActiveLB lattice_switch;

/**
 * @brief Propagate the LB fluid.
 */
void lb_lbfluid_propagate();

/**
 * @brief Event handler for integration start.
 */
void lb_lbfluid_on_integration_start();

/** Perform a full initialization of
 *  the lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values.
 */
void lb_lbfluid_init();

/** (Re-)initialize the fluid. */
void lb_lbfluid_reinit_fluid();

/** (Re-)initialize the derived parameters for the lattice Boltzmann system.
 *  The current state of the fluid is unchanged.
 */
void lb_lbfluid_reinit_parameters();

/**
 * @brief Get the current counter of the Philox RNG.
 */
uint64_t lb_lbfluid_get_rng_state();

/**
 * @brief Set the current counter of the Philox RNG.
 */
void lb_lbfluid_set_rng_state(uint64_t counter);

/**
 * @brief Return the instance of the Lattice within the LB method.
 */
const Lattice &lb_lbfluid_get_lattice();

/**
 * @brief Get the global variable @ref lattice_switch.
 */
ActiveLB lb_lbfluid_get_lattice_switch();

/**
 * @brief Set the global variable @ref lattice_switch.
 */
void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch);

/**
 * @brief Set the LB time step.
 */
void lb_lbfluid_set_tau(double p_tau);

/**
 * @brief Set the global LB density.
 */
void lb_lbfluid_set_density(double p_dens);

/**
 * @brief Set the global LB vicosity.
 */
void lb_lbfluid_set_viscosity(double p_visc);

/**
 * @brief Set the global LB bulk vicosity.
 */
void lb_lbfluid_set_bulk_viscosity(double p_bulk_visc);

/**
 * @brief Set the global LB relaxation parameter for odd modes.
 */
void lb_lbfluid_set_gamma_odd(double p_gamma_odd);

/**
 * @brief Set the global LB relaxation parameter for even modes.
 */
void lb_lbfluid_set_gamma_even(double p_gamma_even);

/**
 * @brief Set the global LB lattice spacing.
 */
void lb_lbfluid_set_agrid(double p_agrid);

/**
 * @brief Set the external force density acting on the LB fluid.
 */
void lb_lbfluid_set_ext_force_density(const Utils::Vector3d &force_density);

/**
 * @brief Set the LB fluid thermal energy.
 */
void lb_lbfluid_set_kT(double kT);

/**
 * @brief Invalidate the particle allocation on the GPU.
 */
void lb_lbfluid_invalidate_particle_allocation();

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
 * @brief Set the LB fluid populations for a single node.
 */
void lb_lbnode_set_pop(const Utils::Vector3i &ind, const Utils::Vector19d &pop);

/**
 * @brief Get the LB time step.
 */
double lb_lbfluid_get_tau();

/**
 * @brief Get the LB grid spacing.
 */
double lb_lbfluid_get_agrid();

/**
 * @brief Get the global LB relaxation parameter for even modes.
 */
double lb_lbfluid_get_gamma_even();

/**
 * @brief Get the global LB bulk viscosity.
 */
double lb_lbfluid_get_bulk_viscosity();

/**
 * @brief Get the global LB viscosity.
 */
double lb_lbfluid_get_viscosity();

/**
 * @brief Get the global LB density.
 */
double lb_lbfluid_get_density();

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
 * @brief Get the LB fluid density for a single node.
 */
double lb_lbnode_get_density(const Utils::Vector3i &ind);

/**
 * @brief Get the LB fluid velocity for a single node.
 */
const Utils::Vector3d lb_lbnode_get_velocity(const Utils::Vector3i &ind);
const Utils::Vector6d lb_lbnode_get_pi(const Utils::Vector3i &ind);
const Utils::Vector6d lb_lbnode_get_pi_neq(const Utils::Vector3i &ind);

/** calculates the average stress of all nodes by iterating
 * over all nodes and deviding by the number_of_nodes.
 */
const Utils::Vector6d lb_lbfluid_get_stress();

/**
 * @brief Get the LB fluid boundary bool for a single node.
 */
int lb_lbnode_get_boundary(const Utils::Vector3i &ind);

/**
 * @brief Get the LB fluid populations for a single node.
 */
const Utils::Vector19d lb_lbnode_get_pop(const Utils::Vector3i &ind);

/* IO routines */
void lb_lbfluid_print_vtk_boundary(const std::string &filename);
void lb_lbfluid_print_vtk_velocity(const std::string &filename,
                                   std::vector<int> = {-1, -1, -1},
                                   std::vector<int> = {-1, -1, -1});

void lb_lbfluid_print_boundary(const std::string &filename);
void lb_lbfluid_print_velocity(const std::string &filename);

void lb_lbfluid_save_checkpoint(const std::string &filename, int binary);
void lb_lbfluid_load_checkpoint(const std::string &filename, int binary);

/**
 * @brief Checks whether the given node index is within the LB lattice.
 */
bool lb_lbnode_is_index_valid(const Utils::Vector3i &ind);

void lb_lbfluid_on_lb_params_change(LBParam field);

Utils::Vector3d lb_lbfluid_calc_fluid_momentum();

#endif
