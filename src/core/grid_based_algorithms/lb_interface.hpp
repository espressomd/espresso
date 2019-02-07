#ifndef CORE_LB_INTERFACE
#define CORE_LB_INTERFACE

#include "config.hpp"
#include "lattice.hpp"
#include "utils/Vector.hpp"

#if defined(LB) || defined(LB_GPU)

void lb_lbfluid_update();
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

#ifdef LB
extern int transfer_momentum;
#endif

#ifdef LB_GPU
extern int transfer_momentum_gpu;
#endif

uint64_t lb_lbfluid_get_rng_state();
void lb_lbfluid_set_rng_state(uint64_t counter);

/** calculates the fluid velocity at a given position of the
 * lattice. Note that it can lead to undefined behaviour if the
 * position is not within the local lattice. */
const Vector3d lb_lbfluid_get_interpolated_velocity(const Vector3d &p);
void lb_lbfluid_add_force_density(const Vector3d &p,
                                  const Vector3d &force_density);
const Lattice &lb_lbfluid_get_lattice();
int lb_lbfluid_get_lattice_switch();

void lb_lbfluid_set_lattice_switch(int local_lattice_switch);
void lb_lbfluid_set_tau(double p_tau);
void lb_lbfluid_set_rho(double p_dens);
void lb_lbfluid_set_visc(double p_visc);
void lb_lbfluid_set_bulk_visc(double p_bulk_visc);
void lb_lbfluid_set_gamma_odd(double p_gamma_odd);
void lb_lbfluid_set_gamma_even(double p_gamma_even);
void lb_lbfluid_set_friction(double p_friction);
void lb_lbfluid_set_couple_flag(int couple_flag);
void lb_lbfluid_set_agrid(double p_agrid);
void lb_lbfluid_set_ext_force_density(int component,
                                      const Vector3d &force_density);
void lb_lbfluid_set_kT(double kT);
void lb_lbfluid_invalidate_particle_allocation();

void lb_lbnode_set_rho(const Vector3i &ind, double rho);
void lb_lbnode_set_u(const Vector3i &ind, const Vector3d &u);
void lb_lbnode_set_pop(const Vector3i &ind, const Vector<19, double> &pop);

double lb_lbfluid_get_tau();
double lb_lbfluid_get_agrid();
int lb_lbfluid_get_couple_flag();
double lb_lbfluid_get_friction();
double lb_lbfluid_get_gamma_even();
double lb_lbfluid_get_bulk_visc();
double lb_lbfluid_get_visc();
double lb_lbfluid_get_rho();
const Vector3d lb_lbfluid_get_ext_force_density();
double lb_lbfluid_get_kT();

double lb_lbnode_get_rho(const Vector3i &ind);
const Vector3d lb_lbnode_get_u(const Vector3i &ind);
const Vector<6, double> lb_lbnode_get_pi(const Vector3i &ind);
const Vector<6, double> lb_lbnode_get_pi_neq(const Vector3i &ind);
int lb_lbnode_get_boundary(const Vector3i &ind);
const Vector<19, double> lb_lbnode_get_pop(const Vector3i &ind);

/* IO routines */
void lb_lbfluid_print_vtk_boundary(const std::string &filename);
void lb_lbfluid_print_vtk_velocity(const std::string &filename,
                                   std::vector<int> = {-1, -1, -1},
                                   std::vector<int> = {-1, -1, -1});

void lb_lbfluid_print_boundary(const std::string &filename);
void lb_lbfluid_print_velocity(const std::string &filename);

void lb_lbfluid_save_checkpoint(const std::string &filename, int binary);
void lb_lbfluid_load_checkpoint(const std::string &filename, int binary);

bool lb_lbnode_is_index_valid(const Vector3i &ind);

/** Calculate the fluid velocity at a given position of the lattice.
 *  Note that it can lead to undefined behaviour if the
 *  position is not within the local lattice. This version of the function
 *  can be called without the position needing to be on the local processor.
 */
int lb_lbfluid_get_interpolated_velocity_global(Vector3d &p, double *v);
void lb_lbfluid_on_lb_params_change(int field);
#endif

#endif
