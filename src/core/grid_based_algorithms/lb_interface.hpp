#ifndef CORE_LB_INTERFACE
#define CORE_LB_INTERFACE

#include "config.hpp"
#include "utils/Vector.hpp"
#include "lattice.hpp"

#if defined(LB) || defined(LB_GPU)

void lb_update();
/**
 * @brief Event handler for integration start.
 */
void lb_on_integration_start();

/** Perform a full initialization of
 *  the lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values.
 */
void lb_init();

/** (Re-)initialize the fluid. */
void lb_reinit_fluid();

/** (Re-)initialize the derived parameters for the lattice Boltzmann system.
 *  The current state of the fluid is unchanged.
 */
void lb_reinit_parameters();

#ifdef LB
extern int transfer_momentum;
#endif

#ifdef LB_GPU
extern int transfer_momentum_gpu;
#endif

uint64_t lb_fluid_rng_state();
void lb_fluid_set_rng_state(uint64_t counter);

/** calculates the fluid velocity at a given position of the
 * lattice. Note that it can lead to undefined behaviour if the
 * position is not within the local lattice. */
Vector3d lb_lbfluid_get_interpolated_velocity(const Vector3d &p);
void lb_lbfluid_add_force_density(const Vector3d &p, const Vector3d &force_density);
const Lattice& lb_lbfluid_get_lattice();

int lb_lbfluid_set_density(double *p_dens);
int lb_lbfluid_get_density(double *p_dens);
int lb_lbfluid_set_visc(double *p_visc);
int lb_lbfluid_get_visc(double *p_visc);
int lb_lbfluid_set_bulk_visc(double *p_bulk_visc);
int lb_lbfluid_get_bulk_visc(double *p_bulk_visc);
int lb_lbfluid_set_gamma_odd(double *p_gamma_odd);
int lb_lbfluid_set_gamma_even(double *p_gamma_even);
int lb_lbfluid_get_gamma_even(double *p_gamma_even);
int lb_lbfluid_set_friction(double *p_friction);
int lb_lbfluid_get_friction(double *p_friction);
int lb_lbfluid_set_couple_flag(int couple_flag);
int lb_lbfluid_get_couple_flag(int *couple_flag);
int lb_lbfluid_set_agrid(double p_agrid);
int lb_lbfluid_get_agrid(double *p_agrid);
int lb_lbfluid_set_ext_force_density(int component, double p_fx, double p_fy,
                                     double p_fz);
int lb_lbfluid_get_ext_force_density(double *p_f);

int lb_lbfluid_set_tau(double p_tau);
int lb_lbfluid_get_tau(double *p_tau);
#ifdef SHANCHEN
int lb_lbfluid_set_remove_momentum(void);
int lb_lbfluid_set_shanchen_coupling(double *p_coupling);
int lb_lbfluid_set_mobility(double *p_mobility);
#endif
int lb_set_lattice_switch(int local_lattice_switch);

/* IO routines */
int lb_lbfluid_print_vtk_boundary(char *filename);
int lb_lbfluid_print_vtk_velocity(char *filename,
                                  std::vector<int> = {-1, -1, -1},
                                  std::vector<int> = {-1, -1, -1});

int lb_lbfluid_print_boundary(char *filename);
int lb_lbfluid_print_velocity(char *filename);

int lb_lbfluid_save_checkpoint(char *filename, int binary);
int lb_lbfluid_load_checkpoint(char *filename, int binary);

bool lb_lbnode_is_index_valid(const Vector3i &ind);
int lb_lbnode_get_rho(const Vector3i &ind, double *p_rho);
int lb_lbnode_get_u(const Vector3i &ind, double *u);
int lb_lbnode_get_pi(const Vector3i &ind, double *pi);
int lb_lbnode_get_pi_neq(const Vector3i &ind, double *pi_neq);
int lb_lbnode_get_boundary(const Vector3i &ind, int *p_boundary);
int lb_lbnode_get_pop(const Vector3i &ind, double *pop);

int lb_lbnode_set_rho(const Vector3i &ind, double *rho);
int lb_lbnode_set_u(const Vector3i &ind, double *u);
int lb_lbnode_set_pop(const Vector3i &ind, double *pop);

/** Calculate the fluid velocity at a given position of the lattice.
 *  Note that it can lead to undefined behaviour if the
 *  position is not within the local lattice. This version of the function
 *  can be called without the position needing to be on the local processor.
 */
int lb_lbfluid_get_interpolated_velocity_global(Vector3d &p, double *v);

#endif

#endif
