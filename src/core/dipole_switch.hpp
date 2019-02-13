#ifndef ESPRESSO_DIPOLE_SWITCH_HPP
#define ESPRESSO_DIPOLE_SWITCH_HPP

#include "statistics.hpp"

void pressure_n_dipolar(int &n_dipolar);
void pressure_calc_long_range_dipole_force(Observable_stat &virials,
                                           Observable_stat &p_tensor);

void nonbonded_interaction_data_dipole_sanity_checks(int &state);
void nonbonded_interaction_data_calc_dipolar_cutoff(double &ret);

void integrate_dipole_sanity_check();

void initialize_on_observable_calc_dipole();
void initialize_on_coulomb_change_dipole();
void initialize_on_boxl_change_dipole();
void initialize_init_dipole();

void forces_inline_calc_pair_dipole_force(Particle *p1, Particle *p2, double *d,
                                          double dist, double dist2,
                                          Vector3d &force);

void forces_calc_long_range_dipole();

void energy_inline_add_pair_dipole_energy(Particle *p1, Particle *p2, double *d,
                                          double dist, double dist2,
                                          Observable_stat &energy);

void energy_calc_long_range_dipole_energy(Observable_stat &energy);
void energy_n_dipolar(int &n_dipolar);

int mdlc_correction_set_dipole_mesh();

#endif // ESPRESSO_DIPOLE_SWITCH_HPP
