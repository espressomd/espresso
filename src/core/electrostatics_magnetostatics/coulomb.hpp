#ifndef ESPRESSO_COULOMB_SWITCH_HPP
#define ESPRESSO_COULOMB_SWITCH_HPP

#include "statistics.hpp" // Observable_stat

#ifdef ELECTROSTATICS

void pressure_inline_add_pair_pressure(Particle *p1, Particle *p2, double *d,
                                       double dist, double dist2,
                                       Observable_stat &virials,
                                       Observable_stat &p_tensor);
void pressure_n_coulomb(int &n_coulomb);
void pressure_calc_long_range_coulomb_force(Observable_stat &virials,
                                            Observable_stat &p_tensor);

void nonbonded_interaction_data_coulomb_sanity_checks(int &state);
void nonbonded_interaction_data_calc_electrostatics_cutoff(double &ret);
void nonbonded_interaction_data_deactivate_coulomb_method();

void integrate_coulomb_sanity_check();

void initialize_on_observable_calc();
void initialize_on_coulomb_change();
void initialize_on_resort_particles();
void initialize_on_boxl_change();
void initialize_init_coulomb();

void forces_inline_calc_pair_coulomb_force(Particle *p1, Particle *p2,
                                           double *d, double dist, double dist2,
                                           Vector3d &force);

void forces_calc_long_range_coulomb_force();

void energy_inline_add_pair_coulomb_energy(Particle *p1, Particle *p2,
                                           double *d, double dist, double dist2,
                                           Observable_stat &energy);

void energy_calc_long_range_coulomb_energy(Observable_stat &energy);
void energy_n_coulomb(int &n_coulomb);

void icc_calc_pair_coulomb_force(Particle *p1, Particle *p2, double *d,
                                 double dist, double dist2, double *force);
void icc_calc_long_range_force_contribution_iccp3m();
int iccp3m_sanity_check();

int elc_switch_error();

void bcast_coulomb_params();

#endif // ELECTROSTATICS
#endif // ESPRESSO_COULOMB_SWITCH_HPP