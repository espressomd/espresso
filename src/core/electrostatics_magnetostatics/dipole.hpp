#ifndef ESPRESSO_DIPOLE_HPP
#define ESPRESSO_DIPOLE_HPP

#include "statistics.hpp"

#ifdef ELECTROSTATICS
#ifdef DIPOLES

/** \name Compounds for Dipole interactions */
/*@{*/

/** field containing the interaction parameters for
 *  the Dipole interaction.  */
struct Dipole_parameters {
  double prefactor;

  DipolarInteraction method;
};
/*@}*/

/** Structure containing the Dipole parameters. */
extern Dipole_parameters dipole;

namespace Dipole {
// pressure
void pressure_n(int &n_dipolar);
void calc_pressure_long_range(Observable_stat &virials,
                              Observable_stat &p_tensor);

// nonbonded_interaction_data
void nonbonded_sanity_check(int &state);
void cutoff(double &ret);

// integrate
void integrate_sanity_check();

// initialize
void on_observable_calc();
void on_coulomb_change();
void on_boxl_change();
void init();

// forces
void calc_long_range_force();

// energy
void calc_energy_long_range(Observable_stat &energy);
void energy_n(int &n_dipolar);

// mdlc_correction
int set_mesh();

// communication
void bcast_params();

/** @brief Set the dipolar prefactor */
int set_Dprefactor(double prefactor);

} // namespace Dipole

#endif // DIPOLES
#endif // ELECTROSTATICS
#endif // ESPRESSO_DIPOLE_HPP
