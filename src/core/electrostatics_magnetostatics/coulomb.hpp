#ifndef ESPRESSO_COULOMB_HPP
#define ESPRESSO_COULOMB_HPP

#include "config.hpp"
#include "utils/Vector.hpp"

extern double coulomb_cutoff;

#ifdef ELECTROSTATICS
#include "Observable_stat.hpp"

/** \name Type codes for the type of Coulomb interaction
    Enumeration of implemented methods for the electrostatic
    interaction.
*/
/************************************************************/
/*@{*/

enum CoulombMethod {
  COULOMB_NONE,      //< Coulomb interaction switched off (NONE)
  COULOMB_DH,        //< Coulomb method is Debye-Hueckel
  COULOMB_P3M,       //< Coulomb method is P3M
  COULOMB_MMM1D,     //< Coulomb method is one-dimensional MMM
  COULOMB_MMM2D,     //< Coulomb method is two-dimensional MMM
  COULOMB_ELC_P3M,   //< Coulomb method is P3M plus ELC
  COULOMB_RF,        //< Coulomb method is Reaction-Field
  COULOMB_P3M_GPU,   //< Coulomb method is P3M with GPU based long range part
                     // calculation
  COULOMB_MMM1D_GPU, //< Coulomb method is one-dimensional MMM running on GPU
  COULOMB_SCAFACOS,  //< Coulomb method is scafacos
};
/*@}*/

/** \name Compounds for Coulomb interactions */
/*@{*/

/** field containing the interaction parameters for
 *  the Coulomb  interaction.  */
struct Coulomb_parameters {
  /** bjerrum length times temperature. */
  double prefactor = 0.;

  double field_induced = 0.;
  double field_applied = 0.;

  /** Method to treat Coulomb interaction. */
  CoulombMethod method = COULOMB_NONE;
};

/** Structure containing the Coulomb parameters. */
extern Coulomb_parameters coulomb;

namespace Coulomb {
void pressure_n(int &n_coulomb);
void calc_pressure_long_range(Observable_stat &virials,
                              Observable_stat &p_tensor);

void sanity_checks(int &state);
double cutoff(const Utils::Vector3d &box_l);
void deactivate();

void integrate_sanity_check();
void on_observable_calc();
void on_coulomb_change();
void on_resort_particles();
void on_boxl_change();
void init();

void calc_long_range_force();

void calc_energy_long_range(Observable_stat &energy);
int energy_n();

int iccp3m_sanity_check();

int elc_sanity_check();

void bcast_coulomb_params();

/** @brief Set the electrostatics prefactor */
int set_prefactor(double prefactor);

/** @brief Deactivates the current Coulomb method
    This was part of coulomb_set_bjerrum()
*/
void deactivate_method();
} // namespace Coulomb
#endif // ELECTROSTATICS
#endif // ESPRESSO_COULOMB_HPP