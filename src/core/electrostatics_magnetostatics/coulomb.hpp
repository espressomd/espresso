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
#ifndef ESPRESSO_COULOMB_HPP
#define ESPRESSO_COULOMB_HPP

#include "config.hpp"

#ifdef ELECTROSTATICS
#include "Observable_stat.hpp"

#include <ParticleRange.hpp>
#include <utils/Vector.hpp>

/** \name Type codes for the type of Coulomb interaction
    Enumeration of implemented methods for the electrostatic
    interaction.
*/
/************************************************************/
/*@{*/

enum CoulombMethod {
  COULOMB_NONE,      ///< %Coulomb interaction switched off (NONE)
  COULOMB_DH,        ///< %Coulomb method is Debye-Hueckel
  COULOMB_P3M,       ///< %Coulomb method is P3M
  COULOMB_P3M_GPU,   ///< %Coulomb method is P3M with GPU-based long-range part
  COULOMB_ELC_P3M,   ///< %Coulomb method is P3M plus ELC
  COULOMB_MMM1D,     ///< %Coulomb method is one-dimensional MMM
  COULOMB_RF,        ///< %Coulomb method is Reaction-Field
  COULOMB_MMM1D_GPU, ///< %Coulomb method is one-dimensional MMM running on GPU
  COULOMB_SCAFACOS,  ///< %Coulomb method is scafacos
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
                              Observable_stat &p_tensor,
                              const ParticleRange &particles);

void sanity_checks(int &state);
double cutoff(const Utils::Vector3d &box_l);
void deactivate();

void on_observable_calc();
void on_coulomb_change();
void on_boxl_change();
void init();

void calc_long_range_force(const ParticleRange &particles);

void calc_energy_long_range(Observable_stat &energy,
                            const ParticleRange &particles);
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

/** @brief Update particles with properties depending on other particles
 *   e.g., charges of ICC particles
 */
void update_dependent_particles();
} // namespace Coulomb
#endif // ELECTROSTATICS
#endif // ESPRESSO_COULOMB_HPP
