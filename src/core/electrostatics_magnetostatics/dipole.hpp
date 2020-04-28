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
#ifndef ESPRESSO_DIPOLE_HPP
#define ESPRESSO_DIPOLE_HPP

#include "config.hpp"

extern double dipolar_cutoff;

#ifdef DIPOLES
#include "Observable_stat.hpp"

#include <utils/Vector.hpp>

#include <ParticleRange.hpp>
#include <boost/mpi/communicator.hpp>

/** \name Compounds for Dipole interactions */
/*@{*/

/** \name Type codes for the type of dipolar interaction
  Enumeration of implemented methods for the magnetostatic
  interaction.
 */
/************************************************************/
/*@{*/
enum DipolarInteraction {
  /** dipolar interaction switched off (NONE). */
  DIPOLAR_NONE = 0,
  /** dipolar method is P3M. */
  DIPOLAR_P3M,
  /** Dipolar method is P3M plus DLC. */
  DIPOLAR_MDLC_P3M,
  /** Dipolar method is all with all and no replicas */
  DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA,
  /** Dipolar method is magnetic dipolar direct sum */
  DIPOLAR_DS,
  /** Dipolar method is direct sum plus DLC. */
  DIPOLAR_MDLC_DS,
  /** Direct summation on gpu */
  DIPOLAR_DS_GPU,
#ifdef DIPOLAR_BARNES_HUT
  /** Direct summation on gpu by Barnes-Hut algorithm */
  DIPOLAR_BH_GPU,
#endif
  /** Scafacos library */
  DIPOLAR_SCAFACOS
};

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
int pressure_n();
void calc_pressure_long_range();

void nonbonded_sanity_check(int &state);
double cutoff(const Utils::Vector3d &box_l);

void on_observable_calc();
void on_coulomb_change();
void on_boxl_change();
void init();

void calc_long_range_force(const ParticleRange &particles);

void calc_energy_long_range(Observable_stat &energy,
                            const ParticleRange &particles);
void energy_n(int &n_dipolar);

int set_mesh();
void bcast_params(const boost::mpi::communicator &comm);

/** @brief Set the dipolar prefactor */
int set_Dprefactor(double prefactor);

void set_method_local(DipolarInteraction method);

} // namespace Dipole

#endif // DIPOLES
#endif // ESPRESSO_DIPOLE_HPP
