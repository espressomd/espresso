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

#include <cstddef>

#ifdef DIPOLES
#include "Observable_stat.hpp"

#include <utils/Vector.hpp>

#include <ParticleRange.hpp>
#include <boost/mpi/communicator.hpp>

/** Type codes for the type of dipolar interaction.
 *  Enumeration of implemented methods for the magnetostatic interaction.
 */
enum DipolarInteraction {
  /** Dipolar interaction switched off. */
  DIPOLAR_NONE = 0,
  /** Dipolar method is P3M. */
  DIPOLAR_P3M,
  /** Dipolar method is P3M plus DLC. */
  DIPOLAR_MDLC_P3M,
  /** Dipolar method is all with all and no replicas. */
  DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA,
  /** Dipolar method is magnetic dipolar direct summation. */
  DIPOLAR_DS,
  /** Dipolar method is direct summation plus DLC. */
  DIPOLAR_MDLC_DS,
  /** Dipolar method is direct summation on GPU. */
  DIPOLAR_DS_GPU,
#ifdef DIPOLAR_BARNES_HUT
  /** Dipolar method is direct summation on GPU by Barnes-Hut algorithm. */
  DIPOLAR_BH_GPU,
#endif
  /** Dipolar method is ScaFaCoS. */
  DIPOLAR_SCAFACOS
};

/** Interaction parameters for the %dipole interaction. */
struct Dipole_parameters {
  double prefactor;

  DipolarInteraction method;
};

/** Structure containing the %dipole parameters. */
extern Dipole_parameters dipole;

namespace Dipole {
/** Number of electrostatic contributions to the system pressure calculation. */
constexpr size_t pressure_n() { return 0; }

/** Number of electrostatic contributions to the system energy calculation.
 *  - slot 0: energies from particle pairs and magnetic field constraints
 *  - slot 1: energies from magnetostatics solvers
 *  - slot 2: energy corrections
 */
inline size_t energy_n() {
  switch (dipole.method) {
  case DIPOLAR_NONE:
    return 1; // because there may be an external magnetic field
  case DIPOLAR_P3M:
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
  case DIPOLAR_DS:
  case DIPOLAR_DS_GPU:
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
#endif
  case DIPOLAR_SCAFACOS:
    return 2;
  case DIPOLAR_MDLC_P3M:
  case DIPOLAR_MDLC_DS:
    return 3;
  default:
    return 0;
  }
}

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

int set_mesh();
void bcast_params(const boost::mpi::communicator &comm);

/** @brief Set the dipolar prefactor */
int set_Dprefactor(double prefactor);

void set_method_local(DipolarInteraction method);
} // namespace Dipole
#else  // DIPOLES
namespace Dipole {
constexpr size_t pressure_n() { return 0; }
constexpr size_t energy_n() { return 0; }
} // namespace Dipole
#endif // DIPOLES
#endif // ESPRESSO_DIPOLE_HPP
