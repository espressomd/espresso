/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *  Pressure calculation. Really similar to energy.hpp.
 */

#ifndef CORE_PRESSURE_HPP
#define CORE_PRESSURE_HPP

#include "statistics.hpp"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat virials, total_pressure, p_tensor, total_p_tensor;
///
extern Observable_stat_non_bonded virials_non_bonded, total_pressure_non_bonded,
    p_tensor_non_bonded, total_p_tensor_non_bonded;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/
void init_virials(Observable_stat *stat);
void init_virials_non_bonded(Observable_stat_non_bonded *stat_nb);
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb);
void init_p_tensor(Observable_stat *stat);
void master_pressure_calc(int v_comp);

/** Calculates the pressure in the system from a virial expansion.
 *  @param[out] result Calculated scalar pressure
 *  @param[out] result_t Calculated stress tensor
 *  @param[out] result_nb Calculated intra- and inter-molecular nonbonded
 *                        contributions to the scalar pressure
 *  @param[out] result_t_nb Calculated intra- and inter-molecular nonbonded
 *                          contributions to the stress tensor
 *  @param[in] v_comp flag which enables (1) compensation of the velocities
 *                    required for deriving a pressure reflecting
 *                    \ref nptiso_struct::p_inst (hence it only works with
 *                    domain decomposition); naturally it therefore doesn't
 *                    make sense to use it without NpT.
 */
void pressure_calc(double *result, double *result_t, double *result_nb,
                   double *result_t_nb, int v_comp);

/** Function to calculate stress tensor for the observables */
int observable_compute_stress_tensor(int v_comp, double *A);

void update_pressure(int v_comp);

/*@}*/

#endif
