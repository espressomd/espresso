/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file pressure.hpp
    Pressure calculation. Really similar to \ref energy.hpp "energy.h".
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
extern Observable_stat_non_bonded virials_non_bonded, total_pressure_non_bonded, p_tensor_non_bonded, total_p_tensor_non_bonded;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/
void init_virials(Observable_stat *stat);
void init_virials_non_bonded(Observable_stat_non_bonded *stat_nb);
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb);
void init_p_tensor(Observable_stat *stat);
void master_pressure_calc(int v_comp);


/** Calculates the pressure in the system from a virial expansion using the terms from \ref calculate_verlet_virials or \ref nsq_calculate_virials dependeing on the used cell system.<BR>
    @param result here the data about the scalar pressure are stored
    @param result_t here the data about the stress tensor are stored
    @param result_nb here the data about the intra- and inter- molecular nonbonded contributions to scalar pressure are stored
    @param result_t_nb here the data about the intra- and inter- molecular nonbonded contributions to stress tensor are stored
    @param v_comp flag which enables (1) compensation of the velocities required
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT.
*/
void pressure_calc(double *result, double *result_t, double *result_nb, double *result_t_nb, int v_comp);

/** implementation of 'analyse local_stress_tensor */
int local_stress_tensor_calc (DoubleList *TensorInBin, int bins[3], int periodic[3], double range_start[3], double range[3]);

/** function to calculate stress tensor for the observables */
int observable_compute_stress_tensor(int v_comp, double *A);

void update_pressure(int v_comp);
void update_stress_tensor(int v_comp);
int analyze_local_stress_tensor(int* periodic, double* range_start, double* range, int* bins, DoubleList* local_stress_tensor);

/*@}*/

#endif
