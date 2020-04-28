/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
/** \file
 *  Implementation of pressure.hpp.
 */

#include "cells.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "integrate.hpp"
#include "npt.hpp"
#include "pressure_inline.hpp"
#include "virtual_sites.hpp"
#include <boost/range/algorithm/copy.hpp>

#include "short_range_loop.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

Observable_stat virials{};
Observable_stat total_pressure{};
Observable_stat p_tensor{};
Observable_stat total_p_tensor{};

/* Observables used in the calculation of intra- and inter- molecular
   non-bonded contributions to pressure and to stress tensor */
Observable_stat_non_bonded virials_non_bonded{};
Observable_stat_non_bonded total_pressure_non_bonded{};
Observable_stat_non_bonded p_tensor_non_bonded{};
Observable_stat_non_bonded total_p_tensor_non_bonded{};

nptiso_struct nptiso = {0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        true,
                        0,
                        {NPTGEOM_XDIR, NPTGEOM_YDIR, NPTGEOM_ZDIR},
                        0,
                        false,
                        0};

/************************************************************/
/* local prototypes                                         */
/************************************************************/

/** Calculate long range virials (P3M, ...). */
void calc_long_range_virials(const ParticleRange &particles);

/** Initializes a virials Observable stat. */
void init_virials(Observable_stat *stat);

/** Initializes a virials Observable stat. */
void init_virials_non_bonded(Observable_stat_non_bonded *stat_nb);

/** on the master node: calc energies only if necessary
 *  @param v_comp flag which enables (1) compensation of the velocities required
 *                for deriving a pressure reflecting \ref nptiso_struct::p_inst
 *                (hence it only works with domain decomposition); naturally it
 *                therefore doesn't make sense to use it without NpT.
 */
void master_pressure_calc(int v_comp);

/** Initializes stat to be used by \ref pressure_calc. */
void init_p_tensor(Observable_stat *stat);

/** Initializes stat_nb to be used by \ref pressure_calc. */
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb);

/*********************************/
/* Scalar and Tensorial Pressure */
/*********************************/
static void add_single_particle_virials(Particle &p) {
  cell_structure.execute_bond_handler(p, add_bonded_stress);
}

void pressure_calc(double *result, double *result_t, double *result_nb,
                   double *result_t_nb, int v_comp) {
  auto const volume = box_geo.volume();

  if (!interactions_sanity_checks())
    return;

  init_virials(&virials);

  init_p_tensor(&p_tensor);

  init_virials_non_bonded(&virials_non_bonded);

  init_p_tensor_non_bonded(&p_tensor_non_bonded);

  on_observable_calc();

  for (auto const &p : cell_structure.local_particles()) {
    add_kinetic_virials(p, v_comp);
  }

  short_range_loop([](Particle &p) { add_single_particle_virials(p); },
                   [](Particle &p1, Particle &p2, Distance const &d) {
                     add_non_bonded_pair_virials(p1, p2, d.vec21,
                                                 sqrt(d.dist2));
                   });

  /* rescale kinetic energy (=ideal contribution) */
  virials.data[0] /= (3.0 * volume * time_step * time_step);

  calc_long_range_virials(cell_structure.local_particles());

#ifdef VIRTUAL_SITES
  {
    auto const vs_stress = virtual_sites()->stress_tensor();

    *virials.virtual_sites += trace(vs_stress);
    boost::copy(flatten(vs_stress), p_tensor.virtual_sites);
  }
#endif

  for (size_t n = 1; n < virials.data.size(); n++)
    virials.data[n] /= 3.0 * volume;

  for (int i = 0; i < 9; i++)
    p_tensor.data[i] /= (volume * time_step * time_step);

  for (size_t i = 9; i < p_tensor.data.size(); i++)
    p_tensor.data[i] /= volume;

  /* Intra- and Inter- part of nonbonded interaction */
  for (auto &value : virials_non_bonded.data_nb) {
    value /= 3.0 * volume;
  }

  for (auto &value : p_tensor_non_bonded.data_nb) {
    value /= volume;
  }

  /* gather data */
  MPI_Reduce(virials.data.data(), result, virials.data.size(), MPI_DOUBLE,
             MPI_SUM, 0, comm_cart);
  MPI_Reduce(p_tensor.data.data(), result_t, p_tensor.data.size(), MPI_DOUBLE,
             MPI_SUM, 0, comm_cart);

  MPI_Reduce(virials_non_bonded.data_nb.data(), result_nb,
             virials_non_bonded.data_nb.size(), MPI_DOUBLE, MPI_SUM, 0,
             comm_cart);
  MPI_Reduce(p_tensor_non_bonded.data_nb.data(), result_t_nb,
             p_tensor_non_bonded.data_nb.size(), MPI_DOUBLE, MPI_SUM, 0,
             comm_cart);
}

/************************************************************/

void calc_long_range_virials(const ParticleRange &particles) {
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  Coulomb::calc_pressure_long_range(virials, p_tensor, particles);
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  Dipole::calc_pressure_long_range();
#endif /*ifdef DIPOLES */
}

/* Initialize the virials used in the calculation of the scalar pressure */
/************************************************************/
void init_virials(Observable_stat *stat) {
  // Determine number of contribution for different interaction types
  // bonded, nonbonded, Coulomb, dipolar, rigid bodies
  int n_pre, n_non_bonded, n_coulomb(0), n_dipolar(0), n_vs(0);

  n_pre = 1;
  n_non_bonded = (max_seen_particle_type * (max_seen_particle_type + 1)) / 2;

#ifdef ELECTROSTATICS
  Coulomb::pressure_n(n_coulomb);
#endif
#ifdef DIPOLES
  Dipole::pressure_n();
#endif
#ifdef VIRTUAL_SITES
  n_vs = 1;
#endif

  // Allocate memory for the data
  obsstat_realloc_and_clear(stat, n_pre, bonded_ia_params.size(), n_non_bonded,
                            n_coulomb, n_dipolar, n_vs, 1);
  stat->init_status = 0;
}

/************************************************************/
void init_virials_non_bonded(Observable_stat_non_bonded *stat_nb) {
  int n_non_bonded;

  n_non_bonded = (max_seen_particle_type * (max_seen_particle_type + 1)) / 2;

  obsstat_realloc_and_clear_non_bonded(stat_nb, n_non_bonded, 1);
}

/* Initialize the p_tensor */
/***************************/
void init_p_tensor(Observable_stat *stat) {
  // Determine number of contribution for different interaction types
  // bonded, nonbonded, Coulomb, dipolar, rigid bodies
  int n_pre, n_non_bonded, n_coulomb(0), n_vs(0);

  n_pre = 1;
  n_non_bonded = (max_seen_particle_type * (max_seen_particle_type + 1)) / 2;

#ifdef ELECTROSTATICS
  Coulomb::pressure_n(n_coulomb);
#endif
#ifdef DIPOLES
  auto const n_dipolar = Dipole::pressure_n();
#else
  auto const n_dipolar = 0;
#endif
#ifdef VIRTUAL_SITES
  n_vs = 1;
#endif

  obsstat_realloc_and_clear(stat, n_pre, bonded_ia_params.size(), n_non_bonded,
                            n_coulomb, n_dipolar, n_vs, 9);
  stat->init_status = 0;
}

/***************************/
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb) {
  int n_nonbonded;
  n_nonbonded = (max_seen_particle_type * (max_seen_particle_type + 1)) / 2;

  obsstat_realloc_and_clear_non_bonded(stat_nb, n_nonbonded, 9);
}

/************************************************************/
void master_pressure_calc(int v_comp) {
  if (v_comp)
    mpi_gather_stats(3, total_pressure.data.data(), total_p_tensor.data.data(),
                     total_pressure_non_bonded.data_nb.data(),
                     total_p_tensor_non_bonded.data_nb.data());
  else
    mpi_gather_stats(2, total_pressure.data.data(), total_p_tensor.data.data(),
                     total_pressure_non_bonded.data_nb.data(),
                     total_p_tensor_non_bonded.data_nb.data());

  total_pressure.init_status = 1 + v_comp;
  total_p_tensor.init_status = 1 + v_comp;
}

/************************************************************/
int observable_compute_stress_tensor(int v_comp, double *A) {
  /* if desired (v_comp==1) replace ideal component with instantaneous one */
  if (total_pressure.init_status != 1 + v_comp) {
    init_virials(&total_pressure);
    init_p_tensor(&total_p_tensor);

    init_virials_non_bonded(&total_pressure_non_bonded);
    init_p_tensor_non_bonded(&total_p_tensor_non_bonded);

    if (v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) &&
        !(nptiso.invalidate_p_vel)) {
      if (total_pressure.init_status == 0)
        master_pressure_calc(0);
      p_tensor.data[0] = 0.0;
      Utils::Vector3d p_vel;
      MPI_Reduce(nptiso.p_vel.data(), p_vel.data(), 3, MPI_DOUBLE, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      for (int i = 0; i < 3; i++)
        if (nptiso.geometry & nptiso.nptgeom_dir[i])
          p_tensor.data[0] += p_vel[i];
      p_tensor.data[0] /= (nptiso.dimension * nptiso.volume);
      total_pressure.init_status = 1 + v_comp;
    } else
      master_pressure_calc(v_comp);
  }

  for (int j = 0; j < 9; j++) {
    double value = total_p_tensor.data[j];
    for (size_t i = 1; i < total_p_tensor.data.size() / 9; i++)
      value += total_p_tensor.data[9 * i + j];
    A[j] = value;
  }
  return 0;
}
