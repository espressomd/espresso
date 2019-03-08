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
 *  Implementation of pressure.hpp.
 */

#include "cells.hpp"
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "npt.hpp"
#include "pressure_inline.hpp"
#include "virtual_sites.hpp"

#include "short_range_loop.hpp"

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
                        0.0,
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        1,
                        0,
                        {NPTGEOM_XDIR, NPTGEOM_YDIR, NPTGEOM_ZDIR},
                        0,
                        0,
                        0};

/************************************************************/
/* callbacks for setmd                                      */
/************************************************************/

/************************************************************/
/* local prototypes                                         */
/************************************************************/

/** Calculate long range virials (P3M, MMM2d...). */
void calc_long_range_virials();

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
inline void add_single_particle_virials(int v_comp, Particle &p) {
  add_kinetic_virials(&p, v_comp);
  add_bonded_virials(&p);
  add_three_body_bonded_stress(&p);
}

void pressure_calc(double *result, double *result_t, double *result_nb,
                   double *result_t_nb, int v_comp) {
  int n, i;
  double volume = box_l[0] * box_l[1] * box_l[2];

  if (!interactions_sanity_checks())
    return;

  init_virials(&virials);

  init_p_tensor(&p_tensor);

  init_virials_non_bonded(&virials_non_bonded);

  init_p_tensor_non_bonded(&p_tensor_non_bonded);

  on_observable_calc();
  // Run short-range loop if max cut >0
  if (max_cut > 0) {
    short_range_loop(
        [&v_comp](Particle &p) { add_single_particle_virials(v_comp, p); },
        [](Particle &p1, Particle &p2, Distance &d) {
          add_non_bonded_pair_virials(&(p1), &(p2), d.vec21.data(),
                                      sqrt(d.dist2), d.dist2);
        });
  } else {
    // Only add single particle virials
    for (auto &p : local_cells.particles()) {
      add_single_particle_virials(v_comp, p);
    }
  }
  /* rescale kinetic energy (=ideal contribution) */
  virials.data.e[0] /= (3.0 * volume * time_step * time_step);

  calc_long_range_virials();

#ifdef VIRTUAL_SITES
  virtual_sites()->pressure_and_stress_tensor_contribution(
      virials.virtual_sites, p_tensor.virtual_sites);
#endif

  for (n = 1; n < virials.data.n; n++)
    virials.data.e[n] /= 3.0 * volume;

  for (i = 0; i < 9; i++)
    p_tensor.data.e[i] /= (volume * time_step * time_step);

  for (i = 9; i < p_tensor.data.n; i++)
    p_tensor.data.e[i] /= volume;

  /* Intra- and Inter- part of nonbonded interaction */
  for (n = 0; n < virials_non_bonded.data_nb.n; n++)
    virials_non_bonded.data_nb.e[n] /= 3.0 * volume;

  for (i = 0; i < p_tensor_non_bonded.data_nb.n; i++)
    p_tensor_non_bonded.data_nb.e[i] /= volume;

  /* gather data */
  MPI_Reduce(virials.data.e, result, virials.data.n, MPI_DOUBLE, MPI_SUM, 0,
             comm_cart);
  MPI_Reduce(p_tensor.data.e, result_t, p_tensor.data.n, MPI_DOUBLE, MPI_SUM, 0,
             comm_cart);

  MPI_Reduce(virials_non_bonded.data_nb.e, result_nb,
             virials_non_bonded.data_nb.n, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
  MPI_Reduce(p_tensor_non_bonded.data_nb.e, result_t_nb,
             p_tensor_non_bonded.data_nb.n, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}

/************************************************************/

void calc_long_range_virials() {
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    fprintf(stderr,
            "WARNING: pressure calculated, but ELC pressure not implemented\n");
    break;
  case COULOMB_P3M_GPU:
    fprintf(
        stderr,
        "WARNING: pressure calculated, but GPU P3M pressure not implemented\n");
    break;
  case COULOMB_P3M: {
    p3m_charge_assign();
    virials.coulomb[1] = p3m_calc_kspace_forces(0, 1);
    p3m_charge_assign();
    p3m_calc_kspace_stress(p_tensor.coulomb + 9);
    break;
  }
#endif
  case COULOMB_MMM2D:
    fprintf(
        stderr,
        "WARNING: pressure calculated, but MMM2D pressure not implemented\n");
    break;
  case COULOMB_MMM1D:
  case COULOMB_MMM1D_GPU:
    fprintf(
        stderr,
        "WARNING: pressure calculated, but MMM1D pressure not implemented\n");
    break;
  default:
    break;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  switch (coulomb.Dmethod) {
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    fprintf(
        stderr,
        "WARNING: pressure calculated, but DAWAANR pressure not implemented\n");
    break;
  case DIPOLAR_MDLC_DS:
    fprintf(stderr,
            "WARNING: pressure calculated, but DLC pressure not implemented\n");
    break;
  case DIPOLAR_DS:
    fprintf(stderr, "WARNING: pressure calculated, but  MAGNETIC DIRECT SUM "
                    "pressure not implemented\n");
    break;

#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    fprintf(stderr,
            "WARNING: pressure calculated, but DLC pressure not implemented\n");
    break;
  case DIPOLAR_P3M: {
    int k;
    dp3m_dipole_assign();
    virials.dipolar[1] = dp3m_calc_kspace_forces(0, 1);

    for (k = 0; k < 3; k++)
      p_tensor.coulomb[9 + k * 3 + k] = virials.dipolar[1] / 3.;
    fprintf(stderr, "WARNING: stress tensor calculated, but dipolar P3M stress "
                    "tensor not implemented\n");
    fprintf(stderr, "WARNING: things have been added to the coulomb virial and "
                    "p_tensor arrays !!!!!!!\n");

    break;
  }
#endif
  default:
    break;
  }
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
  switch (coulomb.method) {
  case COULOMB_NONE:
    n_coulomb = 0;
    break;
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    n_coulomb = 2;
    break;
  default:
    n_coulomb = 1;
  }
#endif
#ifdef DIPOLES
  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE:
    n_dipolar = 0;
    break;
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    n_dipolar = 0;
    break;
  case DIPOLAR_DS:
    n_dipolar = 0;
    break;
  case DIPOLAR_P3M:
    n_dipolar = 2;
    break;
  default:
    n_dipolar = 0;
    break;
  }
#endif
#ifdef VIRTUAL_SITES
  n_vs = virtual_sites()->n_pressure_contribs();
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
  stat_nb->init_status_nb = 0;
}

/* Initialize the p_tensor */
/***************************/
void init_p_tensor(Observable_stat *stat) {
  // Determine number of contribution for different interaction types
  // bonded, nonbonded, Coulomb, dipolar, rigid bodies
  int n_pre, n_non_bonded, n_coulomb(0), n_dipolar(0), n_vs(0);

  n_pre = 1;
  n_non_bonded = (max_seen_particle_type * (max_seen_particle_type + 1)) / 2;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE:
    n_coulomb = 0;
    break;
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    n_coulomb = 2;
    break;
  default:
    n_coulomb = 1;
  }
#endif
#ifdef DIPOLES
  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE:
    n_dipolar = 0;
    break;
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    n_dipolar = 0;
    break;
  case DIPOLAR_DS:
    n_dipolar = 0;
    break;
  case DIPOLAR_P3M:
    n_dipolar = 2;
    break;
  default:
    n_dipolar = 0;
  }
#endif
#ifdef VIRTUAL_SITES
  n_vs = virtual_sites()->n_pressure_contribs();
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
  stat_nb->init_status_nb = 0;
}

/************************************************************/
void master_pressure_calc(int v_comp) {
  if (v_comp)
    mpi_gather_stats(3, total_pressure.data.e, total_p_tensor.data.e,
                     total_pressure_non_bonded.data_nb.e,
                     total_p_tensor_non_bonded.data_nb.e);
  else
    mpi_gather_stats(2, total_pressure.data.e, total_p_tensor.data.e,
                     total_pressure_non_bonded.data_nb.e,
                     total_p_tensor_non_bonded.data_nb.e);

  total_pressure.init_status = 1 + v_comp;
  total_p_tensor.init_status = 1 + v_comp;
  total_pressure_non_bonded.init_status_nb = 1 + v_comp;
  total_p_tensor_non_bonded.init_status_nb = 1 + v_comp;
}

/************************************************************/
int observable_compute_stress_tensor(int v_comp, double *A) {
  int i, j;
  double value;
  double p_vel[3];

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
      p_tensor.data.e[0] = 0.0;
      MPI_Reduce(nptiso.p_vel, p_vel, 3, MPI_DOUBLE, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      for (i = 0; i < 3; i++)
        if (nptiso.geometry & nptiso.nptgeom_dir[i])
          p_tensor.data.e[0] += p_vel[i];
      p_tensor.data.e[0] /= (nptiso.dimension * nptiso.volume);
      total_pressure.init_status = 1 + v_comp;
    } else
      master_pressure_calc(v_comp);
  }

  for (j = 0; j < 9; j++) {
    value = total_p_tensor.data.e[j];
    for (i = 1; i < total_p_tensor.data.n / 9; i++)
      value += total_p_tensor.data.e[9 * i + j];
    A[j] = value;
  }
  return 0;
}
