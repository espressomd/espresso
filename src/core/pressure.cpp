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

namespace {
auto cpn = &Coulomb::pressure_n;
auto dpn = &Dipole::pressure_n;
} // namespace

Observable_stat virials{1, ::cpn, ::dpn};
Observable_stat total_pressure{1, ::cpn, ::dpn};
Observable_stat p_tensor{9, ::cpn, ::dpn};
Observable_stat total_p_tensor{9, ::cpn, ::dpn};

/* Observables used in the calculation of intra- and inter- molecular
 * non-bonded contributions to pressure and to stress tensor */
Observable_stat_non_bonded virials_non_bonded{1};
Observable_stat_non_bonded total_pressure_non_bonded{1};
Observable_stat_non_bonded p_tensor_non_bonded{9};
Observable_stat_non_bonded total_p_tensor_non_bonded{9};

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

/** Calculate long-range virials (P3M, ...). */
void calc_long_range_virials(const ParticleRange &particles) {
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  Coulomb::calc_pressure_long_range(virials, p_tensor, particles);
#endif
#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  Dipole::calc_pressure_long_range();
#endif
}

static void add_single_particle_virials(Particle &p) {
  cell_structure.execute_bond_handler(p, add_bonded_stress);
}

void pressure_calc(Observable_stat *result, Observable_stat *result_t,
                   Observable_stat_non_bonded *result_nb,
                   Observable_stat_non_bonded *result_t_nb, bool v_comp) {
  auto const volume = box_geo.volume();

  if (!interactions_sanity_checks())
    return;

  virials.realloc_and_clear();

  p_tensor.realloc_and_clear();

  virials_non_bonded.realloc_and_clear();

  p_tensor_non_bonded.realloc_and_clear();

  on_observable_calc();

  for (auto const &p : cell_structure.local_particles()) {
    add_kinetic_virials(p, v_comp);
  }

  short_range_loop([](Particle &p) { add_single_particle_virials(p); },
                   [](Particle &p1, Particle &p2, Distance const &d) {
                     add_non_bonded_pair_virials(p1, p2, d.vec21,
                                                 sqrt(d.dist2));
                   });

  calc_long_range_virials(cell_structure.local_particles());

#ifdef VIRTUAL_SITES
  if (!virials.virtual_sites.empty()) {
    auto const vs_stress = virtual_sites()->stress_tensor();

    virials.virtual_sites[0] += trace(vs_stress);
    boost::copy(flatten(vs_stress), p_tensor.virtual_sites.begin());
  }
#endif

  /* rescale kinetic energy (=ideal contribution) */
  virials.rescale(3.0 * volume, time_step);

  p_tensor.rescale(volume, time_step);

  /* Intra- and Inter- part of nonbonded interaction */
  virials_non_bonded.rescale(3.0 * volume);

  p_tensor_non_bonded.rescale(volume);

  /* gather data */
  virials.reduce(result);
  p_tensor.reduce(result_t);
  virials_non_bonded.reduce(result_nb);
  p_tensor_non_bonded.reduce(result_t_nb);
}

/** on the master node: calc energies only if necessary
 *  @param v_comp flag which enables compensation of the velocities required
 *                for deriving a pressure reflecting \ref nptiso_struct::p_inst
 *                (hence it only works with domain decomposition); naturally it
 *                therefore doesn't make sense to use it without NpT.
 */
void master_pressure_calc(bool v_comp) {
  mpi_gather_stats(v_comp ? GatherStats::pressure_v_comp
                          : GatherStats::pressure,
                   reinterpret_cast<void *>(&total_pressure),
                   reinterpret_cast<void *>(&total_p_tensor),
                   reinterpret_cast<void *>(&total_pressure_non_bonded),
                   reinterpret_cast<void *>(&total_p_tensor_non_bonded));

  total_pressure.is_initialized = true;
  total_p_tensor.is_initialized = true;
  total_pressure.v_comp = v_comp;
  total_p_tensor.v_comp = v_comp;
}

/** Calculate the sum of the ideal gas components of the instantaneous
 *  pressure, rescaled by the box dimensions.
 */
double calculate_npt_p_vel_rescaled() {
  Utils::Vector3d p_vel;
  MPI_Reduce(nptiso.p_vel.data(), p_vel.data(), 3, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  double acc = 0.0;
  for (int i = 0; i < 3; i++)
    if (nptiso.geometry & nptiso.nptgeom_dir[i])
      acc += p_vel[i];
  return acc / (nptiso.dimension * nptiso.volume);
}

void update_pressure(bool v_comp) {
  if (!(total_pressure.is_initialized && total_pressure.v_comp == v_comp)) {
    total_pressure.realloc_and_clear();
    total_p_tensor.realloc_and_clear();

    total_pressure_non_bonded.realloc_and_clear();
    total_p_tensor_non_bonded.realloc_and_clear();

    if (v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) &&
        !(nptiso.invalidate_p_vel)) {
      if (!total_pressure.is_initialized)
        master_pressure_calc(false);
      total_pressure.kinetic[0] = calculate_npt_p_vel_rescaled();
      total_pressure.v_comp = true;
    } else {
      master_pressure_calc(v_comp);
    }
  }
}

int observable_compute_stress_tensor(bool v_comp, double *A) {
  if (!(total_pressure.is_initialized && total_pressure.v_comp == v_comp)) {
    total_pressure.realloc_and_clear();
    total_p_tensor.realloc_and_clear();

    total_pressure_non_bonded.realloc_and_clear();
    total_p_tensor_non_bonded.realloc_and_clear();

    if (v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) &&
        !(nptiso.invalidate_p_vel)) {
      if (!total_pressure.is_initialized)
        master_pressure_calc(false);
      p_tensor.kinetic[0] = calculate_npt_p_vel_rescaled();
      total_pressure.v_comp = true;
    } else {
      master_pressure_calc(v_comp);
    }
  }

  for (int j = 0; j < 9; j++) {
    A[j] = total_p_tensor.accumulate(0, j);
  }
  return 0;
}
