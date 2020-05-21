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

#include "Observable_stat.hpp"
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

/** Scalar pressure of the system */
Observable_stat_wrapper obs_scalar_pressure{1};
/** Pressure tensor of the system */
Observable_stat_wrapper obs_pressure_tensor{9};
/** Contribution from the intra- and inter-molecular non-bonded interactions
 *  to the scalar pressure of the system.
 */
Observable_stat_non_bonded_wrapper obs_scalar_pressure_non_bonded{1};
/** Contribution from the intra- and inter-molecular non-bonded interactions
 *  to the pressure tensor of the system.
 */
Observable_stat_non_bonded_wrapper obs_pressure_tensor_non_bonded{9};

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
  Coulomb::calc_pressure_long_range(obs_scalar_pressure.local,
                                    obs_pressure_tensor.local, particles);
#endif
#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  Dipole::calc_pressure_long_range();
#endif
}

static void add_single_particle_virials(Particle &p) {
  cell_structure.execute_bond_handler(p, add_bonded_pressure_tensor);
}

void pressure_calc(bool v_comp) {
  auto const volume = box_geo.volume();

  if (!interactions_sanity_checks())
    return;

  obs_scalar_pressure.local.resize_and_clear();
  obs_scalar_pressure_non_bonded.local.resize_and_clear();
  obs_pressure_tensor.local.resize_and_clear();
  obs_pressure_tensor_non_bonded.local.resize_and_clear();

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
  if (!obs_scalar_pressure.local.virtual_sites.empty()) {
    auto const vs_pressure_tensor = virtual_sites()->pressure_tensor();

    obs_scalar_pressure.local.virtual_sites[0] += trace(vs_pressure_tensor);
    boost::copy(flatten(vs_pressure_tensor),
                obs_pressure_tensor.local.virtual_sites.begin());
  }
#endif

  /* rescale kinetic energy (=ideal contribution) */
  obs_scalar_pressure.local.rescale(3.0 * volume);

  obs_pressure_tensor.local.rescale(volume);

  /* Intra- and Inter- part of nonbonded interaction */
  obs_scalar_pressure_non_bonded.local.rescale(3.0 * volume);

  obs_pressure_tensor_non_bonded.local.rescale(volume);

  /* gather data */
  obs_scalar_pressure.reduce();
  obs_pressure_tensor.reduce();
  obs_scalar_pressure_non_bonded.reduce();
  obs_pressure_tensor_non_bonded.reduce();
}

/** Reduce the system scalar pressure and pressure tensor from all MPI ranks.
 *  @param v_comp flag which enables compensation of the velocities required
 *                for deriving a pressure reflecting \ref nptiso_struct::p_inst
 *                (hence it only works with domain decomposition); naturally it
 *                therefore doesn't make sense to use it without NpT.
 */
void master_pressure_calc(bool v_comp) {
  mpi_gather_stats(v_comp ? GatherStats::pressure_v_comp
                          : GatherStats::pressure);

  obs_scalar_pressure.v_comp = v_comp;
  obs_pressure_tensor.v_comp = v_comp;
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
  obs_scalar_pressure.resize();
  obs_scalar_pressure_non_bonded.resize();
  obs_pressure_tensor.resize();
  obs_pressure_tensor_non_bonded.resize();

  if (v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) &&
      !(nptiso.invalidate_p_vel)) {
    master_pressure_calc(false);
    obs_scalar_pressure.kinetic[0] = calculate_npt_p_vel_rescaled();
    obs_scalar_pressure.v_comp = true;
  } else {
    master_pressure_calc(v_comp);
  }
}

Utils::Vector9d observable_compute_pressure_tensor() {
  update_pressure(true);
  Utils::Vector9d pressure_tensor{};
  for (size_t j = 0; j < 9; j++) {
    pressure_tensor[j] = obs_pressure_tensor.accumulate(0, j);
  }
  return pressure_tensor;
}
