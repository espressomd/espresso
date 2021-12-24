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
 *  Functions to compute the electric field acting on the induced charges,
 *  excluding forces other than the electrostatic ones. Detailed information
 *  about the ICC* method is included in the corresponding header file
 *  \ref icc.hpp.
 */

#include "config.hpp"

#ifdef ELECTROSTATICS

#include "icc.hpp"

#include "CellStructure.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/coulomb_inline.hpp"
#include "errorhandling.hpp"
#include "event.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <tuple>
#include <vector>

icc_struct icc_cfg;

/** Calculate the electrostatic forces between source charges (= real charges)
 *  and wall charges. For each electrostatic method, the proper functions
 *  for short- and long-range parts are called. Long-range parts are calculated
 *  directly, short-range parts need helper functions according to the particle
 *  data organisation. This is a modified version of \ref force_calc.
 */
void force_calc_icc(CellStructure &cell_structure,
                    const ParticleRange &particles,
                    const ParticleRange &ghost_particles);

/** Variant of @ref add_non_bonded_pair_force where only %Coulomb
 *  contributions are calculated
 */
inline void add_non_bonded_pair_force_icc(Particle &p1, Particle &p2,
                                          Utils::Vector3d const &d, double dist,
                                          double dist2) {
  auto forces = Coulomb::pair_force(p1, p2, d, dist);

  p1.f.f += std::get<0>(forces);
  p2.f.f -= std::get<0>(forces);
#ifdef P3M
  p1.f.f += std::get<1>(forces);
  p2.f.f += std::get<2>(forces);
#endif
}

void icc_iteration(CellStructure &cell_structure,
                   const ParticleRange &particles,
                   const ParticleRange &ghost_particles) {
  if (icc_cfg.n_icc == 0)
    return;

  Coulomb::icc_sanity_check();

  auto const pref = 1.0 / (coulomb.prefactor * 2 * Utils::pi());
  icc_cfg.citeration = 0;

  double globalmax = 0.;

  for (int j = 0; j < icc_cfg.num_iteration; j++) {
    double charge_density_max = 0.;

    // calculate electrostatic forces (SR+LR) excluding self-interactions
    force_calc_icc(cell_structure, particles, ghost_particles);
    cell_structure.ghosts_reduce_forces();

    double diff = 0;

    for (auto &p : particles) {
      if (p.p.identity < icc_cfg.n_icc + icc_cfg.first_id &&
          p.p.identity >= icc_cfg.first_id) {
        auto const id = p.p.identity - icc_cfg.first_id;
        /* the dielectric-related prefactor: */
        auto const del_eps =
            (icc_cfg.ein[id] - icc_cfg.eout) / (icc_cfg.ein[id] + icc_cfg.eout);
        /* calculate the electric field at the certain position */
        auto const local_e_field = p.f.f / p.p.q + icc_cfg.ext_field;

        if (local_e_field.norm2() == 0) {
          runtimeErrorMsg()
              << "ICC found zero electric field on a charge. This must "
                 "never happen";
        }

        auto const charge_density_old = p.p.q / icc_cfg.areas[id];

        charge_density_max =
            std::max(charge_density_max, std::abs(charge_density_old));

        auto const charge_density_update =
            del_eps * pref * (local_e_field * icc_cfg.normals[id]) +
            2 * icc_cfg.eout / (icc_cfg.eout + icc_cfg.ein[id]) *
                icc_cfg.sigma[id];
        /* relative variation: never use an estimator which can be negative
         * here */
        auto const charge_density_new =
            (1. - icc_cfg.relax) * charge_density_old +
            (icc_cfg.relax) * charge_density_update;

        /* Take the largest error to check for convergence */
        auto const relative_difference =
            std::abs((charge_density_new - charge_density_old) /
                     (charge_density_max +
                      std::abs(charge_density_new + charge_density_old)));

        diff = std::max(diff, relative_difference);

        p.p.q = charge_density_new * icc_cfg.areas[id];

        /* check if the charge now is more than 1e6, to determine if ICC still
         * leads to reasonable results. This is kind of an arbitrary measure
         * but does a good job spotting divergence! */
        if (std::abs(p.p.q) > 1e6) {
          runtimeErrorMsg()
              << "too big charge assignment in icc! q >1e6 , assigned "
                 "charge= "
              << p.p.q;

          diff = 1e90; /* A very high value is used as error code */
          break;
        }
      }
    } /* cell particles */
    /* Update charges on ghosts. */
    cell_structure.ghosts_update(Cells::DATA_PART_PROPERTIES);

    icc_cfg.citeration++;

    boost::mpi::all_reduce(comm_cart, diff, globalmax,
                           boost::mpi::maximum<double>());

    if (globalmax < icc_cfg.convergence)
      break;
  } /* iteration */

  if (globalmax > icc_cfg.convergence) {
    runtimeErrorMsg()
        << "ICC failed to converge in the given number of maximal steps.";
  }

  on_particle_charge_change();
}

void force_calc_icc(CellStructure &cell_structure,
                    const ParticleRange &particles,
                    const ParticleRange &ghost_particles) {
  // reset forces
  for (auto &p : particles) {
    p.f.f = {};
  }
  for (auto &p : ghost_particles) {
    p.f.f = {};
  }

  // calc ICC forces
  cell_structure.non_bonded_loop(
      [](Particle &p1, Particle &p2, Distance const &d) {
        add_non_bonded_pair_force_icc(p1, p2, d.vec21, sqrt(d.dist2), d.dist2);
      });

  Coulomb::calc_long_range_force(particles);
}

void mpi_icc_init_local(const icc_struct &icc_cfg_) {
  icc_cfg = icc_cfg_;

  on_particle_charge_change();
  check_runtime_errors(comm_cart);
}

REGISTER_CALLBACK(mpi_icc_init_local)

int mpi_icc_init() {
  mpi_call(mpi_icc_init_local, icc_cfg);

  on_particle_charge_change();
  return check_runtime_errors(comm_cart);
}

void icc_set_params(int n_icc, double convergence, double relaxation,
                    Utils::Vector3d const &ext_field, int max_iterations,
                    int first_id, double eps_out, std::vector<double> &areas,
                    std::vector<double> &e_in, std::vector<double> &sigma,
                    std::vector<Utils::Vector3d> &normals) {
  if (n_icc < 0)
    throw std::runtime_error("ICC: invalid number of particles. " +
                             std::to_string(n_icc));
  if (convergence <= 0)
    throw std::runtime_error("ICC: invalid convergence value. " +
                             std::to_string(convergence));
  if (relaxation < 0 or relaxation > 2)
    throw std::runtime_error("ICC: invalid relaxation value. " +
                             std::to_string(relaxation));
  if (max_iterations <= 0)
    throw std::runtime_error("ICC: invalid max_iterations. " +
                             std::to_string(max_iterations));
  if (first_id < 0)
    throw std::runtime_error("ICC: invalid first_id. " +
                             std::to_string(first_id));
  if (eps_out <= 0)
    throw std::runtime_error("ICC: invalid eps_out. " +
                             std::to_string(eps_out));
  if (areas.size() != n_icc)
    throw std::runtime_error("ICC: invalid areas vector.");
  if (e_in.size() != n_icc)
    throw std::runtime_error("ICC: invalid e_in vector.");
  if (sigma.size() != n_icc)
    throw std::runtime_error("ICC: invalid sigma vector.");
  if (normals.size() != n_icc)
    throw std::runtime_error("ICC: invalid normals vector.");

  icc_cfg.n_icc = n_icc;
  icc_cfg.convergence = convergence;
  icc_cfg.relax = relaxation;
  icc_cfg.ext_field = ext_field;
  icc_cfg.num_iteration = max_iterations;
  icc_cfg.first_id = first_id;
  icc_cfg.eout = eps_out;

  icc_cfg.areas = std::move(areas);
  icc_cfg.ein = std::move(e_in);
  icc_cfg.sigma = std::move(sigma);
  icc_cfg.normals = std::move(normals);

  mpi_icc_init();
}

void icc_deactivate() {
  icc_cfg.n_icc = 0;
  icc_cfg.areas.resize(0);
  icc_cfg.ein.resize(0);
  icc_cfg.normals.resize(0);
  icc_cfg.sigma.resize(0);

  mpi_icc_init();
}
#endif
