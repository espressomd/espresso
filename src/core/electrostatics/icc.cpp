/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "icc.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "actor/visitors.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "electrostatics/coulomb.hpp"
#include "electrostatics/coulomb_inline.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "integrate.hpp"

#include <utils/constants.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/operations.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

/** Calculate the electrostatic forces between source charges (= real charges)
 *  and wall charges. For each electrostatic method, the proper functions
 *  for short- and long-range parts are called. Long-range parts are calculated
 *  directly, short-range parts need helper functions according to the particle
 *  data organisation. This is a modified version of \ref force_calc.
 */
static void force_calc_icc(
    CellStructure &cell_structure, ParticleRange const &particles,
    ParticleRange const &ghost_particles,
    Coulomb::ShortRangeForceKernel::result_type const &coulomb_kernel,
    Coulomb::ShortRangeForceCorrectionsKernel::result_type const &elc_kernel) {
  // reset forces
  for (auto &p : particles) {
    p.force() = {};
  }
  for (auto &p : ghost_particles) {
    p.force() = {};
  }

  // calc ICC forces
  cell_structure.non_bonded_loop(
      [coulomb_kernel_ptr = coulomb_kernel.get_ptr(),
       elc_kernel_ptr = elc_kernel.get_ptr()](Particle &p1, Particle &p2,
                                              Distance const &d) {
        auto const q1q2 = p1.q() * p2.q();
        if (q1q2 != 0.) {
          auto force = (*coulomb_kernel_ptr)(q1q2, d.vec21, std::sqrt(d.dist2));
          p1.force() += force;
          p2.force() -= force;
#ifdef P3M
          if (elc_kernel_ptr) {
            (*elc_kernel_ptr)(p1, p2, q1q2);
          }
#endif // P3M
        }
      });

  Coulomb::calc_long_range_force(particles);
}

void ICCStar::iteration(CellStructure &cell_structure,
                        ParticleRange const &particles,
                        ParticleRange const &ghost_particles) {

  try {
    sanity_check();
  } catch (std::runtime_error const &err) {
    runtimeErrorMsg() << err.what();
    return;
  }

  auto const prefactor =
      boost::apply_visitor(GetCoulombPrefactor(), *electrostatics_actor);
  auto const pref = 1. / (prefactor * 2. * Utils::pi());
  auto const kernel = Coulomb::pair_force_kernel();
  auto const elc_kernel = Coulomb::pair_force_elc_kernel();
  icc_cfg.citeration = 0;

  auto global_max_rel_diff = 0.;

  for (int j = 0; j < icc_cfg.max_iterations; j++) {
    auto charge_density_max = 0.;

    // calculate electrostatic forces (SR+LR) excluding self-interactions
    force_calc_icc(cell_structure, particles, ghost_particles, kernel,
                   elc_kernel);
    cell_structure.ghosts_reduce_forces();

    auto max_rel_diff = 0.;

    for (auto &p : particles) {
      auto const pid = p.id();
      if (pid >= icc_cfg.first_id and pid < icc_cfg.n_icc + icc_cfg.first_id) {
        if (p.q() == 0.) {
          runtimeErrorMsg()
              << "ICC found zero electric charge on a particle. This must "
                 "never happen";
          break;
        }
        auto const id = p.id() - icc_cfg.first_id;
        /* the dielectric-related prefactor: */
        auto const eps_in = icc_cfg.epsilons[id];
        auto const eps_out = icc_cfg.eps_out;
        auto const del_eps = (eps_in - eps_out) / (eps_in + eps_out);
        /* calculate the electric field at the certain position */
        auto const local_e_field = p.force() / p.q() + icc_cfg.ext_field;

        if (local_e_field.norm2() == 0.) {
          runtimeErrorMsg()
              << "ICC found zero electric field on a charge. This must "
                 "never happen";
        }

        auto const charge_density_old = p.q() / icc_cfg.areas[id];

        charge_density_max =
            std::max(charge_density_max, std::abs(charge_density_old));

        auto const charge_density_update =
            del_eps * pref * (local_e_field * icc_cfg.normals[id]) +
            2. * icc_cfg.eps_out / (icc_cfg.eps_out + icc_cfg.epsilons[id]) *
                icc_cfg.sigmas[id];
        /* relative variation: never use an estimator which can be negative
         * here */
        auto const charge_density_new =
            (1. - icc_cfg.relaxation) * charge_density_old +
            (icc_cfg.relaxation) * charge_density_update;

        /* Take the largest error to check for convergence */
        auto const relative_difference =
            std::abs((charge_density_new - charge_density_old) /
                     (charge_density_max +
                      std::abs(charge_density_new + charge_density_old)));

        max_rel_diff = std::max(max_rel_diff, relative_difference);

        p.q() = charge_density_new * icc_cfg.areas[id];

        /* check if the charge now is more than 1e6, to determine if ICC still
         * leads to reasonable results. This is kind of an arbitrary measure
         * but does a good job of spotting divergence! */
        if (std::abs(p.q()) > 1e6) {
          runtimeErrorMsg()
              << "Particle with id " << p.id() << " has a charge (q=" << p.q()
              << ") that is too large for the ICC algorithm";

          max_rel_diff = std::numeric_limits<double>::max();
          break;
        }
      }
    }

    /* Update charges on ghosts. */
    cell_structure.ghosts_update(Cells::DATA_PART_PROPERTIES);

    icc_cfg.citeration++;

    boost::mpi::all_reduce(comm_cart, max_rel_diff, global_max_rel_diff,
                           boost::mpi::maximum<double>());

    if (global_max_rel_diff < icc_cfg.convergence)
      break;
  }

  if (global_max_rel_diff > icc_cfg.convergence) {
    runtimeErrorMsg()
        << "ICC failed to converge in the given number of maximal steps.";
  }

  on_particle_charge_change();
}

void icc_data::sanity_checks() const {
  if (convergence <= 0.)
    throw std::domain_error("Parameter 'convergence' must be > 0");
  if (relaxation < 0. or relaxation > 2.)
    throw std::domain_error("Parameter 'relaxation' must be >= 0 and <= 2");
  if (max_iterations <= 0)
    throw std::domain_error("Parameter 'max_iterations' must be > 0");
  if (first_id < 0)
    throw std::domain_error("Parameter 'first_id' must be >= 0");
  if (eps_out <= 0.)
    throw std::domain_error("Parameter 'eps_out' must be > 0");

  assert(n_icc >= 1);
  assert(areas.size() == n_icc);
  assert(epsilons.size() == n_icc);
  assert(sigmas.size() == n_icc);
  assert(normals.size() == n_icc);
}

ICCStar::ICCStar(icc_data data) {
  data.sanity_checks();
  icc_cfg = std::move(data);
}

void ICCStar::on_activation() const {
  sanity_check();
  on_particle_charge_change();
}

struct SanityChecksICC : public boost::static_visitor<void> {
  template <typename T>
  void operator()(std::shared_ptr<T> const &actor) const {}
#ifdef P3M
#ifdef CUDA
  [[noreturn]] void
  operator()(std::shared_ptr<CoulombP3MGPU> const &actor) const {
    throw std::runtime_error("ICC does not work with P3MGPU");
  }
#endif // CUDA
  void
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &actor) const {
    if (actor->elc.dielectric_contrast_on) {
      throw std::runtime_error("ICC conflicts with ELC dielectric contrast");
    }
    boost::apply_visitor(*this, actor->base_solver);
  }
#endif // P3M
  [[noreturn]] void operator()(std::shared_ptr<DebyeHueckel> const &) const {
    throw std::runtime_error("ICC does not work with DebyeHueckel.");
  }
  [[noreturn]] void operator()(std::shared_ptr<ReactionField> const &) const {
    throw std::runtime_error("ICC does not work with ReactionField.");
  }
};

void ICCStar::sanity_check() const {
  sanity_checks_active_solver();
#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    throw std::runtime_error("ICC does not work in the NPT ensemble");
  }
#endif
}

void ICCStar::sanity_checks_active_solver() const {
  if (electrostatics_actor) {
    boost::apply_visitor(SanityChecksICC(), *electrostatics_actor);
  } else {
    throw std::runtime_error("An electrostatics solver is needed by ICC");
  }
}

void update_icc_particles() {
  if (electrostatics_extension) {
    if (auto icc = boost::get<std::shared_ptr<ICCStar>>(
            electrostatics_extension.get_ptr())) {
      (**icc).iteration(cell_structure, cell_structure.local_particles(),
                        cell_structure.ghost_particles());
    }
  }
}

#endif // ELECTROSTATICS
