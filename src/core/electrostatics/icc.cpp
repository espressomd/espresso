/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_ELECTROSTATICS

#include "icc.hpp"

#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "actor/visitors.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "electrostatics/coulomb.hpp"
#include "electrostatics/coulomb_inline.hpp"
#include "electrostatics/p3m.hpp"
#include "errorhandling.hpp"
#include "integrators/Propagation.hpp"
#include "short_range_cabana.hpp"
#include "system/System.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/operations.hpp>

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
#include <Kokkos_Core.hpp>
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <variant>
#include <vector>

/** Calculate the electrostatic forces between source charges (= real charges)
 *  and wall charges. For each electrostatic method, the proper functions
 *  for short- and long-range parts are called. Long-range parts are calculated
 *  directly, short-range parts need helper functions according to the particle
 *  data organisation. This is a modified version of
 *  @ref System::System::calculate_forces.
 */
static void force_calc_icc(
    CellStructure &cell_structure,
    Coulomb::ShortRangeForceKernel::result_type const &coulomb_kernel,
    Coulomb::ShortRangeForceCorrectionsKernel::result_type const &elc_kernel) {
  // reset forces
  auto const reset_kernel = [](Particle &p) { p.force_and_torque() = {}; };
  cell_structure.for_each_local_particle(reset_kernel);
  cell_structure.for_each_ghost_particle(reset_kernel);
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  cell_structure.reset_local_force();
#endif

  // calc ICC forces
  cell_structure.non_bonded_loop(
      [coulomb_kernel_ptr = get_ptr(coulomb_kernel),
       elc_kernel_ptr = get_ptr(elc_kernel)](Particle &p1, Particle &p2,
                                             Distance const &d) {
        auto const q1q2 = p1.q() * p2.q();
        if (q1q2 != 0.) {
          auto force = (*coulomb_kernel_ptr)(q1q2, d.vec21, std::sqrt(d.dist2));
          p1.force() += force;
          p2.force() -= force;
#ifdef ESPRESSO_P3M
          if (elc_kernel_ptr) {
            (*elc_kernel_ptr)(p1.pos(), p2.pos(), p1.force_and_torque().f,
                              p2.force_and_torque().f, q1q2);
          }
#endif // ESPRESSO_P3M
        }
      });
}

void ICCStar::iteration() {
  try {
    sanity_check();
  } catch (std::runtime_error const &err) {
    runtimeErrorMsg() << err.what();
    return;
  }

  auto &system = get_system();
  auto &cell_structure = *system.cell_structure;
  auto const &coulomb = system.coulomb;
  auto const particles = cell_structure.local_particles();
  auto const prefactor = std::visit(
      [](auto const &ptr) { return ptr->prefactor; }, *coulomb.impl->solver);
  auto const pref = 1. / (prefactor * 2. * std::numbers::pi);
  auto const kernel = coulomb.pair_force_kernel();
  auto const elc_kernel = coulomb.pair_force_elc_kernel();
  icc_cfg.citeration = 0;

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  using execution_space = Kokkos::DefaultExecutionSpace;
  auto const &unique_particles = cell_structure.get_unique_particles();
  auto const &local_force = cell_structure.get_local_force();
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM

  auto global_max_rel_diff = 0.;

  for (int j = 0; j < icc_cfg.max_iterations; j++) {
    auto charge_density_max = 0.;

    // calculate electrostatic forces (SR+LR) excluding self-interactions
    force_calc_icc(cell_structure, kernel, elc_kernel);
    system.coulomb.calc_long_range_force();
    cell_structure.ghosts_reduce_forces();
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    // force reduction
    int num_threads = execution_space().concurrency();
    kokkos_parallel_range_for<Kokkos::RangePolicy<execution_space>>(
        "reduction", std::size_t{0}, unique_particles.size(),
        [&local_force, &unique_particles, num_threads](std::size_t const i) {
          auto &force = unique_particles.at(i)->force();
          for (int tid = 0; tid < num_threads; ++tid) {
            force[0] += local_force(i, tid, 0);
            force[1] += local_force(i, tid, 1);
            force[2] += local_force(i, tid, 2);
          }
        });
    Kokkos::fence();
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM

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
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    // refresh local properties
    update_aosoa_charges(cell_structure);
#endif

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

  system.on_particle_charge_change();
}

void icc_data::sanity_checks() const {
  if (n_icc <= 0)
    throw std::domain_error("Parameter 'n_icc' must be >= 1");
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
  if (areas.size() != static_cast<std::size_t>(n_icc))
    throw std::invalid_argument("Parameter 'areas' has incorrect shape");
  if (epsilons.size() != static_cast<std::size_t>(n_icc))
    throw std::invalid_argument("Parameter 'epsilons' has incorrect shape");
  if (sigmas.size() != static_cast<std::size_t>(n_icc))
    throw std::invalid_argument("Parameter 'sigmas' has incorrect shape");
  if (normals.size() != static_cast<std::size_t>(n_icc))
    throw std::invalid_argument("Parameter 'normals' has incorrect shape");
}

ICCStar::ICCStar(icc_data data) {
  data.sanity_checks();
  icc_cfg = std::move(data);
}

void ICCStar::on_activation() const {
  sanity_check();
  auto &system = get_system();
  system.on_particle_charge_change();
}

struct SanityChecksICC {
  template <typename T> void operator()(std::shared_ptr<T> const &) const {}
#ifdef ESPRESSO_P3M
#ifdef ESPRESSO_CUDA
  void operator()(std::shared_ptr<CoulombP3M> const &p) const {
    if (p->is_gpu()) {
      throw std::runtime_error("ICC does not work with P3M on GPU");
    }
  }
#endif // ESPRESSO_CUDA
  void
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &actor) const {
    if (actor->elc.dielectric_contrast_on) {
      throw std::runtime_error("ICC conflicts with ELC dielectric contrast");
    }
    std::visit(*this, actor->base_solver);
  }
#endif // ESPRESSO_P3M
  [[noreturn]] void operator()(std::shared_ptr<DebyeHueckel> const &) const {
    throw std::runtime_error("ICC does not work with DebyeHueckel.");
  }
  [[noreturn]] void operator()(std::shared_ptr<ReactionField> const &) const {
    throw std::runtime_error("ICC does not work with ReactionField.");
  }
};

void ICCStar::sanity_check() const {
  sanity_checks_active_solver();
#ifdef ESPRESSO_NPT
  if (get_system().has_npt_enabled()) {
    throw std::runtime_error("ICC does not work in the NpT ensemble");
  }
#endif
}

void ICCStar::sanity_checks_active_solver() const {
  auto &system = get_system();
  if (system.coulomb.impl->solver) {
    std::visit(SanityChecksICC(), *system.coulomb.impl->solver);
  } else {
    throw std::runtime_error("An electrostatics solver is needed by ICC");
  }
}

bool System::System::has_icc_enabled() const {
  return coulomb.impl->extension and
         std::holds_alternative<std::shared_ptr<ICCStar>>(
             *coulomb.impl->extension);
}

void System::System::update_icc_particles() {
  if (coulomb.impl->extension) {
    if (auto icc = std::get_if<std::shared_ptr<ICCStar>>(
            get_ptr(coulomb.impl->extension))) {
      (**icc).iteration();
    }
  }
}

#endif // ESPRESSO_ELECTROSTATICS
