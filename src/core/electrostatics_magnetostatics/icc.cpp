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
  Detailed Information about the method is included in the corresponding header
  file \ref icc.hpp.
*/

#include "icc.hpp"

#ifdef ELECTROSTATICS

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "electrostatics_magnetostatics/p3m_gpu.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

#include "short_range_loop.hpp"
#include <utils/NoOp.hpp>

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/coulomb_inline.hpp"

iccp3m_struct iccp3m_cfg;

/* functions that are used in icc* to compute the electric field acting on the
 * induced charges, excluding forces other than the electrostatic ones */
void init_forces_iccp3m();

/** Calculation of the electrostatic forces between source charges (= real
 * charges) and wall charges. For each electrostatic method the proper functions
 * for short and long range parts are called. Long Range Parts are calculated
 * directly, short range parts need helper functions according to the particle
 * data organisation. A modified version of \ref force_calc in \ref forces.hpp.
 */
void force_calc_iccp3m();

/** Variant of add_non_bonded_pair_force where only Coulomb
 *  contributions are calculated   */
inline void add_non_bonded_pair_force_iccp3m(Particle *p1, Particle *p2,
                                             double d[3], double dist,
                                             double dist2) {
  Utils::Vector3d force{};
  Coulomb::calc_pair_force(p1, p2, d, dist, force);

  p1->f.f += force;
  p2->f.f -= force;
}

void iccp3m_alloc_lists() {
  auto const n_ic = iccp3m_cfg.n_ic;

  iccp3m_cfg.areas.resize(n_ic);
  iccp3m_cfg.ein.resize(n_ic);
  iccp3m_cfg.normals.resize(n_ic);
  iccp3m_cfg.sigma.resize(n_ic);
}

int iccp3m_iteration() {
  if (iccp3m_cfg.n_ic == 0)
    return 0;

  Coulomb::iccp3m_sanity_check();

  if ((iccp3m_cfg.eout <= 0)) {
    runtimeErrorMsg()
        << "ICCP3M: nonpositive dielectric constant is not allowed.";
  }

  auto const pref = 1.0 / (coulomb.prefactor * 6.283185307);
  iccp3m_cfg.citeration = 0;

  double globalmax = 1e100;

  for (int j = 0; j < iccp3m_cfg.num_iteration; j++) {
    double hmax = 0.;

    force_calc_iccp3m(); /* Calculate electrostatic forces (SR+LR) excluding
                            source source interaction*/
    ghost_communicator(&cell_structure.collect_ghost_force_comm);

    double diff = 0;

    for (auto &p : local_cells.particles()) {
      if (p.p.identity < iccp3m_cfg.n_ic + iccp3m_cfg.first_id &&
          p.p.identity >= iccp3m_cfg.first_id) {
        auto const id = p.p.identity - iccp3m_cfg.first_id;
        /* the dielectric-related prefactor: */
        auto const del_eps = (iccp3m_cfg.ein[id] - iccp3m_cfg.eout) /
                             (iccp3m_cfg.ein[id] + iccp3m_cfg.eout);
        /* calculate the electric field at the certain position */
        auto const E = p.f.f / p.p.q + iccp3m_cfg.ext_field;

        if (E[0] == 0 && E[1] == 0 && E[2] == 0) {
          runtimeErrorMsg()
              << "ICCP3M found zero electric field on a charge. This must "
                 "never happen";
        }

        /* recalculate the old charge density */
        auto const hold = p.p.q / iccp3m_cfg.areas[id];
        /* determine if it is higher than the previously highest charge
         * density */
        hmax = std::max(hmax, std::abs(hold));

        auto const f1 = del_eps * pref * (E * iccp3m_cfg.normals[id]);
        auto const f2 = (not iccp3m_cfg.sigma.empty())
                            ? (2 * iccp3m_cfg.eout) /
                                  (iccp3m_cfg.eout + iccp3m_cfg.ein[id]) *
                                  (iccp3m_cfg.sigma[id])
                            : 0.;
        /* relative variation: never use an estimator which can be negative
         * here */
        auto const hnew =
            (1. - iccp3m_cfg.relax) * hold + (iccp3m_cfg.relax) * (f1 + f2);

        /* Take the largest error to check for convergence */
        auto const relative_difference =
            std::abs(1 * (hnew - hold) / (hmax + std::abs(hnew + hold)));

        diff = std::max(diff, relative_difference);

        p.p.q = hnew * iccp3m_cfg.areas[id];

        /* check if the charge now is more than 1e6, to determine if ICC still
         * leads to reasonable results */
        /* this is kind a arbitrary measure but, does a good job spotting
         * divergence !*/
        if (std::abs(p.p.q) > 1e6) {
          runtimeErrorMsg()
              << "too big charge assignment in iccp3m! q >1e6 , assigned "
                 "charge= "
              << p.p.q;

          diff = 1e90; /* A very high value is used as error code */
          break;
        }
      }
    } /* cell particles */
    /* Update charges on ghosts. */
    ghost_communicator(&cell_structure.exchange_ghosts_comm);

    iccp3m_cfg.citeration++;

    MPI_Allreduce(&diff, &globalmax, 1, MPI_DOUBLE, MPI_MAX, comm_cart);

    if (globalmax < iccp3m_cfg.convergence)
      break;
  } /* iteration */

  if (globalmax > iccp3m_cfg.convergence) {
    runtimeErrorMsg()
        << "ICC failed to converge in the given number of maximal steps.";
  }

  on_particle_charge_change();

  return iccp3m_cfg.citeration;
}

void force_calc_iccp3m() {
  init_forces_iccp3m();

  short_range_loop(Utils::NoOp{}, [](Particle &p1, Particle &p2, Distance &d) {
    /* calc non bonded interactions */
    add_non_bonded_pair_force_iccp3m(&(p1), &(p2), d.vec21.data(),
                                     sqrt(d.dist2), d.dist2);
  });

  Coulomb::calc_long_range_force();
}

void init_forces_iccp3m() {
  for (auto &p : local_cells.particles()) {
    p.f = ParticleForce{};
  }

  for (auto &p : ghost_cells.particles()) {
    p.f = ParticleForce{};
  }
}

#endif
