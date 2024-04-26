/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "ShapeBasedConstraint.hpp"

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "dpd.hpp"
#include "energy_inline.hpp"
#include "errorhandling.hpp"
#include "forces_inline.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>

namespace Constraints {
/** Check if a non-bonded interaction is defined */
static bool is_active(IA_parameters const &data) {
  return data.max_cut != INACTIVE_CUTOFF;
}

Utils::Vector3d ShapeBasedConstraint::total_force() const {
  return all_reduce(comm_cart, m_local_force, std::plus<>());
}

double ShapeBasedConstraint::total_normal_force() const {
  return all_reduce(comm_cart, m_outer_normal_force, std::plus<double>());
}

double ShapeBasedConstraint::min_dist(BoxGeometry const &box_geo,
                                      ParticleRange const &particles) const {
  auto global_mindist = std::numeric_limits<double>::infinity();

  auto const local_mindist = std::accumulate(
      particles.begin(), particles.end(),
      std::numeric_limits<double>::infinity(),
      [this, &box_geo](double min, Particle const &p) {
        auto const &ia_params = get_ia_param(p.type());
        if (is_active(ia_params)) {
          double dist;
          Utils::Vector3d vec;
          m_shape->calculate_dist(box_geo.folded_position(p.pos()), dist, vec);
          return std::min(min, dist);
        }
        return min;
      });
  boost::mpi::reduce(comm_cart, local_mindist, global_mindist,
                     boost::mpi::minimum<double>(), 0);
  return global_mindist;
}

ParticleForce ShapeBasedConstraint::force(Particle const &p,
                                          Utils::Vector3d const &folded_pos,
                                          double) {
  ParticleForce pf{};
  auto const &ia_params = get_ia_param(p.type());

  if (is_active(ia_params)) {
    double dist = 0.;
    Utils::Vector3d dist_vec;
    m_shape->calculate_dist(folded_pos, dist, dist_vec);
    auto const coulomb_kernel = m_system.coulomb.pair_force_kernel();

#ifdef DPD
    Utils::Vector3d dpd_force{};
#endif
    Utils::Vector3d outer_normal_vec{};

    if (dist > 0) {
      outer_normal_vec = -dist_vec / dist;
      pf = calc_central_radial_force(p, part_rep, ia_params, dist_vec, dist) +
           calc_central_radial_charge_force(p, part_rep, ia_params, dist_vec,
                                            dist, get_ptr(coulomb_kernel)) +
           calc_non_central_force(p, part_rep, ia_params, dist_vec, dist);

#ifdef DPD
      if (m_system.thermostat->thermo_switch & THERMO_DPD) {
        dpd_force = dpd_pair_force(p, part_rep, *m_system.thermostat->dpd,
                                   *m_system.box_geo, ia_params, dist_vec, dist,
                                   dist * dist);
        // Additional use of DPD here requires counter increase
        m_system.thermostat->dpd->rng_increment();
      }
#endif
    } else if (m_penetrable && (dist <= 0)) {
      if ((!m_only_positive) && (dist < 0)) {
        pf =
            calc_central_radial_force(p, part_rep, ia_params, dist_vec, -dist) +
            calc_central_radial_charge_force(p, part_rep, ia_params, dist_vec,
                                             -dist, get_ptr(coulomb_kernel)) +
            calc_non_central_force(p, part_rep, ia_params, dist_vec, -dist);

#ifdef DPD
        if (m_system.thermostat->thermo_switch & THERMO_DPD) {
          dpd_force = dpd_pair_force(p, part_rep, *m_system.thermostat->dpd,
                                     *m_system.box_geo, ia_params, dist_vec,
                                     dist, dist * dist);
          // Additional use of DPD here requires counter increase
          m_system.thermostat->dpd->rng_increment();
        }
#endif
      }
    } else {
      runtimeErrorMsg() << "Constraint violated by particle " << p.id()
                        << " dist " << dist;
    }

#ifdef ROTATION
    part_rep.torque() += calc_opposing_force(pf, dist_vec).torque;
#endif
#ifdef DPD
    pf.f += dpd_force;
#endif
    m_local_force -= pf.f;
    m_outer_normal_force -= outer_normal_vec * pf.f;
  }
  return pf;
}

void ShapeBasedConstraint::add_energy(const Particle &p,
                                      const Utils::Vector3d &folded_pos, double,
                                      Observable_stat &obs_energy) const {
  double energy = 0.0;
  auto const &ia_params = get_ia_param(p.type());

  if (is_active(ia_params)) {
    auto const coulomb_kernel = m_system.coulomb.pair_energy_kernel();
    double dist = 0.0;
    Utils::Vector3d vec;
    m_shape->calculate_dist(folded_pos, dist, vec);
    if (dist > 0) {
      energy = calc_non_bonded_pair_energy(p, part_rep, ia_params, vec, dist,
                                           get_ptr(coulomb_kernel));
    } else if ((dist <= 0) && m_penetrable) {
      if (!m_only_positive && (dist < 0)) {
        energy = calc_non_bonded_pair_energy(p, part_rep, ia_params, vec, -dist,
                                             get_ptr(coulomb_kernel));
      }
    } else {
      runtimeErrorMsg() << "Constraint violated by particle " << p.id();
    }
  }
  // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
  if (part_rep.type() >= 0) {
    obs_energy.add_non_bonded_contribution(
        p.type(), part_rep.type(), p.mol_id(), part_rep.mol_id(), energy);
  }
}
} // namespace Constraints
