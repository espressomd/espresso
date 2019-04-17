/*
Copyright (C) 2010-2018 The ESPResSo project

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
#include <boost/mpi/collectives.hpp>

#include "ShapeBasedConstraint.hpp"
#include "communication.hpp"
#include "energy_inline.hpp"
#include "errorhandling.hpp"
#include "forces_inline.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

namespace Constraints {
Utils::Vector3d ShapeBasedConstraint::total_force() const {
  return all_reduce(comm_cart, m_local_force, std::plus<>());
}

double ShapeBasedConstraint::total_normal_force() const {
  return all_reduce(comm_cart, m_outer_normal_force, std::plus<double>());
}

double ShapeBasedConstraint::min_dist() {
  double global_mindist = std::numeric_limits<double>::infinity();
  auto parts = local_cells.particles();

  auto const local_mindist = std::accumulate(
      parts.begin(), parts.end(), std::numeric_limits<double>::infinity(),
      [this](double min, Particle const &p) {
        IA_parameters *ia_params;
        ia_params = get_ia_param(p.p.type, part_rep.p.type);
        if (checkIfInteraction(ia_params)) {
          double vec[3], dist;
          m_shape->calculate_dist(folded_position(p), &dist, vec);
          return std::min(min, dist);
        }
        return min;
      });
  boost::mpi::reduce(comm_cart, local_mindist, global_mindist,
                     boost::mpi::minimum<double>(), 0);
  return global_mindist;
}

ParticleForce ShapeBasedConstraint::force(const Particle &p,
                                          const Utils::Vector3d &folded_pos,
                                          double t) {

  double dist = 0.;
  Utils::Vector3d dist_vec, force, torque1, torque2, outer_normal_vec;

  IA_parameters *ia_params = get_ia_param(p.p.type, part_rep.p.type);

  for (int j = 0; j < 3; j++) {
    force[j] = 0;
#ifdef ROTATION
    torque1[j] = torque2[j] = 0;
#endif
  }

  if (checkIfInteraction(ia_params)) {
    m_shape->calculate_dist(folded_pos, &dist, dist_vec.data());

    if (dist > 0) {
      outer_normal_vec = -dist_vec / dist;
      auto const dist2 = dist * dist;
      calc_non_bonded_pair_force(&p, &part_rep, ia_params, dist_vec.data(),
                                 dist, dist2, force.data(), torque1.data(),
                                 torque2.data());
#ifdef DPD
      if (thermo_switch & THERMO_DPD) {
        force += dpd_pair_force(&p, &part_rep, ia_params, dist_vec.data(), dist,
                                dist2);
      }
#endif
    } else if (m_penetrable && (dist <= 0)) {
      if ((!m_only_positive) && (dist < 0)) {
        auto const dist2 = dist * dist;
        calc_non_bonded_pair_force(&p, &part_rep, ia_params, dist_vec.data(),
                                   -1.0 * dist, dist2, force.data(),
                                   torque1.data(), torque2.data());
#ifdef DPD
        if (thermo_switch & THERMO_DPD) {
          force += dpd_pair_force(&p, &part_rep, ia_params, dist_vec.data(),
                                  dist, dist2);
        }
#endif
      }
    } else {
      runtimeErrorMsg() << "Constraint"
                        << " violated by particle " << p.p.identity << " dist "
                        << dist;
    }
  }

  m_local_force -= force;
  m_outer_normal_force -= outer_normal_vec * force;

#ifdef ROTATION
  part_rep.f.torque += torque2;
  return {force, torque1};
#else
  return force;
#endif
}

void ShapeBasedConstraint::add_energy(const Particle &p,
                                      const Utils::Vector3d &folded_pos,
                                      double t, Observable_stat &energy) const {
  double dist;
  IA_parameters *ia_params;
  double nonbonded_en = 0.0;

  ia_params = get_ia_param(p.p.type, part_rep.p.type);

  dist = 0.;
  if (checkIfInteraction(ia_params)) {
    double vec[3];
    m_shape->calculate_dist(folded_pos, &dist, vec);
    if (dist > 0) {
      nonbonded_en = calc_non_bonded_pair_energy(&p, &part_rep, ia_params, vec,
                                                 dist, dist * dist);
    } else if ((dist <= 0) && m_penetrable) {
      if (!m_only_positive && (dist < 0)) {
        nonbonded_en = calc_non_bonded_pair_energy(
            &p, &part_rep, ia_params, vec, -1.0 * dist, dist * dist);
      }
    } else {
      runtimeErrorMsg() << "Constraint "
                        << " violated by particle " << p.p.identity;
    }
  }
  if (part_rep.p.type >= 0)
    *obsstat_nonbonded(&energy, p.p.type, part_rep.p.type) += nonbonded_en;
}
} // namespace Constraints
