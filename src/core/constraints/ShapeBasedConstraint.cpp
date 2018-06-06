#include <boost/mpi/collectives.hpp>

#include "ShapeBasedConstraint.hpp"
#include "communication.hpp"
#include "energy_inline.hpp"
#include "errorhandling.hpp"
#include "forces_inline.hpp"
#include "interaction_data.hpp"

namespace Constraints {
Vector3d ShapeBasedConstraint::total_force() const {
  Vector3d total_force;
  boost::mpi::all_reduce(comm_cart, m_local_force, total_force,
                         std::plus<Vector3d>());

  return total_force;
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
          m_shape->calculate_dist(folded_position(p).data(), &dist, vec);
          return std::min(min, dist);
        } else
          return min;
      });
  boost::mpi::reduce(comm_cart, local_mindist, global_mindist,
                     boost::mpi::minimum<double>(), 0);
  return global_mindist;
}

    void ShapeBasedConstraint::add_force(Particle *p, const Vector3d &folded_pos) {
  double dist, vec[3], force[3], torque1[3], torque2[3];

  IA_parameters *ia_params = get_ia_param(p->p.type, part_rep.p.type);

  dist = 0.;
  for (int j = 0; j < 3; j++) {
    force[j] = 0;
#ifdef ROTATION
    torque1[j] = torque2[j] = 0;
#endif
  }

  if (checkIfInteraction(ia_params)) {
    m_shape->calculate_dist(folded_pos.data(), &dist, vec);

    if (dist > 0) {
      auto const dist2 = dist * dist;
      calc_non_bonded_pair_force(p, &part_rep, ia_params, vec, dist, dist2,
                                 force, torque1, torque2);
#ifdef DPD
      if (thermo_switch & THERMO_DPD) {
          dpd_pair_force(p, &part_rep, ia_params, vec, dist, dist2);
      }
#endif
    } else if (m_penetrable && (dist <= 0)) {
      if ((!m_only_positive) && (dist < 0)) {
        auto const dist2 = dist * dist;
        calc_non_bonded_pair_force(p, &part_rep, ia_params, vec, -1.0 * dist,
                                   dist * dist, force, torque1, torque2);
#ifdef DPD
        if (thermo_switch & THERMO_DPD) {
            dpd_pair_force(p, &part_rep, ia_params, vec, dist, dist2);
        }
#endif
      }
    } else {
        runtimeErrorMsg() << "Constraint"
                          << " violated by particle " << p->p.identity
                          << " dist " << dist;
    }
  }
  for (int j = 0; j < 3; j++) {
    p->f.f[j] += force[j];
    m_local_force[j] -= force[j];
#ifdef ROTATION
    p->f.torque[j] += torque1[j];
    part_rep.f.torque[j] += torque2[j];
#endif
  }
}

void ShapeBasedConstraint::add_energy(const Particle *p, const Vector3d &folded_pos,
                                      Observable_stat &energy) const {
  double dist, vec[3];
  IA_parameters *ia_params;
  double nonbonded_en = 0.0;

  ia_params = get_ia_param(p->p.type, part_rep.p.type);

  dist = 0.;
  if (checkIfInteraction(ia_params)) {
    m_shape->calculate_dist(folded_pos.data(), &dist, vec);
    if (dist > 0) {
      nonbonded_en = calc_non_bonded_pair_energy(p, &part_rep, ia_params, vec,
                                                 dist, dist * dist);
    } else if ((dist <= 0) && m_penetrable) {
      if (!m_only_positive && (dist < 0)) {
        nonbonded_en = calc_non_bonded_pair_energy(p, &part_rep, ia_params, vec,
                                                   -1.0 * dist, dist * dist);
      }
    } else {
      runtimeErrorMsg() << "Constraint "
                        << " violated by particle " << p->p.identity;
    }
  }
  if (part_rep.p.type >= 0)
    *obsstat_nonbonded(&energy, p->p.type, part_rep.p.type) += nonbonded_en;
}
}
