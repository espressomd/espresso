/*
  Copyright (C) 2017 The ESPResSo project

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

#ifndef CORE_CONSTAINTS_SLIP_PLANE_HPP
#define CORE_CONSTAINTS_SLIP_PLANE_HPP

#include "Constraint.hpp"
#include "core/shapes/Wall.hpp"
#include "random.hpp"
#include "integrate.hpp"

#include <cmath>
#include <vector>

namespace Constraints {
class SlipPlane : public Constraint {
  Vector3d m_v;
  double m_gamma, m_kT, m_rcut;
  Shapes::Wall m_wall;
  std::vector<int> m_types;

public:
  SlipPlane() = default;

  Vector3d &v() { return m_v; }
  double &gamma() { return m_gamma; }
  double &kT() { return m_kT; }
  double &r_cut() { return m_rcut; }
  Shapes::Wall &wall() { return m_wall; }
  std::vector<int> &types() { return m_types; }

  void add_force(Particle *p, double *folded_pos) override {
    if (std::find(m_types.begin(), m_types.end(), p->p.type) == m_types.end()) {
      return;
    }

    double dist;
    double d[3];
    m_wall.calculate_dist(folded_pos, &dist, d);

    if (dist < m_rcut) {
      auto const omega = 1 - (dist / m_rcut);
      auto const omega_sqrt = std::sqrt(omega);

      auto const pre_diss = m_gamma / time_step;
      auto const pre_rand = std::sqrt(24. * m_kT * m_gamma / time_step);
      auto const v_scaled = time_step * m_v;
      auto const rand =
          Vector3d{d_random() - 0.5, d_random() - 0.5, d_random() - 0.5};

      auto const force = pre_diss * omega * (Vector3d{p->m.v} - v_scaled) +
                         pre_rand * omega_sqrt * rand;

      for (int i = 0; i < 3; i++) {
        p->f.f[i] += force[i];
      }
    }
  }
};
}

#endif
