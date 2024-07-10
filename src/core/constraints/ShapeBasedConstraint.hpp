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

#pragma once

#include "BoxGeometry.hpp"
#include "Constraint.hpp"
#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <shapes/NoWhere.hpp>
#include <shapes/Shape.hpp>

#include <utils/Vector.hpp>

#include <memory>
#include <utility>

namespace System {
class System;
}

namespace Constraints {

class ShapeBasedConstraint : public Constraint {
public:
  ShapeBasedConstraint()
      : part_rep{}, m_shape{std::make_shared<Shapes::NoWhere>()},
        m_penetrable{false}, m_only_positive{false}, m_local_force{},
        m_outer_normal_force{}, m_system{} {}

  void add_energy(const Particle &p, const Utils::Vector3d &folded_pos,
                  double time, Observable_stat &energy) const override;

  ParticleForce force(const Particle &p, const Utils::Vector3d &folded_pos,
                      double time) override;

  bool fits_in_box(Utils::Vector3d const &) const override { return true; }

  /** @brief Calculate the minimum distance between all particle pairs. */
  double min_dist(BoxGeometry const &box_geo,
                  ParticleRange const &particles) const;

  /** @brief Calculate distance from the constraint */
  void calc_dist(Utils::Vector3d const &pos, double &dist,
                 Utils::Vector3d &vec) const {
    m_shape->calculate_dist(pos, dist, vec);
  }

  void set_shape(std::shared_ptr<Shapes::Shape> const &shape) {
    m_shape = shape;
  }

  Shapes::Shape const &shape() const { return *m_shape; }

  void reset_force() override {
    m_local_force = Utils::Vector3d{0, 0, 0};
    m_outer_normal_force = 0.0;
  }

  bool &only_positive() { return m_only_positive; }
  bool &penetrable() { return m_penetrable; }
  int &type() { return part_rep.type(); }
  Utils::Vector3d &velocity() { return part_rep.v(); }

  void set_type(int type);

  Utils::Vector3d total_force() const;
  double total_normal_force() const;
  void bind_system(std::shared_ptr<System::System const> const &system) {
    m_system = system;
  }

private:
  Particle part_rep;
  std::shared_ptr<Shapes::Shape> m_shape;
  bool m_penetrable;
  bool m_only_positive;
  Utils::Vector3d m_local_force;
  double m_outer_normal_force;
  std::weak_ptr<System::System const> m_system;

  IA_parameters const &get_ia_param(int type) const;
};

} // namespace Constraints
