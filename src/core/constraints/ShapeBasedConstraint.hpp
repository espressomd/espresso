/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
#include "particle_store/ParticleStore.hpp"

#include <shapes/NoWhere.hpp>
#include <shapes/Shape.hpp>

#include <utils/Vector.hpp>

#include <memory>
#include <utility>

// forward declaration
struct KokkosHandle;

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

  ~ShapeBasedConstraint() override {
    // Release part_rep's ParticleStore columns while Kokkos is still alive
    // (see the equivalent handling in ~CellStructure). No-op if never attached.
    m_part_rep_store.release_columns();
  }

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
  int &type() {
    ensure_part_rep_attached();
    return part_rep.type();
  }
  // v() returns a write-through proxy (not an lvalue reference) once velocity
  // moves into the ParticleStore columns, so expose value get/set instead of a
  // bound reference (mirrors the ParticleHandle proxy accessors).
  Utils::Vector3d velocity() const {
    ensure_part_rep_attached();
    return part_rep.v();
  }
  void set_velocity(Utils::Vector3d const &velocity) {
    ensure_part_rep_attached();
    part_rep.v() = velocity;
  }

  void set_type(int type);

  Utils::Vector3d total_force() const;
  double total_normal_force() const;
  void bind_system(std::shared_ptr<System::System const> const &system) {
    m_system = system;
  }

private:
  /** @brief Attach @c part_rep to its single-row store on first use.
   *  A @ref Particle is a view -- EVERY accessor reads the store column, so
   *  @c part_rep must be bound before its type/velocity/force is touched.
   *  Idempotent. Attached lazily (not in the constructor) because Kokkos may
   *  not be initialized then; every entry point that touches @c part_rep
   *  (set_type/type/velocity/get_ia_param/force/add_energy) calls this first,
   *  and all of them run after the System (hence Kokkos) exists. The mutated
   *  members are @c mutable so the const accessors can attach. Defined in the
   *  .cpp (needs @c ::kokkos_handle from communication.hpp). */
  void ensure_part_rep_attached() const;

  mutable Particle part_rep;
  /** Standalone store backing @c part_rep's columns.
   *  @c part_rep is a representative wall particle owned by the constraint,
   *  not by any cell structure, so it needs its own single-row store for the
   *  view accessors to work. Attached lazily by @c ensure_part_rep_attached.
   */
  mutable ParticleStore m_part_rep_store;
  /** Co-ownership of the Kokkos runtime (mirrors @ref CellStructure's
   *  @c m_kokkos_handle). @c m_part_rep_store holds Kokkos Views that must be
   *  destroyed before @c Kokkos::finalize(); if this constraint outlives the
   *  last CellStructure, holding a handle keeps the runtime alive until the
   *  destructor releases the columns. Captured in @c ensure_part_rep_attached,
   *  where Kokkos is guaranteed initialized. */
  mutable std::shared_ptr<KokkosHandle> m_kokkos_handle;
  std::shared_ptr<Shapes::Shape> m_shape;
  bool m_penetrable;
  bool m_only_positive;
  Utils::Vector3d m_local_force;
  double m_outer_normal_force;
  std::weak_ptr<System::System const> m_system;

  IA_parameters const &get_ia_param(int type) const;
};

} // namespace Constraints
