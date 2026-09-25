/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#include "BondList.hpp"
#include "cell_system/CellStructure.hpp"
#include "particle_node.hpp"
#include "particle_store/VectorReference.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <optional>
#include <vector>

/** Fixture to create particles during a test and remove them at the end. */
struct ParticleFactory {
  ParticleFactory() = default;

  ~ParticleFactory() { clear_particles(); }

  void create_particle(Utils::Vector3d const &pos, int p_id, int type) {
    ::make_new_particle(p_id, pos);
    set_particle_property(p_id, &Particle::type, type);
    particle_cache.emplace_back(p_id);
  }

  void set_particle_type(int p_id, int type) const {
    set_particle_property(p_id, &Particle::type, type);
  }

  void set_particle_v(int p_id, Utils::Vector3d const &vel) const {
    set_particle_property(p_id, &Particle::v, vel);
  }

  void insert_particle_bond(int p_id, int bond_id,
                            std::vector<int> const &partner_ids) const {
    auto &system = System::get_system();
    // get_local_particle resolves an id to a store row and needs a
    // synchronized store; a preceding create/move left the store dirty.
    system.cell_structure->ensure_particle_store_synchronized();
    auto p = system.cell_structure->get_local_particle(p_id);
    if (p and not p->is_ghost()) {
      p->bonds().insert(BondView(bond_id, partner_ids));
    }
    system.on_particle_change();
  }

  template <typename T>
  void set_particle_property(int p_id, T &(Particle::*getter)(),
                             T const &value) const {
    auto &system = System::get_system();
    // get_local_particle resolves an id to a store row and needs a
    // synchronized store; a preceding create/move left the store dirty.
    system.cell_structure->ensure_particle_store_synchronized();
    auto p = system.cell_structure->get_local_particle(p_id);
    if (p) {
      ((*p).*getter)() = value;
    }
  }

  /** Overload for accessors returning a write-through proxy by value
   *  (force/torque live in the ParticleStore columns). */
  void set_particle_property(int p_id, VectorReference (Particle::*getter)(),
                             Utils::Vector3d const &value) const {
    auto &system = System::get_system();
    system.cell_structure->ensure_particle_store_synchronized();
    auto p = system.cell_structure->get_local_particle(p_id);
    if (p) {
      ((*p).*getter)() = value;
    }
  }

  /** Overload for the integer-column proxy (image box). */
  void set_particle_property(int p_id,
                             IntegerVectorReference (Particle::*getter)(),
                             Utils::Vector3i const &value) const {
    auto &system = System::get_system();
    system.cell_structure->ensure_particle_store_synchronized();
    auto p = system.cell_structure->get_local_particle(p_id);
    if (p) {
      ((*p).*getter)() = value;
    }
  }

  template <typename T>
  auto get_particle_property(int p_id, T &(Particle::*getter)()) const {
    auto &system = System::get_system();
    // get_local_particle resolves an id to a store row and needs a
    // synchronized store; a preceding create/move left the store dirty.
    system.cell_structure->ensure_particle_store_synchronized();
    auto p = system.cell_structure->get_local_particle(p_id);
    std::optional<T> opt = std::nullopt;
    if (p) {
      opt = ((*p).*getter)();
    }
    return opt;
  }

  /** Overload for accessors returning a write-through proxy by value
   *  (force/torque live in the ParticleStore columns). Ensures the store is
   *  synchronized so the returned proxy references a valid row, and converts
   *  to a plain @ref Utils::Vector3d. */
  auto get_particle_property(int p_id,
                             VectorReference (Particle::*getter)()) const {
    auto &system = System::get_system();
    system.cell_structure->ensure_particle_store_synchronized();
    auto p = system.cell_structure->get_local_particle(p_id);
    std::optional<Utils::Vector3d> opt = std::nullopt;
    if (p) {
      opt = Utils::Vector3d(((*p).*getter)());
    }
    return opt;
  }

  void clear_particles() {
    for (auto pid : particle_cache) {
      remove_particle(pid);
    }
    particle_cache.clear();
  }

private:
  std::vector<int> particle_cache;
};
