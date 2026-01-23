/*
 * Copyright (C) 2022 The ESPResSo project
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

#include "cell_system/CellStructure.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/interactions/BondedInteractions.hpp"

#include "core/Particle.hpp"
#include "core/exclusions.hpp"
#include "core/system/System.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>

#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace CellSystem {
class CellSystem;
}
namespace Particles {

class ParticleModifier;

inline auto error_msg(std::string const &name, std::string const &reason) {
  std::stringstream msg;
  msg << "attribute '" << name << "' of 'ParticleHandle' " << reason;
  return msg.str();
}

inline auto get_real_particle(boost::mpi::communicator const &comm, int p_id,
                              ::CellStructure &cell_structure) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }
  auto ptr = cell_structure.get_local_particle(p_id);
  if (ptr != nullptr and ptr->is_ghost()) {
    ptr = nullptr;
  }
  auto const n_found = boost::mpi::all_reduce(
      comm, static_cast<int>(ptr != nullptr), std::plus<>());
  if (n_found == 0) {
    throw std::runtime_error("Particle with id " + std::to_string(p_id) +
                             " not found");
  }
  return ptr;
}

inline void particle_checks(int p_id, Utils::Vector3d const &pos) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }
#ifndef __FAST_MATH__
  if (std::isnan(pos[0]) or std::isnan(pos[1]) or std::isnan(pos[2]) or
      std::isinf(pos[0]) or std::isinf(pos[1]) or std::isinf(pos[2])) {
    throw std::domain_error("Particle position must be finite");
  }
#endif // __FAST_MATH__
}

#ifdef ESPRESSO_EXCLUSIONS
/**
 * @brief Locally add an exclusion to a particle.
 * @param pid1 the identity of the first exclusion partner
 * @param pid2 the identity of the second exclusion partner
 * @param cell_structure the cell structure
 */
inline void local_add_exclusion(int pid1, int pid2,
                                ::CellStructure &cell_structure) {
  if (auto p1 = cell_structure.get_local_particle(pid1)) {
    add_exclusion(*p1, pid2);
  }
  if (auto p2 = cell_structure.get_local_particle(pid2)) {
    add_exclusion(*p2, pid1);
  }
}

/**
 * @brief Locally remove an exclusion to a particle.
 * @param pid1 the identity of the first exclusion partner
 * @param pid2 the identity of the second exclusion partner
 * @param cell_structure the cell structure
 */
inline void local_remove_exclusion(int pid1, int pid2,
                                   ::CellStructure &cell_structure) {
  if (auto p1 = cell_structure.get_local_particle(pid1)) {
    delete_exclusion(*p1, pid2);
  }
  if (auto p2 = cell_structure.get_local_particle(pid2)) {
    delete_exclusion(*p2, pid1);
  }
}

inline void particle_exclusion_sanity_checks(int pid1, int pid2,
                                             ::CellStructure &cell_structure,
                                             auto const &comm) {
  if (pid1 == pid2) {
    throw std::runtime_error("Particles cannot exclude themselves (id " +
                             std::to_string(pid1) + ")");
  }
  std::ignore = get_real_particle(comm, pid1, cell_structure);
  std::ignore = get_real_particle(comm, pid2, cell_structure);
}
#endif // ESPRESSO_EXCLUSIONS

class ParticleHandle : public AutoParameters<ParticleHandle> {
  std::function<Variant(VariantMap const &)> cb_get_bond;
  int m_pid;
  mutable std::weak_ptr<CellSystem::CellSystem> m_cell_structure;
  mutable std::weak_ptr<Interactions::BondedInteractions> m_bonded_ias;
  mutable std::weak_ptr<::System::System> m_system;
  auto get_cell_structure() const {
    auto ptr = m_cell_structure.lock();
    assert(ptr != nullptr);
    return ptr;
  }
  auto get_bonded_ias() const {
    auto ptr = m_bonded_ias.lock();
    assert(ptr != nullptr);
    return ptr;
  }
  auto get_system() const {
    auto ptr = m_system.lock();
    assert(ptr != nullptr);
    return ptr;
  }

  template <typename T>
  T get_particle_property(T const &(Particle::*getter)() const) const;
  template <typename T, class F> T get_particle_property(F const &fun) const;

  template <typename T>
  void set_particle_property(T &(Particle::*setter)(),
                             Variant const &value) const;

  template <class F> void set_particle_property(F const &fun) const;

  std::size_t setup_hidden_args(VariantMap const &params);

public:
  ParticleHandle();
  friend class ParticleModifier;

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  void do_construct(VariantMap const &params) override;

  void attach(std::weak_ptr<::System::System> system) {
    assert(m_system.expired());
    m_system = system;
  }
};

/** @brief Thin wrapper to read and write particle attributes in batches. */
class ParticleModifier : public ParticleHandle {
public:
  void do_construct(VariantMap const &params) override {
    auto const n_extra_args = setup_hidden_args(params);
    if (n_extra_args) {
      throw std::logic_error(
          "ParticleModifier is not meant to create new particles");
    }
    m_pid = params.contains("id") ? get_value<int>(params, "id") : -1;
  }

  void set_pid(int pid) {
    assert(pid >= 0);
    m_pid = pid;
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "remove_particle") {
      throw std::logic_error("ParticleModifier has no method '" + name + "'");
    }
    return ParticleHandle::do_call_method(name, params);
  }
};

} // namespace Particles
} // namespace ScriptInterface
