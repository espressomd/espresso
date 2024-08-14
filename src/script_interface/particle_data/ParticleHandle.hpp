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

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/interactions/BondedInteractions.hpp"

#include "core/Particle.hpp"
#include "core/system/System.hpp"

#include <cassert>
#include <functional>
#include <memory>
#include <string>

namespace ScriptInterface {
namespace CellSystem {
class CellSystem;
}
namespace Particles {

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
#ifdef EXCLUSIONS
  void particle_exclusion_sanity_checks(int pid1, int pid2) const;
#endif // EXCLUSIONS

public:
  ParticleHandle();

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  void do_construct(VariantMap const &params) override;

  void attach(std::weak_ptr<::System::System> system) {
    assert(m_system.expired());
    m_system = system;
  }
};

} // namespace Particles
} // namespace ScriptInterface
