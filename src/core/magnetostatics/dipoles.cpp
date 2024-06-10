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

#include "config/config.hpp"

#include "magnetostatics/solver.hpp"

#ifdef DIPOLES

#include "magnetostatics/dipoles.hpp"

#include "ParticleRange.hpp"
#include "actor/traits.hpp"
#include "actor/visit_try_catch.hpp"
#include "actor/visitors.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/demangle.hpp>

#include <cassert>
#include <optional>
#include <stdexcept>

namespace Dipoles {

Solver::Solver() {
  impl = std::make_unique<Implementation>();
  reinit_on_observable_calc = false;
}

Solver const &get_dipoles() { return System::get_system().dipoles; }

void Solver::sanity_checks() const {
  if (impl->solver) {
    std::visit([](auto &ptr) { ptr->sanity_checks(); }, *impl->solver);
  }
}

void Solver::on_dipoles_change() {
  reinit_on_observable_calc = true;
  if (impl->solver) {
    visit_try_catch([](auto &ptr) { ptr->init(); }, *impl->solver);
  }
}

void Solver::on_boxl_change() {
  if (impl->solver) {
    visit_try_catch([](auto &ptr) { ptr->on_boxl_change(); }, *impl->solver);
  }
}

void Solver::on_node_grid_change() {
  if (impl->solver) {
    std::visit([](auto &ptr) { ptr->on_node_grid_change(); }, *impl->solver);
  }
}

void Solver::on_periodicity_change() {
  if (impl->solver) {
    visit_try_catch([](auto &ptr) { ptr->on_periodicity_change(); },
                    *impl->solver);
  }
}

void Solver::on_cell_structure_change() {
  if (impl->solver) {
    visit_try_catch([](auto &ptr) { ptr->on_cell_structure_change(); },
                    *impl->solver);
  }
}

double Solver::cutoff() const {
#ifdef DP3M
  if (impl->solver) {
    if (auto dp3m = get_actor_by_type<DipolarP3M>(impl->solver)) {
      return dp3m->dp3m.params.r_cut;
    }
  }
#endif
  return -1.;
}

void Solver::on_observable_calc() {
  if (reinit_on_observable_calc) {
#ifdef DP3M
    if (impl->solver) {
      if (auto dp3m = get_actor_by_type<DipolarP3M>(impl->solver)) {
        dp3m->count_magnetic_particles();
      }
    }
#endif
    reinit_on_observable_calc = false;
  }
}

struct LongRangeForce {
  ParticleRange const &m_particles;
  explicit LongRangeForce(ParticleRange const &particles)
      : m_particles(particles) {}

#ifdef DP3M
  void operator()(std::shared_ptr<DipolarP3M> const &actor) const {
    actor->add_long_range_forces(m_particles);
  }
#endif // DP3M
  void operator()(std::shared_ptr<DipolarLayerCorrection> const &actor) const {
    actor->add_force_corrections(m_particles);
    std::visit(*this, actor->base_solver);
  }
  void operator()(std::shared_ptr<DipolarDirectSum> const &actor) const {
    actor->add_long_range_forces(m_particles);
  }
#ifdef DIPOLAR_DIRECT_SUM
  void operator()(std::shared_ptr<DipolarDirectSumGpu> const &actor) const {
    actor->add_long_range_forces();
  }
#endif
#ifdef SCAFACOS_DIPOLES
  void operator()(std::shared_ptr<DipolarScafacos> const &actor) const {
    actor->add_long_range_forces();
  }
#endif
};

struct LongRangeEnergy {
  ParticleRange const &m_particles;
  explicit LongRangeEnergy(ParticleRange const &particles)
      : m_particles(particles) {}

#ifdef DP3M
  double operator()(std::shared_ptr<DipolarP3M> const &actor) const {
    return actor->long_range_energy(m_particles);
  }
#endif // DP3M
  double
  operator()(std::shared_ptr<DipolarLayerCorrection> const &actor) const {
    auto energy = std::visit(*this, actor->base_solver);
    return energy + actor->energy_correction(m_particles);
  }
  double operator()(std::shared_ptr<DipolarDirectSum> const &actor) const {
    return actor->long_range_energy(m_particles);
  }
#ifdef DIPOLAR_DIRECT_SUM
  double operator()(std::shared_ptr<DipolarDirectSumGpu> const &actor) const {
    actor->long_range_energy();
    return 0.;
  }
#endif
#ifdef SCAFACOS_DIPOLES
  double operator()(std::shared_ptr<DipolarScafacos> const &actor) const {
    return actor->long_range_energy();
  }
#endif
};

#ifdef DIPOLE_FIELD_TRACKING
struct LongRangeField {
  ParticleRange const &m_particles;
  explicit LongRangeField(ParticleRange const &particles)
      : m_particles(particles) {}

  void operator()(std::shared_ptr<DipolarDirectSum> const &actor) const {
    actor->dipole_field_at_part(m_particles);
  }

  template <typename T,
            std::enable_if_t<!traits::has_dipole_fields<T>::value> * = nullptr>
  void operator()(std::shared_ptr<T> const &) const {
    runtimeErrorMsg() << "Dipoles field calculation not implemented by "
                      << "dipolar method " << Utils::demangle<T>();
  }
};
#endif

void Solver::calc_pressure_long_range() const {
  if (impl->solver) {
    runtimeWarningMsg() << "pressure calculated, but pressure not implemented.";
  }
}

void Solver::calc_long_range_force(ParticleRange const &particles) const {
  if (impl->solver) {
    std::visit(LongRangeForce(particles), *impl->solver);
  }
}

double Solver::calc_energy_long_range(ParticleRange const &particles) const {
  if (impl->solver) {
    return std::visit(LongRangeEnergy(particles), *impl->solver);
  }
  return 0.;
}

#ifdef DIPOLE_FIELD_TRACKING
void Solver::calc_long_range_field(ParticleRange const &particles) const {
  if (impl->solver) {
    std::visit(LongRangeField(particles), *impl->solver);
  }
}
#endif

} // namespace Dipoles
#endif // DIPOLES
