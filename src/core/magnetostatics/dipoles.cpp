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

#include <config/config.hpp>

#ifdef ESPRESSO_DIPOLES

#include "magnetostatics/dipoles.hpp"

#include "actor/traits.hpp"
#include "actor/visit_try_catch.hpp"
#include "actor/visitors.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "magnetostatics/solver.hpp"

#include <utils/demangle.hpp>

#include <cassert>
#include <optional>
#include <stdexcept>

namespace Dipoles {

Solver::Solver() {
  impl = std::make_unique<Implementation>();
  reinit_on_observable_calc = false;
}

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

struct LongRangePressure {
  auto operator()(std::shared_ptr<DipolarDirectSum> const &actor) const {
    return actor->long_range_pressure();
  }
#ifdef ESPRESSO_DP3M
  auto operator()(std::shared_ptr<DipolarP3M> const &actor) const {
    return actor->long_range_pressure();
  }
#endif
  template <typename T>
    requires(not traits::has_pressure<T>::value)
  auto operator()(std::shared_ptr<T> const &) const {
    runtimeWarningMsg() << "Pressure not implemented for this dipolar method.";
    return Utils::Vector9d{};
  }
};

double Solver::cutoff() const {
#ifdef ESPRESSO_DP3M
  if (impl->solver) {
    if (auto dp3m = get_actor_by_type<DipolarP3M>(impl->solver)) {
      return dp3m->dp3m_params.r_cut;
    }
  }
#endif
  return inactive_cutoff;
}

void Solver::on_observable_calc() {
  if (reinit_on_observable_calc) {
#ifdef ESPRESSO_DP3M
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
  template <class Solver>
  void operator()(std::shared_ptr<Solver> const &actor) const {
    actor->add_long_range_forces();
  }
};

struct LongRangeEnergy {
  template <class Solver>
  double operator()(std::shared_ptr<Solver> const &actor) const {
    return actor->long_range_energy();
  }
};

Utils::Vector9d Solver::calc_pressure_long_range() const {
  if (impl->solver) {
    return std::visit(LongRangePressure{}, *impl->solver);
  }
  return {};
}

void Solver::calc_long_range_force() const {
  if (impl->solver) {
    std::visit(LongRangeForce{}, *impl->solver);
  }
}

double Solver::calc_energy_long_range() const {
  if (impl->solver) {
    return std::visit(LongRangeEnergy{}, *impl->solver);
  }
  return 0.;
}

} // namespace Dipoles
#endif // ESPRESSO_DIPOLES
