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
#ifndef ESPRESSO_SRC_CORE_MAGNETOSTATICS_DIPOLE_REGISTRATION_HPP
#define ESPRESSO_SRC_CORE_MAGNETOSTATICS_DIPOLE_REGISTRATION_HPP

#include "config/config.hpp"

#ifdef DIPOLES

#include "magnetostatics/dipoles.hpp"

#include "actor/registration.hpp"
#include "actor/visitors.hpp"

#include "event.hpp"
#include "system/System.hpp"

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <memory>
#include <stdexcept>
#include <type_traits>

namespace Dipoles {

template <typename T, std::enable_if_t<traits::is_solver<T>::value> * = nullptr>
void add_actor(std::shared_ptr<T> const &actor) {
  auto &system = System::get_system();
  if (system.dipoles.solver) {
    auto const name = get_actor_name(*system.dipoles.solver);
    throw std::runtime_error("A magnetostatics solver is already active (" +
                             name + ")");
  }
  add_actor(system.dipoles.solver, actor, ::on_dipoles_change,
            detail::flag_all_reduce);
}

template <typename T, std::enable_if_t<traits::is_solver<T>::value> * = nullptr>
void remove_actor(std::shared_ptr<T> const &actor) {
  auto &system = System::get_system();
  if (not is_already_stored(actor, system.dipoles.solver)) {
    throw std::runtime_error(
        "The given magnetostatics solver is not currently active");
  }
  remove_actor(system.dipoles.solver, actor, ::on_dipoles_change);
}

} // namespace Dipoles

#endif // DIPOLES
#endif
