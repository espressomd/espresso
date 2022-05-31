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
#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_COULOMB_REGISTRATION_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_COULOMB_REGISTRATION_HPP

#include "config.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics/coulomb.hpp"

#include "actor/registration.hpp"
#include "actor/visitors.hpp"

#include "event.hpp"

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <memory>
#include <stdexcept>
#include <type_traits>

namespace Coulomb {

template <typename T, std::enable_if_t<traits::is_solver<T>::value> * = nullptr>
void add_actor(std::shared_ptr<T> const &actor) {
  if (::electrostatics_actor) {
    auto const name = get_actor_name(*::electrostatics_actor);
    throw std::runtime_error("An electrostatics solver is already active (" +
                             name + ")");
  }
  add_actor(::electrostatics_actor, actor, ::on_coulomb_change,
            detail::flag_all_reduce);
}

template <typename T, std::enable_if_t<traits::is_solver<T>::value> * = nullptr>
void remove_actor(std::shared_ptr<T> const &actor) {
  if (not is_already_stored(actor, electrostatics_actor)) {
    throw std::runtime_error(
        "The given electrostatics solver is not currently active");
  }
  remove_actor(::electrostatics_actor, actor, ::on_coulomb_change);
}

template <typename T,
          std::enable_if_t<traits::is_extension<T>::value> * = nullptr>
void add_actor(std::shared_ptr<T> const &actor) {
  if (::electrostatics_extension) {
    auto const name = get_actor_name(*::electrostatics_extension);
    throw std::runtime_error("An electrostatics extension is already active (" +
                             name + ")");
  }
  add_actor(::electrostatics_extension, actor, ::on_coulomb_change,
            detail::flag_all_reduce);
}

template <typename T,
          std::enable_if_t<traits::is_extension<T>::value> * = nullptr>
void remove_actor(std::shared_ptr<T> const &actor) {
  if (not is_already_stored(actor, ::electrostatics_extension)) {
    throw std::runtime_error(
        "The given electrostatics extension is not currently active");
  }
  remove_actor(::electrostatics_extension, actor, ::on_coulomb_change);
}

} // namespace Coulomb

#endif // ELECTROSTATICS
#endif
