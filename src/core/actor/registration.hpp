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

#ifndef ESPRESSO_SRC_CORE_ACTOR_REGISTRATION_HPP
#define ESPRESSO_SRC_CORE_ACTOR_REGISTRATION_HPP

#include "actor/visitors.hpp"

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <memory>

/** @brief Register an actor in a thread-safe manner. */
template <typename Variant, typename T>
void add_actor(boost::optional<Variant> &active_actor,
               std::shared_ptr<T> const &actor, void (&on_actor_change)(),
               bool (&flag_all_reduce)(bool)) {
  auto const cleanup_if_any_rank_failed = [&](bool this_failed) {
    auto const any_failed = flag_all_reduce(this_failed);
    if (any_failed) {
      active_actor = boost::none;
      on_actor_change();
    }
  };
  try {
    active_actor = actor;
    actor->on_activation();
    on_actor_change();
    cleanup_if_any_rank_failed(false);
  } catch (...) {
    cleanup_if_any_rank_failed(true);
    throw;
  }
}

template <typename Variant, typename T>
void remove_actor(boost::optional<Variant> &active_actor,
                  std::shared_ptr<T> const &actor, void (&on_actor_change)()) {
  active_actor = boost::none;
  on_actor_change();
}

#endif
