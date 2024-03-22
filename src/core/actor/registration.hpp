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

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <functional>
#include <memory>
#include <optional>

template <typename Variant, typename T>
void add_actor(boost::mpi::communicator const &comm,
               std::optional<Variant> &active_actor,
               std::shared_ptr<T> const &actor, void (&on_actor_change)()) {
  std::optional<Variant> other = actor;
  auto const cleanup_if_any_rank_failed = [&](bool failed) {
    if (boost::mpi::all_reduce(comm, failed, std::logical_or<>())) {
      active_actor.swap(other);
      on_actor_change();
    }
  };
  try {
    active_actor.swap(other);
    actor->on_activation();
    on_actor_change();
    cleanup_if_any_rank_failed(false);
  } catch (...) {
    cleanup_if_any_rank_failed(true);
    throw;
  }
}

#endif
