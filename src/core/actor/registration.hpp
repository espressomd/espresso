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

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <functional>
#include <memory>
#include <optional>
#include <variant>

namespace System {
class System;
}

template <typename Variant, typename T, class F>
void add_actor(boost::mpi::communicator const &comm,
               std::shared_ptr<System::System> const &system,
               std::optional<Variant> &active_actor,
               std::shared_ptr<T> const &actor, F &&on_actor_change) {
  std::optional<Variant> other = actor;
  auto const activate = [&system](auto &leaf) {
    leaf->bind_system(system);
    leaf->on_activation();
  };
  auto const deactivate = [&system](auto &leaf) {
    leaf->detach_system(system);
  };
  auto const cleanup_if_any_rank_failed = [&](bool failed) {
    if (boost::mpi::all_reduce(comm, failed, std::logical_or<>())) {
      deactivate(actor);
      active_actor.swap(other);
      if (active_actor) {
        std::visit([&](auto &leaf) { activate(leaf); }, *active_actor);
      }
      on_actor_change();
    }
  };
  try {
    active_actor.swap(other);
    if (other) {
      std::visit([&](auto &leaf) { deactivate(leaf); }, *other);
    }
    activate(actor);
    on_actor_change();
    cleanup_if_any_rank_failed(false);
  } catch (...) {
    cleanup_if_any_rank_failed(true);
    throw;
  }
}
