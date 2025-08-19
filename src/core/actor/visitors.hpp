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

#include "actor/traits.hpp"

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <variant>

/** @brief Get an actor of a specific type, recursively. */
template <typename Actor> struct GetActorByType {
  template <typename T> auto operator()(std::shared_ptr<T> const &obj) const {
    if constexpr (std::is_same_v<T, Actor>) {
      return obj;
    }
    if constexpr (traits::is_layer_correction<T>::value) {
      return std::visit(*this, obj->base_solver);
    }
    return std::shared_ptr<Actor>{nullptr};
  }
};

/** @brief Get an active actor of a specific type, recursively. */
template <typename Actor, typename Variant>
std::shared_ptr<Actor> get_actor_by_type(Variant const &variant) {
  return std::visit(GetActorByType<Actor>(), variant);
}

template <typename Actor, typename Variant>
std::shared_ptr<Actor>
get_actor_by_type(std::optional<Variant> const &optional) {
  return (optional) ? get_actor_by_type<Actor>(*optional) : nullptr;
}

/** @brief Check if an actor of a specific type is active, recursively. */
template <typename Actor> struct HasActorOfType {
  template <typename T> auto operator()(std::shared_ptr<T> const &obj) const {
    if constexpr (std::is_same_v<T, Actor>) {
      return true;
    }
    if constexpr (traits::is_layer_correction<T>::value) {
      return std::visit(*this, obj->base_solver);
    }
    return false;
  }
};

/** @brief Check if an actor of a specific type is active, recursively. */
template <typename Actor, typename Variant>
auto has_actor_of_type(Variant const &variant) {
  return std::visit(HasActorOfType<Actor>(), variant);
}

template <typename Actor, typename Variant>
auto has_actor_of_type(std::optional<Variant> const &optional) {
  return (optional) ? has_actor_of_type<Actor>(*optional) : false;
}
