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

#ifndef ESPRESSO_SRC_CORE_ACTOR_VISITORS_HPP
#define ESPRESSO_SRC_CORE_ACTOR_VISITORS_HPP

#include "actor/traits.hpp"

#include <utils/demangle.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>

struct GetActorName : public boost::static_visitor<std::string> {
  template <typename T> auto operator()(std::shared_ptr<T> const &) const {
    return Utils::demangle<T>();
  }
};

/** @brief Get the symbol name of an actor. */
template <class Variant> auto get_actor_name(Variant const &variant) {
  return boost::apply_visitor(GetActorName(), variant);
}

/** @brief Get an actor of a specific type, recursively. */
template <typename Actor>
struct GetActorByType : public boost::static_visitor<std::shared_ptr<Actor>> {
private:
  template <typename T>
  static constexpr bool is_exact_match_v = std::is_same_v<T, Actor>;
  template <typename T>
  static constexpr bool is_layer_correction_v =
      traits::is_layer_correction<T>::value;

public:
  template <typename T, std::enable_if_t<is_exact_match_v<T>> * = nullptr>
  auto operator()(std::shared_ptr<T> const &obj) const {
    return obj;
  }

  template <typename T, std::enable_if_t<not is_exact_match_v<T> and
                                         is_layer_correction_v<T>> * = nullptr>
  auto operator()(std::shared_ptr<T> const &obj) const {
    return boost::apply_visitor(*this, obj->base_solver);
  }

  template <typename T,
            std::enable_if_t<not is_exact_match_v<T> and
                             not is_layer_correction_v<T>> * = nullptr>
  auto operator()(std::shared_ptr<T> const &) const {
    return std::shared_ptr<Actor>{nullptr};
  }
};

/** @brief Get an active actor of a specific type, recursively. */
template <typename Actor, typename Variant>
std::shared_ptr<Actor> get_actor_by_type(Variant const &variant) {
  return boost::apply_visitor(GetActorByType<Actor>(), variant);
}

template <typename Actor, typename Variant>
std::shared_ptr<Actor>
get_actor_by_type(boost::optional<Variant> const &optional) {
  return (optional) ? get_actor_by_type<Actor>(*optional) : nullptr;
}

/** @brief Check if an actor of a specific type is active, recursively. */
template <typename Actor>
struct HasActorOfType : public boost::static_visitor<bool> {
private:
  template <typename T>
  static constexpr bool is_exact_match_v = std::is_same_v<T, Actor>;
  template <typename T>
  static constexpr bool is_layer_correction_v =
      traits::is_layer_correction<T>::value;

public:
  template <typename T, std::enable_if_t<is_exact_match_v<T>> * = nullptr>
  auto operator()(std::shared_ptr<T> const &) const {
    return true;
  }

  template <typename T, std::enable_if_t<not is_exact_match_v<T> and
                                         is_layer_correction_v<T>> * = nullptr>
  auto operator()(std::shared_ptr<T> const &obj) const {
    return boost::apply_visitor(*this, obj->base_solver);
  }

  template <typename T,
            std::enable_if_t<not is_exact_match_v<T> and
                             not is_layer_correction_v<T>> * = nullptr>
  auto operator()(std::shared_ptr<T> const &) const {
    return false;
  }
};

/** @brief Check if an actor of a specific type is active, recursively. */
template <typename Actor, typename Variant>
auto has_actor_of_type(Variant const &variant) {
  return boost::apply_visitor(HasActorOfType<Actor>(), variant);
}

template <typename Actor, typename Variant>
auto has_actor_of_type(boost::optional<Variant> const &optional) {
  return (optional) ? has_actor_of_type<Actor>(*optional) : false;
}

/** Check whether two actors are identical by pointer. */
struct ActorEquality : public boost::static_visitor<bool> {
  template <
      typename T, typename U,
      typename std::enable_if_t<std::is_same_v<T, U>, std::nullptr_t> = nullptr>
  bool operator()(std::shared_ptr<T> const &lhs,
                  std::shared_ptr<U> const &rhs) const {
    return lhs == rhs;
  }

  template <typename T, typename U,
            typename std::enable_if_t<!std::is_same_v<T, U>, std::nullptr_t> =
                nullptr>
  bool operator()(std::shared_ptr<T> const &,
                  std::shared_ptr<U> const &) const {
    return false;
  }
};

/** @brief Check if an actor is already stored in an optional. */
template <typename T, class Variant>
bool is_already_stored(std::shared_ptr<T> const &actor,
                       boost::optional<Variant> const &active_actor) {
  auto const visitor = std::bind(ActorEquality(), actor, std::placeholders::_1);
  return active_actor and boost::apply_visitor(visitor, *active_actor);
}

#endif
