/*
 * Copyright (C) 2025 The ESPResSo project
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

#include <boost/serialization/split_free.hpp>

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <variant>

namespace boost::serialization {

namespace detail {
template <std::size_t I, class Archive, class Variant>
void load_impl(Archive &ar, std::size_t index, Variant &obj) {
  if (index == I) {
    std::variant_alternative_t<I, Variant> opt{};
    ar >> opt;
    obj.template emplace<I>(std::move(opt));
  } else if constexpr (I + 1 < std::variant_size_v<Variant>) {
    load_impl<I + 1>(ar, index, obj);
  }
}
} // namespace detail

template <class Archive, class... Ts>
void save(Archive &ar, std::variant<Ts...> const &obj, unsigned const) {
  ar << obj.index();
  std::visit([&](const auto &value) { ar << value; }, obj);
}

template <class Archive, class... Ts>
void load(Archive &ar, std::variant<Ts...> &obj, unsigned const) {
  std::size_t index = 0;
  ar >> index;
  if (index >= std::variant_size_v<std::variant<Ts...>>) {
    throw std::domain_error("std::variant cannot be reloaded (type mismatch)");
  }
  detail::load_impl<0>(ar, index, obj);
}

template <class Archive, class... Ts>
void serialize(Archive &ar, std::variant<Ts...> &obj, unsigned const version) {
  split_free(ar, obj, version);
}

} // namespace boost::serialization
