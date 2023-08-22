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

#pragma once

#include <boost/serialization/split_free.hpp>
#include <boost/serialization/utility.hpp>

#include <unordered_map>

namespace boost::serialization {

template <typename Archive, typename K, typename V>
void load(Archive &ar, std::unordered_map<K, V> &map, unsigned const int) {
  using value_type = typename std::unordered_map<K, V>::value_type;
  using size_type = typename std::unordered_map<K, V>::size_type;

  size_type count;
  ar >> count;
  map.reserve(count);

  value_type pair{};
  for (size_type i = 0; i < count; i++) {
    ar >> pair;
    map.emplace_hint(map.end(), pair);
  }
}

template <typename Archive, typename K, typename V>
void save(Archive &ar, std::unordered_map<K, V> const &map,
          unsigned const int) {
  auto const count = map.size();
  ar << count;

  for (auto const &pair : map) {
    ar << pair;
  }
}

template <typename Archive, typename K, typename V>
void serialize(Archive &ar, std::unordered_map<K, V> &map,
               unsigned int const version) {
  split_free(ar, map, version);
}

} // namespace boost::serialization
