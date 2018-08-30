/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CORE_UTILS_SERIALIZATION_FLAT_SET_HPP
#define CORE_UTILS_SERIALIZATION_FLAT_SET_HPP

#include <boost/container/flat_set.hpp>

namespace boost {
namespace serialization {
template <typename Archive, typename V, typename C>
void load(Archive &ar, boost::container::flat_set<V, C> &v,
          const unsigned int) {
  using value_type = typename boost::container::flat_set<V>::value_type;
  using size_type = typename boost::container::flat_set<V>::size_type;
  size_type count;

  ar >> count;
  v.reserve(count);

  for (; count > 0; --count) {
    value_type e;

    ar >> e;
    v.emplace_hint(v.end(), e);
  }
}

template <typename Archive, typename V, typename C>
void save(Archive &ar, boost::container::flat_set<V, C> const &v,
          const unsigned int) {
  typename boost::container::flat_set<V>::size_type count(v.size());

  ar << count;

  for (auto const &e : v) {
    ar << e;
  }
}

template <typename Archive, typename V, typename C>
void serialize(Archive &ar, boost::container::flat_set<V, C> &v,
               const unsigned int version) {
  split_free(ar, v, version);
}
} // namespace serialization
} // namespace boost

#endif
