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
#ifndef CORE_UTILS_SERIALIZATION_LIST_HPP
#define CORE_UTILS_SERIALIZATION_LIST_HPP

#include <boost/serialization/array.hpp>
#if BOOST_VERSION >= 106400 && BOOST_VERSION < 106500
#include <boost/serialization/array_wrapper.hpp>
#endif
#include <boost/serialization/split_free.hpp>

#include "utils/List.hpp"

namespace boost {
namespace serialization {
template <typename T, class Archive>
void load(Archive &ar, Utils::List<T> &v, const unsigned int /*file_version*/) {
  typename Utils::List<T>::size_type n;
  ar >> n;
  v.resize(n);

  ar >> make_array(v.data(), n);
}

template <typename T, class Archive>
void save(Archive &ar, Utils::List<T> const &v,
          const unsigned int /*file_version*/) {
  typename Utils::List<T>::size_type n = v.size();
  ar << n;
  ar << make_array(v.data(), v.size());
}

template <typename T, class Archive>
void serialize(Archive &ar, Utils::List<T> &v,
               const unsigned int file_version) {
  split_free(ar, v, file_version);
}
} // namespace serialization
} // namespace boost

#endif
