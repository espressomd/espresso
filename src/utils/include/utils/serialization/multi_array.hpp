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
#ifndef UTILS_SERIALIZATION_MULTI_ARRAY_HPP
#define UTILS_SERIALIZATION_MULTI_ARRAY_HPP

#include <boost/multi_array.hpp>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

namespace boost {
namespace serialization {

template <typename Archive, class T, std::size_t N, class Allocator>
void load(Archive &ar, boost::multi_array<T, N, Allocator> &marray, unsigned) {
  boost::array<std::size_t, N> shape;
  ar &make_array(shape.data(), N);

  marray.resize(shape);

  ar &make_array(marray.data(), marray.num_elements());
}

template <typename Archive, class T, std::size_t N, class Allocator>
void save(Archive &ar, const boost::multi_array<T, N, Allocator> &marray,
          unsigned) {
  ar &make_array(marray.shape(), marray.num_dimensions());

  ar &make_array(marray.data(), marray.num_elements());
}

template <typename Archive, class T, std::size_t N, class Allocator>
void serialize(Archive &ar, boost::multi_array<T, N, Allocator> &v,
               const unsigned int version) {
  split_free(ar, v, version);
}
} // namespace serialization
} // namespace boost

#endif
