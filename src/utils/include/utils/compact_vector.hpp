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

#ifndef UTILS_INCLUDE_UTILS_COMPACT_VECTOR_HPP
#define UTILS_INCLUDE_UTILS_COMPACT_VECTOR_HPP

#include "serialization/array.hpp"

#include <boost/container/vector.hpp>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/serialization/split_member.hpp>

#include <cstdint>

namespace Utils {
namespace detail {
using container_size_type = std::uint16_t;
using container_options = boost::container::vector_options<
    boost::container::stored_size<detail::container_size_type>>::type;
} // namespace detail

/**
 * Custom vector container optimized for size.
 * Allocate only 16 bits for the number of elements, and take only
 * <tt>16 + N * sizeof(T)</tt> bits during serialization in a binary archive.
 */
template <typename T>
class compact_vector // NOLINT(bugprone-exception-escape)
    : public boost::container::vector<T, boost::container::new_allocator<T>,
                                      detail::container_options> {
public:
  using boost::container::vector<T, boost::container::new_allocator<T>,
                                 detail::container_options>::vector;

  template <class Archive> void save(Archive &oa, unsigned int const) const {
    auto const size = static_cast<detail::container_size_type>(this->size());
    oa << size << boost::serialization::make_array(this->data(), this->size());
  }

  template <class Archive> void load(Archive &ia, unsigned int const) {
    auto size = static_cast<detail::container_size_type>(0u);
    ia >> size;
    this->resize(size);
    ia >> boost::serialization::make_array(this->data(), size);
  }

  template <class Archive>
  void serialize(Archive &ar, unsigned int const file_version) {
    boost::serialization::split_member(ar, *this, file_version);
  }
};

} // namespace Utils

UTILS_ARRAY_BOOST_MPI_T(Utils::compact_vector, 0)
UTILS_ARRAY_BOOST_CLASS(Utils::compact_vector, 0, object_serializable)
UTILS_ARRAY_BOOST_TRACK(Utils::compact_vector, 0, track_never)

#endif
