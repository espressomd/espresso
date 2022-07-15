/*
 * Copyright (C) 2020-2022 The ESPResSo project
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
#ifndef UTILS_SERIALIZATION_PACK_HPP
#define UTILS_SERIALIZATION_PACK_HPP

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include <sstream>
#include <string>

namespace Utils {
/**
 * @brief Pack a serialize type into a string.
 *
 * @tparam T Serializable type
 * @param v Value to serialize
 * @return String representation
 */
template <class T> std::string pack(T const &v) {
  std::stringstream ss;
  boost::archive::binary_oarchive(ss) << v;

  return ss.str();
}

/**
 * @brief Unpack a serialize type into a string.
 *
 * @tparam T Serializable type
 * @param state String to construct the value from, as returned by @ref
 * Utils::pack.
 * @return Unpacked value
 */
template <class T> T unpack(std::string const &state) {
  namespace iostreams = boost::iostreams;

  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);

  T val;
  boost::archive::binary_iarchive(ss) >> val;

  return val;
}
} // namespace Utils

#endif // UTILS_SERIALIZATION_PACK_HPP
