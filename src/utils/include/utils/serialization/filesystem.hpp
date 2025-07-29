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
#include <boost/serialization/utility.hpp>

#include <filesystem>
#include <string>

namespace boost::serialization {

template <typename Archive>
void load(Archive &ar, std::filesystem::path &path, unsigned int const) {
  std::string source;
  ar >> source;
  path = std::filesystem::path(source, std::filesystem::path::generic_format);
}

template <typename Archive>
void save(Archive &ar, std::filesystem::path const &path, unsigned int const) {
  ar << path.generic_string();
}

template <typename Archive>
void serialize(Archive &ar, std::filesystem::path &path,
               unsigned int const version) {
  split_free(ar, path, version);
}

} // namespace boost::serialization
