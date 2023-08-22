/*
 * Copyright (C) 2023 The ESPResSo project
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

#include <optional>

namespace boost::serialization {

template <class Archive, class T>
void save(Archive &ar, std::optional<T> const &opt, unsigned int const) {
  bool has_value = opt.has_value();
  ar << has_value;
  if (has_value) {
    ar << *opt;
  }
}

template <class Archive, class T>
void load(Archive &ar, std::optional<T> &opt, unsigned int const) {
  bool has_value;
  ar >> has_value;
  if (has_value) {
    if (not opt.has_value()) {
      opt = T{}; // optional must be initialized with something
    }
    ar >> *opt;
  } else {
    opt.reset();
  }
}

template <class Archive, class T>
void serialize(Archive &ar, std::optional<T> &opt, unsigned int const version) {
  split_free(ar, opt, version);
}

} // namespace boost::serialization
