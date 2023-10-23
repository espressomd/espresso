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

#pragma once

#include "packed_variant.hpp"

#include <boost/serialization/utility.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <utils/serialization/unordered_map.hpp>

#include <string>
#include <utility>
#include <vector>

namespace ScriptInterface {
/**
 * @brief State of an object ready for serialization.
 *
 * This specifies the internal serialization format and
 * should not be used outside of the class.
 */
struct ObjectState {
  std::string name;
  PackedMap params;
  std::vector<std::pair<ObjectId, std::string>> objects;
  std::string internal_state;

  template <class Archive> void serialize(Archive &ar, long int) {
    ar &name &params &objects &internal_state;
  }
};
} // namespace ScriptInterface
