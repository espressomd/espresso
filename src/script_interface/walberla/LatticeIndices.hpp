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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_LATTICE_INDICES_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_LATTICE_INDICES_HPP

#include "script_interface/ScriptInterface.hpp"

#include <utils/Vector.hpp>

#include <sstream>
#include <stdexcept>

namespace ScriptInterface {

/** @brief Interface to carry out simple operations on lattice indices. */
class LatticeIndices : public ObjectHandle {
protected:
  [[nodiscard]] bool is_index_valid(Utils::Vector3i const &index,
                                    Utils::Vector3i const &shape) const {
    return index < shape and index >= Utils::Vector3i{};
  }

  void throw_invalid_index(Utils::Vector3i const &index,
                           Utils::Vector3i const &shape) const {
    if (context()->is_head_node()) {
      auto constexpr formatter = Utils::Vector3i::formatter(", ");
      std::stringstream ss;
      ss << "provided index [" << formatter << index << "] is out of range "
         << "for shape [" << formatter << shape << "]";
      throw std::out_of_range(ss.str());
    }
    throw Exception("");
  }

  [[nodiscard]] Utils::Vector3i
  get_mapped_index(Utils::Vector3i const &index,
                   Utils::Vector3i const &shape) const {
    auto output = index;
    for (auto i : {0, 1, 2}) {
      if (output[i] < 0) {
        output[i] += shape[i];
      }
    }
    if (not is_index_valid(output, shape)) {
      throw_invalid_index(index, shape);
    }
    return output;
  }
};

} // namespace ScriptInterface

#endif
