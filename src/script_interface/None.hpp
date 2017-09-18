/*
  Copyright (C) 2017 The ESPResSo project

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

#ifndef SCRIPT_INTERFACE_NONE_HPP
#define SCRIPT_INTERFACE_NONE_HPP

#include <boost/serialization/access.hpp>

namespace ScriptInterface {

/**
 * @brief Type to indicate no value in Variant.
 */
class None {
public:
  constexpr None() = default;
  constexpr None(std::nullptr_t) {}

  constexpr bool operator==(None const &) const { return true; }
  constexpr bool operator!=(None const &) const { return false; }

  constexpr operator bool() const { return false; }
private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &, long int /* version */) const {}
};

}

#endif
