/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef SCRIPT_INTERFACE_VARIANT_HPP
#define SCRIPT_INTERFACE_VARIANT_HPP

#include <string>
#include <vector>

#include <boost/variant.hpp>

#include "Vector.hpp"
#include "utils/ObjectId.hpp"

namespace ScriptInterface {
class ScriptInterfaceBase;
using ObjectId = Utils::ObjectId<ScriptInterfaceBase>;

/**
 * @brief Possible types for parameters.
 */
typedef boost::make_recursive_variant<
    bool, int, double, std::string, std::vector<int>, std::vector<double>,
    Vector2d, Vector3d, ObjectId, std::vector<boost::recursive_variant_>>::type
    Variant;

enum class VariantType {
  BOOL = 0,
  INT,
  DOUBLE,
  STRING,
  INT_VECTOR,
  DOUBLE_VECTOR,
  VECTOR2D,
  VECTOR3D,
  OBJECTID,
  VECTOR
};

inline bool is_objectid(Variant const &v) {
  return (v.which() == static_cast<int>(VariantType::OBJECTID));
}

} /* namespace ScriptInterface */

#endif
