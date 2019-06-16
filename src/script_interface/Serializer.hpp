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
#ifndef SCRIPT_INTERFACE_SERIALIZER_HPP
#define SCRIPT_INTERFACE_SERIALIZER_HPP

#include "ObjectHandle.hpp"

namespace ScriptInterface {
/**
 * @brief Serialize a Variant into a Variant with type info.
 *
 * ObjectId values are flattened by the get_state function of
 * the ScriptObject they refer to.
 */
class Serializer : public recursive_visitor<Serializer, Variant, Variant> {
public:
  template <typename T> Variant operator()(T const &val) const {
    return std::vector<Variant>{{val}};
  }

  Variant operator()(const ObjectHandle * so_ptr) const {
    if (so_ptr) {
      return std::vector<Variant>{{so_ptr->name(),
                                      static_cast<int>(so_ptr->policy()),
                                      so_ptr->get_state()}};
    }
    return std::vector<Variant>{None{}};
  }

  Variant operator()(ObjectRef const &so_ptr) const {
    return this->operator()(so_ptr.get());
  }
};

/**
 * @brief Serialize a Variant into a Variant with type info.
 *
 * ObjectId values are flattened by the get_state function of
 * the ScriptObject they refer to.
 */
class UnSerializer : public recursive_visitor<UnSerializer, Variant, Variant> {
public:
  template <typename T> Variant operator()(T const & /* val */) {
    throw std::runtime_error("Invalid format.");
  }

  Variant operator()(std::vector<Variant> const &val) {
    using boost::get;
    switch (val.size()) {
    case 1: /* Normal value */
      return val[0];
    case 3: /* Object value */
    {
      return ObjectHandle::make_shared(
          get<std::string>(val[0]),
          ObjectHandle::CreationPolicy(get<int>(val[1])));
      // , val[2]
    }
    default: /* Error */
      throw std::runtime_error("Invalid format.");
    }
  }
};
} // namespace ScriptInterface
#endif
