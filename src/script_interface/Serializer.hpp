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
#include "PackedVariant.hpp"

namespace ScriptInterface {
/**
 * @brief Serialize a Variant into a Variant with type info.
 *
 * ObjectId values are flattened by the get_state function of
 * the ScriptObject they refer to.
 */
class Serializer : public recursive_visitor<Serializer, Variant, PackedVariant> {
public:
  using recursive_visitor<Serializer, Variant, PackedVariant>::operator();

  template <typename T> PackedVariant operator()(T const &val) const {
    return std::vector<PackedVariant>{{val}};
  }

  PackedVariant operator()(const ObjectHandle *so_ptr) const {
    if (so_ptr) {
      return std::vector<PackedVariant>{{so_ptr->name(), so_ptr->get_state()}};
    }
    return std::vector<PackedVariant>{None{}};
  }

  PackedVariant operator()(ObjectRef const &so_ptr) const {
    return this->operator()(so_ptr.get());
  }
};

/**
 * @brief Serialize a Variant into a Variant with type info.
 *
 * ObjectId values are flattened by the get_state function of
 * the ScriptObject they refer to.
 */
class UnSerializer : public boost::static_visitor<Variant> {
public:
  template <typename T> Variant operator()(const T& val) const {
    return val;
  }

  Variant operator()(std::vector<PackedVariant> const &val) const {
    using boost::get;
    switch (val.size()) {
    case 1: /* Normal value */
      return operator()(val[0]);
    case 2: /* Object value */
    {
      auto so_ptr = ObjectHandle::make_shared(get<std::string>(val.at(0)));
      so_ptr->set_state(val.at(1));

      return so_ptr;
    }
    default: /* Error */
      throw std::runtime_error("Invalid format.");
    }
  }
};
} // namespace ScriptInterface
#endif
