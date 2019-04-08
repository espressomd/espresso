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

#include "get_value.hpp"

#include <boost/variant/static_visitor.hpp>

namespace ScriptInterface {
/**
 * @brief Serialize a Variant into a Variant with type info.
 *
 * ObjectId values are flattened by the get_state function of
 * the ScriptObject they refer to.
 */
class Serializer : public boost::static_visitor<Variant> {
public:
  template <typename T> Variant operator()(T const &val) const {
    return std::vector<Variant>{{val}};
  }

  Variant operator()(ObjectId const &oid) const {
    auto so_ptr = get_value<std::shared_ptr<ScriptInterfaceBase>>(oid);
    if (so_ptr) {
      return std::vector<Variant>{{so_ptr->name(),
                                   static_cast<int>(so_ptr->policy()),
                                   so_ptr->get_state()}};
    }
    return std::vector<Variant>{None{}};
  }
};

/**
 * @brief Serialize a Variant into a Variant with type info.
 *
 * ObjectId values are flattened by the get_state function of
 * the ScriptObject they refer to.
 */
class UnSerializer : public boost::static_visitor<Variant> {
  std::vector<std::shared_ptr<ScriptInterfaceBase>> m_created_objects;

public:
  std::vector<std::shared_ptr<ScriptInterfaceBase>> const &
  created_objects() const {
    return m_created_objects;
  }

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
      auto so_ptr = ScriptInterfaceBase::make_shared(
          get<std::string>(val[0]),
          ScriptInterfaceBase::CreationPolicy(get<int>(val[1])), val[2]);
      /* Store a copy to keep the so alive. */
      m_created_objects.push_back(so_ptr);

      return so_ptr->id();
    }
    default: /* Error */
      throw std::runtime_error("Invalid format.");
    }
  }
};
} // namespace ScriptInterface
#endif
