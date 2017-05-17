/*
  Copyright (C) 2016,2017 The ESPResSo project

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

#ifndef SCRIPT_INTERFACE_GET_VALUE_HPP
#define SCRIPT_INTERFACE_GET_VALUE_HPP

#include "ScriptInterfaceBase.hpp"
#include "Variant.hpp"

namespace ScriptInterface {

namespace detail {
/**
 * @brief Implementation of get_value.
 *
 * Helper struct is needed because partial specialization of functions
 * is not allowed.
 */
template <typename T, typename = void> struct get_value_helper {
  T operator()(Variant const &v) const { return boost::get<T>(v); }
};

/* Vector case */
template <size_t N, typename T> struct get_value_helper<Vector<N, T>, void> {
  Vector<N, T> operator()(Variant const &v) const {
    return Vector<N, T>(boost::get<std::vector<double>>(v));
  }
};

/* This allows direct retrieval of a shared_ptr to the object from
   an ObjectId variant. If the type is a derived type, the type is
   also checked.

   We do a couple of checks: First we chack if the id is actualy the
   empty id, which means None and is a valid value, represented by
   an empty ptr.
   If the id is not empty, we try to retieve an instance for that id.
   If it does not exist we throw, this means the caller supplied an id
   with no corresponding instance.
   If we can find an instance, we check if it actualy has the right
   type, and if so, we return it, therwise we throw.
*/
template <typename T>
struct get_value_helper<
    std::shared_ptr<T>,
    typename std::enable_if<std::is_base_of<ScriptInterfaceBase, T>::value,
                            void>::type> {
  std::shared_ptr<T> operator()(Variant const &v) const {
    auto const object_id = boost::get<ObjectId>(v);
    if (object_id == ObjectId()) {
      return nullptr;
    } else {
      auto so_ptr = ScriptInterfaceBase::get_instance(object_id).lock();
      if (!so_ptr) {
        throw std::runtime_error("Unknown Object.");
      }

      auto t_ptr = std::dynamic_pointer_cast<T>(so_ptr);

      if (t_ptr) {
        return t_ptr;
      } else {
        throw std::runtime_error("Wrong type.");
      }
    }
  }
};
}

/**
 * @brief Extract value of specific type T from a Variant.
 *
 * This is a wrapper around boost::get that allows us to
 * customize the behavior for different types. This is
 * needed e.g. to deal with Vector types that are not
 * explicitly contained in Variant, but can easily
 * be converted.
 */
template <typename T> T get_value(Variant const &v) {
  return detail::get_value_helper<T>{}(v);
}

} /* namespace ScriptInterface */

#endif
