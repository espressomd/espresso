/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef SCRIPT_INTERFACE_OBJECT_MAP_HPP
#define SCRIPT_INTERFACE_OBJECT_MAP_HPP

#include "script_interface/ObjectContainer.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace ScriptInterface {
/**
 * @brief Owning map of ObjectHandles
 * @tparam ManagedType Type of the managed objects, needs to be
 *         derived from @ref ObjectHandle
 */
template <typename ManagedType, class BaseType = ObjectHandle,
          class KeyType = int>
class ObjectMap : public ObjectContainer<ObjectMap, ManagedType, BaseType> {
public:
  using Base = ObjectContainer<ObjectMap, ManagedType, BaseType>;
  using Base::add_parameters;

private:
  std::unordered_map<KeyType, std::shared_ptr<ManagedType>> m_elements;

  virtual KeyType
  insert_in_core(std::shared_ptr<ManagedType> const &obj_ptr) = 0;
  virtual void insert_in_core(KeyType const &key,
                              std::shared_ptr<ManagedType> const &obj_ptr) = 0;
  virtual void erase_in_core(KeyType const &key) = 0;

public:
  ObjectMap() {
    add_parameters({
        {"_objects", AutoParameter::read_only,
         [this]() { return make_unordered_map_of_variants(m_elements); }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_elements = get_value_or<decltype(m_elements)>(params, "_objects", {});
    for (auto const &[key, element] : m_elements) {
      insert_in_core(key, element);
    }
  }

  /**
   * @brief Add an element to the map.
   *
   * @param key Identifier of the element to add.
   * @param element The element to add.
   */
  void insert(KeyType const &key, std::shared_ptr<ManagedType> const &element) {
    insert_in_core(key, element);
    m_elements[key] = element;
  }

  /**
   * @brief Add an element to the map.
   * A free key is generated automatically.
   *
   * @param element The element to add.
   */
  KeyType insert(std::shared_ptr<ManagedType> const &element) {
    auto const key = insert_in_core(element);
    m_elements[key] = element;
    return key;
  }

  /**
   * @brief Removes all occurrences of an element from the map.
   *
   * @param key Identifier of the element to remove.
   */
  void erase(KeyType const &key) {
    erase_in_core(key);
    m_elements.erase(key);
  }

  /**
   * @brief Map elements.
   */
  auto const &elements() const { return m_elements; }

  /**
   * @brief Clear the map.
   */
  void clear() {
    for (auto const &kv : m_elements) {
      erase_in_core(kv.first);
    }

    m_elements.clear();
  }

protected:
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {

    if (method == "insert") {
      auto obj_ptr =
          get_value<std::shared_ptr<ManagedType>>(parameters.at("object"));

      if (parameters.count("key")) {
        auto const key = get_key(parameters.at("key"));
        insert(key, obj_ptr);
        return none;
      }
      return insert(obj_ptr);
    }

    if (method == "erase") {
      auto const key = get_key(parameters.at("key"));
      erase(key);
      return none;
    }

    if (method == "get") {
      auto const key = get_key(parameters.at("key"));
      return Variant{m_elements.at(key)};
    }

    if (method == "get_map") {
      return make_unordered_map_of_variants(m_elements);
    }

    if (method == "keys") {
      std::vector<Variant> res;
      for (auto const &kv : m_elements) {
        res.push_back(kv.first);
      }
      return res;
    }

    if (method == "clear") {
      clear();
      return none;
    }

    if (method == "size") {
      return static_cast<int>(m_elements.size());
    }

    if (method == "empty") {
      return m_elements.empty();
    }

    if (method == "contains") {
      return m_elements.find(get_key(parameters.at("key"))) != m_elements.end();
    }

    return Base::do_call_method(method, parameters);
  }

  KeyType get_key(Variant const &key) const {
    try {
      return get_value<KeyType>(key);
    } catch (...) {
      using namespace detail::demangle;
      auto const actual = simplify_symbol_variant(key);
      auto const target = simplify_symbol(static_cast<KeyType *>(nullptr));
      if (Base::context()->is_head_node()) {
        throw std::invalid_argument("Key has to be of type '" + target +
                                    "', got type '" + actual + "'");
      }
      throw;
    }
  }
};
} // Namespace ScriptInterface
#endif
