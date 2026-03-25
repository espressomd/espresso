/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include "script_interface/ObjectContainer.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

#include <memory>
#include <ranges>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace ScriptInterface {
/**
 * @brief Owning map of object handles.
 *
 * Mapped elements are cleared from the core during destruction.
 * Due to how dynamic dispatch works, derived types must mark
 * @ref erase_in_core as @c final` and call @ref do_destruct in
 * their virtual destructor. This is to ensure that pure virtual
 * functions called by @ref clear cannot be executed in a type
 * derived from the type currently being destructed, since derived
 * types no longer exist at this point of the destruction sequence.
 *
 * @tparam ManagedType Type of the managed objects, needs to be
 *         derived from @ref ObjectHandle
 */
template <typename ManagedType, class BaseType = ObjectHandle,
          class KeyType = int>
class ObjectMap : public ObjectContainer<ObjectMap, ManagedType, BaseType> {
public:
  using Base = ObjectContainer<ObjectMap, ManagedType, BaseType>;
  using Base::add_parameters;
  using key_type = KeyType;
  using mapped_type = std::shared_ptr<ManagedType>;
  using container_type = std::unordered_map<key_type, mapped_type>;

private:
  container_type m_elements;
  bool dtor_sequence_initiated = false;

  virtual key_type insert_in_core(mapped_type const &obj_ptr) = 0;
  virtual void insert_in_core(key_type const &key,
                              mapped_type const &obj_ptr) = 0;
  virtual void erase_in_core(key_type const &key) = 0;

protected:
  void do_destruct() {
    assert(not dtor_sequence_initiated);
    clear();
    dtor_sequence_initiated = true;
  }

public:
  ObjectMap() {
    add_parameters({
        {"_objects", AutoParameter::read_only,
         [this]() { return make_unordered_map_of_variants(m_elements); }},
    });
  }

  ~ObjectMap() override { assert(dtor_sequence_initiated); }

  // prevent inheritance of the default implementation
  void do_construct(VariantMap const &params) override = 0;

  /**
   * @brief Add an element to the map.
   *
   * @param key Identifier of the element to add.
   * @param element The element to add.
   */
  void insert(key_type const &key, mapped_type const &element) {
    insert_in_core(key, element);
    m_elements[key] = element;
  }

  /**
   * @brief Add an element to the map.
   * A free key is generated automatically.
   *
   * @param element The element to add.
   */
  key_type insert(mapped_type const &element) {
    auto const key = insert_in_core(element);
    m_elements[key] = element;
    return key;
  }

  /**
   * @brief Removes all occurrences of an element from the map.
   *
   * @param key Identifier of the element to remove.
   */
  void erase(key_type const &key) {
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
    for (auto const &key : std::views::elements<0>(m_elements)) {
      erase_in_core(key);
    }

    m_elements.clear();
  }

protected:
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {

    if (method == "insert") {
      auto obj_ptr = get_value<mapped_type>(parameters.at("object"));

      if (parameters.contains("key")) {
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
      auto const view = std::views::elements<0>(m_elements);
      return std::vector<KeyType>{view.begin(), view.end()};
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
      auto const key = get_key(parameters.at("key"));
      return m_elements.contains(key);
    }

    return Base::do_call_method(method, parameters);
  }

  void restore_from_checkpoint(VariantMap const &params) {
    m_elements = get_value_or<decltype(m_elements)>(params, "_objects", {});
    for (auto const &[key, element] : m_elements) {
      insert_in_core(key, element);
    }
  }

  key_type get_key(Variant const &key) const {
    try {
      return get_value<key_type>(key);
    } catch (...) {
      using namespace detail::demangle;
      auto const actual = simplify_symbol_variant(key);
      auto const target = simplify_symbol(static_cast<key_type *>(nullptr));
      if (Base::context()->is_head_node()) {
        throw std::invalid_argument("Key has to be of type '" + target +
                                    "', got type '" + actual + "'");
      }
      throw;
    }
  }
};
} // Namespace ScriptInterface
