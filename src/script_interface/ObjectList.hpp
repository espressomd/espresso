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

#pragma once

#include "script_interface/ObjectContainer.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/get_value.hpp"

#include <utils/serialization/pack.hpp>

#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ScriptInterface {

/**
 * @brief Owning list of object handles.
 *
 * Stored elements are cleared from the core during destruction.
 * Due to how dynamic dispatch works, derived types must mark
 * @ref remove_in_core as @c final` and call @ref do_destruct in
 * their virtual destructor. This is to ensure that pure virtual
 * functions called by @ref clear cannot be executed in a type
 * derived from the type currently being destructed, since derived
 * types no longer exist at this point of the destruction sequence.
 *
 * @tparam ManagedType Type of the managed objects, needs to be
 *         derived from @ref ObjectHandle
 */
template <typename ManagedType, class BaseType = ObjectHandle>
class ObjectList : public ObjectContainer<ObjectList, ManagedType, BaseType> {
public:
  using Base = ObjectContainer<ObjectList, ManagedType, BaseType>;
  using Base::add_parameters;
  using Base::context;
  using value_type = std::shared_ptr<ManagedType>;

private:
  std::vector<value_type> m_elements;
  bool dtor_sequence_initiated = false;

  virtual void add_in_core(value_type const &obj_ptr) = 0;
  virtual void remove_in_core(value_type const &obj_ptr) = 0;
  virtual bool has_in_core(value_type const &obj_ptr) const = 0;

protected:
  void do_destruct() {
    assert(not dtor_sequence_initiated);
    clear();
    dtor_sequence_initiated = true;
  }

public:
  ObjectList() {
    add_parameters({
        {"_objects", AutoParameter::read_only,
         [this]() { return make_vector_of_variants(m_elements); }},
    });
  }

  ~ObjectList() override { assert(dtor_sequence_initiated); }

  void do_construct(VariantMap const &params) override {
    m_elements = get_value_or<decltype(m_elements)>(params, "_objects", {});
    for (auto const &object : m_elements) {
      add_in_core(object);
    }
  }

  /**
   * @brief Add an element to the list.
   *
   * @param element The element to add.
   */
  void add(value_type const &element) {
    if (has_in_core(element)) {
      if (Base::context()->is_head_node()) {
        throw std::runtime_error("This object is already present in the list");
      }
      throw Exception("");
    }
    context()->parallel_try_catch([&]() { add_in_core(element); });
    m_elements.push_back(element);
  }

  /**
   * @brief Removes all occurrences of an element from the list.
   *
   * @param element The element to remove.
   */
  void remove(value_type const &element) {
    if (not has_in_core(element)) {
      if (Base::context()->is_head_node()) {
        throw std::runtime_error("This object is absent from the list");
      }
      throw Exception("");
    }
    remove_in_core(element);
    std::erase(m_elements, element);
  }

  /**
   * @brief List elements.
   */
  auto const &elements() const { return m_elements; }

  /**
   * @brief Clear the list.
   */
  void clear() {
    for (auto const &e : m_elements) {
      remove_in_core(e);
    }

    m_elements.clear();
  }

protected:
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {

    if (method == "add") {
      auto obj_ptr = get_value<value_type>(parameters.at("object"));

      add(obj_ptr);
      return none;
    }

    if (method == "remove") {
      auto obj_ptr = get_value<value_type>(parameters.at("object"));

      remove(obj_ptr);
      return none;
    }

    if (method == "get_elements") {
      return make_vector_of_variants(m_elements);
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

    return Base::do_call_method(method, parameters);
  }
};
} // Namespace ScriptInterface
