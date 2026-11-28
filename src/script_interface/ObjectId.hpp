/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "ObjectHandle.hpp"

#include <boost/serialization/access.hpp>

#include <cstdint>
#include <functional>

namespace ScriptInterface {

/**
 * @brief Strongly typed integer type to hold a unique identifier
 * for a @ref ObjectHandle object, using its memory address.
 * Therefore, this object can only be constructed from a
 * @ref ObjectHandle object on the head node.
 */
struct ObjectId {
#ifdef UINTPTR_MAX
  using value_type = std::uintptr_t;
#else
  using value_type = std::size_t;
  static_assert(sizeof(void *) <= sizeof(value_type));
#endif

  ObjectId() = default;
  ObjectId(ObjectHandle const *p) : m_id{reinterpret_cast<value_type>(p)} {}

  constexpr bool operator==(ObjectId const &) const = default;
  constexpr bool operator!=(ObjectId const &) const = default;

  value_type m_id;

private:
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, unsigned const /*version*/) {
    ar & m_id;
  }
};

} // namespace ScriptInterface

namespace std {
template <> struct hash<ScriptInterface::ObjectId> {
  std::size_t operator()(ScriptInterface::ObjectId const &oid) const noexcept {
    return std::hash<ScriptInterface::ObjectId::value_type>{}(oid.m_id);
  }
};
} // namespace std
