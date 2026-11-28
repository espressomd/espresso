/*
 * Copyright (C) 2020-2026 The ESPResSo project
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

#include "ObjectId.hpp"
#include "Variant.hpp"

#include <utils/serialization/pack.hpp>
#include <utils/serialization/unordered_map.hpp>
#include <utils/serialization/variant.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace ScriptInterface {
/**
 * @brief Packed version of @ref Variant.
 *
 * When packing variants by @ref PackVisitor, objects of type
 * @ref ObjectRef are packed as @ref ObjectId. Other than that,
 * all other types allowed in @ref Variant also appear here.
 */
using PackedVariant = make_recursive_variant<ObjectId>;

using PackedMap = std::vector<std::pair<std::string, PackedVariant>>;

/**
 * @brief Visitor that converts a Variant to a @ref PackedVariant.
 *
 * While packing, keeps track of all the @ref ObjectRef values that
 * were encountered and stores them. This also keeps the
 * referees alive if there are no other owners.
 */
struct PackVisitor {
private:
  mutable std::unordered_map<ObjectId, ObjectRef> m_objects;

public:
  /** @brief Map of objects whose references were replaced by ids. */
  auto const &objects() const { return m_objects; }

  /* For the vector, we recurse into each element. */
  PackedVariant operator()(const std::vector<Variant> &vec) const {
    std::vector<PackedVariant> ret(vec.size());

    std::ranges::transform(vec, ret.begin(), [this](const Variant &v) {
      return std::visit(*this, v);
    });

    return ret;
  }

  /* For the map, we recurse into each element. */
  template <typename K>
  PackedVariant operator()(const std::unordered_map<K, Variant> &map) const {
    std::unordered_map<K, PackedVariant> ret{};

    for (auto const &[key, variant] : map) {
      ret.emplace(key, std::visit(*this, variant));
    }

    return ret;
  }

  /* For object references we store the object reference, and
   * replace it by just an id. */
  PackedVariant operator()(const ObjectRef &so_ptr) const {
    auto const oid = ObjectId(so_ptr.get());
    m_objects[oid] = so_ptr;

    return oid;
  }

  /* Regular value are just verbatim copied into the result. */
  template <class T> PackedVariant operator()(T &&val) const {
    return std::forward<T>(val);
  }
};

/**
 * @brief Visitor that converts a @ref PackedVariant to a @ref Variant.
 *
 * ObjectId are replaced according to the provided object map.
 */
struct UnpackVisitor {
  std::unordered_map<ObjectId, ObjectRef> const &objects;

  explicit UnpackVisitor(std::unordered_map<ObjectId, ObjectRef> const &objects)
      : objects(objects) {}

  /* For the vector, we recurse into each element. */
  Variant operator()(const std::vector<PackedVariant> &vec) const {
    std::vector<Variant> ret(vec.size());

    std::ranges::transform(vec, ret.begin(), [this](const PackedVariant &v) {
      return std::visit(*this, v);
    });

    return ret;
  }

  /* For the map, we recurse into each element. */
  template <typename K>
  Variant operator()(const std::unordered_map<K, PackedVariant> &map) const {
    std::unordered_map<K, Variant> ret{};

    for (auto const &[key, packed_variant] : map) {
      ret.emplace(key, std::visit(*this, packed_variant));
    }

    return ret;
  }

  /* Regular value are just verbatim copied into the result. */
  template <class T> Variant operator()(T &&val) const {
    return std::forward<T>(val);
  }

  /* For object id's they are replaced by references according to the map. */
  Variant operator()(const ObjectId &id) const { return objects.at(id); }
};

/**
 * @brief Transform a Variant to a @ref PackedVariant
 *
 * Applies @ref PackVisitor to a @ref Variant.
 *
 * @param v Input @ref Variant
 * @return Packed variant.
 */
inline PackedVariant pack(const Variant &v) {
  return std::visit(PackVisitor(), v);
}

/**
 * @brief Unpack a @ref PackedVariant.
 *
 * Applies @ref UnpackVisitor to a @ref Variant.
 *
 * @param v Packed @ref Variant.
 * @param objects Map of ids to reference.
 * @return Transformed variant.
 */
inline Variant unpack(const PackedVariant &v,
                      std::unordered_map<ObjectId, ObjectRef> const &objects) {
  return std::visit(UnpackVisitor(objects), v);
}

/**
 * @brief Pack a @ref VariantMap.
 *
 * Applies @ref pack to every value in the
 * input map.
 */
inline PackedMap pack(const VariantMap &v) {
  PackedMap ret(v.size());

  std::ranges::transform(v, ret.begin(), [](auto const &kv) {
    return std::pair<std::string, PackedVariant>{kv.first, pack(kv.second)};
  });

  return ret;
}

/**
 * @brief Unpack a @ref PackedMap.
 *
 * Applies @ref unpack to every value in the
 * input map.
 */
inline VariantMap
unpack(const PackedMap &v,
       std::unordered_map<ObjectId, ObjectRef> const &objects) {
  VariantMap ret;

  std::ranges::transform(
      v, std::inserter(ret, ret.end()),
      [&objects](auto const &kv) -> std::pair<std::string, Variant> {
        return {kv.first, unpack(kv.second, objects)};
      });

  return ret;
}
} // namespace ScriptInterface
