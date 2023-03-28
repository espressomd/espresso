/*
 * Copyright (C) 2020-2022 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_PACKED_VARIANT_HPP
#define SCRIPT_INTERFACE_PACKED_VARIANT_HPP

#include "Variant.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
using ObjectId = std::size_t;

/**
 * @brief Id for object.
 *
 * This assigns every ObjectHandle a unique id.
 */
inline ObjectId object_id(const ObjectHandle *p) {
  // NOLINTNEXTLINE(bugprone-sizeof-expression)
  static_assert(sizeof(const ObjectHandle *) <= sizeof(ObjectId));
  // Use the pointer value as the unique identifier.
  // This function is only called on the head node.
  return reinterpret_cast<ObjectId>(p);
}

/**
 * @brief Packed version of @ref Variant.
 *
 * When packing variants by @ref PackVisitor, objects of type
 * @ref ObjectRef are packed as @ref ObjectId. Other than that,
 * all other types allowed in @ref Variant must appear here.
 */
using PackedVariant = boost::make_recursive_variant<
    None, bool, int, std::size_t, double, std::string, ObjectId,
    Utils::Vector3b, Utils::Vector3i, Utils::Vector2d, Utils::Vector3d,
    Utils::Vector4d, std::vector<int>, std::vector<double>,
    std::vector<boost::recursive_variant_>,
    std::unordered_map<int, boost::recursive_variant_>,
    std::unordered_map<std::string, boost::recursive_variant_>>::type;

using PackedMap = std::vector<std::pair<std::string, PackedVariant>>;

/**
 * @brief Visitor that converts a Variant to a PackedVariant.
 *
 * While packing, keeps track of all the ObjectRef values that
 * were encountered and stores them. This also keeps the
 * referees alive if there are no other owners.
 */
struct PackVisitor : boost::static_visitor<PackedVariant> {
private:
  mutable std::unordered_map<ObjectId, ObjectRef> m_objects;

public:
  /** @brief Map of objects whose references were replaced by ids. */
  auto const &objects() const { return m_objects; }

  /* For the vector, we recurse into each element. */
  auto operator()(const std::vector<Variant> &vec) const {
    std::vector<PackedVariant> ret(vec.size());

    boost::transform(vec, ret.begin(), [this](const Variant &v) {
      return boost::apply_visitor(*this, v);
    });

    return ret;
  }

  /* For the map, we recurse into each element. */
  template <typename K>
  auto operator()(const std::unordered_map<K, Variant> &map) const {
    std::unordered_map<K, PackedVariant> ret{};

    for (auto const &it : map) {
      ret.insert({it.first, boost::apply_visitor(*this, it.second)});
    }

    return ret;
  }

  /* For object references we store the object reference, and
   * replace it by just an id. */
  PackedVariant operator()(const ObjectRef &so_ptr) const {
    auto const oid = object_id(so_ptr.get());
    m_objects[oid] = so_ptr;

    return oid;
  }

  /* Regular value are just verbatim copied into the result. */
  template <class T> PackedVariant operator()(T &&val) const {
    return std::forward<T>(val);
  }
};

/**
 * @brief Visitor that converts a PackedVariant to a Variant.
 *
 * ObjectId are replaced according to the provided object map.
 */
struct UnpackVisitor : boost::static_visitor<Variant> {
  std::unordered_map<ObjectId, ObjectRef> const &objects;

  explicit UnpackVisitor(std::unordered_map<ObjectId, ObjectRef> const &objects)
      : objects(objects) {}

  /* For the vector, we recurse into each element. */
  auto operator()(const std::vector<PackedVariant> &vec) const {
    std::vector<Variant> ret(vec.size());

    boost::transform(vec, ret.begin(), [this](const PackedVariant &v) {
      return boost::apply_visitor(*this, v);
    });

    return ret;
  }

  /* For the map, we recurse into each element. */
  template <typename K>
  auto operator()(const std::unordered_map<K, PackedVariant> &map) const {
    std::unordered_map<K, Variant> ret{};

    for (auto const &it : map) {
      ret.insert({it.first, boost::apply_visitor(*this, it.second)});
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
 * @brief Transform a Variant to a PackedVariant
 *
 * Applies @ref PackVisitor to a @ref Variant.
 *
 * @param v Input Variant
 * @return Packed variant.
 */
inline PackedVariant pack(const Variant &v) {
  return boost::apply_visitor(PackVisitor(), v);
}

/**
 * @brief Unpack a PackedVariant.
 *
 * Applies @ref UnpackVisitor to a @ref Variant.
 *
 * @param v Packed Variant.
 * @param objects Map of ids to reference.
 * @return Transformed variant.
 */
inline Variant unpack(const PackedVariant &v,
                      std::unordered_map<ObjectId, ObjectRef> const &objects) {
  return boost::apply_visitor(UnpackVisitor(objects), v);
}

/**
 * @brief Pack a VariantMap.
 *
 * Applies @ref pack to every value in the
 * input map.
 */
inline PackedMap pack(const VariantMap &v) {
  PackedMap ret(v.size());

  boost::transform(v, ret.begin(), [](auto const &kv) {
    return std::pair<std::string, PackedVariant>{kv.first, pack(kv.second)};
  });

  return ret;
}

/**
 * @brief Unpack a PackedMap.
 *
 * Applies @ref unpack to every value in the
 * input map.
 */
inline VariantMap
unpack(const PackedMap &v,
       std::unordered_map<ObjectId, ObjectRef> const &objects) {
  VariantMap ret;

  boost::transform(
      v, std::inserter(ret, ret.end()),
      [&objects](auto const &kv) -> std::pair<std::string, Variant> {
        return {kv.first, unpack(kv.second, objects)};
      });

  return ret;
}
} // namespace ScriptInterface

#endif
