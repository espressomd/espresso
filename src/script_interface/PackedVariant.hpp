#ifndef SCRIPT_INTERFACE_PACKED_VARIANT_HPP
#define SCRIPT_INTERFACE_PACKED_VARIANT_HPP

#include "Variant.hpp"

#include <functional>
#include <unordered_map>
#include <utility>

namespace ScriptInterface {
using ObjectId = std::size_t;

inline ObjectId object_id(const ObjectHandle *p) {
  return std::hash<const ObjectHandle *>{}(p);
}
inline ObjectId object_id(ObjectRef const &p) { return object_id(p.get()); }

using PackedVariant = boost::make_recursive_variant<
    None, bool, int, double, std::string, std::vector<int>, std::vector<double>,
    ObjectId, std::vector<boost::recursive_variant_>, Utils::Vector2d,
    Utils::Vector3d, Utils::Vector4d>::type;

using PackedMap = std::vector<std::pair<std::string, PackedVariant>>;

struct PackVisitor
    : recursive_visitor<PackVisitor, Variant, PackedVariant> {
  using recursive_visitor<PackVisitor, Variant, PackedVariant>::
  operator();

  template <class T> PackedVariant operator()(T &&val) const {
    return std::forward<T>(val);
  }

  PackedVariant operator()(const ObjectRef &so_ptr) const {
    return object_id(so_ptr);
  }
};

struct UnpackVisitor
    : recursive_visitor<UnpackVisitor, PackedVariant, Variant> {
  std::unordered_map<ObjectId, ObjectRef> const &local_objects;

  explicit UnpackVisitor(
      std::unordered_map<ObjectId, ObjectRef> const &local_objects)
      : local_objects(local_objects) {}

  using recursive_visitor<UnpackVisitor, PackedVariant, Variant>::
  operator();

  template <class T> Variant operator()(T &&val) const {
    return std::forward<T>(val);
  }

  Variant operator()(const ObjectId &id) const { return local_objects.at(id); }
};

inline PackedVariant pack(const Variant &v) {
  return boost::apply_visitor(PackVisitor{}, v);
}

inline Variant unpack(const PackedVariant &v,
               std::unordered_map<ObjectId, ObjectRef> const &local_objects) {
  return boost::apply_visitor(UnpackVisitor{local_objects}, v);
}

inline PackedMap pack(const VariantMap &v) {
  std::vector<std::pair<std::string, PackedVariant>> ret(v.size());

  boost::transform(v, ret.begin(), [](auto const &kv) {
    return std::pair<std::string, PackedVariant>{kv.first, pack(kv.second)};
  });

  return ret;
}

inline VariantMap unpack(const PackedMap &v,
                  std::unordered_map<ObjectId, ObjectRef> const&local_objects) {
  VariantMap ret;

  boost::transform(v, std::inserter(ret, ret.end()), [&local_objects](auto const &kv) {
    return std::pair<std::string, Variant>{kv.first,
                                           unpack(kv.second, local_objects)};
  });

  return ret;
}
} // namespace ScriptInterface

#endif