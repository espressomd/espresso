#ifndef SCRIPT_INTERFACE_PACKED_VARIANT_HPP
#define SCRIPT_INTERFACE_PACKED_VARIANT_HPP

#include "Variant.hpp"

#include <functional>
#include <unordered_map>
#include <utility>

namespace ScriptInterface {
    using ObjectId = std::size_t;

        ObjectId object_id(const ObjectHandle *p) {
            return std::hash<const ObjectHandle *>{}(p);
        }
        ObjectId object_id(ObjectRef const &p) { return object_id(p.get()); }

        std::unordered_map<ObjectId, std::shared_ptr<ObjectHandle>> local_objects;

        using PackedVariant = boost::make_recursive_variant<
                None, bool, int, double, std::string, std::vector<int>, std::vector<double>,
        ObjectId, std::vector<boost::recursive_variant_>, Utils::Vector2d,
        Utils::Vector3d, Utils::Vector4d>::type;

        using PackedMap = std::vector<std::pair<std::string, PackedVariant>>;

        struct VariantToTransport
                : recursive_visitor<VariantToTransport, Variant, PackedVariant> {
            using recursive_visitor<VariantToTransport, Variant, PackedVariant>::
            operator();

            template <class T> PackedVariant operator()(T &&val) const {
                return std::forward<T>(val);
            }

            PackedVariant operator()(const ObjectRef &so_ptr) const {
                return object_id(so_ptr);
            }
        };

        struct TransportToVariant
                : recursive_visitor<TransportToVariant, PackedVariant, Variant> {
            using recursive_visitor<TransportToVariant, PackedVariant, Variant>::
            operator();

            template <class T> Variant operator()(T &&val) const {
                return std::forward<T>(val);
            }

            Variant operator()(const ObjectId &id) const { return local_objects.at(id); }
        };

        PackedVariant pack(const Variant &v) {
            return boost::apply_visitor(VariantToTransport{}, v);
        }

        Variant unpack(const PackedVariant &v) {
            return boost::apply_visitor(TransportToVariant{}, v);
        }

        PackedMap pack(const VariantMap &v) {
            std::vector<std::pair<std::string, PackedVariant>> ret(v.size());

            boost::transform(v, ret.begin(), [](auto const &kv) {
                return std::pair<std::string, PackedVariant>{kv.first, pack(kv.second)};
            });

            return ret;
        }

        VariantMap unpack(const PackedMap &v) {
            VariantMap ret;

            boost::transform(v, std::inserter(ret, ret.end()), [](auto const &kv) {
                return std::pair<std::string, Variant>{kv.first, unpack(kv.second)};
            });

            return ret;
        }
}

#endif