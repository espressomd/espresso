//
// Created by florian on 14.05.19.
//

#ifndef ESPRESSO_PACK_HPP
#define ESPRESSO_PACK_HPP

#include "Variant.hpp"
#include "get_value.hpp"

#include <unordered_map>
#include <utility>

namespace ScriptInterface {
    template <typename T, typename U>
    std::vector<Variant> pack_pair(const std::pair<T, U> &pair) {
        return {{pair.first, pair.second}};
    }

    template <typename T, typename U>
    const std::pair<T, U> unpack_pair(const std::vector<Variant> &v) {
        return {boost::get<T>(v.at(0)), boost::get<U>(v.at(1))};
    }

/**
 * @brief Pack a map into a vector of Variants
 *        by serializing the key-value pairs.
 *
 */
    template <typename K, typename V>
    std::vector<Variant> pack_map(const std::unordered_map<K, V> &map) {
        std::vector<Variant> ret(map.size());

        std::transform(map.begin(), map.end(), ret.begin(),
                       [](const std::pair<K, V> &p) { return pack_pair(p); });

        return ret;
    }

    template <typename K, typename V>
    std::unordered_map<K, V> unpack_map(const std::vector<Variant> &v) {
        std::unordered_map<K, V> ret;

        for (auto const &pair : v) {
            ret.insert(unpack_pair<K, V>(boost::get<const std::vector<Variant>>(pair)));
        }

        return ret;
    }
}
#endif //ESPRESSO_PACK_HPP
