//
// Created by florian on 14.05.19.
//

#ifndef ESPRESSO_PACK_HPP
#define ESPRESSO_PACK_HPP
#include "get_value.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include <sstream>
#include <unordered_map>
#include <utility>

namespace ScriptInterface {
template <class T> std::string pack(T const &v) {
  std::stringstream ss;
  boost::archive::binary_oarchive(ss) << v;

  return ss.str();
}

template <class T> T unpack(std::string const &state) {
  namespace iostreams = boost::iostreams;

  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);

  T val;
  boost::archive::binary_iarchive(ss) >> val;

  return val;
}

template <typename T, typename U>
std::vector<Variant> pack_pair(const std::pair<T, U> &pair) {
  return {{pair.first, pair.second}};
}

template <typename T, typename U>
const std::pair<T, U> unpack_pair(const std::vector<Variant> &v) {
  return {get_value<T>(v.at(0)), get_value<U>(v.at(1))};
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
    ret.insert(unpack_pair<K, V>(get_value<const std::vector<Variant>>(pair)));
  }

  return ret;
}
} // namespace ScriptInterface
#endif // ESPRESSO_PACK_HPP
