#ifndef CORE_UTILS_SERIALIZATION_FLAT_MAP_HPP
#define CORE_UTILS_SERIALIZATION_FLAT_MAP_HPP

#include <boost/container/flat_map.hpp>

namespace boost {
namespace serialization {

template <typename Archive, typename K, typename V>
void load(Archive &ar, boost::container::flat_map<K, V> &v,
          const unsigned int) {
  using value_type = typename boost::container::flat_map<K, V>::value_type;
  typename boost::container::flat_map<K, V>::size_type count;

  ar &count;
  v.reserve(count);

  for (int i = 0; i < count; i++) {
    value_type e;

    ar >> e;
    v.emplace_hint(v.end(), e);
  }
}

template <typename Archive, typename K, typename V>
void save(Archive &ar, boost::container::flat_map<K, V> const &v,
          const unsigned int) {
  typename boost::container::flat_map<K, V>::size_type count(v.size());

  ar << count;

  for (auto const &e : v) {
    ar << e;
  }
}

template <typename Archive, typename K, typename V>
void serialize(Archive &ar, boost::container::flat_map<K, V> &v,
               const unsigned int version) {
  split_free(ar, v, version);
}
} // namespace serialization
} // namespace boost

#endif
