#ifndef CORE_UTILS_SERIALIZATION_FLAT_SET_HPP
#define CORE_UTILS_SERIALIZATION_FLAT_SET_HPP

#include <boost/container/flat_set.hpp>

namespace boost {
namespace serialization {
template <typename Archive, typename V, typename C>
void load(Archive &ar, boost::container::flat_set<V, C> &v,
          const unsigned int) {
  using value_type = typename boost::container::flat_set<V>::value_type;
  using size_type = typename boost::container::flat_set<V>::size_type;
  size_type count;

  ar >> count;
  v.reserve(count);

  for (; count > 0; --count) {
    value_type e;

    ar >> e;
    v.emplace_hint(v.end(), e);
  }
}

template <typename Archive, typename V, typename C>
void save(Archive &ar, boost::container::flat_set<V, C> const &v,
          const unsigned int) {
  typename boost::container::flat_set<V>::size_type count(v.size());

  ar << count;

  for (auto const &e : v) {
    ar << e;
  }
}

template <typename Archive, typename V, typename C>
void serialize(Archive &ar, boost::container::flat_set<V, C> &v,
               const unsigned int version) {
  split_free(ar, v, version);
}
}
}

#endif
