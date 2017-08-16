#ifndef CORE_UTILS_SERIALIZATION_LIST_HPP
#define CORE_UTILS_SERIALIZATION_LIST_HPP

#include <boost/serialization/split_free.hpp>

#include "core/utils/List.hpp"

namespace boost {
namespace serialization {
template <typename T, class Archive>
void load(Archive &ar, Utils::List<T> &v, const unsigned int file_version) {
  typename Utils::List<T>::size_type n;
  ar >> n;
  v.resize(n);

  ar >> make_array(v.data(), n);
}

template <typename T, class Archive>
void save(Archive &ar, Utils::List<T> const &v,
          const unsigned int file_version) {
  typename Utils::List<T>::size_type n = v.size();
  ar << n;
  ar << make_array(v.data(), v.size());
}

template <typename T, class Archive>
void serialize(Archive &ar, Utils::List<T> &v,
               const unsigned int file_version) {
  split_free(ar, v, file_version);
}
}
}

#endif
