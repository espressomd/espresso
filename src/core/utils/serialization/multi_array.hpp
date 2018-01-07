#ifndef UTILS_SERIALIZATION_MULTI_ARRAY_HPP
#define UTILS_SERIALIZATION_MULTI_ARRAY_HPP

#include <boost/multi_array.hpp>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

namespace boost {
namespace serialization {

template <typename Archive, class T, std::size_t N, class Allocator>
void load(Archive &ar, boost::multi_array<T, N, Allocator> &marray, unsigned) {
  boost::array<std::size_t, N> shape;
  ar &make_array(shape.data(), N);

  marray.resize(shape);

  ar &make_array(marray.data(), marray.num_elements());
}

template <typename Archive, class T, std::size_t N, class Allocator>
void save(Archive &ar, const boost::multi_array<T, N, Allocator> &marray,
          unsigned) {
  ar &make_array(marray.shape(), marray.num_dimensions());

  ar &make_array(marray.data(), marray.num_elements());
}

template <typename Archive, class T, std::size_t N, class Allocator>
void serialize(Archive &ar, boost::multi_array<T, N, Allocator> &v,
               const unsigned int version) {
  split_free(ar, v, version);
}
}
}

#endif
