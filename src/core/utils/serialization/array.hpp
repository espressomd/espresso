#ifndef UTILS_SERIALIZATION_ARRAY_HPP
#define UTILS_SERIALIZATION_ARRAY_HPP

#include <boost/version.hpp>

/* Newer versions of boost already contain this function */
#if BOOST_VERSION < 105600
#include <array>
namespace boost {
namespace serialization {
template <typename Archive, typename T, std::size_t N>
void serialize(Archive &ar, std::array<T, N> &a, const unsigned int) {
  ar &*static_cast<T(*)[N]>(static_cast<void *>(a.data()));
}
} // namespace serialization
} // namespace boost
#else
#include <boost/serialization/array.hpp>
#endif
#endif
