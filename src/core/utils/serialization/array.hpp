#ifndef UTILS_SERIALIZATION_ARRAY_HPP
#define UTILS_SERIALIZATION_ARRAY_HPP

#include <array>

#include <boost/version.hpp>

/* New versions of boost alrady containt this
 * function
 * */
#if BOOST_VERSION < 105600

namespace boost {
namespace serialization {
template <typename Archive, typename T, std::size_t N>
void serialize(Archive &ar, std::array<T, N> &a, const unsigned int) {
  ar &*static_cast<T(*)[N]>(static_cast<void *>(a.data()));
}
}
}
#endif
#endif
