#ifndef UTILS_SERIALIZATION_ARRAY_HPP
#define UTILS_SERIALIZATION_ARRAY_HPP

#include <boost/version.hpp>
#include <boost/mpi/datatype.hpp>

/* New versions of boost alrady containt this
 * function
 * */
#if BOOST_VERSION < 105600
#include <array>
#include <boost/serialization/serialization.hpp>

namespace boost {
namespace serialization {
template <typename Archive, typename T, std::size_t N>
void serialize(Archive &ar, std::array<T, N> &a, const unsigned int) {
  ar &*static_cast<T(*)[N]>(static_cast<void *>(a.data()));
}
}
}

#else
#include <array>
#include <boost/serialization/serialization.hpp>
#endif

namespace boost {
namespace mpi {
template <> struct is_mpi_datatype<std::array<double, 3>> : mpl::true_ {};
}
}

#endif
