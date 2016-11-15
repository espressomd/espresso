#ifndef UTILS_SERIALIZATION_ARRAY_HPP
#define UTILS_SERIALIZATION_ARRAY_HPP

#include <boost/version.hpp>

/* New versions of boost already containt this
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

/* Will be included in boost from verison 1.63 */
#if BOOST_VERSION <= 106200
#include <boost/mpi/datatype_fwd.hpp>
namespace boost {
namespace mpi {
template <class T, std::size_t N>
struct is_mpi_datatype<std::array<T, N>> : public is_mpi_datatype<T> {};
} /* namespace mpi */
} /* namespace boost */
#endif
#endif
