#ifndef ESPRESSO_UTILS_MPI_DATATYPE_HPP
#define ESPRESSO_UTILS_MPI_DATATYPE_HPP

#include <boost/mpi/datatype_fwd.hpp>
#include <boost/optional/optional_fwd.hpp>

namespace boost {
namespace mpi {
template <class T>
struct is_mpi_datatype<boost::optional<T>> : is_mpi_datatype<T> {};
} // namespace mpi
} // namespace boost

#endif
