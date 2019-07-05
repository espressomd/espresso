#ifndef UTILS_MPI_CART_COMM_HPP
#define UTILS_MPI_CART_COMM_HPP

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/exception.hpp>

#include <mpi.h>

#include <utility>

namespace Utils {
namespace Mpi {

/**
 * @brief Wrapper around MPI_Dims_create.
 *
 * @tparam dim Number of dimensions
 */
template <size_t dim> Vector<int, dim> dims_create(int nodes) {
  Vector<int, dim> dims{};
  BOOST_MPI_CHECK_RESULT(MPI_Dims_create,
                         (nodes, static_cast<int>(dim), dims.data()))

  return dims;
}

/**
 * @brief Wrapper around MPI_Cart_create.
 *
 * @tparam dim Number of dimensions
 */
template <size_t dim>
boost::mpi::communicator cart_create(
    boost::mpi::communicator const &comm, Vector<int, dim> const &dims,
    bool reorder = true,
    Vector<int, dim> const &periodicity = Vector<int, dim>::broadcast(1)) {
  MPI_Comm temp_comm;
  BOOST_MPI_CHECK_RESULT(MPI_Cart_create,
                         (comm, dim, dims.data(), periodicity.data(),
                          static_cast<int>(reorder), &temp_comm))

  return boost::mpi::communicator(temp_comm, boost::mpi::comm_take_ownership);
}

/**
 * @brief Wrapper around MPI_Cart_coords.
 *
 * @tparam dim Number of dimensions
 */
template <size_t dims>
Vector3i cart_coords(boost::mpi::communicator const &comm, int rank) {
  Vector3i pos;
  BOOST_MPI_CHECK_RESULT(MPI_Cart_coords, (comm, rank, dims, pos.data()))
  return pos;
}

/**
 * @brief Wrapper around MPI_Cart_rank.
 *
 * @tparam dim Number of dimensions
 */
template <size_t dims>
int cart_rank(boost::mpi::communicator const &comm,
              const Vector<int, dims> &pos) {
  int rank;
  BOOST_MPI_CHECK_RESULT(MPI_Cart_rank, (comm, pos.data(), &rank))
  return rank;
}

/**
 * @brief Wrapper around MPI_Cart_shift.
 *
 * @return pair of source and destination rank.
 */
inline std::pair<int, int> cart_shift(boost::mpi::communicator const &comm,
                                      int direction, int displacement) {
  int src = -1, dst = -1;
  BOOST_MPI_CHECK_RESULT(MPI_Cart_shift,
                         (comm, direction, displacement, &src, &dst))

  return {src, dst};
}

/**
 * @brief Calculate face neighbors in cartesian communicator.
 *
 * This calculates the ranks of neighbors with shift
 * (+/- 1, 0, 0), (0, +/-1, 0), ... in a cartesian communicator.
 * This is also the order in which they are returned.
 *
 * @tparam dim Number of dimensions
 * @param comm Cartesian communicator
 * @return Array of neighbor ranks
 */
template <size_t dim>
Utils::Vector<int, 2 * dim>
calc_face_neighbors(const boost::mpi::communicator &comm) {
  using std::get;

  auto ret = Utils::Vector<int, 2 * dim>::broadcast(-1);

  for (int i = 0; i < static_cast<int>(dim); i++) {
    ret[2 * i + 0] = get<1>(cart_shift(comm, i, -1));
    ret[2 * i + 1] = get<1>(cart_shift(comm, i, +1));
  }

  return ret;
}

} // namespace Mpi
} // namespace Utils

#endif
