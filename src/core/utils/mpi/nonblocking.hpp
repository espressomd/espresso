//
// Created by fweik on 2/14/19.
//

#ifndef ESPRESSO_NONBLOCKING_HPP
#define ESPRESSO_NONBLOCKING_HPP

#include <boost/mpi/nonblocking.hpp>

namespace Utils {
namespace Mpi {
template <typename ReqRange> void wait_all(ReqRange rng) {
  using std::begin;
  using std::end;

  boost::mpi::wait_all(begin(rng), end(rng));
}
} // namespace Mpi
} // namespace Utils

#endif // ESPRESSO_NONBLOCKING_HPP
