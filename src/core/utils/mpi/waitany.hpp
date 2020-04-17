/*
 * Copyright (C) 2010-2020 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef WAITANY_HPP_INCLUDED
#define WAITANY_HPP_INCLUDED

#include <boost/mpi/request.hpp>
#include <boost/mpi/status.hpp>
#include <mpi.h>
#include <utility> // std::pair

namespace Utils {
namespace Mpi {

/** A pair of boost::mpi::requests that supports test().
 */
class RequestPair {
  std::pair<boost::mpi::request, boost::mpi::request> m_reqs;

public:
  const boost::mpi::request &first() const { return m_reqs.first; }
  const boost::mpi::request &second() const { return m_reqs.second; }

  boost::mpi::request &first() { return m_reqs.first; }
  boost::mpi::request &second() { return m_reqs.second; }

  /** Returns true if both requests are completed.
   * In accordance with boost::mpi, this function should not be called
   * after it returned true. (Will always return false in this case.)
   */
  bool test();

  RequestPair(boost::mpi::request &&r1, boost::mpi::request &&r2)
      : m_reqs(std::move(r1), std::move(r2)) {}
  RequestPair(RequestPair &&) = default;
  RequestPair &operator=(RequestPair &&) = default;
  RequestPair() = default;
};

namespace __details {
/** Checks if a requests still has pending communication.
 */
bool is_valid(const boost::mpi::request &r) {
  return r.m_requests[0] != MPI_REQUEST_NULL ||
         r.m_requests[1] != MPI_REQUEST_NULL;
}

/** Checks if any of the two requests in a pair still has pending communication.
 */
bool is_valid(const RequestPair &rp) {
  return is_valid(rp.first()) || is_valid(rp.second());
}
} // namespace __details

bool RequestPair::test() {
  // If a request is *not* valid anymore, it has been completed before.
  bool completed1 = !__details::is_valid(first());
  bool completed2 = !__details::is_valid(second());

  if (!completed1 && first().test()) {
    if (completed2)
      return true;
    // Update "completed1" to allow the following ifs to directly return true.
    completed1 = true;
  }

  if (!completed2 && second().test()) {
    if (completed1)
      return true;
    // No need to set "completed2" here. Function will exit afterwards.
  }

  return false;
}

/** Status return policy for wait_any.
 */
enum class Status { Return, Ignore };

template <Status return_policy, typename It, typename Ret>
struct wait_any_return_type;

template <typename It, typename Ret>
struct wait_any_return_type<Status::Ignore, It, Ret> {
  /** Iterator to the ready element.
   */
  It iterator;

  wait_any_return_type(It _iterator, Ret _ret) : iterator(_iterator) {}
  wait_any_return_type(wait_any_return_type &&) = default;
};

template <typename It, typename Ret>
struct wait_any_return_type<Status::Return, It, Ret> {
  /** Iterator to the ready element
   */
  It iterator;
  /** Return value of the test() call.
   * In case of a RequestPair, this will simply be a boolean that is true.
   * In case of a raw boost::mpi::request, this will be a boost::optional
   * with an associated boost::mpi::status.
   */
  Ret return_value;

  wait_any_return_type(It _iterator, Ret _ret)
      : iterator(_iterator), return_value(_ret) {}
  wait_any_return_type(wait_any_return_type &&) = default;
};

/** Wait_any
 *
 * 1. This implements a fixed version of boost::mpi::wait_any. This function
 * is necessary until the bugfixed version of the boost::mpi function is widely
 * available.
 * 2. This function also works on iterators to RequestPairs.
 *
 * As first template argument, choose if the mpi::status in case 1 should be
 * returned or not.
 *
 */
template <Status return_policy, typename It>
auto wait_any(It first, It last)
    -> wait_any_return_type<return_policy, It, decltype(first->test())> {
  using RetType =
      wait_any_return_type<return_policy, It, decltype(first->test())>;

  bool exists_non_null;
  decltype(first->test()) ret;

  while (true) {
    exists_non_null = false;
    for (It el = first; el != last; el++) {

      // Don't call test() on an already completed boost::mpi::request
      if (!__details::is_valid(*el))
        continue;

      exists_non_null = true;
      if ((ret = el->test()))
        return RetType{el, ret};
    }

    // Prevent infinite loop
    if (!exists_non_null)
      return RetType{last, {}};
  }
}

} // namespace Mpi
} // namespace Utils

#endif
