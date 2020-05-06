/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE waitany_test
#define BOOST_TEST_DYN_LINK
#include <boost/mpi.hpp>
#include <boost/test/unit_test.hpp>

#include "utils/mpi/waitany.hpp"

#include <string>

namespace mpi = boost::mpi;

BOOST_AUTO_TEST_CASE(single_request_wait_any) {
  // Manually implemented all-to-all. Payload is rank number.
  // Checks for each completed request, if the associated has been received
  // correctly.
  mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();

  std::vector<int> to_send(size, rank), received(size, -1);
  std::vector<mpi::request> s_requests, r_requests;

  for (int i = 0; i < size; ++i) {
    r_requests.push_back(world.irecv(i, 0, received[i]));
    s_requests.push_back(world.isend(i, 0, to_send[i]));
  }

  for (int i = 0; i < size; ++i) {
    const auto el =
        Utils::Mpi::wait_any(r_requests.begin(), r_requests.end()).iterator;
    // Check validity of return value
    BOOST_CHECK(el != r_requests.end());
    const int source = std::distance(r_requests.begin(), el);
    // Check if data is valid and request is internally marked as completed
    BOOST_CHECK(received[source] == source);
    BOOST_CHECK(!Utils::Mpi::__details::is_active(r_requests[source]));
  }

  // All receives done. wait_any must return end iterator.
  {
    const auto p = Utils::Mpi::wait_any(r_requests.begin(), r_requests.end());
    BOOST_CHECK(p.iterator == r_requests.end());
    BOOST_CHECK(!p.return_value);
  }

  mpi::wait_all(s_requests.begin(), s_requests.end());
}

BOOST_AUTO_TEST_CASE(pair_request_wait_any) {
  // Manually implemented all-to-all. Payload is rank number.
  // Checks for each completed request, if the associated has been received
  // correctly.
  mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();

  std::vector<int> to_send1(size, rank), to_send2(size, rank * 2),
      received1(size, -1), received2(size, -1);
  std::vector<Utils::Mpi::RequestPair> s_requests, r_requests;

  for (int i = 0; i < size; ++i) {
    r_requests.emplace_back(world.irecv(i, 0, received1[i]),
                            world.irecv(i, 1, received2[i]));
    s_requests.emplace_back(world.isend(i, 0, to_send1[i]),
                            world.isend(i, 1, to_send2[i]));
  }

  for (int i = 0; i < size; ++i) {
    const auto el =
        Utils::Mpi::wait_any(r_requests.begin(), r_requests.end()).iterator;
    // Check validity of return value
    BOOST_CHECK(el != r_requests.end());
    const int source = std::distance(r_requests.begin(), el);
    // Check if data is valid
    BOOST_CHECK(received1[source] == source);
    BOOST_CHECK(received2[source] == 2 * source);
  }

  // All receives done. wait_any must return end iterator.
  {
    const auto p = Utils::Mpi::wait_any(r_requests.begin(), r_requests.end());
    BOOST_CHECK(p.iterator == r_requests.end());
    BOOST_CHECK(!p.return_value);
  }

  for (int i = 0; i < size; ++i) {
    const auto p = Utils::Mpi::wait_any<Utils::Mpi::Status::Ignore>(
        s_requests.begin(), s_requests.end());
    BOOST_CHECK(p.iterator != s_requests.end());
  }

  // Check if all requests are internally marked as invalid now.
  for (int i = 0; i < size; ++i) {
    BOOST_CHECK(!Utils::Mpi::__details::is_active(r_requests[i]));
    BOOST_CHECK(!Utils::Mpi::__details::is_active(s_requests[i]));
  }
}

int main(int argc, char **argv) {
  mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
