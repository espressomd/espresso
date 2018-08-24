/*
  Copyright (C) 2017 The ESPResSo project
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <random>
#include <string>
#include <vector>

#include <boost/mpi.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE gather_buffer test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/mpi/gather_buffer.hpp"
using Utils::Mpi::gather_buffer;
namespace mpi = boost::mpi;

void check_pointer(const mpi::communicator & comm, int root) {
  if (comm.rank() == root) {
    auto const n = comm.size();
    const int total_size = n * (n + 1) / 2;

    std::vector<int> buf(total_size, comm.rank() + 1);
    auto const ret_size =
        gather_buffer(buf.data(), comm.rank() + 1, comm, root);

    BOOST_CHECK(ret_size == total_size);

    /* Check order in result */
    BOOST_CHECK(std::is_sorted(buf.begin(), buf.end()));

    /* Check values */
    for (int i = 1; i <= n; i++) {
      std::vector<int>::iterator lower, upper;
      std::tie(lower, upper) = std::equal_range(buf.begin(), buf.end(), i);

      BOOST_CHECK(i == std::distance(lower, upper));
    }
  } else {
    std::vector<int> buf(comm.rank() + 1, comm.rank() + 1);
    gather_buffer(buf.data(), buf.size(), comm, root);

    /* Check that buffer is unchanged */
    BOOST_CHECK(buf.size() == comm.rank() + 1);
    for (auto const &i : buf) {
      BOOST_CHECK(i == comm.rank() + 1);
    }
  }
}

void check_vector(const mpi::communicator & comm, int root) {
  std::vector<int> buf(comm.rank() + 1, comm.rank() + 1);

  gather_buffer(buf, comm, root);

  if (comm.rank() == root) {
    auto const n = comm.size();
    const int total_size = n * (n + 1) / 2;

    BOOST_CHECK(buf.size() == total_size);

    /* Check order in result */
    BOOST_CHECK(std::is_sorted(buf.begin(), buf.end()));

    /* Check values */
    for (int i = 1; i <= n; i++) {
      std::vector<int>::iterator lower, upper;
      std::tie(lower, upper) = std::equal_range(buf.begin(), buf.end(), i);

      BOOST_CHECK(i == std::distance(lower, upper));
    }
  } else {
    /* Check that buffer is unchanged */
    BOOST_CHECK(buf.size() == comm.rank() + 1);
    for (auto const &i : buf) {
      BOOST_CHECK(i == comm.rank() + 1);
    }
  }
}

void check_vector_empty(const mpi::communicator & comm, int empty) {
  std::vector<int> buf((comm.rank() == empty) ? 0 : 11, comm.rank());
  gather_buffer(buf, comm);

  if (comm.rank() == 0) {
    BOOST_CHECK(buf.size() == (comm.size() - 1) * 11);

    for (int i = 0; i < comm.size(); i++) {
      std::vector<int>::iterator lower, upper;
      std::tie(lower, upper) = std::equal_range(buf.begin(), buf.end(), i);

      if (i == empty) {
        BOOST_CHECK(0 == std::distance(lower, upper));
      } else {
        BOOST_CHECK(11 == std::distance(lower, upper));
      }
    }
  }
}

void check_pointer_empty(const mpi::communicator & comm, int empty) {
  auto const n_elem = (comm.rank() == empty) ? 0 : 11;
  std::vector<int> buf(n_elem, comm.rank());

  if (comm.rank() == 0) {
    buf.resize((comm.size() - 1) * 11);
  }

  gather_buffer(buf.data(), n_elem, comm);

  if (comm.rank() == 0) {
    for (int i = 0; i < comm.size(); i++) {
      std::vector<int>::iterator lower, upper;
      std::tie(lower, upper) = std::equal_range(buf.begin(), buf.end(), i);

      if (i == empty) {
        BOOST_CHECK(0 == std::distance(lower, upper));
      } else {
        BOOST_CHECK(11 == std::distance(lower, upper));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(pointer) {
  mpi::communicator world;
  check_pointer(world, 0);
}

BOOST_AUTO_TEST_CASE(pointer_root) {
  mpi::communicator world;

  auto root = (world.size() >= 3) ? world.size() - 2 : world.size() - 1;
  check_pointer(world, root);
}

BOOST_AUTO_TEST_CASE(vector) {
  mpi::communicator world;
  check_pointer(world, 0);
}

BOOST_AUTO_TEST_CASE(vector_root) {
  mpi::communicator world;

  auto empty = (world.size() >= 3) ? world.size() - 2 : world.size() - 1;
  check_pointer(world, empty);
}

BOOST_AUTO_TEST_CASE(vector_empty) {
  mpi::communicator world;

  check_vector_empty(world, 0);
}

BOOST_AUTO_TEST_CASE(vector_empty_root) {
  mpi::communicator world;
  auto root = (world.size() >= 3) ? world.size() - 2 : world.size() - 1;

  check_vector_empty(world, root);
}

BOOST_AUTO_TEST_CASE(pointer_empty) {
  mpi::communicator world;

  check_pointer_empty(world, 0);
}

BOOST_AUTO_TEST_CASE(pointer_empty_root) {
  mpi::communicator world;
  auto root = (world.size() >= 3) ? world.size() - 2 : world.size() - 1;

  check_pointer_empty(world, root);
}

BOOST_AUTO_TEST_CASE(non_trivial_type) {
  mpi::communicator world;

  std::string s(
      "A long string that is too long for short-string optimization.");
  BOOST_CHECK(s.size() > sizeof(std::string));

  std::vector<std::string> buf(world.rank() + 1, s);

  gather_buffer(buf, world);

  if(world.rank() == 0) {
    auto const n = world.size();
    BOOST_CHECK(buf.size() == (n * (n+1) / 2));

    for(auto const& e : buf) {
      BOOST_CHECK(e == s);
    }
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
