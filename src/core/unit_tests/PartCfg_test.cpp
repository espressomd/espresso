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

/** \file MpiCallbacks_test.cpp Unit tests for the MpiCallbacks class.
 *
*/

#include <boost/mpi.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE PartCfg test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "mock/Particle.hpp"
using Testing::Particle;

#include "core/PartCfg.hpp"

struct Particles {
  std::vector<Particle> parts;

  std::vector<Particle> &particles() { return parts; }
};

BOOST_AUTO_TEST_CASE(update) {
  Particles local_parts;

  auto const rank = Communication::mpiCallbacks().comm().rank();
  auto const size = Communication::mpiCallbacks().comm().size();
  auto const n_part = 10000;

  local_parts.parts.reserve(n_part);

  for (int i = 0; i < n_part; i++) {
    local_parts.parts.emplace_back(rank * n_part + i);
  }

  PartCfg<Particles> part_cfg{local_parts};

  if (rank == 0) {
    part_cfg.update(0);
    Communication::mpiCallbacks().abort_loop();

    BOOST_CHECK(part_cfg.size() == size * n_part);

    for (int i = 0; i < size; i++) {
      BOOST_CHECK(i == part_cfg[i].identity());
    }
  } else
    Communication::mpiCallbacks().loop();
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
