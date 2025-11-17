/*
 * Copyright (C) 2025 The ESPResSo project
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

#pragma once

#include "communication.hpp"

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/test/unit_test.hpp>

#include <memory>

/**
 * @brief Global variables manager for the unit tests.
 *
 * This fixture provides a fallback strategy to safely initialize and destruct
 * MPI-related global variables from the ESPResSo core without race conditions.
 */
struct EspressoCoreGlobalConfig {
  EspressoCoreGlobalConfig() {
    auto mpi_env = std::make_shared<boost::mpi::environment>(
        boost::unit_test::framework::master_test_suite().argc,
        boost::unit_test::framework::master_test_suite().argv,
        boost::mpi::threading::multiple);
    ::communication_environment =
        std::make_unique<CommunicationEnvironment>(mpi_env);
  }
  ~EspressoCoreGlobalConfig() { ::communication_environment.reset(); }
};
