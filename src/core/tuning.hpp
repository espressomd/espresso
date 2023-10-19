/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "system/System.hpp"

#include <stdexcept>
#include <string>

class TuningFailed : public std::runtime_error {
  std::string get_first_error() const;

public:
  TuningFailed() : std::runtime_error{get_first_error()} {}
};

/**
 * @brief Benchmark the integration loop.
 * Call @ref integrate() several times and measure the elapsed time
 * without propagating the system. It therefore doesn't include e.g.
 * Verlet list updates.
 * @param int_steps   Number of integration steps.
 * @return Average time per integration loop in milliseconds.
 */
double benchmark_integration_step(System::System &system, int int_steps);
