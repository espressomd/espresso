/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
/** \file
 *  Force capping calculation.
 */

#include "forcecap.hpp"

#include "Particle.hpp"
#include "communication.hpp"
#include "event.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>

static double force_cap = 0.0;

double forcecap_get() { return force_cap; }

void forcecap_cap(ParticleRange const &particles) {
  if (force_cap <= 0) {
    return;
  }

  auto const force_cap_sq = Utils::sqr(force_cap);

  for (auto &p : particles) {
    auto const force_sq = p.force().norm2();
    if (force_sq > force_cap_sq) {
      p.force() *= force_cap / std::sqrt(force_sq);
    }
  }
}

void mpi_set_forcecap_local(double force_cap) {
  ::force_cap = force_cap;
  on_forcecap_change();
}

REGISTER_CALLBACK(mpi_set_forcecap_local)

void mpi_set_forcecap(double force_cap) {
  mpi_call_all(mpi_set_forcecap_local, force_cap);
}
