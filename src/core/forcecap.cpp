/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016,2017 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
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
/** \file forcecap.cpp force cap calculation.
 *
 *  For more information see \ref forcecap.hpp "forcecap.h".
*/

#include "forcecap.hpp"
#include "utils.hpp"
#include "global.hpp"

double force_cap = 0.0;

void forcecap_set(double forcecap) {
  force_cap = forcecap;
  mpi_bcast_parameter(FIELD_FORCE_CAP);
}

double forcecap_get() {
  return force_cap;
}

void forcecap_cap(ParticleRange particles) {
  if (force_cap <= 0) {
    return;
  }

  auto const fc2 = force_cap * force_cap;

  for (auto &p : particles) {
    auto const f2 = sqrlen(p.f.f);
    if (f2 > fc2) {
      auto const scale = force_cap / std::sqrt(f2);

      for (int i = 0; i < 3; i++) {
        p.f.f[i] *= scale;
      }
    }
  }
}
