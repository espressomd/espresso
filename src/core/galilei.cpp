/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file galilei.cpp
 *
 */

#include "galilei.hpp"
#include "cells.hpp"
#include "forces.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "utils.hpp"

galilei_struct gal;

/* Stop the particle motion by setting the
   velocity of each particle to zero */
void local_kill_particle_motion(int omega) {
  for (auto &p : local_cells.particles()) {
    p.m.v[0] = 0.0;
    p.m.v[1] = 0.0;
    p.m.v[2] = 0.0;

    if (omega != 0) {
#ifdef ROTATION
      p.m.omega[0] = 0.0;
      p.m.omega[1] = 0.0;
      p.m.omega[2] = 0.0;
#endif
    }
  }
}

/* Set all the forces acting on the particles
   to zero */
void local_kill_particle_forces(int torque) {
  for (auto &p : local_cells.particles()) {
    p.f.f[0] = 0.0;
    p.f.f[1] = 0.0;
    p.f.f[2] = 0.0;

    if (torque != 0) {
#ifdef ROTATION
      p.f.torque[0] = 0.0;
      p.f.torque[1] = 0.0;
      p.f.torque[2] = 0.0;
#endif
    }
  }
}

/* Calculate the CMS of the system */
void local_system_CMS(double *sdata) {
  double x = 0.0, y = 0.0, z = 0.0;
  double mass = 0.0;

  for (auto const &p : local_cells.particles()) {
    double M = p.p.mass;
    mass += M;

    Vector3d ppos = unfolded_position(p);

    x += M * ppos[0];
    y += M * ppos[1];
    z += M * ppos[2];
  }

  sdata[0] = x;
  sdata[1] = y;
  sdata[2] = z;
  sdata[3] = mass;
}

/* Calculate the CMS velocity of the system */
void local_system_CMS_velocity(double *sdata) {
  double x = 0.0, y = 0.0, z = 0.0;
  double mass = 0.0;

  for (auto const &p : local_cells.particles()) {
    double M = p.p.mass;
    mass += M;

    x += M * p.m.v[0];
    y += M * p.m.v[1];
    z += M * p.m.v[2];
  }

  sdata[0] = x;
  sdata[1] = y;
  sdata[2] = z;
  sdata[3] = mass;
}

/* Remove the CMS velocity */
void local_galilei_transform(double *sdata) {
  for (auto &p : local_cells.particles()) {
    p.m.v[0] -= sdata[0];
    p.m.v[1] -= sdata[1];
    p.m.v[2] -= sdata[2];
  }
}
