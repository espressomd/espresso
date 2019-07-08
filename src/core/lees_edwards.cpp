/*
  Copyright (C) 2018 The ESPResSo project

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
/** \file lees_edwards.cpp
 */

#include "lees_edwards.hpp"
#include "integrate.hpp"
#include <cmath>

lees_edwards_protocol_struct lees_edwards_protocol = {
    LEES_EDWARDS_PROTOCOL_OFF, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 1};

/* Functions to determine the current offset and shear rate with respect to the
 * chosen protocol */

double lees_edwards_get_offset(double time) {

  if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_OFF) {
    return 0.;
  } else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_STEP) {
    return lees_edwards_protocol.offset;
  } else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_STEADY_SHEAR) {
    return lees_edwards_protocol.velocity *
           (time - lees_edwards_protocol.time0);
  } else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_OSC_SHEAR) {
    return lees_edwards_protocol.amplitude *
           std::sin(lees_edwards_protocol.frequency *
                    (time - lees_edwards_protocol.time0));
  } else {
    return 0.0;
  }
}

double lees_edwards_get_velocity(double time) {
  if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_OFF) {
    return 0.;
  } else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_STEP) {
    return 0.;
  } else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_STEADY_SHEAR) {
    return lees_edwards_protocol.velocity;
  } else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_OSC_SHEAR) {
    return lees_edwards_protocol.frequency * lees_edwards_protocol.amplitude *
           std::cos(lees_edwards_protocol.frequency *
                    (time - lees_edwards_protocol.time0));
  } else {
    return 0.0;
  }
}

bool less_edwards_supports_verlet_list() { return false; }

#ifdef LEES_EDWARDS
void local_lees_edwards_image_reset() {
  for (auto &p : local_cells.particles()) {
    p.l.i[0] = 0;
    p.l.i[1] = 0;
    p.l.i[2] = 0;
    p.p.lees_edwards_offset = 0;
  }
}
#endif
