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

#include <cmath>
#include "lees_edwards.hpp"
#include "integrate.hpp"

lees_edwards_protocol_struct lees_edwards_protocol;
extern double box_l[3];
double const box_l_y = box_l[1];

/* Functions to determine the current offset and shear rate with respect to the chosen protocol */

void setup_lees_edwards_protocol() {
  if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_OFF) {
    lees_edwards_protocol.offset = 0.0;
    lees_edwards_protocol.velocity = 0.0;
    lees_edwards_protocol.amplitude = 0.0;
    lees_edwards_protocol.frequency = 0.0;
  } 


  else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_STEP) {
    lees_edwards_protocol.velocity = 0.0;
    lees_edwards_protocol.amplitude = 0.0;
    lees_edwards_protocol.frequency = 0.0;
  } 

  else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_STEADY_SHEAR) {
    lees_edwards_protocol.offset = lees_edwards_protocol.velocity * (sim_time-lees_edwards_protocol.time0);
    lees_edwards_protocol.amplitude = 0.0;
    lees_edwards_protocol.frequency = 0.0;
  } 
  
  else if (lees_edwards_protocol.type == LEES_EDWARDS_PROTOCOL_OSC_SHEAR) {
    lees_edwards_protocol.offset = lees_edwards_protocol.amplitude * std::sin(lees_edwards_protocol.frequency*(sim_time-lees_edwards_protocol.time0));
    lees_edwards_protocol.velocity = lees_edwards_protocol.frequency * lees_edwards_protocol.amplitude * std::cos(lees_edwards_protocol.frequency*(sim_time-lees_edwards_protocol.time0)); 
  } 

  else {
  lees_edwards_protocol.offset = 0.0;
  lees_edwards_protocol.velocity = 0.0;
  }

}
