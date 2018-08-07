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
 *
 */

#include <cmath>
#include "lees_edwards.hpp"
#include "integrate.hpp"

lees_edwards_protocol_struct lees_edwards_protocol;

/* Functions to determine the current offset and shear rate with respect to the chosen protocol */

void setup_lees_edwards_protocol(lees_edwards_protocol_struct lees_edwards_protocol, double sim_time, double box_l_y) {
  if (lees_edwards_protocol.type == 0 ) {
    //lees_edwards_protocol.offset = Current value;
    lees_edwards_protocol.velocity = 0.0;
    lees_edwards_protocol.amplitude = 0.0;
    lees_edwards_protocol.frequency = 0.0;
  } 


  if (lees_edwards_protocol.type == 1 ) {
    lees_edwards_protocol.velocity = 0.0;
    lees_edwards_protocol.amplitude = 0.0;
    lees_edwards_protocol.frequency = 0.0;
  } 

  if (lees_edwards_protocol.type == 2 ) {
    lees_edwards_protocol.offset = (lees_edwards_protocol.velocity / box_l_y) * sim_time;
    lees_edwards_protocol.amplitude = 0.0;
    lees_edwards_protocol.frequency = 0.0;
  } 
  
  if (lees_edwards_protocol.type == 3 ) {
    lees_edwards_protocol.offset = (lees_edwards_protocol.amplitude / box_l_y) * std::sin(lees_edwards_protocol.frequency*sim_time);
    lees_edwards_protocol.velocity = (lees_edwards_protocol.amplitude / box_l_y) * std::cos(lees_edwards_protocol.frequency*sim_time); 
  } 
}
