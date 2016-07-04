/*
  Copyright (C) 2010,2012,2016 The ESPResSo project
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
/** \file lees_edwards.h
    Data and methods for Lees-Edwards periodic boundary conditions.  The gist of LE-PBCs is to impose shear flow
    by constantly moving the PBC wrap such that:  $x_{unfolded} == x_{folded} + x_{img} \times L  + (y_{img} * \gamma * t)$,
    and $vx_{unfolded} =  vx_{folded} + (y_{img} * \gamma)$.
*/
#ifndef LEES_EDWARDS_H
#define LEES_EDWARDS_H

#include "config.hpp"

extern double lees_edwards_offset, lees_edwards_rate;
extern int    lees_edwards_count;

#ifdef LEES_EDWARDS

void lees_edwards_step_boundaries();

#endif //LEES_EDWARDS

#endif //LEES_EDWARDS_H
