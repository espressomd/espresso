/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file tuning.hpp
    This contains a timing loop for the force calculation. Via the global variable timings you can specify how many
    force evaluations are sampled. Via \ref markTime and \ref diffTime you can also easily time anything other than
    the force evaluation.
*/

#ifndef TUNING_H
#define TUNING_H

/** if positive, the number of samples for timing */
extern int timing_samples;

/** returns the time for some force calculations.
    Actually performs \ref mpi_integrate (0)
    @param default_samples the number of samples to take if
    \ref timing_samples is not set. */
double time_force_calc(int default_samples);

/** set a time marker. \ref diffTime always gives the time in ms between
    the last two calls to markTime. */
void markTime();

/** calculate milliseconds between last two calls to \ref markTime. */
double diffTime();

#endif
