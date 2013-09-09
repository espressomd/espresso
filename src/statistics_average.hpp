/*
  Copyright (C) 2012,2013 The ESPResSo project
  
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
#ifndef _STATISTICS_AVERAGE_H
#define _STATISTICS_AVERAGE_H

#include "statistics_observable.hpp"

typedef struct {
  observable* obs;
  double* average;
  double* variance;
  int n_obs;
  int n_sweeps; 
  int ready;
  int with_variance;
} average;

/* The actual containers! */
extern average* averages;
extern int n_averages;

/* The Interface */
int averages_add_empty_average();
int average_init(average* self, observable* obs, int with_variance=1);
int average_update(average* self);
int average_print(average* self);
int average_print_formatted(average* self, Tcl_Interp* interp);


#endif

