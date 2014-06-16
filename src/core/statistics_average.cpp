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
#include "statistics_average.hpp"


average* averages;
int n_averages;

int averages_add_empty_average() {
  n_averages++;
  realloc(averages, n_averages*sizeof(Average*));
  averages[n_averages-1]=(Average*) malloc(sizeof(Average));
  return 0;
}

int average_init(Average* self, observable* obs, int with_variance) {
  self->obs=obs;
  self->n_obs=obs->n;
  self->with_variance=with_variance;
  self->n_sweeps=0;
  self->average=(double*) malloc(self->n_obs*sizeof(double));
  if (self->with_variance)
    self->variance=(double*) malloc(self->n_obs*sizeof(double)); 
  self->ready=self->obs->ready;
}

int 


