/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#include <limits>
#include <algorithm>
#include "minimize_energy.hpp"

static double steepest_descent_step(double gamma) {
  Cell *cell;
  Particle *p;
  int c, i, j, np;
  double f_max = -std::numeric_limits<double>::max();
  double f;

  force_calc();

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      f = 0.0;
#ifdef VIRTUAL_SITES
      if (ifParticleIsVirtual(&p[i])) continue;
#endif
      for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
        if (!(p[i].p.ext_flag & COORD_FIXED(j)))
#endif
          {
            f += SQR(p[i].f.f[j]);
            p[i].r.p[j] += gamma * p[i].f.f[j];
          }
      }
      f_max = std::max(f_max, f);
    }
  }
  return f_max;
}

struct MinimizeEnergyParameters {
  double f_max;
  double gamma;
  int max_steps;
};

static MinimizeEnergyParameters *params = 0;

void minimize_energy_init(const double f_max, const double gamma, const int max_steps) {
  if(!params)
    params = new MinimizeEnergyParameters;

  params->f_max = f_max;
  params->gamma = gamma;
  params->max_steps = max_steps;
}

bool minimize_energy(void) {
  MPI_Bcast(params, sizeof(MinimizeEnergyParameters), MPI_BYTE, 0, comm_cart);
  double f_max_local = 2*params->f_max;
  double f_max_global;
  for(int i = 0; i < params->max_steps; ++i) {
    f_max_local = steepest_descent_step(params->gamma);
    MPI_Allreduce(&f_max_local, &f_max_global, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
    if(this_node == 0)
      printf("minimize_energy(): step %d, f_max %e\n", i, f_max_global);
    if(f_max_global < params->f_max)
      return true;
  }
  return false;

}

