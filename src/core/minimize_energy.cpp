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
#include "integrate.hpp"
#include "initialize.hpp"

#ifdef MINIMIZE_ENERGY_DEBUG
#define MINIMIZE_ENERGY_TRACE(A) A
#else
#define MINIMIZE_ENERGY_TRACE(A)
#endif

struct MinimizeEnergyParameters {
  double f_max;
  double gamma;
  int max_steps;
  double max_displacement;
};

static MinimizeEnergyParameters *params = 0;

/* Signum of val */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

bool steepest_descent_step(void) {
  Cell *cell;
  Particle *p;
  int c, i, j, np;
  double f_max = -std::numeric_limits<double>::max();
  double f;
  double dp, dp2, dp2_max = -std::numeric_limits<double>::max();
  const double skin2 = SQR(skin);
  double f_max_global;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      f = 0.0;
      dp2 = 0.0;
#ifdef VIRTUAL_SITES
      if (ifParticleIsVirtual(&p[i])) continue;
#endif
      for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
        if (!(p[i].p.ext_flag & COORD_FIXED(j)))
#endif
          {
            f += SQR(p[i].f.f[j]);	    	    
	    dp = params->gamma * p[i].f.f[j];
	    if(fabs(dp) > params->max_displacement)
	      dp = sgn<double>(dp)*params->max_displacement;
	    dp2 += SQR(dp);
            p[i].r.p[j] += dp;
	    MINIMIZE_ENERGY_TRACE(printf("part %d dim %d dp %e gamma*f %e\n", i, j, dp, params->gamma * p[i].f.f[j]));
          }
      }
      f_max = std::max(f_max, f);
      dp2_max = std::max(dp2_max, dp2);
      resort_particles = 1;
    }
  }
  MINIMIZE_ENERGY_TRACE(printf("dp_max %e f_max %e resort_particles %d\n", sqrt(dp2_max), f_max, resort_particles));
  announce_resort_particles();
  MPI_Allreduce(&f_max, &f_max_global, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  return (sqrt(f_max_global) < params->f_max);
}

void minimize_energy_init(const double f_max, const double gamma, const int max_steps, const double max_displacement) {
  if(!params)
    params = new MinimizeEnergyParameters;

  params->f_max = f_max;
  params->gamma = gamma;
  params->max_steps = max_steps;
  params->max_displacement = max_displacement;
}

bool minimize_energy(void) {
  MPI_Bcast(params, sizeof(MinimizeEnergyParameters), MPI_BYTE, 0, comm_cart);
  int integ_switch_old = integ_switch;
  integ_switch = INTEG_METHOD_STEEPEST_DESCENT;
  integrate_vv(params->max_steps, -1);
  integ_switch = integ_switch_old;
  
  return true;
}

