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
  /* Verlet list criterion */
  const double skin2 = SQR(0.5*skin);
  double f_max_global;
  double dx[3], dx2;
  const double max_dx2 = SQR(params->max_displacement);

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      f = 0.0;
      dx2 = 0.0;
#ifdef VIRTUAL_SITES
      if (ifParticleIsVirtual(&p[i])) continue;
#endif
      for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
        if (!(p[i].p.ext_flag & COORD_FIXED(j)))
#endif
          {
            f += SQR(p[i].f.f[j]);	    	    
	    dx[j] = params->gamma * p[i].f.f[j];	    
	    dx2 += SQR(dx[j]);
	    MINIMIZE_ENERGY_TRACE(printf("part %d dim %d dx %e gamma*f %e\n", i, j, dx[j], params->gamma * p[i].f.f[j]));
	  }
#ifdef EXTERNAL_FORCES
	else {
	  dx[j] = 0.0;	  
	}
#endif
      }

      if(dx2 <= max_dx2) {
	p[i].r.p[0] += dx[0];
	p[i].r.p[1] += dx[1];
	p[i].r.p[2] += dx[2];
      } else {
	const double c = params->max_displacement/std::sqrt(dx2);
	p[i].r.p[0] += c*dx[0];
	p[i].r.p[1] += c*dx[1];
	p[i].r.p[2] += c*dx[2];
      }
      if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) resort_particles = 1;
      f_max = std::max(f_max, f);
    }
  }
  MINIMIZE_ENERGY_TRACE(printf("f_max %e resort_particles %d\n", f_max, resort_particles));
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
  if(!params)
    params = new MinimizeEnergyParameters;

  MPI_Bcast(params, sizeof(MinimizeEnergyParameters), MPI_BYTE, 0, comm_cart);
  int integ_switch_old = integ_switch;
  integ_switch = INTEG_METHOD_STEEPEST_DESCENT;
  integrate_vv(params->max_steps, -1);
  integ_switch = integ_switch_old;
  
  return true;
}

