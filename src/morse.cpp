/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file morse.cpp
 *
 *  Implementation of \ref morse.hpp
 */
#include "morse.hpp"
#include "communication.hpp"

#ifdef MORSE

int morse_set_params(int part_type_a, int part_type_b,
		     double eps, double alpha,
		     double rmin, double cut, double cap_radius)
{
  double add1, add2;
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  data->MORSE_eps   = eps;
  data->MORSE_alpha = alpha;
  data->MORSE_rmin  = rmin;
  data->MORSE_cut   = cut;

  /* calculate dependent parameter */
  add1 = exp(-2.0*data->MORSE_alpha*(data->MORSE_cut - data->MORSE_rmin));
  add2 = 2.0*exp(-data->MORSE_alpha*(data->MORSE_cut - data->MORSE_rmin));
  data->MORSE_rest = data->MORSE_eps * (add1 - add2); 
 
  if (cap_radius > 0) {
    data->MORSE_capradius = cap_radius;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  mpi_cap_forces(force_cap);

  return ES_OK;
}

void calc_morse_cap_radii()
{
  /* do not compute cap radii if force capping is "individual" */
  if( force_cap != -1.0){
    int i,j,cnt=0;
    IA_parameters *params;
    double force=0.0, rad=0.0, step, add1, add2;

    for(i=0; i<n_particle_types; i++) {
      for(j=0; j<n_particle_types; j++) {
        params = get_ia_param(i,j);
        if(force_cap > 0.0 && params->MORSE_eps > 0.0) {

    /* I think we have to solve this numerically... and very crude as well */

    cnt=0;
    rad = params->MORSE_rmin - 0.69314719 / params->MORSE_alpha;
    step = -0.1 * rad;
    force=0.0;
    
    while(step != 0) {
      add1 = exp(-2.0 * params->MORSE_alpha * (rad - params->MORSE_rmin));
      add2 = exp( -params->MORSE_alpha * (rad - params->MORSE_rmin));
      force   = -params->MORSE_eps * 2.0 * params->MORSE_alpha * (add2 - add1) / rad;

      if((step < 0 && force_cap < force) || (step > 0 && force_cap > force)) {
        step = - (step/2.0); 
      }
      if(fabs(force-force_cap) < 1.0e-6) step=0;
      rad += step; cnt++;
    } 
          params->MORSE_capradius = rad;
        }
        else {
    params->MORSE_capradius = 0.0; 
        }
        FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
          this_node,i,j,rad,force,cnt));
      }
    }
  }
}

#endif

