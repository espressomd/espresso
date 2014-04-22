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
#include "utils.hpp"

#ifdef LENNARD_JONES
#include "lj.hpp"
#include "mol_cut.hpp"
#include "communication.hpp"

/** set the force cap for the LJ interaction.
    @param ljforcecap the maximal force, 0 to disable, -1 for individual cutoff
    for each of the interactions.
*/
int ljforcecap_set_params(double ljforcecap)
{
  mpi_cap_forces(ljforcecap);
  
  return ES_OK;
}

int lennard_jones_set_params(int part_type_a, int part_type_b,
				      double eps, double sig, double cut,
				      double shift, double offset,
				      double cap_radius, double min)
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  data->LJ_eps    = eps;
  data->LJ_sig    = sig;
  data->LJ_cut    = cut;
  data->LJ_shift  = shift;
  data->LJ_offset = offset;
  if (cap_radius > 0) {
    data->LJ_capradius = cap_radius;
  }
  if (min > 0) {
    data->LJ_min = min;
  }
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  mpi_cap_forces(force_cap);

  return ES_OK;
}

/** calculate lj_capradius from force_cap */
void calc_lj_cap_radii()
{
  /* do not compute cap radii if force capping is "individual" */
  if( force_cap != -1.0 ){
    int i,j,cnt=0;
    IA_parameters *params;
    double force=0.0, rad=0.0, step, frac2, frac6;

    for(i=0; i<n_particle_types; i++) {
      for(j=0; j<n_particle_types; j++) {
        params = get_ia_param(i,j);
        if(force_cap > 0.0 && params->LJ_eps > 0.0) {
    /* I think we have to solve this numerically... and very crude as well */
    cnt=0;
    rad = params->LJ_sig;
    step = -0.1 * params->LJ_sig;
    force=0.0;
    
    while(step != 0) {
      frac2 = SQR(params->LJ_sig/rad);
      frac6 = frac2*frac2*frac2;
      force = 48.0 * params->LJ_eps * frac6*(frac6 - 0.5)/rad;
      if((step < 0 && force_cap < force) || (step > 0 && force_cap > force)) {
        step = - (step/2.0); 
      }
      if(fabs(force-force_cap) < 1.0e-6) step=0;
      rad += step; cnt++;
    } 
          params->LJ_capradius = rad;
        }
        else {
    params->LJ_capradius = 0.0; 
        }
        FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
          this_node,i,j,rad,force,cnt));
      }
    }
  }
}

#endif /* ifdef LENNARD_JONES */
