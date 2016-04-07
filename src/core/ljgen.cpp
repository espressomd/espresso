/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

/** \file ljgen.cpp Routines to calculate the generalized lennard jones
 *  energy and/or force for a particle pair. "Generalized" here means
 *  that the LJ energy is of the form
 *
 *  eps * [ b1 * (sigma/(r-r_offset))^a1 - b2 * (sigma/(r-r_offset))^a2 + shift]
 *
 *  \ref forces.cpp
*/

#include "config.hpp"

#ifdef LENNARD_JONES_GENERIC

// we use their force cap
#include "lj.hpp"
#include "ljgen.hpp"
#include "communication.hpp"

int ljgen_set_params(int part_type_a, int part_type_b,
			       double eps, double sig, double cut,
			       double shift, double offset,
			       int a1, int a2, double b1, double b2,
			       double cap_radius
#ifdef LJGEN_SOFTCORE
             , double lambda, double softrad
#endif
             )
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);
  
  if (!data) return ES_ERROR;

  data->LJGEN_eps    = eps;
  data->LJGEN_sig    = sig;
  data->LJGEN_cut    = cut;
  data->LJGEN_shift  = shift;
  data->LJGEN_offset = offset;
  data->LJGEN_a1     = a1;
  data->LJGEN_a2     = a2;
  data->LJGEN_b1     = b1;
  data->LJGEN_b2     = b2;
  if (cap_radius > 0) 
    data->LJGEN_capradius = cap_radius;
#ifdef LJGEN_SOFTCORE
  if (lambda >= 0.0 && lambda <= 1.0)
    data->LJGEN_lambda  = lambda;
  if (softrad >=0.0)
    data->LJGEN_softrad = softrad;  
#endif

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  mpi_cap_forces(force_cap);

  return ES_OK;
}

/** calculate lj_capradius from force_cap */
void calc_ljgen_cap_radii()
{
  /* do not compute cap radii if force capping is "individual" */
  if( force_cap != -1.0){
    int i,j,cnt=0;
    IA_parameters *params;
    
    for(i=0; i<n_particle_types; i++) {
      for(j=0; j<n_particle_types; j++) {
        params = get_ia_param(i,j);
        FORCE_TRACE(fprintf(stderr, "calc_ljgen_cap_radii(), [%d, %d]: LJGEN_sig = %e, LJGEN_eps = %e, force_cap = %e, exponents (%d, %d)\n",
                           i, j, params->LJGEN_sig, params->LJGEN_eps, force_cap, params->LJGEN_a1, params->LJGEN_a2));
            
        if(force_cap > 0.0 && params->LJGEN_eps > 0.0) {
          const double dx = params->LJGEN_cut / 100000.;
          
          for(double rad = params->LJGEN_cut; rad > 0.; rad -=dx) {                  
            const double frac = params->LJGEN_sig/rad;
            const double force =  params->LJGEN_eps 
#ifdef LJGEN_SOFTCORE
                * params->LJGEN_lambda
#endif
                * (params->LJGEN_b1 * params->LJGEN_a1 * pow(frac, params->LJGEN_a1) 
                   - params->LJGEN_b2 * params->LJGEN_a2 * pow(frac, params->LJGEN_a2))/rad;

            if(fabs(force) < force_cap) {
              params->LJGEN_capradius = rad;
              break;
            }
          }
        } else {
          params->LJGEN_capradius = 0.0; 
        }
        FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f\n",
                            this_node,i,j,params->LJGEN_capradius));
      }
    }
  }
}

#endif
