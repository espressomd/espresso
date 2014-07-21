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
/** \file buckingham.cpp
 *
 *  Implementation of \ref buckingham.hpp
 */
#include "buckingham.hpp"

#ifdef BUCKINGHAM
#include "communication.hpp"

int buckingham_set_params(int part_type_a, int part_type_b,
			  double A, double B, double C, double D, double cut,
			  double discont, double shift, double cap_radius,
			  double F1, double F2)
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  data->BUCK_A       = A;
  data->BUCK_B       = B;
  data->BUCK_C       = C;
  data->BUCK_D       = D;
  data->BUCK_cut     = cut;
  data->BUCK_discont = discont;
  data->BUCK_shift   = shift;
  if (cap_radius > 0.0) {
    data->BUCK_capradius = cap_radius;
  }

  /* Replace the buckingham potential for interatomic dist. less
     than or equal to discontinuity by a straight line (F1+F2*r) */
  F1 = buck_energy_r(A, B, C, D, shift, discont) +
    discont*buck_force_r(A, B, C, D, discont);
  F2 = -buck_force_r(A, B, C, D, discont);

  data->BUCK_F1 = F1;
  data->BUCK_F2 = F2;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  mpi_cap_forces(force_cap);

  return ES_OK;
}

/** calculate buck_capradius from force_cap */
void calc_buck_cap_radii()
{
  /* do not compute cap radii if force capping is individual */
  if( force_cap != -1.0 ){
    int i,j,cnt=0;
    IA_parameters *params = NULL;
    double force=0.0, frac2, frac6, frac8;
    double r0,r1 = 0.0,diff,exp_term,C_R,D_R;
    for(i=0; i<n_particle_types; i++) {
      for(j=0; j<n_particle_types; j++) {
        params = get_ia_param(i,j);
        if(force_cap>0.0 ) {
          force = -params->BUCK_F2;
          if (force_cap<force)
    {
      /* Solving numerically using Newton Raphson Technique */
      force = 0.0;
      cnt = 0;
      r1 = (params->BUCK_cut+params->BUCK_discont)/2.0; //First guess value
      r0 = 0.0 ;
      while(fabs(r1 - r0)>1.0e-10) {
         r0 = r1;
         exp_term = params->BUCK_A*params->BUCK_B*exp(-params->BUCK_B*r0);
         frac2 = SQR(r0);
         frac6 = frac2*frac2*frac2;
         frac8 = frac6*frac2;
         C_R = 6.0*params->BUCK_C/frac8;
         D_R = 4.0*params->BUCK_D/frac6;
         diff = (exp_term - C_R*r0 - D_R*r0 - force_cap)/(-params->BUCK_B*exp_term + 7.0*C_R + 5.0*D_R);
         r1 = r0 - diff;
         if(r1>params->BUCK_discont)
            r1=0.5*(params->BUCK_discont+r0);
         cnt++;
         if(cnt>500)
         {
           fprintf(stderr,"%d: ERROR@buckingham.h: Failed to converge while determining Buckingham cap radius!!",this_node);
           fprintf(stderr,"%d: tolerance = %f",this_node, diff);
                 exit (0);
         }
            }
      frac2 = SQR(r1);
      frac6 = frac2*frac2*frac2;
      force = params->BUCK_A*params->BUCK_B*exp(-params->BUCK_B*r1) - 6.0*params->BUCK_C/(r1*frac6) - 4.0*params->BUCK_D*r1/(frac6);
            params->BUCK_capradius = r1;
    }
    else
       params->BUCK_capradius = params->BUCK_discont;

        }
        else
    params->BUCK_capradius = 0.0;

        FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
          this_node,i,j,r1,force,cnt));
      }
    }
    BUCK_TRACE(fprintf(stderr,"%d: BUCK: Buckingham force cap imposed %f, Calculated force %f and Cap radius %f after %d iterations\n",this_node,force_cap,force,params->BUCK_capradius,cnt);
    );
  }
}

#endif
