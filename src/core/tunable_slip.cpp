/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2009,2010 
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
#include "tunable_slip.hpp"
#include "random.hpp"
#include "communication.hpp"

#ifdef TUNABLE_SLIP

int tunable_slip_set_params(int part_type_a, int part_type_b,
			    double temp, double gamma, double r_cut,
			    double time, double vx, double vy, double vz)
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);
  
  if (!data) return ES_ERROR;
  
  /* TUNABLE SLIP should be symmetrical! */
  data->TUNABLE_SLIP_temp     = temp;
  data->TUNABLE_SLIP_gamma    = gamma;
  data->TUNABLE_SLIP_r_cut    = r_cut;
  data->TUNABLE_SLIP_time     = time;
  data->TUNABLE_SLIP_vx       = vx;
  data->TUNABLE_SLIP_vy       = vy;
  data->TUNABLE_SLIP_vz       = vz;
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  
  return ES_OK;
}


void add_tunable_slip_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params, double d[3], double dist, double force[3])
{
  double gamma_t=0.0;
  double gamma_t_sqrt=0.0;
  double pre_diss=0.0;
  double pre_rand=0.0;
  double scale_vx=0.0;
  double scale_vy=0.0;
  double scale_vz=0.0; 

  if(dist < ia_params->TUNABLE_SLIP_r_cut) {
    gamma_t += 1-(dist/(ia_params->TUNABLE_SLIP_r_cut));
   }else{
    gamma_t += 0.0;
   }
  gamma_t_sqrt += sqrt(gamma_t);
  pre_diss += ((ia_params->TUNABLE_SLIP_gamma)/ia_params->TUNABLE_SLIP_time);
  pre_rand += sqrt((24.0*ia_params->TUNABLE_SLIP_temp*ia_params->TUNABLE_SLIP_gamma)/ia_params->TUNABLE_SLIP_time);

  /* Rescaling of reference velocity */ 
  scale_vx += ia_params->TUNABLE_SLIP_vx*ia_params->TUNABLE_SLIP_time;
  scale_vy += ia_params->TUNABLE_SLIP_vy*ia_params->TUNABLE_SLIP_time;
  scale_vz += ia_params->TUNABLE_SLIP_vz*ia_params->TUNABLE_SLIP_time;

  force[0] += -(pre_diss*gamma_t*(p1->m.v[0]-scale_vx))+(pre_rand*gamma_t_sqrt*(d_random()-0.5));
  force[1] += -(pre_diss*gamma_t*(p1->m.v[1]-scale_vy))+(pre_rand*gamma_t_sqrt*(d_random()-0.5));
  force[2] += -(pre_diss*gamma_t*(p1->m.v[2]-scale_vz))+(pre_rand*gamma_t_sqrt*(d_random()-0.5));
    
  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: TUNABLE_SLIP   f = (%.3e,%.3e,%.3e) with part id=%d\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: TUNABLE_SLIP   f = (%.3e,%.3e,%.3e) with part id=%d\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity));
 
}

#endif
