/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file nemd.cpp

    For more information see \ref nemd.hpp
 */
#include "nemd.hpp"
#include <cstdio>
#include "integrate.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

/************************************************************/


int nemd_method = NEMD_METHOD_OFF;

#ifdef NEMD
Nemd nemddata = { -1, 0, 0, 0.0, 1.0, NULL, 0, 0, NULL, 0.0, 0.0, 0.0, 0};

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Free all associated memory. */
int nemd_free(void) 
{
  INTEG_TRACE(fprintf(stderr,"%d: nemd_free\n",this_node));
  if(nemddata.n_exchange > 0) {
    free(nemddata.slab[nemddata.mid_slab].fastest);
    free(nemddata.slab[nemddata.top_slab].fastest);
  }
  free(nemddata.velocity_profile);
  free(nemddata.slab);
  nemddata.n_slabs          = -1;
  nemddata.mid_slab         = 0;
  nemddata.top_slab         = 0;
  nemddata.thickness        = 0.0;
  nemddata.invthickness     = 1.0;
  nemddata.velocity_profile = NULL;
  nemddata.profile_norm     = 0;
  nemddata.n_exchange       = 0;
  nemddata.slab             = NULL;
  nemddata.shear_rate       = 0.0;
  nemddata.slab_vel         = 0.0;
  return ES_OK;
}

/** Initialize all data structures for nemd. */
void nemd_init(int n_slabs, int n_exchange, double shear_rate) 
{
  int i;
  
  INTEG_TRACE(fprintf(stderr,"%d: nemd_init: n_slabs=%d n_exchange=%d\n",this_node, n_slabs, n_exchange));

  /* check node grid */
  if( n_nodes > 1 ) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{037 NEMD is a single node feature} ");
    return;
  }

  /* first free old structures befor initializing new ones */
  if(nemddata.n_slabs > -1) nemd_free();
  /* exit nemd integration */
  if(n_slabs == 0) return;

  /* fill nemd structure */
  nemddata.n_slabs          = n_slabs;
  nemddata.top_slab         = 0;
  nemddata.mid_slab         = n_slabs/2;

  nemddata.thickness        = box_l[2]/(double)nemddata.n_slabs;
  nemddata.invthickness     = 1.0 / nemddata.thickness;

  nemddata.shear_rate       = shear_rate;
  nemddata.slab_vel         = time_step*shear_rate*box_l[2]/4.0;

  nemddata.n_exchange       = n_exchange;
 
  nemddata.slab             = (Slab *)malloc(nemddata.n_slabs*sizeof(Slab));
  nemddata.velocity_profile = (double *)malloc(nemddata.n_slabs*sizeof(double));

  nemddata.momentum = 0.0;
  nemddata.momentum_norm = 0;

  /* initialize slabs and velocity profile */
  for(i=0;i<nemddata.n_slabs;i++) {
    nemddata.velocity_profile[i]     = 0.0;

    nemddata.slab[i].v_mean          = 0.0;
    nemddata.slab[i].n_parts_in_slab = 0;
    nemddata.slab[i].v_min           = 0.0;
    nemddata.slab[i].ind_min         = 0;
    nemddata.slab[i].fastest         = NULL;
    nemddata.slab[i].n_fastest       = 0;
    nemddata.slab[i].vel_diff        = 0.0;
  }
  /* allocate arrays for indices of fastest particles in slab */
  if(nemddata.n_exchange > 0) {
    nemddata.slab[nemddata.top_slab].fastest = (int *)malloc(nemddata.n_exchange*sizeof(int));
    nemddata.slab[nemddata.mid_slab].fastest = (int *)malloc(nemddata.n_exchange*sizeof(int));
  }
  for(i=0;i<nemddata.n_exchange;i++) {
    nemddata.slab[nemddata.top_slab].fastest[i] = -1;
    nemddata.slab[nemddata.mid_slab].fastest[i] = -1;
  }
  nemddata.slab[nemddata.top_slab].v_min   = -1e10;
  nemddata.slab[nemddata.mid_slab].v_min   = +1e10;
}



/*@}*/


void nemd_change_momentum() 
{
  int i;
  double tmp_v0;
  Slab *mid_slab, *top_slab;

  if(nemd_method == NEMD_METHOD_OFF) return;

  INTEG_TRACE(fprintf(stderr,"%d: nemd_change_momentum: Method %d\n",this_node,nemd_method));

  mid_slab = &nemddata.slab[nemddata.mid_slab];
  top_slab = &nemddata.slab[nemddata.top_slab];

  if(nemd_method == NEMD_METHOD_EXCHANGE ) {
    /* exit if there are not enough particles */
    INTEG_TRACE(fprintf(stderr,"%d: parts_in_slabs: top %d mid %d\n",this_node,top_slab->n_parts_in_slab,mid_slab->n_parts_in_slab));
    if(mid_slab->n_fastest != nemddata.n_exchange || 
       top_slab->n_fastest != nemddata.n_exchange) {
      char *errtxt = runtime_error(128 + ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{038 nemd_exchange_momentum: Not enough particles in slab!} ");
      /* cannot continue */
      return;
    }

    /* perform momentum exchange */
    for(i=0;i<nemddata.n_exchange;i++) {
      /* store momentum change */
      nemddata.momentum += local_particles[top_slab->fastest[i]]->m.v[0]; 
      nemddata.momentum -= local_particles[mid_slab->fastest[i]]->m.v[0]; 
      tmp_v0 = local_particles[mid_slab->fastest[i]]->m.v[0];
      local_particles[mid_slab->fastest[i]]->m.v[0] = local_particles[top_slab->fastest[i]]->m.v[0];
      local_particles[top_slab->fastest[i]]->m.v[0] = tmp_v0;
    }

    /* prepare next round */
    top_slab->n_fastest = 0;
    top_slab->v_min     = -1e10;
    mid_slab->n_fastest = 0;
    mid_slab->v_min     = 1e10;
  }
  else if (nemd_method ==  NEMD_METHOD_SHEARRATE) {
    double vel_mean; 
    /* calculate velocity difference vel_diff = vel_required - vel_actual */
    vel_mean = top_slab->v_mean/(double)top_slab->n_parts_in_slab;
    top_slab->vel_diff = nemddata.slab_vel - vel_mean;
    vel_mean = mid_slab->v_mean/(double)mid_slab->n_parts_in_slab;
    mid_slab->vel_diff = - (nemddata.slab_vel + vel_mean);
    /* store momentum change */
    nemddata.momentum += top_slab->n_parts_in_slab*top_slab->vel_diff;
    nemddata.momentum -= mid_slab->n_parts_in_slab*mid_slab->vel_diff;
  }
  nemddata.momentum_norm ++;
}

/************************************************************/

void nemd_store_velocity_profile() 
{
  int i;

  if(nemddata.n_slabs == -1) return;
  INTEG_TRACE(fprintf(stderr,"%d: nemd_store_velocity_profile:\n",this_node));

  for(i=0;i<nemddata.n_slabs;i++) {
    if(nemddata.slab[i].n_parts_in_slab == 0) fprintf(stderr,"Zero particles in slab %d!!!\n",i);
    nemddata.velocity_profile[i]     += nemddata.slab[i].v_mean/(double)nemddata.slab[i].n_parts_in_slab;
    nemddata.slab[i].v_mean          = 0.0;
    nemddata.slab[i].n_parts_in_slab = 0;
  }
  nemddata.profile_norm++;
}

/************************************************************/

#endif

