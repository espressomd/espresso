// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file nemd.c

    For more information see \ref nemd.h

    <b>Responsible:</b>
    <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 */
#include "nemd.h"

/************************************************************/

#ifdef NEMD
Nemd nemddata = { -1, 0, 0, 0.0, 1.0, NULL, 0, NEMD_METHOD_NOTSET, 0, NULL, 0.0, 0.0};
#endif

/************************************************************/

int nemd(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
#ifdef NEMD
  int n_slabs, n_exchange=0;
  double shearrate=0.0;
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];

  INTEG_TRACE(fprintf(stderr,"%d: nemd:\n",this_node));
  Tcl_ResetResult(interp);

 /* print nemd status */
  if(argc == 1) {
    sprintf(buffer, "%d", nemddata.n_slabs);
    Tcl_AppendResult(interp, "nemd ", buffer," ", (char *)NULL);
    if(nemddata.method == NEMD_METHOD_EXCHANGE) {
      sprintf(buffer, "%d", nemddata.n_exchange);
      Tcl_AppendResult(interp, "exchange ", buffer, (char *)NULL);
    }
    else if (nemddata.method == NEMD_METHOD_SHEARRATE) {
      sprintf(buffer, "%f", nemddata.shear_rate);
      Tcl_AppendResult(interp, "shearrate ", buffer, (char *)NULL);
    }
    return (TCL_OK);
  }

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
		     argv[0], " <n_slabs> <method> <value>\" or \n\"",
		     argv[0], " profile\"", (char *)NULL);
    return (TCL_ERROR);
  }

  /* print velocity profile */
  if (ARG1_IS_S("profile")) {
    return nemd_print_profile(interp);
  }  
  /* set nemd parameters */
  else {
    /* check number of arguments */
    if (argc < 3) {
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
		       argv[0], " <n_slabs> <method> <value> \"", (char *)NULL);
      return (TCL_ERROR);
    }
    /* copy parameters */
    if ( ! ARG_IS_I(1, n_slabs) ) {
      Tcl_AppendResult(interp, "first parameter of nemd is <n_slabs> of type",
		       " <INT> or \"profile\"", (char *)NULL);
      return (TCL_ERROR);
    }
    if ( ARG_IS_S(2,"exchange") ) {
      if ( ! ARG_IS_I(3, n_exchange) ) {
	Tcl_AppendResult(interp, "third parameter of nemd exchange is n_exchange",
			 " of type <INT> ", (char *)NULL);
	return (TCL_ERROR);
      }
      nemddata.method = NEMD_METHOD_EXCHANGE;
    }
    else if ( ARG_IS_S(2,"shearrate") ) {
      if ( ! ARG_IS_D(3, shearrate) ) {
	Tcl_AppendResult(interp, "third parameter of nemd shearrate is shearrate",
			 " of type <DOUBLE> ", (char *)NULL);
	return (TCL_ERROR);
      }
      nemddata.method = NEMD_METHOD_SHEARRATE;
    }

    /* parameter sanity */
    if ( n_slabs<0 || n_slabs%2!=0 ) {
      Tcl_AppendResult(interp, "nemd <n_slabs> must be non negative and even!",(char *)NULL);
      return (TCL_ERROR);
    }  
    if ( n_slabs > 0 && n_exchange < 0 ) {
      Tcl_AppendResult(interp, "nemd <n_exchange> must be positive!",(char *)NULL);
      return (TCL_ERROR);
    } 
    
    /* communicat eparameters here */
    
    nemd_init(n_slabs, n_exchange, shearrate);
  }
#endif
  return (TCL_OK);
}

/************************************************************/

void nemd_init(int n_slabs, int n_exchange, double shear_rate) 
{
#ifdef NEMD
  int i;
  
  INTEG_TRACE(fprintf(stderr,"%d: nemd_init: n_slabs=%d n_exchange=%d\n",this_node, n_slabs, n_exchange));

  /* check node grid */
  if( this_node > 0 ) {
    fprintf(stderr,"%d: NEMD is a single node feature. Exiting.\n",this_node);
    errexit();
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
#endif
}

/************************************************************/

void nemd_free() 
{
#ifdef NEMD
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
  nemddata.n_exchange       = 0;
  nemddata.slab             = NULL;
  nemddata.velocity_profile = NULL;
  nemddata.profile_norm     = 0;
#endif
}

/************************************************************/

void nemd_change_momentum() 
{
#ifdef NEMD
  int i;
  double tmp_v0;
  Slab *mid_slab, *top_slab;

  if(nemddata.n_slabs == -1) return;

  INTEG_TRACE(fprintf(stderr,"%d: nemd_change_momentum:\n",this_node));

  mid_slab = &nemddata.slab[nemddata.mid_slab];
  top_slab = &nemddata.slab[nemddata.top_slab];

  if(nemddata.method == NEMD_METHOD_EXCHANGE ) {

    if(mid_slab->n_fastest != nemddata.n_exchange || 
       top_slab->n_fastest != nemddata.n_exchange) {
      fprintf(stderr,"%d: nemd_exchange_momentum: Not enough particles in slab!\n",this_node);
      errexit();
    }

    INTEG_TRACE(fprintf(stderr,"%d: parts_in_slabs: top %d mid %d\n",this_node,top_slab->n_parts_in_slab,mid_slab->n_parts_in_slab));
    
    for(i=0;i<nemddata.n_exchange;i++) {

      INTEG_TRACE(fprintf(stderr,"%d: Exchange part[%d].m.v[0]=%f and part[%d].m.v[0]=%f\n",this_node,local_particles[mid_slab->fastest[i]]->p.identity,local_particles[mid_slab->fastest[i]]->m.v[0],local_particles[top_slab->fastest[i]]->p.identity,local_particles[top_slab->fastest[i]]->m.v[0]));

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
  else if (nemddata.method ==  NEMD_METHOD_SHEARRATE) {
    double vel_mean; 
    vel_mean = top_slab->v_mean/(double)top_slab->n_parts_in_slab;
    top_slab->vel_diff = nemddata.slab_vel - vel_mean;
    vel_mean = mid_slab->v_mean/(double)mid_slab->n_parts_in_slab;
    mid_slab->vel_diff = - nemddata.slab_vel - vel_mean;
  }
#endif
}

/************************************************************/

void nemd_store_velocity_profile() 
{
#ifdef NEMD
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
#endif
}

/************************************************************/

int nemd_print_profile(Tcl_Interp *interp)
{
#ifdef NEMD
  int i;
  double val;
  char buffer[50 + TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  
  INTEG_TRACE(fprintf(stderr,"%d: nemd_print_profile:\n",this_node));
  
  if(nemddata.n_slabs == -1) {
    Tcl_AppendResult(interp, "{ nemd is off }", (char *)NULL);
    return (TCL_OK);
  }
  
  if(nemddata.profile_norm==0) return (TCL_OK);

  for(i=0;i<nemddata.n_slabs;i++) {
    /* note: output velocities as usual have to be resacled by 1/time_step! */
    val = nemddata.velocity_profile[i]/(nemddata.profile_norm*time_step);
    Tcl_PrintDouble(interp, val, buffer);
    Tcl_AppendResult(interp,buffer, " ", (char *)NULL);
    
    nemddata.velocity_profile[i] = 0.0;
  }
  
  nemddata.profile_norm = 0;
#endif
  return (TCL_OK);
}
