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
Nemd nemddata = { -1, 0, 0, 0.0, 1.0, 0, NULL, NULL, 0};
#endif

/************************************************************/

int nemd(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
#ifdef NEMD
  int n_slabs, n_exchange;
  INTEG_TRACE(fprintf(stderr,"%d: nemd:\n",this_node));

  //  Tcl_ResetResult(interp);

  if (argc < 1) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <n_slabs> <n_exchange> \" or \"profile\"", (char *)NULL);
    return (TCL_ERROR);
  }

  /* print velocity profile */
  if (ARG1_IS_S("profile")) {
    return nemd_print_profile(interp);
  }  
  /* set nemd parameters */
  else {
    /* check number of arguments */
    if (argc < 2) {
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
		       argv[0], " <n_slabs> <n_exchange> \"", (char *)NULL);
      return (TCL_ERROR);
    }
    /* copy parameters */
    if ( (! ARG_IS_I(1, n_slabs)) || (! ARG_IS_I(2, n_exchange)) ) {
      Tcl_AppendResult(interp, " needs 2 INTEGER parameters: "
		       "<n_slabs> <n_exchange>",(char *)NULL);
      return (TCL_ERROR);
    }
    /* parameter sanity */
    if ( n_slabs<0 || n_slabs%2!=0 ) {
      Tcl_AppendResult(interp, "nemd <n_slabs> must be non negative and even!",(char *)NULL);
      return (TCL_ERROR);
    }  
    if ( n_slabs > 0 && n_exchange < 1 ) {
      Tcl_AppendResult(interp, "nemd <n_exchange> must be positive!",(char *)NULL);
      return (TCL_ERROR);
    } 
    
    /* communicat eparameters here */
    
    nemd_init(n_slabs, n_exchange);
  }
#endif
  return (TCL_OK);
}

/************************************************************/

void nemd_init(int n_slabs, int n_exchange) 
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
  }
  /* allocate arrays for indices of fastest particles in slab */
  nemddata.slab[nemddata.top_slab].fastest = (int *)malloc(nemddata.n_exchange*sizeof(int));
  nemddata.slab[nemddata.mid_slab].fastest = (int *)malloc(nemddata.n_exchange*sizeof(int));
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
  free(nemddata.slab[nemddata.mid_slab].fastest);
  free(nemddata.slab[nemddata.top_slab].fastest);
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

void nemd_exchange_momentum() 
{
#ifdef NEMD
  int i;
  double tmp_v0;
  Slab *mid_slab, *top_slab;

  if(nemddata.n_slabs == -1) return;
  INTEG_TRACE(fprintf(stderr,"%d: nemd_exchange_momentum\n",this_node));

  mid_slab = &nemddata.slab[nemddata.mid_slab];
  top_slab = &nemddata.slab[nemddata.top_slab];

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
  
  //Tcl_AppendResult(interp, "{", (char *)NULL);
  
  for(i=0;i<nemddata.n_slabs;i++) {
    /* note: output velocities as usual have to be resacled by 1/time_step! */
    val = nemddata.velocity_profile[i]/(nemddata.profile_norm*time_step);
    Tcl_PrintDouble(interp, val, buffer);
    Tcl_AppendResult(interp," ", buffer, (char *)NULL);
    
    nemddata.velocity_profile[i] = 0.0;
  }
  
  //Tcl_AppendResult(interp, " }", (char *)NULL);
  
  nemddata.profile_norm = 0;
#endif
  return (TCL_OK);
}
