// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file nemd.c

    For more information see \ref nemd.h

    <b>Responsible:</b>
    <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 */
#include "nemd.h"
#include "integrate.h"
#include "communication.h"

/************************************************************/


int nemd_method = NEMD_METHOD_OFF;

#ifdef NEMD
Nemd nemddata = { -1, 0, 0, 0.0, 1.0, NULL, 0, 0, NULL, 0.0, 0.0, 0.0, 0};
#endif

/************************************************************/

/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Hand over nemd usage information to tcl interpreter. */
int nemd_usage(Tcl_Interp *interp) 
{
#ifdef NEMD
  Tcl_AppendResult(interp, "Usage of tcl command nemd:\n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd\" for returning the status or \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd off\" \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd exchange <INT n_slabs> <INT n_exchange>\" \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd shearrate <INT n_slabs> <DOUBLE shearrate>\" \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd profile\" for returning the velocity profile \n", (char *)NULL);
  Tcl_AppendResult(interp, "\"nemd viscosity\" for returning the viscosity \n", (char *)NULL);
  return (TCL_ERROR);
#else
  Tcl_AppendResult(interp, "nemd not compiled in!", (char *)NULL);
  return (TCL_ERROR);
#endif
}

#ifdef NEMD
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
  return (TCL_OK);
}

/** Initialize all data structures for nemd. */
void nemd_init(int n_slabs, int n_exchange, double shear_rate) 
{
  int i;
  
  INTEG_TRACE(fprintf(stderr,"%d: nemd_init: n_slabs=%d n_exchange=%d\n",this_node, n_slabs, n_exchange));

  /* check node grid */
  if( n_nodes > 0 ) {
    char *errtxt = runtime_error(128);
    sprintf(errtxt, "{NEMD is a single node feature} ");
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

/** Hand over nemd status information to tcl interpreter. */
int nemd_print_status(Tcl_Interp *interp) 
{
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];
  switch (nemd_method) {
  case NEMD_METHOD_OFF:
    Tcl_AppendResult(interp, "off", (char *)NULL);
    return (TCL_OK);
    break;
  case NEMD_METHOD_EXCHANGE:
    sprintf(buffer, "%d", nemddata.n_slabs);
    Tcl_AppendResult(interp, "exchange ",buffer, (char *)NULL);
    sprintf(buffer, "%d", nemddata.n_exchange);
    Tcl_AppendResult(interp, " ",buffer, (char *)NULL);
    return (TCL_OK);
    break;
  case NEMD_METHOD_SHEARRATE:
    sprintf(buffer, "%d", nemddata.n_slabs);
    Tcl_AppendResult(interp, "shearrate ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nemddata.shear_rate, buffer);
    Tcl_AppendResult(interp, " ",buffer, (char *)NULL);
    return (TCL_OK);
    break;
  default:
    return (TCL_ERROR);
  }
  return (TCL_ERROR);
}

/** Set nemd method to exchange and set nemd parameters */
int nemd_set_exchange(Tcl_Interp *interp, int argc, char **argv) 
{
  int n_slabs, n_exchange;
  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args:  ", (char *)NULL);
    return nemd_usage(interp);
  }
  if ( !ARG_IS_I(2, n_slabs) || !ARG_IS_I(3, n_exchange) ) {
    Tcl_AppendResult(interp, "wrong argument type:  ", (char *)NULL);
    return nemd_usage(interp);
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

  nemd_method = NEMD_METHOD_EXCHANGE;
  nemd_init(n_slabs, n_exchange, 0.0);
  return (TCL_OK);
}

/** Set nemd method to shearrate and set nemd parameters */
int nemd_set_shearrate(Tcl_Interp *interp, int argc, char **argv) 
{
  int n_slabs;
  double shearrate;
  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args:  ", (char *)NULL);
    return nemd_usage(interp);
  }
  if ( !ARG_IS_I(2, n_slabs) || !ARG_IS_D(3, shearrate) ) {
    Tcl_AppendResult(interp, "wrong argument type:  ", (char *)NULL);
    return nemd_usage(interp);
  }
 /* parameter sanity */
  if ( n_slabs<0 || n_slabs%2!=0 ) {
    Tcl_AppendResult(interp, "nemd <n_slabs> must be non negative and even!",(char *)NULL);
    return (TCL_ERROR);
  }  
  nemd_method = NEMD_METHOD_SHEARRATE;
  nemd_init(n_slabs, 0, shearrate);
  return (TCL_OK);
}

/** Hand over velocity profile to tcl interpreter */
int nemd_print_profile(Tcl_Interp *interp)
{
  int i;
  double val;
  char buffer[TCL_DOUBLE_SPACE];
  
  INTEG_TRACE(fprintf(stderr,"%d: nemd_print_profile:\n",this_node));
  if(nemd_method == NEMD_METHOD_OFF) {
    Tcl_AppendResult(interp, "nemd is off", (char *)NULL);
    return (TCL_OK);
  }
  
  for(i=0;i<nemddata.n_slabs;i++) {
    /* note: output velocities as usual have to be resacled by 1/time_step! */
    val = nemddata.velocity_profile[i]/(nemddata.profile_norm*time_step);
    Tcl_PrintDouble(interp, val, buffer);
    Tcl_AppendResult(interp," ", buffer, (char *)NULL);
    
    nemddata.velocity_profile[i] = 0.0;
  }
  
  nemddata.profile_norm = 0;
  return (TCL_OK);
}

int nemd_print_viscosity(Tcl_Interp *interp)
{
  double shear_rate=0.0, mean_force, viscosity;
  char buffer[TCL_DOUBLE_SPACE];
  INTEG_TRACE(fprintf(stderr,"%d: nemd_print_viscosity:\n",this_node));

  /* calculate shear_rate */
  switch (nemd_method) {
  case NEMD_METHOD_OFF:
    Tcl_AppendResult(interp, "nemd is off", (char *)NULL);
    return (TCL_OK);
  case NEMD_METHOD_EXCHANGE:
    shear_rate = 1.0;
  case NEMD_METHOD_SHEARRATE:
    shear_rate = nemddata.shear_rate;
  }
  /* rescale momentum exchange (vel_internal = time_step * vel) */
  nemddata.momentum /= time_step;
  /* Calculate average Force := Momentum transfer per time unit */
  mean_force = nemddata.momentum / (nemddata.momentum_norm*time_step);
  /* Calculate viscosity := mean_force / (shearrate * area) */
  viscosity = mean_force / (shear_rate*4.0*box_l[0]*box_l[1]);
  Tcl_PrintDouble(interp, viscosity, buffer);
  Tcl_AppendResult(interp,buffer, (char *)NULL);

  nemddata.momentum = 0.0;
  nemddata.momentum_norm = 0;
  return (TCL_OK);
}
#endif

/*@}*/

/************************************************************
 *            Exported Functions                            *
 ************************************************************/

int nemd(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
#ifdef NEMD
  int status = TCL_OK;

  INTEG_TRACE(fprintf(stderr,"%d: nemd:\n",this_node));
  Tcl_ResetResult(interp);

  /* print nemd status */
  if(argc == 1) {
    status = nemd_print_status(interp) ;
  }
  else if (ARG1_IS_S("off")) {
    nemd_method = NEMD_METHOD_OFF;
    status = nemd_free();
  }  
  else if (ARG1_IS_S("exchange")) {
    status = nemd_set_exchange(interp,argc,argv);
  } 
  else if (ARG1_IS_S("shearrate")) {
    status = nemd_set_shearrate(interp,argc,argv);
  } 
  else if (ARG1_IS_S("profile")) {
    status = nemd_print_profile(interp);
  } 
  else if (ARG1_IS_S("viscosity")) {
    status = nemd_print_viscosity(interp);
  } 
  else {
    Tcl_AppendResult(interp, "Unkwnown keyword: \n", (char *)NULL);
    return nemd_usage(interp);
  }

  return mpi_gather_runtime_errors(interp, status);

#endif
  INTEG_TRACE(fprintf(stderr,"%d: call to nemd but not compiled in!\n",this_node));
  return nemd_usage(interp);
}

/************************************************************/
#ifdef NEMD

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
      char *errtxt = runtime_error(128);
      sprintf(errtxt,"%d: nemd_exchange_momentum: Not enough particles in slab!\n",this_node);
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

