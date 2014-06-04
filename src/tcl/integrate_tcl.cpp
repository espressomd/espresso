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

/** \file integrate.cpp   Molecular dynamics integrator.
 *
 *  For more information about the integrator 
 *  see \ref integrate.hpp "integrate.h".
*/

#include "integrate.hpp"
#include "npt.hpp"
#include "interaction_data.hpp"
#include "lb.hpp"
#include "pressure_tcl.hpp"
#include "communication.hpp"
#include "parser.hpp"
#include "statistics_correlation.hpp"
#include "statistics_observable.hpp"

int tclcommand_invalidate_system(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  Tcl_AppendResult(interp, "invalidate_system is no longer supported ",
                   "as it's behavior is actually undefined.\n",
                   "For consistency after reading a checkpoint, use sort_particles.\n",
                   "For timing purposes, use the recalc_forces flag of integrate\n", (char *)NULL);
  return TCL_ERROR;
}

/**  Hand over integrate usage information to tcl interpreter. */
int tclcommand_integrate_print_usage(Tcl_Interp *interp) 
{
  Tcl_AppendResult(interp, "Usage of tcl-command integrate:\n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate [reuse_forces] <INT n steps>' for integrating n steps and reusing unconditionally the given forces for the first step \n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate set' for printing integrator status \n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate set nvt' for enabling NVT integration or \n" , (char *)NULL);
#ifdef NPT
  Tcl_AppendResult(interp, "'integrate set npt_isotropic <DOUBLE p_ext> [<DOUBLE piston>] [<INT, INT, INT system_geometry>] [-cubic_box]' for enabling isotropic NPT integration \n" , (char *)NULL);
#endif
  return (TCL_ERROR);
}

/** Hand over integrate status information to tcl interpreter. */
int tclcommand_integrate_print_status(Tcl_Interp *interp) 
{
  int i;
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];
  switch (integ_switch) {
  case INTEG_METHOD_NVT:
    Tcl_AppendResult(interp, "{ set nvt }", (char *)NULL);
    return (TCL_OK);
  case INTEG_METHOD_NPT_ISO:
    Tcl_PrintDouble(interp, nptiso.p_ext, buffer);
    Tcl_AppendResult(interp, "{ set npt_isotropic ", buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nptiso.piston, buffer);
    Tcl_AppendResult(interp, " ",buffer, (char *)NULL);
    for ( i = 0 ; i < 3 ; i++){
      if ( nptiso.geometry & nptiso.nptgeom_dir[i] ) {
	sprintf(buffer, " %d", 1 );
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      } else {
	sprintf(buffer, " %d", 0 );
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
    }
    if ( nptiso.cubic_box ) {
      Tcl_AppendResult(interp, " -cubic_box", (char *)NULL);
    }
    Tcl_AppendResult(interp, " } ", (char *)NULL);

    return (TCL_OK);
  }
  return (TCL_ERROR);
}

/** Parse integrate nvt command */
int tclcommand_integrate_set_nvt(Tcl_Interp *interp, int argc, char **argv)
{
  integ_switch = INTEG_METHOD_NVT;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
  return (TCL_OK);
}

/** Parse integrate npt_isotropic command */
int tclcommand_integrate_set_npt_isotropic(Tcl_Interp *interp, int argc, char **argv)
{
  int xdir, ydir, zdir;
  xdir = ydir = zdir = nptiso.cubic_box = 0;

  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args: \n", (char *)NULL);
    return tclcommand_integrate_print_usage(interp);
  }  
  /* set parameters p_ext and piston */
  if ( !ARG_IS_D(3, nptiso.p_ext) )  return tclcommand_integrate_print_usage(interp);
  tclcallback_p_ext(interp, &nptiso.p_ext);
  if ( argc > 4 ) { 
    if(!ARG_IS_D(4, nptiso.piston) ) return tclcommand_integrate_print_usage(interp);
    tclcallback_npt_piston(interp, &nptiso.piston); }
  else if ( nptiso.piston <= 0.0 ) {
    Tcl_AppendResult(interp, "You must set <piston> as well before you can use this integrator! \n", (char *)NULL);
    return tclcommand_integrate_print_usage(interp);
  }

  if ( argc > 5 ) {
    if (!ARG_IS_I(5,xdir) || !ARG_IS_I(6,ydir) || !ARG_IS_I(7,zdir) ) {
      return tclcommand_integrate_print_usage(interp);}
    else {
      /* set the geometry to include rescaling specified directions only*/
      nptiso.geometry = 0; nptiso.dimension = 0; nptiso.non_const_dim = -1;
      if ( xdir ) { 
	nptiso.geometry = ( nptiso.geometry | NPTGEOM_XDIR ); 
	nptiso.dimension += 1;
	nptiso.non_const_dim = 0;
      }
      if ( ydir ) { 
	nptiso.geometry = ( nptiso.geometry | NPTGEOM_YDIR );
	nptiso.dimension += 1;
	nptiso.non_const_dim = 1;
      }
      if ( zdir ) { 
	nptiso.geometry = ( nptiso.geometry | NPTGEOM_ZDIR );
	nptiso.dimension += 1;
	nptiso.non_const_dim = 2;
      }
    }
  } else {
    /* set the geometry to include rescaling in all directions; the default*/
    nptiso.geometry = 0;
    nptiso.geometry = ( nptiso.geometry | NPTGEOM_XDIR );
    nptiso.geometry = ( nptiso.geometry | NPTGEOM_YDIR );
    nptiso.geometry = ( nptiso.geometry | NPTGEOM_ZDIR );
    nptiso.dimension = 3; nptiso.non_const_dim = 2;
  }

  if ( argc > 8 ) {
    /* enable if the volume fluctuations should also apply to dimensions which are switched off by the above flags
       and which do not contribute to the pressure (3D) / tension (2D, 1D) */
    if (!ARG_IS_S(8,"-cubic_box")) {
      return tclcommand_integrate_print_usage(interp);
    } else {
      nptiso.cubic_box = 1;
    }
  }

  /* Sanity Checks */
#ifdef ELECTROSTATICS      
  if ( nptiso.dimension < 3 && !nptiso.cubic_box && coulomb.bjerrum > 0 ){
    fprintf(stderr,"WARNING: If electrostatics is being used you must use the -cubic_box option!\n");
    fprintf(stderr,"Automatically reverting to a cubic box for npt integration.\n");
    fprintf(stderr,"Be aware though that all of the coulombic pressure is added to the x-direction only!\n");
    nptiso.cubic_box = 1;
  }
#endif

#ifdef DIPOLES     
  if ( nptiso.dimension < 3 && !nptiso.cubic_box && coulomb.Dbjerrum > 0 ){
    fprintf(stderr,"WARNING: If magnetostatics is being used you must use the -cubic_box option!\n");
    fprintf(stderr,"Automatically reverting to a cubic box for npt integration.\n");
    fprintf(stderr,"Be aware though that all of the magnetostatic pressure is added to the x-direction only!\n");
    nptiso.cubic_box = 1;
  }
#endif


  if( nptiso.dimension == 0 || nptiso.non_const_dim == -1) {
    Tcl_AppendResult(interp, "You must enable at least one of the x y z components as fluctuating dimension(s) for box length motion!", (char *)NULL);
    Tcl_AppendResult(interp, "Cannot proceed with npt_isotropic, reverting to nvt integration... \n", (char *)NULL);
    integ_switch = INTEG_METHOD_NVT;
    mpi_bcast_parameter(FIELD_INTEG_SWITCH);
    return (TCL_ERROR);
  }

  /* set integrator switch */
  integ_switch = INTEG_METHOD_NPT_ISO;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);

  /* broadcast npt geometry information to all nodes */
  mpi_bcast_nptiso_geom();
  return (TCL_OK);
}

int tclcommand_integrate(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int  n_steps, reuse_forces = 0;
  
  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if (argc < 1) {
    Tcl_AppendResult(interp, "wrong # args: \n\"", (char *) NULL);
    return tclcommand_integrate_print_usage(interp);  }
  else if (argc < 2) {                    return tclcommand_integrate_print_status(interp); }

  if (ARG1_IS_S("set")) {
    if      (argc < 3)                    return tclcommand_integrate_print_status(interp);
    if      (ARG_IS_S(2,"nvt"))           return tclcommand_integrate_set_nvt(interp, argc, argv);
#ifdef NPT
    else if (ARG_IS_S(2,"npt_isotropic")) return tclcommand_integrate_set_npt_isotropic(interp, argc, argv);
#endif
    else {
      Tcl_AppendResult(interp, "unknown integrator method:\n", (char *)NULL);
      return tclcommand_integrate_print_usage(interp);
    }
  } else {
    if ( !ARG_IS_I(1,n_steps) ) return tclcommand_integrate_print_usage(interp);

    // actual integration
    if ((argc == 3) && ARG_IS_S(2, "reuse_forces")) {
      reuse_forces = 1;
    }
    else if ((argc == 3) && ARG_IS_S(2, "recalc_forces")) {
      reuse_forces = -1;
    }
    else if (argc != 2) return tclcommand_integrate_print_usage(interp);
  }
  /* go on with integrate <n_steps> */
  if(n_steps < 0) {
    Tcl_AppendResult(interp, "illegal number of steps (must be >0) \n", (char *) NULL);
    return tclcommand_integrate_print_usage(interp);;
  }
  /* perform integration */
  if (!correlations_autoupdate && !observables_autoupdate) {
    if (mpi_integrate(n_steps, reuse_forces))
      return gather_runtime_errors(interp, TCL_OK);
  } else  {
    for (int i=0; i<n_steps; i++) {
      if (mpi_integrate(1, reuse_forces))
        return gather_runtime_errors(interp, TCL_OK);
      autoupdate_observables();
      autoupdate_correlations();
    }
  }
  return TCL_OK;
}

/************************************************************/
/* Callback functions */
 
int tclcallback_skin(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "skin must be positive.", (char *) NULL);
    return (TCL_ERROR);
  }
  skin = data;
  mpi_bcast_parameter(FIELD_SKIN);
  return (TCL_OK);
}

int tclcallback_time_step(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
#ifdef LB_GPU
  float ts = (float)data;
#endif
  if (data < 0.0) {
    Tcl_AppendResult(interp, "time step must be positive.", (char *) NULL);
    return (TCL_ERROR);
  }
#ifdef LB
  else if ((lbpar.tau >= 0.0) && (data > lbpar.tau)) {
    Tcl_AppendResult(interp, "MD time step must be smaller than LB time step.", (char *)NULL);
    return (TCL_ERROR);
  }
#endif
#ifdef LB_GPU	
  else if ((lbpar_gpu.tau >= 0.0) && (ts > lbpar_gpu.tau)) { 
    Tcl_AppendResult(interp, "MD time step must be smaller than LB time step.", (char *)NULL);
    return (TCL_ERROR);
  }
#endif
  mpi_set_time_step(data);

  return (TCL_OK);
}

int tclcallback_time(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  sim_time = data;
  mpi_bcast_parameter(FIELD_SIMTIME);
  return (TCL_OK);
}
