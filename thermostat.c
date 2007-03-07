// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file thermostat.c
    Implementation of \ref thermostat.h "thermostat.h"
 */
#include "thermostat.h"
#include "communication.h"
#include "lattice.h"

/* thermostat switch */
int thermo_switch = THERMO_OFF;
/** Temperature */
double temperature = -1.0;

/* LANGEVIN THERMOSTAT */
/* Langevin friction coefficient gamma. */
double langevin_gamma = 0.0;
/* Friction coefficient gamma for rotation */
double langevin_gamma_rotation;

/* DPD THERMOSTAT */
/* DPD longitudinal friction coefficient gamma. */
double dpd_gamma = 0.0;
/* DPD thermostat cutoff */
double dpd_r_cut = 0.0;
/* inverse off DPD thermostat cutoff */
double dpd_r_cut_inv = 0.0;
/* DPD transversal friction coefficient gamma. */
double dpd_tgamma = 0.0;

/* NPT ISOTROPIC THERMOSTAT */
// INSERT COMMENT
double nptiso_gamma0 = 0.0;
// INSERT COMMENT
double nptiso_gammav = 0.0;

double langevin_pref1, langevin_pref2, langevin_pref2_rotation;
/** buffers for the work around for the correlated random values which cool the system,
    and require a magical heat up whenever reentering the integrator. */
static double langevin_pref2_buffer, langevin_pref2_rotation_buffer;

#ifdef DPD
double dpd_pref1;
double dpd_pref2;
double dpd_pref3;
double dpd_pref4;
static double dpd_pref2_buffer, dpd_pref4_buffer;

#endif
#ifdef NPT
double nptiso_pref1;
double nptiso_pref2;
double nptiso_pref3;
double nptiso_pref4;
#endif


int thermo_ro_callback(Tcl_Interp *interp, void *_data)
{
  Tcl_AppendResult(interp, "variable is readonly: use the thermostat command to set thermostat parameters.", (char *) NULL);
  return (TCL_ERROR);
}
int thermo_parse_off(Tcl_Interp *interp, int argc, char **argv) 
{
  /* set temperature to zero */
  temperature = 0;
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  /* langevin thermostat */
  langevin_gamma = 0;
  mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA);
#ifdef DPD
  /* dpd thermostat */  
  dpd_gamma = 0;
  mpi_bcast_parameter(FIELD_DPD_GAMMA);
  dpd_r_cut = 0;
  mpi_bcast_parameter(FIELD_DPD_RCUT);
  dpd_tgamma = 0;
  mpi_bcast_parameter(FIELD_DPD_TGAMMA);
#endif
#ifdef NPT
  /* npt isotropic thermostat */
  nptiso_gamma0 = 0;
  mpi_bcast_parameter(FIELD_NPTISO_G0);
  nptiso_gammav = 0;
  mpi_bcast_parameter(FIELD_NPTISO_GV);
#endif
  /* switch thermostat off */
  thermo_switch = THERMO_OFF;
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  return (TCL_OK);
}

int thermo_parse_langevin(Tcl_Interp *interp, int argc, char **argv) 
{
  double temp, gamma;

  /* check number of arguments */
  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
		     argv[0]," ",argv[1]," <temp> <gamma>\"", (char *)NULL);
    return (TCL_ERROR);
  }

  /* check argument types */
  if ( !ARG_IS_D(2, temp) || !ARG_IS_D(3, gamma)) {
    Tcl_AppendResult(interp, argv[0]," ",argv[1]," needs two DOUBLES", (char *)NULL);
    return (TCL_ERROR);
  }

  if (temp < 0 || gamma < 0) {
    Tcl_AppendResult(interp, "temperature and gamma must be positive", (char *)NULL);
    return (TCL_ERROR);
  }

  /* broadcast parameters */
  temperature = temp;
  langevin_gamma = gamma;
  thermo_switch = ( thermo_switch | THERMO_LANGEVIN );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA);
  return (TCL_OK);
}
 
#ifdef DPD
int thermo_parse_dpd(Tcl_Interp *interp, int argc, char **argv) 
{
  double temp, gamma, r_cut,tgamma;

#ifdef ROTATION
    fprintf(stderr,"WARNING: Do not use DPD with ROTATION compiled in\n");
    fprintf(stderr,"         You should first check if a combiantion of a DPD thermostat\n");
    fprintf(stderr,"         for the translational degrees of freedom and a LANGEVIN thermostat\n");
    fprintf(stderr,"         for the rotational ones yields correct physics!\n");
    fprintf(stderr,"         After this you may remove these lines (thermostat.c::thermo_parse_dpd)!\n");
#endif

  /* check number of arguments */
  if (argc < 5) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
		     argv[0]," ",argv[1]," <temp> <gamma> <r_cut> {<tgamma>}\"", (char *)NULL);
    return (TCL_ERROR);
  }

  /* check argument types */
  if ( !ARG_IS_D(2, temp) || !ARG_IS_D(3, gamma) || !ARG_IS_D(4, r_cut)) {
    Tcl_AppendResult(interp, argv[0]," ",argv[1]," needs at least three DOUBLES", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if (argc == 6){
    if ( !ARG_IS_D(5,tgamma) ){
      Tcl_AppendResult(interp, argv[0]," ",argv[1],": tgamma should be double",(char *)NULL);
      return (TCL_ERROR);
    }
  }
  else{
    tgamma=0.0;
  }

  /* broadcast parameters */
  temperature = temp;
  dpd_gamma = gamma;
  dpd_r_cut = r_cut;
  dpd_tgamma = tgamma;
  thermo_switch = ( thermo_switch | THERMO_DPD );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_parameter(FIELD_DPD_GAMMA);
  mpi_bcast_parameter(FIELD_DPD_RCUT);
  mpi_bcast_parameter(FIELD_DPD_TGAMMA);
  return (TCL_OK);
}
#endif

#ifdef NPT
int thermo_parse_nptiso(Tcl_Interp *interp, int argc, char **argv) 
{
  double temp, gamma0, gammav;
  /* check number of arguments */
  if (argc < 5) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
		     argv[0]," set ",argv[1]," <temp> <gamma0> <gammav>\"", (char *)NULL);
    return (TCL_ERROR);
  }
  /* check argument types */
  if ( !ARG_IS_D(2, temp) || !ARG_IS_D(3, gamma0) || !ARG_IS_D(4, gammav) ) {
    Tcl_AppendResult(interp, argv[0]," ",argv[1]," needs four DOUBLES", (char *)NULL);
    return (TCL_ERROR);
  }
  /* broadcast parameters */
  temperature = temp;
  nptiso_gamma0 = gamma0;
  nptiso_gammav = gammav;

  thermo_switch = ( thermo_switch | THERMO_NPT_ISO );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_parameter(FIELD_NPTISO_G0);
  mpi_bcast_parameter(FIELD_NPTISO_GV);
  return (TCL_OK);
}
#endif

int thermo_print(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  /* thermostat not initialized */
  if(temperature == -1.0) {
    Tcl_AppendResult(interp,"{ not initialized } ", (char *)NULL);
    return (TCL_OK);
  }

  /* no thermostat on */
  if(thermo_switch == THERMO_OFF) {
    Tcl_AppendResult(interp,"{ off } ", (char *)NULL);
    return (TCL_OK);
  }

  /* langevin */
  if(thermo_switch & THERMO_LANGEVIN ) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ langevin ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, langevin_gamma, buffer);
    Tcl_AppendResult(interp," ",buffer," } ", (char *)NULL);
  }
    
#ifdef DPD
 /* dpd */
  if(thermo_switch & THERMO_DPD) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ dpd ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, dpd_gamma, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, dpd_r_cut, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, dpd_tgamma, buffer);
    Tcl_AppendResult(interp," ",buffer, " } ", (char *)NULL);
  }
#endif

#ifdef NPT
  /* npt_isotropic */
  if(thermo_switch & THERMO_NPT_ISO) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ npt_isotropic ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nptiso_gamma0, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nptiso_gammav, buffer);
    Tcl_AppendResult(interp," ",buffer, " } ", (char *)NULL);
  }
#endif

#ifdef LB
 /* lb */
  if(thermo_switch & THERMO_LB) {
    Tcl_PrintDouble(interp, temperature, buffer);
    Tcl_AppendResult(interp,"{ lb ",buffer, " } ", (char *)NULL);
  }
#endif

  return (TCL_OK);
}

int thermo_usage(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "Usage of tcl-command thermostat:\n", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], "' for status return or \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " off' to deactivate it (=> NVE-ensemble) \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " set langevin <temp> <gamma>' or \n ", (char *)NULL);
#ifdef DPD
  Tcl_AppendResult(interp, "'", argv[0], " set dpd <temp> <gamma> <r_cut>' or \n ", (char *)NULL);
#endif
#ifdef NPT
  Tcl_AppendResult(interp, "'", argv[0], " set npt_isotropic <temp> <gamma0> <gammav>' ", (char *)NULL);
#endif
#ifdef LB
  Tcl_AppendResult(interp, "'", argv[0], " set lb <temperature>" , (char *)NULL);
#endif
  return (TCL_ERROR);
}

int thermostat(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int err = TCL_OK;
  THERMO_TRACE(fprintf(stderr,"%d: thermostat:\n",this_node));

  /* print thermostat status */
  if(argc == 1) return thermo_print(interp);
  
  if ( ARG1_IS_S("set") )          {
    argc--;
    argv++;

    if (argc == 1) {
      Tcl_AppendResult(interp, "wrong # args: \n", (char *)NULL);
      return thermo_usage(interp, argc, argv);
    }
  }
  if ( ARG1_IS_S("off") )
    err = thermo_parse_off(interp, argc, argv);
  else if ( ARG1_IS_S("langevin"))
    err = thermo_parse_langevin(interp, argc, argv);
#ifdef DPD
  else if ( ARG1_IS_S("dpd") )
    err = thermo_parse_dpd(interp, argc, argv);
#endif
#ifdef NPT
  else if ( ARG1_IS_S("npt_isotropic") )
    err = thermo_parse_nptiso(interp, argc, argv);
#endif
#ifdef LB
  else if ( ARG1_IS_S("lb") )
    err = thermo_parse_lb(interp, argc-1, argv+1);
#endif
  else {
    Tcl_AppendResult(interp, "Unknown thermostat ", argv[1], "\n", (char *)NULL);
    return thermo_usage(interp, argc, argv);
  }
  return mpi_gather_runtime_errors(interp, err);
}



void thermo_init_langevin() 
{
  langevin_pref1 = -langevin_gamma/time_step;
  langevin_pref2 = sqrt(24.0*temperature*langevin_gamma/time_step);

#ifdef ROTATION 
  langevin_gamma_rotation = langevin_gamma/3;
  langevin_pref2_rotation = sqrt(24.0*temperature*langevin_gamma_rotation/time_step);
  THERMO_TRACE(fprintf(stderr,"%d: thermo_init_langevin: langevin_gamma_rotation=%f, langevin_pref2_rotation=%f",langevin_gamma_rotation,langevin_pref2_rotation));
#endif
  THERMO_TRACE(fprintf(stderr,"%d: thermo_init_langevin: langevin_pref1=%f, langevin_pref2=%f",this_node,langevin_pref1,langevin_pref2));  
}

#ifdef DPD
void thermo_init_dpd()
{
  /* prefactor friction force */
  dpd_pref1 = dpd_gamma/time_step;
  /* prefactor random force */
  dpd_pref2 = sqrt(24.0*temperature*dpd_gamma/time_step);
  dpd_r_cut_inv = 1.0/dpd_r_cut;
  dpd_pref3 = dpd_tgamma/time_step;
  dpd_pref4 = sqrt(24.0*temperature*dpd_tgamma/time_step);
  THERMO_TRACE(fprintf(stderr,"%d: thermo_init_dpd: dpd_pref1=%f, dpd_pref2=%f, dpd_r_cut_inv=%f,dpd_pref1=%f, dpd_pref2=%f\n",
		       this_node,dpd_pref1,dpd_pref2,dpd_r_cut_inv,dpd_pref3,dpd_pref4));
}
#endif

#ifdef NPT
void thermo_init_npt_isotropic()
{
  if (nptiso.piston != 0.0) {
    nptiso_pref1 = -nptiso_gamma0*0.5 * time_step;
    nptiso_pref2 = sqrt(12.0*temperature*nptiso_gamma0*time_step) * time_step;
    nptiso_pref3 = -nptiso_gammav*(1.0/nptiso.piston)*0.5*time_step;
    nptiso_pref4 = sqrt(12.0*temperature*nptiso_gammav*time_step);
    THERMO_TRACE(fprintf(stderr,"%d: thermo_init_npt_isotropic: nptiso_pref1=%f, nptiso_pref2=%f, nptiso_pref3=%f, nptiso_pref4=%f \n",this_node,nptiso_pref1,nptiso_pref2,nptiso_pref3,nptiso_pref4));
  }
  else {
    thermo_switch = ( thermo_switch ^ THERMO_NPT_ISO );
    THERMO_TRACE(fprintf(stderr,"%d: thermo_init_npt_isotropic: switched off nptiso (piston=%f; thermo_switch=%d) \n",this_node,nptiso.piston,thermo_switch));
  }
}
#endif

void thermo_init()
{
  if(thermo_switch == THERMO_OFF)      return;
  if(thermo_switch & THERMO_LANGEVIN ) thermo_init_langevin();
#ifdef DPD
  if(thermo_switch & THERMO_DPD)       thermo_init_dpd();
#endif
#ifdef NPT
  if(thermo_switch & THERMO_NPT_ISO)   thermo_init_npt_isotropic();
#endif
}

void thermo_heat_up()
{
  if(thermo_switch & THERMO_LANGEVIN) {
    langevin_pref2_buffer          = langevin_pref2;
    langevin_pref2_rotation_buffer = langevin_pref2_rotation;
    langevin_pref2 *= sqrt(3);
    langevin_pref2_rotation *= sqrt(3);
  }
#ifdef DPD
  else if (thermo_switch & THERMO_DPD){
      dpd_pref2_buffer = dpd_pref2;
      dpd_pref2 *= sqrt(3);
      dpd_pref4_buffer = dpd_pref4;
      dpd_pref4 *= sqrt(3);
  }
#endif
}

void thermo_cool_down()
{
  if(thermo_switch & THERMO_LANGEVIN) {
    langevin_pref2          = langevin_pref2_buffer;
    langevin_pref2_rotation = langevin_pref2_rotation_buffer;
  }
#ifdef DPD
  else if (thermo_switch & THERMO_DPD){
      dpd_pref2 = dpd_pref2_buffer;
      dpd_pref4 = dpd_pref4_buffer;
  }
#endif
}

int thermo_parse_lb(Tcl_Interp *interp, int argc, char ** argv)
{
#ifdef LB
  double temp;

  /* get lb interaction type */
  if (argc < 1) {
    Tcl_AppendResult(interp, "lattice-Boltzmann needs 1 parameter: "
		     "<temperature>",
		     (char *) NULL);
    return TCL_ERROR;
  }

  /* copy lattice-boltzmann parameters */
  if (! ARG_IS_D(1, temp)) { return TCL_ERROR; }

  if ( temp < 0.0 ) {
    Tcl_AppendResult(interp, "temperature must be non-negative", (char *) NULL);
    return TCL_ERROR;
  }
  temperature = temp;
  thermo_switch = ( thermo_switch | THERMO_LB );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  
#endif
  
  return TCL_OK;
}
