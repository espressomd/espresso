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
/** \file dpd.cpp
    Implementation of \ref dpd_tcl.hpp "dpd_tcl.h"
 */
#include "dpd_tcl.hpp"

#include "utils.hpp"
#include "parser.hpp"
#include "thermostat.hpp"
#include "interaction_data.hpp"
#include "virtual_sites.hpp"

int tclcommand_thermostat_parse_dpd(Tcl_Interp *interp, int argc, char **argv) 
{
  extern double dpd_gamma,dpd_r_cut;
  extern int dpd_wf;
#ifdef TRANS_DPD
  extern double dpd_tgamma,dpd_tr_cut;
  extern int dpd_twf;
#endif
  double temp, gamma, r_cut;
  int wf=0;
#ifdef TRANS_DPD
  double tgamma=0.0,tr_cut;
  int twf;
  int set_tgamma=0;
#endif

#ifdef ROTATION
    fprintf(stderr,"WARNING: Do not use DPD with ROTATION compiled in\n");
    fprintf(stderr,"         You should first check if a combination of a DPD thermostat\n");
    fprintf(stderr,"         for the translational degrees of freedom and a LANGEVIN thermostat\n");
    fprintf(stderr,"         for the rotational ones yields correct physics!\n");
    fprintf(stderr,"         After this you may remove these lines (thermostat.c::tclcommand_thermostat_parse_dpd)!\n");
#endif

  /* check number of arguments */
  if (argc < 5) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
		     argv[0]," ",argv[1]," <temp> <gamma> <r_cut>", (char *)NULL);
#ifdef TRANS_DPD
    Tcl_AppendResult(interp,"[<tgamma>] [<tR_cut>]", (char *)NULL);
#endif
    Tcl_AppendResult(interp," [WF <wf>]", (char *)NULL);
#ifdef TRANS_DPD
    Tcl_AppendResult(interp," [TWF <twf>]", (char *)NULL);
#endif
    Tcl_AppendResult(interp,"\"", (char *)NULL);
    return (TCL_ERROR);
  }

  /* check argument types */
  if ( !ARG_IS_D(2, temp) || !ARG_IS_D(3, gamma) || !ARG_IS_D(4, r_cut)) {
    Tcl_AppendResult(interp, argv[0]," ",argv[1]," needs at least three DOUBLES", (char *)NULL);
    return (TCL_ERROR);
  }
  argc-=5;
  argv+=5;

#ifdef TRANS_DPD
  tgamma=0;
  tr_cut=r_cut;
  twf=wf;
  if ( (argc>0) && (!ARG0_IS_S("WF")) ) {
     if (!ARG0_IS_D(tgamma)) {
        Tcl_AppendResult(interp," thermostat dpd:  tgamma should be double",(char *)NULL);
        return (TCL_ERROR);
     }
     else{
        argc--;
        argv++;
        set_tgamma++;
     }
  }
#endif
//try for WF
  if ( (argc>0) && (ARG0_IS_S("WF")) ){
    if (!ARG1_IS_I(wf)){
      Tcl_AppendResult(interp," thermostat dpd:  wf should be int",(char *)NULL);
      return (TCL_ERROR);
    }
    else{
      argc-=2;argv+=2;
#ifdef TRANS_DPD
      twf=wf;
#endif
    }
  }
#ifdef TRANS_DPD
  if ( (set_tgamma==0) && (argc>0) && (!ARG0_IS_S("TWF")) ) {
     if (!ARG0_IS_D(tgamma)) {
        Tcl_AppendResult(interp," thermostat dpd:  tgamma should be double",(char *)NULL);
        return (TCL_ERROR);
     }
     else{
        argc--;
        argv++;
        set_tgamma++;
     }
  }

  if ( (argc>0) && (!ARG0_IS_S("TWF")) ) {
    if (set_tgamma!=0) {
      if (!ARG0_IS_D(tr_cut)) {
        Tcl_AppendResult(interp," thermostat dpd:  tr_cut should be double",(char *)NULL);
        return (TCL_ERROR);
      }
      else{
        argc--;
        argv++;
      }
    }
    else{
       Tcl_AppendResult(interp," thermostat dpd: tgamma must be set before twf",(char *)NULL);
       return (TCL_ERROR);
    }
  }

  if ( (argc>0) && (ARG0_IS_S("TWF")) ) {
     if (set_tgamma!=0) {
       if (!ARG1_IS_I(wf)) {
          Tcl_AppendResult(interp," thermostat dpd: twf should be int",(char *)NULL);
          return (TCL_ERROR);
       }
       else{
          argc-=2;argv+=2;
       }
     }
     else{
       Tcl_AppendResult(interp," thermostat dpd: tgamma must be set before twf",(char *)NULL);
       return (TCL_ERROR);
     }
  }
#endif

  if (argc > 0){
       Tcl_AppendResult(interp," thermostat dpd: too many arguments - don't know how to parse them!!!",(char *)NULL);
       return (TCL_ERROR);
  }

  /* broadcast parameters */
  temperature = temp;
  dpd_gamma = gamma;
  dpd_r_cut = r_cut;
  dpd_wf = wf;
#ifdef TRANS_DPD
  dpd_tgamma = tgamma;
  dpd_tr_cut= tr_cut;
  dpd_twf=twf;
#endif
  thermo_switch = ( thermo_switch | THERMO_DPD );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_parameter(FIELD_DPD_GAMMA);
  mpi_bcast_parameter(FIELD_DPD_RCUT);
  mpi_bcast_parameter(FIELD_DPD_WF);
#ifdef TRANS_DPD
  mpi_bcast_parameter(FIELD_DPD_TGAMMA);
  mpi_bcast_parameter(FIELD_DPD_TRCUT);
  mpi_bcast_parameter(FIELD_DPD_TWF);
#endif
  return (TCL_OK);
}


void tclcommand_thermostat_parse_and_print_dpd(Tcl_Interp *interp)
{
  extern double dpd_gamma,dpd_r_cut;
  extern int dpd_wf;
#ifdef TRANS_DPD
  extern double dpd_tgamma,dpd_tr_cut;
  extern int dpd_twf;
#endif
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, temperature, buffer);
  Tcl_AppendResult(interp,"{ dpd ",buffer, (char *)NULL);
  Tcl_PrintDouble(interp, dpd_gamma, buffer);
  Tcl_AppendResult(interp," ",buffer, (char *)NULL);
  Tcl_PrintDouble(interp, dpd_r_cut, buffer);
  Tcl_AppendResult(interp," ",buffer, (char *)NULL);
  sprintf(buffer,"%i",dpd_wf);
  Tcl_AppendResult(interp," WF ",buffer, (char *)NULL);
#ifdef TRANS_DPD
  Tcl_PrintDouble(interp, dpd_tgamma, buffer);
  Tcl_AppendResult(interp," ",buffer, (char *)NULL);
  Tcl_PrintDouble(interp, dpd_tr_cut, buffer);
  Tcl_AppendResult(interp," ",buffer, (char *)NULL);
  sprintf(buffer,"%i",dpd_twf);
  Tcl_AppendResult(interp," TWF ",buffer, (char *)NULL);
#endif
  Tcl_AppendResult(interp," } ", (char *)NULL);
}


void tclcommand_thermostat_print_usage_dpd(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "'", argv[0], " set dpd <temp> <gamma> <r_cut> [WF <wf>]", (char *)NULL);
#ifdef TRANS_DPD
  Tcl_AppendResult(interp, " [<tgamma>] [TWF <twf>]", (char *)NULL);
#endif
  Tcl_AppendResult(interp, " ' or\n ", (char *)NULL);
}


#ifdef INTER_DPD
int tclprint_to_result_inter_dpdIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->dpd_gamma, buffer);
  Tcl_AppendResult(interp, "inter_dpd ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->dpd_r_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%i", data->dpd_wf);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
/*  Tcl_PrintDouble(interp, data->dpd_pref2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  */
  Tcl_PrintDouble(interp, data->dpd_tgamma, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->dpd_tr_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%i", data->dpd_twf);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
//   Tcl_PrintDouble(interp, data->dpd_pref4, buffer);
//   Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  

  return TCL_OK;
}


int tclcommand_thermostat_parse_inter_dpd(Tcl_Interp *interp, int argc, char ** argv)
{
  double temp;

  if (argc < 2) {
    Tcl_AppendResult(interp, "thermostat needs 1 parameter: "
		     "<temperature>",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if (argc>2 && ARG_IS_S(2, "ignore_fixed_particles")) {
    if (argc == 3)
      dpd_ignore_fixed_particles=1;
    else if (argc!= 4 || (!ARG_IS_I(3, dpd_ignore_fixed_particles))) 
      return TCL_ERROR;
    mpi_bcast_parameter(FIELD_DPD_IGNORE_FIXED_PARTICLES);
    return TCL_OK;
  }

  /* copy lattice-boltzmann parameters */
  if (! ARG_IS_D(2, temp)) { return TCL_ERROR; }

  if ( temp < 0.0 ) {
    Tcl_AppendResult(interp, "temperature must be non-negative", (char *) NULL);
    return TCL_ERROR;
  }
  temperature = temp;
  thermo_switch = ( thermo_switch | THERMO_INTER_DPD );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  return (TCL_OK);
}

int tclcommand_inter_parse_inter_dpd(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for LJ */
  extern double temperature;
  double gamma,r_c,tgamma,tr_c;
  int wf,twf;
  int change;

  /* get inter_dpd interaction type */
  if (argc < 7) {
    Tcl_AppendResult(interp, "inter_dpd needs 6 parameters: "
		     "<gamma> <r_cut> <wf> <tgamma> <tr_cut> <twf>",
		     (char *) NULL);
    return 0;
  }
  if (temperature == -1) {
    Tcl_AppendResult(interp, "Please set temperature first:  temperature inter_dpd temp",(char *) NULL);
    return 0;
  }

  /* copy lennard-jones parameters */
  if ((! ARG_IS_D(1, gamma))  ||
      (! ARG_IS_D(2, r_c))    ||
      (! ARG_IS_I(3, wf))     ||
      (! ARG_IS_D(4, tgamma)) ||
      (! ARG_IS_D(5, tr_c))    ||
      (! ARG_IS_I(6, twf))    ) {
    Tcl_AppendResult(interp, "inter_dpd needs 6 parameters: "
		     "<gamma> <r_cut> <wf> <tgamma> <tr_cut> <twf> ",
		     (char *) NULL);
    return 0;
  }
  change = 7;
	
  if (inter_dpd_set_params(part_type_a, part_type_b,
			       gamma,r_c,wf,tgamma,tr_c,twf) == ES_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  inter_dpd_init();

  return change;

}
#endif
