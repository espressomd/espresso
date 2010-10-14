/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
/** \file dpd.c
    Implementation of \ref dpd.h "dpd.h"
 */
#include "dpd.h"

/* DPD THERMOSTAT */
/* DPD longitudinal friction coefficient gamma. */
double dpd_gamma = 0.0;
/* DPD thermostat cutoff */
double dpd_r_cut = 0.0;
/* DPD weightfunction */
int dpd_wf = 0;

/* DPD transversal friction coefficient gamma. */
double dpd_tgamma = 0.0;
/* DPD thermostat trans cutoff */
double dpd_tr_cut = 0.0;
/* trans DPD weightfunction */
int dpd_twf = 0;

#ifdef DPD
/* inverse off DPD thermostat cutoff */
double dpd_r_cut_inv = 0.0;
double dpd_pref1;
double dpd_pref2;
static double dpd_pref2_buffer;

#ifdef TRANS_DPD 
/* inverse off trans DPD thermostat cutoff */
double dpd_tr_cut_inv = 0.0;
double dpd_pref3;
double dpd_pref4;
static double dpd_pref4_buffer;
#endif

void dpd_parse_off(Tcl_Interp *interp, int argc, char **argv)
{
  extern double dpd_gamma,dpd_r_cut;
  extern int dpd_wf;
#ifdef TRANS_DPD
  extern double dpd_tgamma,dpd_tr_cut;
  extern int dpd_twf;
#endif
  dpd_gamma = 0;
  mpi_bcast_parameter(FIELD_DPD_GAMMA);
  dpd_r_cut = 0;
  mpi_bcast_parameter(FIELD_DPD_RCUT);
  dpd_wf=0;
  mpi_bcast_parameter(FIELD_DPD_WF);
#ifdef TRANS_DPD
  dpd_tgamma = 0;
  mpi_bcast_parameter(FIELD_DPD_TGAMMA);
  dpd_tr_cut=0;
  mpi_bcast_parameter(FIELD_DPD_TRCUT);
  dpd_twf=0;
  mpi_bcast_parameter(FIELD_DPD_TWF);
#endif
}


void thermo_init_dpd()
{
  extern double dpd_gamma,dpd_r_cut,dpd_pref1,dpd_pref2;
  /*extern int dpd_wf;*/
#ifdef TRANS_DPD
  extern double dpd_tgamma,dpd_tr_cut,dpd_pref3,dpd_pref4;
  /*extern int dpd_twf;*/
#endif
  /* prefactor friction force */
  /* NOTE: velocities are scaled with time_step, so divide by time_step here*/
  dpd_pref1 = dpd_gamma/time_step;  
  /* prefactor random force */
  /*NOTE random force is propto sqrt(time_step)*/
  dpd_pref2 = sqrt(24.0*temperature*dpd_gamma/time_step);
  dpd_r_cut_inv = 1.0/dpd_r_cut;
#ifdef TRANS_DPD
  /* NOTE: velocities are scaled with time_step, so divide by time_step here*/
  dpd_pref3 = dpd_tgamma/time_step;
  /*NOTE random force is propto sqrt(time_step)*/
  dpd_pref4 = sqrt(24.0*temperature*dpd_tgamma/time_step);
  dpd_tr_cut_inv = 1.0/dpd_tr_cut;
#endif
  THERMO_TRACE(fprintf(stderr,"%d: thermo_init_dpd: dpd_pref1=%f, dpd_pref2=%f",
		       this_node,dpd_pref1,dpd_pref2));
#ifdef TRANS_DPD
  THERMO_TRACE(fprintf(stderr,",dpd_pref3=%f, dpd_pref4=%f\n",dpd_pref3,dpd_pref4));
#endif
  THERMO_TRACE(fprintf(stderr,"\n"));
}


int thermo_parse_dpd(Tcl_Interp *interp, int argc, char **argv) 
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
    fprintf(stderr,"         You should first check if a combiantion of a DPD thermostat\n");
    fprintf(stderr,"         for the translational degrees of freedom and a LANGEVIN thermostat\n");
    fprintf(stderr,"         for the rotational ones yields correct physics!\n");
    fprintf(stderr,"         After this you may remove these lines (thermostat.c::thermo_parse_dpd)!\n");
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


void dpd_print(Tcl_Interp *interp)
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


void dpd_usage(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "'", argv[0], " set dpd <temp> <gamma> <r_cut> [WF <wf>]", (char *)NULL);
#ifdef TRANS_DPD
  Tcl_AppendResult(interp, " [<tgamma>] [TWF <twf>]", (char *)NULL);
#endif
  Tcl_AppendResult(interp, " ' or\n ", (char *)NULL);
}


void dpd_heat_up()
{
   extern double dpd_pref2;
   extern double dpd_pref2_buffer;
#ifdef TRANS_DPD
   extern double dpd_pref4;
   extern double dpd_pref4_buffer;
#endif
      dpd_pref2_buffer = dpd_pref2;
      dpd_pref2 *= sqrt(3);
#ifdef TRANS_DPD
      dpd_pref4_buffer = dpd_pref4;
      dpd_pref4 *= sqrt(3);
#endif
}


void dpd_cool_down()
{
   extern double dpd_pref2;
   extern double dpd_pref2_buffer;
#ifdef TRNAS_DPD
   extern double dpd_pref4;
   extern double dpd_pref4_buffer;
#endif
      dpd_pref2 = dpd_pref2_buffer;
#ifdef TRANS_DPD
      dpd_pref4 = dpd_pref4_buffer;
#endif
}
#endif

#ifdef INTER_DPD
void interdpd_heat_up()
{
	double pref_scale=sqrt(3);
	interdpd_update_params(pref_scale);
}


void interdpd_cool_down()
{
	double pref_scale=1.0/sqrt(3);
	interdpd_update_params(pref_scale);
}

int printinterdpdIAToResult(Tcl_Interp *interp, int i, int j)
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


int interdpd_set_params(int part_type_a, int part_type_b,
				      double gamma, double r_c, int wf,
				      double tgamma, double tr_c,
				      int twf)
{
  extern double temperature;
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);

  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* inter_dpd should be symmetrically */
  data->dpd_gamma  = data_sym->dpd_gamma  = gamma;
  data->dpd_r_cut  = data_sym->dpd_r_cut  = r_c;
  data->dpd_wf     = data_sym->dpd_wf     = wf;
  data->dpd_pref1  = data_sym->dpd_pref1  = gamma/time_step;
  data->dpd_pref2  = data_sym->dpd_pref2  = sqrt(24.0*temperature*gamma/time_step);
  data->dpd_tgamma = data_sym->dpd_tgamma = tgamma;
  data->dpd_tr_cut = data_sym->dpd_tr_cut = tr_c;
  data->dpd_twf    = data_sym->dpd_twf    = twf;
  data->dpd_pref3  = data_sym->dpd_pref3  = tgamma/time_step;
  data->dpd_pref4  = data_sym->dpd_pref4  = sqrt(24.0*temperature*tgamma/time_step);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

int thermo_parse_interdpd(Tcl_Interp *interp, int argc, char ** argv)
{
  double temp;

  if (argc < 2) {
    Tcl_AppendResult(interp, "thermostat needs 1 parameter: "
		     "<temperature>",
		     (char *) NULL);
    return TCL_ERROR;
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

int interdpd_parser(Tcl_Interp * interp,
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
    return TCL_ERROR;
  }
  change = 7;
	
  if (interdpd_set_params(part_type_a, part_type_b,
			       gamma,r_c,wf,tgamma,tr_c,twf) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  interdpd_init();

  return change;

}

void interdpd_init(){
   extern double temperature;
   int type_a,type_b;
   IA_parameters *data;

   for (type_a=0;type_a<n_particle_types;type_a++){
      for (type_b=0;type_b<n_particle_types;type_b++){
         data=get_ia_param(type_a,type_b);
         if ( (data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0) ) {
            data->dpd_pref1=data->dpd_gamma/time_step;
            data->dpd_pref2=sqrt(24.0*temperature*data->dpd_gamma/time_step);
            data->dpd_pref3=data->dpd_tgamma/time_step;
            data->dpd_pref4=sqrt(24.0*temperature*data->dpd_tgamma/time_step);
         }
      }
   }
}

void interdpd_parse_off(){
   int type_a,type_b;
   IA_parameters *data;
   for (type_a=0;type_a<n_particle_types;type_a++){
      for (type_b=0;type_b<n_particle_types;type_b++){
         data=get_ia_param(type_a,type_b);
         data->dpd_gamma  = data->dpd_r_cut  = data->dpd_wf =
         data->dpd_pref1  = data->dpd_pref2  = data->dpd_tgamma =
         data->dpd_tr_cut = data->dpd_twf    = data->dpd_pref3  =
         data->dpd_pref4  = 0.0;
      }
   }
}

void interdpd_update_params(double pref_scale)
{
   int type_a,type_b;
   IA_parameters *data;

   for (type_a=0;type_a<n_particle_types;type_a++){
      for (type_b=0;type_b<n_particle_types;type_b++){
         data=get_ia_param(type_a,type_b);
         if ( (data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0) ) {
            data->dpd_pref2*=pref_scale;
            data->dpd_pref4*=pref_scale;
         }
      }
   }
}
#endif
