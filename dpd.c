// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
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


void dpd_parse_off(Tcl_Interp *interp, int argc, char **argv)
{
  extern double dpd_gamma,dpd_r_cut;
  extern int dpd_wf;
  dpd_gamma = 0;
  mpi_bcast_parameter(FIELD_DPD_GAMMA);
  dpd_r_cut = 0;
  mpi_bcast_parameter(FIELD_DPD_RCUT);
  dpd_wf=0;
  mpi_bcast_parameter(FIELD_DPD_WF);
}


void thermo_init_dpd()
{
  extern double dpd_gamma,dpd_r_cut,dpd_pref1,dpd_pref2;
  /*extern int dpd_wf;*/
  /* prefactor friction force */
  /* NOTE: velocities are scaled with time_step, so divide by time_step here*/
  dpd_pref1 = dpd_gamma/time_step;  
  /* prefactor random force */
  /*NOTE random force is propto sqrt(time_step)*/
  dpd_pref2 = sqrt(24.0*temperature*dpd_gamma/time_step);
  dpd_r_cut_inv = 1.0/dpd_r_cut;
  THERMO_TRACE(fprintf(stderr,"%d: thermo_init_dpd: dpd_pref1=%f, dpd_pref2=%f",
		       this_node,dpd_pref1,dpd_pref2));
  THERMO_TRACE(fprintf(stderr,"\n"));
}


int thermo_parse_dpd(Tcl_Interp *interp, int argc, char **argv) 
{
  extern double dpd_gamma,dpd_r_cut;
  extern int dpd_wf;
  double temp, gamma, r_cut;
  int wf=0;

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
    Tcl_AppendResult(interp," [WF <wf>]", (char *)NULL);
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

//try for WF
  if ( (argc>0) && (ARG0_IS_S("WF")) ){
    if (!ARG1_IS_I(wf)){
      Tcl_AppendResult(interp," thermostat dpd:  wf should be int",(char *)NULL);
      return (TCL_ERROR);
    }
    else{
      argc-=2;argv+=2;
    }
  }

  if (argc > 0){
       Tcl_AppendResult(interp," thermostat dpd: too many arguments - don't know how to parse them!!!",(char *)NULL);
       return (TCL_ERROR);
  }

  /* broadcast parameters */
  temperature = temp;
  dpd_gamma = gamma;
  dpd_r_cut = r_cut;
  dpd_wf = wf;
  thermo_switch = ( thermo_switch | THERMO_DPD );
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);
  mpi_bcast_parameter(FIELD_TEMPERATURE);
  mpi_bcast_parameter(FIELD_DPD_GAMMA);
  mpi_bcast_parameter(FIELD_DPD_RCUT);
  mpi_bcast_parameter(FIELD_DPD_WF);
  return (TCL_OK);
}


void dpd_print(Tcl_Interp *interp)
{
  extern double dpd_gamma,dpd_r_cut;
  extern int dpd_wf;
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, temperature, buffer);
  Tcl_AppendResult(interp,"{ dpd ",buffer, (char *)NULL);
  Tcl_PrintDouble(interp, dpd_gamma, buffer);
  Tcl_AppendResult(interp," ",buffer, (char *)NULL);
  Tcl_PrintDouble(interp, dpd_r_cut, buffer);
  Tcl_AppendResult(interp," ",buffer, (char *)NULL);
  sprintf(buffer,"%i",dpd_wf);
  Tcl_AppendResult(interp," WF ",buffer, (char *)NULL);
  Tcl_AppendResult(interp," } ", (char *)NULL);
}


void dpd_usage(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "'", argv[0], " set dpd <temp> <gamma> <r_cut> [WF <wf>]", (char *)NULL);
  Tcl_AppendResult(interp, " ' or\n ", (char *)NULL);
}


void dpd_heat_up()
{
   extern double dpd_pref2;
   extern double dpd_pref2_buffer;
      dpd_pref2_buffer = dpd_pref2;
      dpd_pref2 *= sqrt(3);
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

  Tcl_PrintDouble(interp, data->dpd_temp, buffer);
  Tcl_AppendResult(interp, "inter_dpd ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->dpd_gamma, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
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


int interdpd_set_params(int part_type_a, int part_type_b,double temp,
				      double gamma, double r_c, int wf,
				      double tgamma, double tr_c,
				      int twf)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);

  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* inter_dpd should be symmetrically */
  data->dpd_temp   = data_sym->dpd_temp   = temp;
  data->dpd_gamma  = data_sym->dpd_gamma  = gamma;
  data->dpd_r_cut  = data_sym->dpd_r_cut  = r_c;
  data->dpd_wf     = data_sym->dpd_wf     = wf;
  data->dpd_pref1  = data_sym->dpd_pref1  = gamma/time_step;
  data->dpd_pref2  = data_sym->dpd_pref2  = sqrt(24.0*temp*gamma/time_step);
  data->dpd_tgamma = data_sym->dpd_tgamma = tgamma;
  data->dpd_tr_cut = data_sym->dpd_tr_cut = tr_c;
  data->dpd_twf    = data_sym->dpd_twf    = twf;
  data->dpd_pref3  = data_sym->dpd_pref3  = tgamma/time_step;
  data->dpd_pref4  = data_sym->dpd_pref4  = sqrt(24.0*temp*tgamma);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}


int interdpd_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for LJ */
  double gamma,r_c,tgamma,tr_c,temp;
  int wf,twf;
  int change;

  /* get inter_dpd interaction type */
  if (argc < 8) {
    Tcl_AppendResult(interp, "inter_dpd needs 7 parameters: "
		     "<temp> <gamma> <r_cut> <wf> <tgamma> <tr_cut> <twf>",
		     (char *) NULL);
    return 0;
  }

  /* copy lennard-jones parameters */
  if ((! ARG_IS_D(1, temp))  ||
      (! ARG_IS_D(2, gamma))    ||
      (! ARG_IS_D(3, r_c))    ||
      (! ARG_IS_I(4, wf))     ||
      (! ARG_IS_D(5, tgamma)) ||
      (! ARG_IS_D(6, tr_c))    ||
      (! ARG_IS_I(7, twf))    ) {
    Tcl_AppendResult(interp, "inter_dpd needs 7 parameters: "
		     "<temp> <gamma> <r_cut> <wf> <tgamma> <tr_cut> <twf> ",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 8;
	
  if (interdpd_set_params(part_type_a, part_type_b,temp,
			       gamma,r_c,wf,tgamma,tr_c,twf) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;

}

void interdpd_init(){
   int type_a,type_b;
   IA_parameters *data;

   for (type_a=0;type_a<n_particle_types;type_a++){
      for (type_b=0;type_b<n_particle_types;type_b++){
         data=get_ia_param(type_a,type_b);
         data->dpd_pref1=data->dpd_gamma/time_step;
         data->dpd_pref2=sqrt(24.0*data->dpd_temp*data->dpd_gamma/time_step);
         data->dpd_pref3=data->dpd_tgamma/time_step;
         data->dpd_pref4=sqrt(24.0*data->dpd_temp*data->dpd_tgamma/time_step);
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
         data->dpd_pref2*=pref_scale;
         data->dpd_pref4*=pref_scale;
      }
   }
}
#endif
