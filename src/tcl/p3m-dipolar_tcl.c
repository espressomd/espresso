/*
  Copyright (C) 2010,2011 The ESPResSo project
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
/** \file p3m-dipolar.c  P3M algorithm for long range magnetic dipole-dipole interaction.
 *
 NB: In general the magnetic dipole-dipole functions bear the same
     name than the charge-charge but, adding in front of the name a D
     and replacing where "charge" appears by "dipole". In this way one
     can recognize the similarity of the functions but avoiding nasty
     confusions in their use.

 PS: By default the magnetic epsilon is metallic = 0.  
*/

#include "p3m-dipolar_tcl.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "integrate.h"
#include "global.h"
#include "grid.h"
#include "domain_decomposition.h"
#include "particle_data.h"
#include "communication.h"
#include "thermostat.h"
#include "cells.h"
#include "tuning.h"
#include "p3m-dipolar.h"
#include "fft-dipolar.h"

#ifdef DP3M

/************************************************
 * DEFINES
 ************************************************/

/* MPI tags for the charge-charge p3m communications: */
///** Tag for communication in P3M_init() -> send_calc_mesh(). */
//#define REQ_P3M_INIT_D   2001
///** Tag for communication in p3m_gather_fft_grid(). */
//#define REQ_P3M_GATHER_D 2011
///** Tag for communication in p3m_spread_force_grid(). */
//#define REQ_P3M_SPREAD_D 2021
//

extern dp3m_data_struct dp3m;


/** \name Private Functions */
/************************************************************/
/*@{*/


static int tclcommand_inter_magnetic_dp3m_print_adaptive_tune_parameters(Tcl_Interp *interp);
static int tclcommand_inter_magnetic_dp3m_print_tune_parameters(Tcl_Interp *interp);

/*@}*/


/* Compute the dipolar surface terms */

/** \name P3M Tuning Functions (private)*/
/************************************************************/
/*@{*/



/******************  functions related to the parsing&tuning  of the dipolar parameters **********/
			 

int tclcommand_inter_magnetic_parse_dp3m_tune_params(Tcl_Interp * interp, int argc, char ** argv, int adaptive)
{
  int mesh = -1, cao = -1, n_interpol = -1;
  double r_cut = -1, accuracy = -1;

  while(argc > 0) {
    if(ARG0_IS_S("r_cut")) {
      if (! (argc > 1 && ARG1_IS_D(r_cut) && r_cut >= -1)) {
	Tcl_AppendResult(interp, "r_cut expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("mesh")) {
      if(! (argc > 1 && ARG1_IS_I(mesh) && mesh >= -1)) {
	Tcl_AppendResult(interp, "mesh expects an integer >= -1",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("cao")) {
      if(! (argc > 1 && ARG1_IS_I(cao) && cao >= -1 && cao < 7)) {
	Tcl_AppendResult(interp, "cao expects an integer between -1 and 7",
			 (char *) NULL);
	return TCL_ERROR;
      } 

    } else if(ARG0_IS_S("accuracy")) {
      if(! (argc > 1 && ARG1_IS_D(accuracy) && accuracy > 0)) {
	Tcl_AppendResult(interp, "accuracy expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }

    } else if (ARG0_IS_S("n_interpol")) {
      if (! (argc > 1 && ARG1_IS_I(n_interpol) && n_interpol >= 0)) {
	Tcl_AppendResult(interp, "n_interpol expects an nonnegative integer",
			 (char *) NULL);
	return TCL_ERROR;
      }
    }
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }
  
  dp3m_set_tune_params(r_cut, mesh, cao, -1.0, accuracy, n_interpol);

  /* check for optional parameters */
  if (argc > 0) {
    if (tclcommand_inter_magnetic_parse_dp3m_opt_params(interp, argc, argv) == TCL_ERROR)
      return TCL_ERROR;
  }

  if (adaptive) {
    if(tclcommand_inter_magnetic_dp3m_print_adaptive_tune_parameters(interp) == TCL_ERROR) 
      return TCL_ERROR;
  }
  else {
    if(tclcommand_inter_magnetic_dp3m_print_tune_parameters(interp) == TCL_ERROR) 
      return TCL_ERROR;
  }

  return TCL_OK;
}



/*****************************************************************************/



int tclcommand_inter_magnetic_parse_dp3m(Tcl_Interp * interp, int argc, char ** argv)
{
  double r_cut, alpha, accuracy = -1.0;
  int mesh, cao, i;

  if (coulomb.Dmethod != DIPOLAR_P3M && coulomb.Dmethod != DIPOLAR_MDLC_P3M)
    coulomb.Dmethod = DIPOLAR_P3M;
    
#ifdef PARTIAL_PERIODIC
  if(PERIODIC(0) == 0 ||
     PERIODIC(1) == 0 ||
     PERIODIC(2) == 0)
    {
      Tcl_AppendResult(interp, "Need periodicity (1,1,1) with dipolar P3M",
		       (char *) NULL);
      return TCL_ERROR;  
    }
#endif

  if (argc < 1) {
    Tcl_AppendResult(interp, "expected: inter dipolar <bjerrum> p3m tune | <r_cut> <mesh> <cao> [<alpha> [<accuracy>]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    Tcl_AppendResult(interp, "Node grid not suited for dipolar P3M. Node grid must be sorted, largest first.", (char *) NULL);
    return TCL_ERROR;  
  }

  if (ARG0_IS_S("tune"))
    return tclcommand_inter_magnetic_parse_dp3m_tune_params(interp, argc-1, argv+1, 0);

  if (ARG0_IS_S("tunev2"))
    return tclcommand_inter_magnetic_parse_dp3m_tune_params(interp, argc-1, argv+1, 1);
      
  if(! ARG0_IS_D(r_cut))
    return TCL_ERROR;  

  if(argc < 3 || argc > 5) {
    Tcl_AppendResult(interp, "wrong # arguments: inter dipolar <bjerrum> p3m <r_cut> <mesh> <cao> [<alpha> [<accuracy>]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if((! ARG_IS_I(1, mesh)) || (! ARG_IS_I(2, cao))) {
    Tcl_AppendResult(interp, "integer expected", (char *) NULL);
    return TCL_ERROR;
  }
	
  if(argc > 3) {
    if(! ARG_IS_D(3, alpha))
      return TCL_ERROR;
  }
  else {
    Tcl_AppendResult(interp, "Automatic p3m tuning not implemented.",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(argc > 4) {
    if(! ARG_IS_D(4, accuracy)) {
      Tcl_AppendResult(interp, "double expected", (char *) NULL);
      return TCL_ERROR;
    }
  }

  if ((i = dp3m_set_params(r_cut, mesh, cao, alpha, accuracy)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "r_cut must be positive", (char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "mesh must be positive", (char *) NULL);
      break;
    case -3:
      Tcl_AppendResult(interp, "cao must be between 1 and 7 and less than mesh",
		       (char *) NULL);
      break;
    case -4:
      Tcl_AppendResult(interp, "alpha must be positive", (char *) NULL);
      break;
    case -5:
      Tcl_AppendResult(interp, "accuracy must be positive", (char *) NULL);
      break;
    default:;
      Tcl_AppendResult(interp, "unspecified error", (char *) NULL);
    }

    return TCL_ERROR;

  }

  return TCL_OK;
}


/*****************************************************************************/



int tclcommand_inter_magnetic_parse_dp3m_opt_params(Tcl_Interp * interp, int argc, char ** argv)
{
  int i; double d1, d2, d3;

  Tcl_ResetResult(interp);

  while (argc > 0) {
    /* p3m parameter: inter */
    if (ARG0_IS_S("n_interpol")) {
      
      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG1_IS_I(i)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 INTEGER parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (dp3m_set_ninterpol(i) == TCL_ERROR) {
	Tcl_AppendResult(interp, argv[0], " argument must be positive",
			 (char *) NULL);
	return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;
    }
    
    /* p3m parameter: mesh_off */
    else if (ARG0_IS_S("mesh_off")) {
      
      if(argc < 4) {
	Tcl_AppendResult(interp, argv[0], " needs 3 parameters",
			 (char *) NULL);
	return TCL_ERROR;
      }
	
      if ((! ARG_IS_D(1, d1)) ||
	  (! ARG_IS_D(2, d2)) ||
	  (! ARG_IS_D(3, d3)))
	{
	  Tcl_AppendResult(interp, argv[0], " needs 3 DOUBLE parameters",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      if (dp3m_set_mesh_offset(d1, d2 ,d3) == TCL_ERROR)
	{
	  Tcl_AppendResult(interp, argv[0], " parameters have to be between 0.0 an 1.0",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      argc -= 4;
      argv += 4;
    }
    
    /* p3m parameter: epsilon */
    else if(ARG0_IS_S( "epsilon")) {

      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }

      if (ARG1_IS_S("metallic")) {
	d1 = P3M_EPSILON_METALLIC;
      }
      else if (! ARG1_IS_D(d1)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 DOUBLE parameter or \"metallic\"",
	                 (char *) NULL);
	return TCL_ERROR;
      }
	
      if (dp3m_set_eps(d1) == TCL_ERROR) {
        Tcl_AppendResult(interp, argv[0], " There is no error msg yet!",
                         (char *) NULL);
        return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;	    
    }
    else {
      Tcl_AppendResult(interp, "Unknown coulomb p3m parameter: \"",argv[0],"\"",(char *) NULL);
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}

int tclprint_to_result_dp3m(Tcl_Interp *interp)
{
#ifdef MAGNETOSTATICS
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, dp3m.params.r_cut, buffer);
  Tcl_AppendResult(interp, "p3m ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",dp3m.params.mesh[0]);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",dp3m.params.cao);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, dp3m.params.alpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, dp3m.params.accuracy, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  Tcl_AppendResult(interp, "} {magnetic epsilon ", (char *) NULL);
  if (dp3m.params.epsilon == P3M_EPSILON_METALLIC)
    Tcl_AppendResult(interp, "metallic ", (char *) NULL);
  else {
    Tcl_PrintDouble(interp, dp3m.params.epsilon, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  }
  sprintf(buffer,"%d",dp3m.params.inter);
  Tcl_AppendResult(interp, "n_interpol ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, dp3m.params.mesh_off[0], buffer);
  Tcl_AppendResult(interp, "mesh_off ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, dp3m.params.mesh_off[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, dp3m.params.mesh_off[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
#endif 

  return TCL_OK;
}




/*****************************************************************************/


/************************************************
 * Functions for dipoloar P3M Parameter tuning
 * This tuning is based on the P3M tunning of the charges
 which in turn is based on the P3M_tune by M. Deserno
 ************************************************/

#define P3M_TUNE_MAX_CUTS 50
/** Tune dipolar P3M parameters to desired accuracy.

    Usage:
    \verbatim inter dipolar <bjerrum> p3m tune accuracy <value> [r_cut <value> mesh <value> cao <value>] \endverbatim

    The parameters are tuned to obtain the desired accuracy in best
    time, by running mpi_integrate(0) for several parameter sets.

    The function utilizes the analytic expression of the error estimate 
    for the dipolar P3M method see JCP,2008 paper by J.J.Cerda et al in 
    order to obtain the rms error in the force for a system of N randomly 
    distributed particles in a cubic box.
    For the real space error the estimate of Kolafa/Perram is used. 

    Parameter range if not given explicit values: For \ref dp3m_struct::r_cut_iL
    the function uses the values (\ref min_local_box_l -\ref #skin) /
    (n * \ref box_l), n being an integer (this implies the assumption that \ref
    dp3m_struct::r_cut_iL is the largest cutoff in the system!). For \ref
    dp3m_struct::mesh the function uses the two values which matches best the
    equation: number of mesh point = number of magnetic dipolar particles. For
    \ref dp3m_struct::cao the function considers all possible values.

    For each setting \ref dp3m_struct::alpha_L is calculated assuming that the
    error contributions of real and reciprocal space should be equal.

    After checking if the total error fulfils the accuracy goal the
    time needed for one force calculation (including verlet list
    update) is measured via \ref mpi_integrate (0).

    The function returns a log of the performed tuning.

    The function is based on routines for charges.
 */


int tclcommand_inter_magnetic_dp3m_print_tune_parameters(Tcl_Interp *interp)
{
  int i,ind, try=0, best_try=0, n_cuts;
  double r_cut_iL, r_cut_iL_min  , r_cut_iL_max, r_cut_iL_best=0, cuts[P3M_TUNE_MAX_CUTS], cut_start;
  int    mesh    , mesh_min      , mesh_max    , mesh_best=0;
  int    cao     , cao_min       , cao_max     , cao_best=0;
  double alpha_L , alpha_L_best=0, accuracy    , accuracy_best=0;
  double mesh_size, k_cut;
  double rs_err, rs_err_best=0, ks_err, ks_err_best=0;
  double int_time=0, min_time=1e20, int_num;
  char b1[TCL_DOUBLE_SPACE + 12],b2[TCL_DOUBLE_SPACE + 12],b3[TCL_DOUBLE_SPACE + 12];
 
  P3M_TRACE(fprintf(stderr,"%d: tclcommand_inter_magnetic_dp3m_print_tune_parameters\n",this_node));
  
  /* preparation */
  mpi_bcast_event(P3M_COUNT_DIPOLES);

  /* calculate r_cut_iL tune range */
  if(dp3m.params.r_cut_iL == 0.0) { 
    n_cuts = P3M_TUNE_MAX_CUTS;
    for(i=0;i<n_cuts;i++) {
      if(min_local_box_l == min_box_l)
	cuts[i] = min_local_box_l/(i+2.0)-(skin);
      else 
	cuts[i] = min_local_box_l/(i+1.0)-(skin);
      cuts[i]*=box_l_i[0];
      if( cuts[i] <= 0.0 ) {
	n_cuts = i; 
	break;
      } 
    }
    r_cut_iL_max = cuts[0];
    r_cut_iL_min = cuts[n_cuts-1];
  }
  else { 
    n_cuts = 1;
    r_cut_iL_min = r_cut_iL_max = dp3m.params.r_cut_iL; 
    cuts[0] = dp3m.params.r_cut_iL;
  }
  /* calculate mesh tune range */
  if(dp3m.params.mesh[0] == 0 ) {
    double expo;
     expo = log(pow((double)dp3m.sum_dip_part,(1.0/3.0)))/log(2.0);    
    mesh_min = (int)(pow(2.0,(double)((int)expo))+0.1);
    mesh_max = mesh_min*4;
    if(mesh_min < 8) { mesh_min = 8; mesh_max = 16; }
  }
  else { mesh_min = mesh_max = dp3m.params.mesh[0]; }
  /* calculate cao tune range */
  if(dp3m.params.cao == 0) { cao_min = 1; cao_max = 7; }
  else             { cao_min = cao_max = dp3m.params.cao; }

  /* Print Status */
  sprintf(b1,"%.5e",dp3m.params.accuracy);
  Tcl_AppendResult(interp, "P3M tune parameters: Accuracy goal = ",b1,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, box_l[0], b1);

      sprintf(b2,"%d",dp3m.sum_dip_part);   

  Tcl_PrintDouble(interp, dp3m.sum_mu2, b3);
  Tcl_AppendResult(interp, "System: box_l = ",b1,", # charged part = ",b2," Sum[q_i^2] = ",b3,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, r_cut_iL_min, b1);  Tcl_PrintDouble(interp, r_cut_iL_max, b2);
  Tcl_AppendResult(interp, "Range for dp3m.params.r_cut_iL: [",b1,"-",b2,"]","\n", (char *) NULL);
  sprintf(b1,"%d",mesh_min);  sprintf(b2,"%d",mesh_max);
  Tcl_AppendResult(interp, "Range for dp3m.params.mesh:     [",b1,"-",b2,"]","\n", (char *) NULL);
  sprintf(b1,"%d",cao_min);  sprintf(b2,"%d",cao_max);
  Tcl_AppendResult(interp, "Range for dp3m.params.cao:      [",b1,"-",b2,"]","\n\n", (char *) NULL);
  Tcl_AppendResult(interp, "set mesh cao r_cut_iL     alpha_L      err          ks_err     rs_err     time [ms]\n", (char *) NULL);

  /* Tuning Loops */
  for(mesh = mesh_min; mesh <= mesh_max; mesh*=2) { /* mesh loop */
    cut_start = box_l[0] * box_l_i[0];

       if(mesh <= 32 || dp3m.sum_dip_part > 2000) int_num=5; else int_num=1;  

    for(cao = cao_min; cao <= cao_max; cao++) {     /* cao loop */
      mesh_size = box_l[0]/(double)mesh;
      k_cut =  mesh_size*cao/2.0;
      if(cao < mesh && k_cut < dmin(min_box_l,min_local_box_l)-skin) {
	ind=0;
	for(i=0;i< n_cuts -1;i++) {
	  if(cut_start <= cuts[i]) ind=i+1;
	}
	while (ind < n_cuts) {           /* r_cut_iL loop */
	  r_cut_iL = cuts[ind];
	  /* calc maximal real space error for setting */
	  //Beginning of JJCP modification 15/5/2006 
	  
	     if(r_cut_iL *box_l[0] < 1.0) break;   // It has little sense checking so little rcuts ...  
	  
	         
		 
	     //Alpha cannot be zero in the dipolar case because real_space formula breaks down	 
	     
	     // Here we follow a method different from Coulombic case because the formula for the real space error
	     // is a trascendental equation for alpha in difference to the coulomb case 
	     
	    
	     
	    rs_err=P3M_DIPOLAR_real_space_error(box_l[0],coulomb.Dprefactor,r_cut_iL,dp3m.sum_dip_part,dp3m.sum_mu2,0.001);
	     

	  if(sqrt(2.0)*rs_err > dp3m.params.accuracy) {
	    /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
	    alpha_L = sqrt(log(sqrt(2.0)*rs_err/dp3m.params.accuracy)) / r_cut_iL;
	    /* calculate real space and k space error for this alpha_L */

	    alpha_L=dp3m_rtbisection(box_l[0],coulomb.Dprefactor,r_cut_iL,dp3m.sum_dip_part,dp3m.sum_mu2,
                    0.0001*box_l[0],5.0*box_l[0],0.0001,dp3m.params.accuracy);
	    
	    
	    rs_err = P3M_DIPOLAR_real_space_error(box_l[0],coulomb.Dprefactor,r_cut_iL,dp3m.sum_dip_part,dp3m.sum_mu2,alpha_L);
	    ks_err = dp3m_k_space_error(box_l[0],coulomb.Dprefactor,mesh,cao,dp3m.sum_dip_part,dp3m.sum_mu2,alpha_L);

	     accuracy = sqrt(SQR(rs_err)+SQR(ks_err));
	    
	    /* check if this matches the accuracy goal */
	    if(accuracy <= dp3m.params.accuracy) {
	      cut_start = cuts[ind];
	      /* broadcast p3m parameters for test run */
	      dp3m.params.r_cut_iL = r_cut_iL;
	      dp3m.params.mesh[0]  = dp3m.params.mesh[1] = dp3m.params.mesh[2] = mesh;
	      dp3m.params.cao      = cao;
	      dp3m.params.alpha_L  = alpha_L;
	      dp3m_scaleby_box_l();
	      /* initialize p3m structures */
	      mpi_bcast_coulomb_params();
	      /* perform force calculation test */
	      int_time = time_force_calc(int_num);
	      if (int_time == -1)
		return TCL_ERROR;
	      try++;
	      P3M_TRACE(fprintf(stderr,"%d ",try));
	      /* print result */
	      sprintf(b1,"%-3d",try); sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
	      Tcl_AppendResult(interp, b1," ", b2," ", b3," ", (char *) NULL);
	      sprintf(b1,"%.5e",r_cut_iL); sprintf(b2,"%.5e",alpha_L); sprintf(b3,"%.5e",accuracy);
	      Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);
	      sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err); sprintf(b3,"%-8d",(int)int_time);
	      Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"\n", (char *) NULL);
	      if(int_time <= min_time  && r_cut_iL > 0) {
		min_time      = int_time;
		r_cut_iL_best = r_cut_iL;
		mesh_best     = mesh;
		cao_best      = cao;
		alpha_L_best  = alpha_L;
		accuracy_best = sqrt(SQR(rs_err)+SQR(ks_err));
		rs_err_best   = rs_err;
		ks_err_best   = ks_err;
		best_try      = try;
	      }
	    }
	  }
	  ind++;
	}
      }
    }
  }
  P3M_TRACE(fprintf(stderr,"\n"));
  if(try==0) {
    Tcl_AppendResult(interp, "\nFailed to tune P3M parameters to required accuracy ! \n", (char *) NULL);
    return (TCL_ERROR);
  }

  /* set tuned p3m parameters */
  dp3m.params.r_cut_iL = r_cut_iL_best;
  dp3m.params.mesh[0]  = dp3m.params.mesh[1] = dp3m.params.mesh[2] = mesh_best;
  dp3m.params.cao      = cao_best;
  dp3m.params.alpha_L  = alpha_L_best;
  dp3m.params.accuracy = accuracy_best;
  dp3m_scaleby_box_l();
  /* broadcast tuned p3m parameters */
  mpi_bcast_coulomb_params();
  /* Tell the user about the outcome */
  sprintf(b1,"%d",try);
  Tcl_AppendResult(interp, "\nTune results of ",b1," trials:\n", (char *) NULL);
  sprintf(b1,"%-3d",best_try); sprintf(b2,"%-4d",mesh_best); sprintf(b3,"%-3d",cao_best);
  Tcl_AppendResult(interp, b1," ", b2," ", b3," ", (char *) NULL);
  sprintf(b1,"%.5e",r_cut_iL_best); sprintf(b2,"%.5e",alpha_L_best); sprintf(b3,"%.5e",accuracy_best);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);
  sprintf(b1,"%.3e",rs_err_best); sprintf(b2,"%.3e",ks_err_best); sprintf(b3,"%-8d",(int)min_time);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"\n", (char *) NULL);
  sprintf(b1,"%g",coulomb.Dbjerrum); sprintf(b2,"%g",dp3m.params.r_cut); sprintf(b3,"%d",mesh_best); 
  Tcl_AppendResult(interp, "=> inter coulomb ", b1, " p3m ", b2, " ", b3, (char *) NULL);
  sprintf(b1,"%d",cao_best); sprintf(b2,"%g",dp3m.params.alpha); sprintf(b3,"%g",accuracy_best);
  Tcl_AppendResult(interp, " ", b1," ", b2," ", b3," \n", (char *) NULL);

  return (TCL_OK);
}

/*****************************************************************************/


/*****************************************************************************/
/** get the optimal alpha and the corresponding computation time for fixed mesh, cao. The r_cut is determined via
    a simple bisection. Returns -1 if the force evaluation does not work, -2 if there is no valid r_cut, and -3 if
    the charge assigment order is to large for this grid */
static double tclcommand_inter_magnetic_p3m_print_mc_time(Tcl_Interp *interp, int mesh, int cao,
			  double r_cut_iL_min, double r_cut_iL_max, double *_r_cut_iL,
			  double *_alpha_L, double *_accuracy)
{
  double int_time;
  double r_cut_iL;
  double rs_err, ks_err, mesh_size, k_cut;
  int i, n_cells;
  char b1[TCL_DOUBLE_SPACE + 12],b2[TCL_DOUBLE_SPACE + 12],b3[TCL_DOUBLE_SPACE + 12];
  /* initial checks. */
  mesh_size = box_l[0]/(double)mesh;
  k_cut =  mesh_size*cao/2.0;
  P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_mc_time: mesh=%d, cao=%d, rmin=%f, rmax=%f\n",
		    mesh, cao, r_cut_iL_min, r_cut_iL_max));
  if(cao >= mesh || k_cut >= dmin(min_box_l,min_local_box_l) - skin) {
    /* print result */
    sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
    Tcl_AppendResult(interp, b2," ", b3," cao too large for this mesh\n", (char *) NULL);
    return -3;
  }

  /* Either low and high boundary are equal (for fixed cut), or the low border is initially 0 and therefore
     has infinite error estimate, as required. Therefore if the high boundary fails, there is no possible r_cut */
  if ((*_accuracy = dp3m_get_accuracy(mesh, cao, r_cut_iL_max, _alpha_L, &rs_err, &ks_err)) > dp3m.params.accuracy) {
    /* print result */
    sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
    Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
    sprintf(b1,"%.5e",r_cut_iL_max); sprintf(b2,"%.5e",*_alpha_L); sprintf(b3,"%.5e",*_accuracy);
    Tcl_AppendResult(interp, b1,"  ", b2,"  ",b3," ", (char *) NULL);
    sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err);
    Tcl_AppendResult(interp, b1,"  ", b2,"  accuracy not achieved\n", (char *) NULL);
    return -2;
  }

  for (;;) {
    P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_mc_time: interval [%f,%f]\n", r_cut_iL_min, r_cut_iL_max));
    r_cut_iL = 0.5*(r_cut_iL_min + r_cut_iL_max);

    if (r_cut_iL_max - r_cut_iL_min < P3M_RCUT_PREC)
      break;

    /* bisection */
    if (dp3m_get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err) > dp3m.params.accuracy)
      r_cut_iL_min = r_cut_iL;
    else
      r_cut_iL_max = r_cut_iL;
  }
  /* final result is always the upper interval boundary, since only there
     we know that the desired minimal accuracy is obtained */
  *_r_cut_iL = r_cut_iL = r_cut_iL_max;

  /* check whether we are running P3M+DLC, and whether we leave a reasonable gap space */
  if (coulomb.Dmethod == DIPOLAR_MDLC_P3M)   fprintf(stderr, "tunning when dlc needs to be fixed, p3m-dipoles.c \n");
  
  /*
  needs to be fixed
  if (coulomb.method == DIPOLAR_MDLC_P3M && elc_params.gap_size <= 1.1*r_cut_iL*box_l[0]) {
    P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_mc_time: mesh %d cao %d r_cut %f reject r_cut %f > gap %f\n", mesh, cao, r_cut_iL,
		      2*r_cut_iL*box_l[0], elc_params.gap_size));
    // print result 
    sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
    Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
    sprintf(b1,"%.5e",r_cut_iL_max); sprintf(b2,"%.5e",*_alpha_L); sprintf(b3,"%.5e",*_accuracy);
    Tcl_AppendResult(interp, b1,"  ", b2,"  ",b3," ", (char *) NULL);
    sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err);
    Tcl_AppendResult(interp, b1,"  ", b2,"  conflict with ELC\n", (char *) NULL);
    return -2;
  }
  */

  /* check whether this radius is too large, so that we would use less cells than allowed */
  n_cells = 1;
  for (i = 0; i < 3; i++)
    n_cells *= (int)(floor(local_box_l[i]/(r_cut_iL*box_l[0] + skin)));
  if (n_cells < min_num_cells) {
    P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_mc_time: mesh %d cao %d r_cut %f reject n_cells %d\n", mesh, cao, r_cut_iL, n_cells));
    /* print result */
    sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
    Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
    sprintf(b1,"%.5e",r_cut_iL_max); sprintf(b2,"%.5e",*_alpha_L); sprintf(b3,"%.5e",*_accuracy);
    Tcl_AppendResult(interp, b1,"  ", b2,"  ",b3," ", (char *) NULL);
    sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err);
    Tcl_AppendResult(interp, b1,"  ", b2,"  radius dangerously high\n", (char *) NULL);
    return -2;
  }

  int_time = dp3m_mcr_time(mesh, cao, r_cut_iL, *_alpha_L);
  if (int_time == -1) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "tuning failed, test integration not possible", (char *)NULL);
    return -1;
  }

  *_accuracy = dp3m_get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err);

  P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_mc_time: mesh %d cao %d r_cut %f time %f\n", mesh, cao, r_cut_iL, int_time));
  /* print result */
  sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
  Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
  sprintf(b1,"%.5e",r_cut_iL); sprintf(b2,"%.5e",*_alpha_L); sprintf(b3,"%.5e",*_accuracy);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ",b3," ", (char *) NULL);
  sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err); sprintf(b3,"%-8d",(int)int_time);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"\n", (char *) NULL);

  return int_time;
}


/*****************************************************************************/

/** get the optimal alpha and the corresponding computation time for fixed mesh. *cao
    should contain an initial guess, which is then adapted by stepping up and down. Returns the time
    upon completion, -1 if the force evaluation does not work, and -2 if the accuracy cannot be met */
static double tclcommand_inter_magnetic_p3m_print_m_time(Tcl_Interp *interp, int mesh,
			 int cao_min, int cao_max, int *_cao,
			 double r_cut_iL_min, double r_cut_iL_max, double *_r_cut_iL,
			 double *_alpha_L, double *_accuracy)
{
  double best_time = -1, tmp_time, tmp_r_cut_iL, tmp_alpha_L=0.0, tmp_accuracy=0.0;
  /* in which direction improvement is possible. Initially, we dont know it yet. */
  int final_dir = 0;
  int cao = *_cao;

  P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_m_time: Dmesh=%d, Dcao_min=%d, Dcao_max=%d, Drmin=%f, Drmax=%f\n",
		    mesh, cao_min, cao_max, r_cut_iL_min, r_cut_iL_max));
  /* the initial step sets a timing mark. If there is no valid r_cut, we can only try
     to increase cao to increase the obtainable precision of the far formula. */
  do {
    tmp_time = tclcommand_inter_magnetic_p3m_print_mc_time(interp, mesh, cao,  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out if the force evaluation is not working */
    if (tmp_time == -1) return -1;
    /* cao is too large for this grid, but still the accuracy cannot be achieved, give up */
    if (tmp_time == -3) {
      P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_m_time: no possible cao found\n"));
      return -2;
    }
    /* we have a valid time, start optimising from there */
    if (tmp_time >= 0) {
      best_time  = tmp_time;
      *_r_cut_iL = tmp_r_cut_iL;
      *_alpha_L  = tmp_alpha_L;
      *_accuracy = tmp_accuracy;
      *_cao      = cao;
      break;
    }
    /* the required accuracy could not be obtained, try higher caos. Therefore optimisation can only be
       obtained with even higher caos, but not lower ones */
    P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_m_time: doesn't give precision, step up\n"));
    cao++;
    final_dir = 1;
  }
  while (cao <= cao_max);
  /* with this mesh, the required accuracy cannot be obtained. */
  if (cao > cao_max) return -2;

  /* at the boundaries, only the opposite direction can be used for optimisation */
  if (cao == cao_min)      final_dir = 1;
  else if (cao == cao_max) final_dir = -1;

  P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_m_time: final constraints dir %d\n", final_dir));

  if (final_dir == 0) {
    /* check in which direction we can optimise. Both directions are possible */
    double dir_times[3];
    for (final_dir = -1; final_dir <= 1; final_dir += 2) {
      dir_times[final_dir + 1] = tmp_time =
	tclcommand_inter_magnetic_p3m_print_mc_time(interp, mesh, cao + final_dir,  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
      /* bail out on errors, as usual */
      if (tmp_time == -1) return -1;
      /* in this direction, we cannot optimise, since we get into precision trouble */
      if (tmp_time < 0) continue;

      if (tmp_time < best_time) {
	best_time  = tmp_time;
	*_r_cut_iL = tmp_r_cut_iL;
	*_alpha_L  = tmp_alpha_L;
	*_accuracy = tmp_accuracy;
	*_cao      = cao + final_dir;
      }
    }
    /* choose the direction which was optimal, if any of the two */
    if      (dir_times[0] == best_time) { final_dir = -1; }
    else if (dir_times[2] == best_time) { final_dir = 1; }
    else {
      /* no improvement in either direction, however if one is only marginally worse, we can still try*/
      /* down is possible and not much worse, while up is either illegal or even worse */
      if ((dir_times[0] >= 0 && dir_times[0] < best_time + P3M_TIME_GRAN) &&
	  (dir_times[2] < 0 || dir_times[2] > dir_times[0]))
	final_dir = -1;
      /* same for up */
      else if ((dir_times[2] >= 0 && dir_times[2] < best_time + P3M_TIME_GRAN) &&
	       (dir_times[0] < 0 || dir_times[0] > dir_times[2]))
	final_dir = 1;
      else {
	/* really no chance for optimisation */
	P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_m_time: Dmesh=%d final Dcao=%d time=%f\n",mesh, cao, best_time));
	return best_time;
      }
    }
    /* we already checked the initial cao and its neighbor */
    cao += 2*final_dir;
  }
  else {
    /* here some constraint is active, and we only checked the initial cao itself */
    cao += final_dir;
  }

  P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_m_time: optimise in direction %d\n", final_dir));

  /* move cao into the optimisation direction until we do not gain anymore. */
  for (; cao >= cao_min && cao <= cao_max; cao += final_dir) {
    tmp_time = tclcommand_inter_magnetic_p3m_print_mc_time(interp, mesh, cao,  r_cut_iL_min, r_cut_iL_max,
			   &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out on errors, as usual */
    if (tmp_time == -1) return -1;
    /* if we cannot meet the precision anymore, give up */
    if (tmp_time < 0) break;

    if (tmp_time < best_time) {
      best_time  = tmp_time;
      *_r_cut_iL = tmp_r_cut_iL;
      *_alpha_L  = tmp_alpha_L;
      *_accuracy = tmp_accuracy;
      *_cao      = cao;
    }
    /* no hope of further optimisation */
    else if (tmp_time > best_time + P3M_TIME_GRAN)
      break;
  }
  P3M_TRACE(fprintf(stderr, "tclcommand_inter_magnetic_p3m_print_m_time: Dmesh=%d final Dcao=%d Dr_cut=%f time=%f\n",mesh, *_cao, *_r_cut_iL, best_time));
  return best_time;
}


/** a probably faster adaptive tuning method. Uses the same error estimates and parameters as
    \ref tclcommand_inter_magnetic_dp3m_print_adaptive_tune_parameters, but a different strategy for finding the optimum. The algorithm
    basically determines the mesh, cao and then the real space cutoff, in this nested order.

    For each mesh, the cao optimal for the mesh tested previously is used as an initial guess,
    and the algorithm tries whether increasing or decreasing it leads to a better solution. This
    is efficient, since the optimal cao only changes little with the meshes in general.

    The real space cutoff for a given mesh and cao is determined via a bisection on the error estimate,
    which determines where the error estimate equals the required accuracy. Therefore the smallest 
    possible, i.e. fastest real space cutoff is determined.

    Both the search over mesh and cao stop to search in a specific direction once the computation time is
    significantly higher than the currently known optimum.

    Compared to \ref tclcommand_inter_magnetic_dp3m_print_tune_parameters, this function will test more parameters sets for efficiency, but
    the error estimate is calculated less often. In general this should be faster and give better results.
 */

static int tclcommand_inter_magnetic_dp3m_print_adaptive_tune_parameters(Tcl_Interp *interp)
{
  int    mesh_max,                   mesh     = -1, tmp_mesh;
  double r_cut_iL_min, r_cut_iL_max, r_cut_iL = -1, tmp_r_cut_iL=0.0;
  int    cao_min, cao_max,           cao      = -1, tmp_cao;

  double                             alpha_L  = -1, tmp_alpha_L=0.0;
  double                             accuracy = -1, tmp_accuracy=0.0;
  double                            time_best=1e20, tmp_time;
  char
    b1[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 12],
    b2[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 12],
    b3[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 17];
 
  P3M_TRACE(fprintf(stderr,"%d: dipolar P3M_adaptive_tune_parameters\n",this_node));

  if (skin == -1) {
    Tcl_AppendResult(interp, "dipolar p3m cannot be tuned, since the skin is not yet set", (char *) NULL);
    return TCL_ERROR;
  }

  /* preparation */
  mpi_bcast_event(P3M_COUNT_DIPOLES);

  /* Print Status */
  sprintf(b1,"%.5e",dp3m.params.accuracy);
  Tcl_AppendResult(interp, "dipolar P3M tune parameters: Accuracy goal = ",b1,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, box_l[0], b1);

  sprintf(b2,"%d",dp3m.sum_dip_part); 

  Tcl_PrintDouble(interp, dp3m.sum_mu2, b3);
  Tcl_AppendResult(interp, "System: box_l = ",b1,", # charged part = ",b2," Sum[q_i^2] = ",b3,"\n", (char *) NULL);

  if (dp3m.sum_dip_part == 0) {
    Tcl_AppendResult(interp, "no dipolar particles in the system, cannot tune dipolar P3M", (char *) NULL);
    return (TCL_ERROR);
  }
  
  
  /* parameter ranges */
  if (dp3m.params.mesh[0] == 0 ) {
    double expo;
    expo = log(pow((double)dp3m.sum_dip_part,(1.0/3.0)))/log(2.0);  

    tmp_mesh = (int)(pow(2.0,(double)((int)expo))+0.1);
    /* this limits the tried meshes if the accuracy cannot
       be obtained with smaller meshes, but normally not all these
       meshes have to be tested */
    mesh_max = tmp_mesh * 256;
    /* avoid using more than 1 GB of FFT arrays (per default, see config.h) */
    if (mesh_max > P3M_MAX_MESH)
      mesh_max = P3M_MAX_MESH;
  }
  else {
    sprintf(b1, "%d", dp3m.params.mesh[0]);
    Tcl_AppendResult(interp, "fixed mesh ", b1, "\n", (char *)NULL);
    tmp_mesh = mesh_max = dp3m.params.mesh[0];
  }

  if(dp3m.params.r_cut_iL == 0.0) {
    r_cut_iL_min = 0;
    r_cut_iL_max = min_local_box_l/2 - skin;
    r_cut_iL_min *= box_l_i[0];
    r_cut_iL_max *= box_l_i[0];
  }
  else {
    sprintf(b1, "%f", dp3m.params.r_cut_iL);
    Tcl_AppendResult(interp, "fixed r_cut_iL ", b1, "\n", (char *)NULL);
    r_cut_iL_min = r_cut_iL_max = dp3m.params.r_cut_iL;
  }

  if(dp3m.params.cao == 0) {
    cao_min = 1;
    cao_max = 7;
    cao = 3;
  }
  else {
    sprintf(b1, "%d", dp3m.params.cao);
    Tcl_AppendResult(interp, "fixed cao ", b1, "\n", (char *)NULL);
    cao_min = cao_max = cao = dp3m.params.cao;
  }

  Tcl_AppendResult(interp, "Dmesh cao Dr_cut_iL   Dalpha_L     Derr         Drs_err    Dks_err    time [ms]\n", (char *) NULL);

  /* mesh loop */
  for (;tmp_mesh <= mesh_max; tmp_mesh *= 2) {
    tmp_cao = cao;
    tmp_time = tclcommand_inter_magnetic_p3m_print_m_time(interp, tmp_mesh,
			  cao_min, cao_max, &tmp_cao,
			  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL,
			  &tmp_alpha_L, &tmp_accuracy);
    /* some error occured during the tuning force evaluation */
    if (tmp_time == -1) return TCL_ERROR;
    /* this mesh does not work at all */
    if (tmp_time < 0) continue;

    /* the optimum r_cut for this mesh is the upper limit for higher meshes,
       everything else is slower */
    r_cut_iL_max = tmp_r_cut_iL;

    /* new optimum */
    if (tmp_time < time_best) {
      time_best = tmp_time;
      mesh      = tmp_mesh;
      cao       = tmp_cao;
      r_cut_iL  = tmp_r_cut_iL;
      alpha_L   = tmp_alpha_L;
      accuracy  = tmp_accuracy;
    }
    /* no hope of further optimisation */
    else if (tmp_time > time_best + P3M_TIME_GRAN)
      break;
  }
  
  P3M_TRACE(fprintf(stderr,"finshed tuning\n"));
  if(time_best == 1e20) {
    Tcl_AppendResult(interp, "failed to tune dipolar P3M parameters to required accuracy", (char *) NULL);
    return (TCL_ERROR);
  }

  /* set tuned p3m parameters */
  dp3m.params.r_cut_iL = r_cut_iL;
  dp3m.params.mesh[0]  = dp3m.params.mesh[1] = dp3m.params.mesh[2] = mesh;
  dp3m.params.cao      = cao;
  dp3m.params.alpha_L  = alpha_L;
  dp3m.params.accuracy = accuracy;
  dp3m_scaleby_box_l();
  /* broadcast tuned p3m parameters */
  mpi_bcast_coulomb_params();
  /* Tell the user about the outcome */
  Tcl_AppendResult(interp, "\nresulting parameters:\n", (char *) NULL);
  sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
  Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
  sprintf(b1,"%.5e",r_cut_iL); sprintf(b2,"%.5e",alpha_L); sprintf(b3,"%.5e",accuracy);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);
  sprintf(b3,"                 %-8d",(int)time_best);
  Tcl_AppendResult(interp, b3, (char *) NULL);
  return (TCL_OK);
}

#endif /* DP3M */

