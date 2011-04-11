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
#include "fft.h"
#include "p3m.h"
#include "thermostat.h"
#include "cells.h"
#include "tuning.h"
#include "elc.h"

#ifdef ELP3M

/* defined only within this header file, for checking that system
   (electro/magnetostatic) specific files are only included from
   here. */
#define P3M_C_CURRENT

/************************************************
 * variables
 ************************************************/

p3m_struct p3m = { 
#ifdef ELECTROSTATICS 
  0.0, 0.0, 
  {0,0,0}, {P3M_MESHOFF, P3M_MESHOFF, P3M_MESHOFF}, 
  0, P3M_N_INTERPOL, 0.0, P3M_EPSILON, 
  {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}, 0.0, 0.0, 0, 0, {0, 0, 0},
#endif
		   
#ifdef MAGNETOSTATICS
  0.0, 0.0, 
  {0,0,0}, {P3M_MESHOFF, P3M_MESHOFF, P3M_MESHOFF}, 
  0, P3M_N_INTERPOL, 0.0, P3M_EPSILON_MAGNETIC, 
  {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}, 0.0, 0.0, 0, 0, {0, 0, 0},
#endif
};


/*********************** miscelanea of functions *************************************/

int tclprint_to_result_ChargeP3M(Tcl_Interp *interp)
{
#ifdef ELECTROSTATICS
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, p3m.r_cut, buffer);
  Tcl_AppendResult(interp, "p3m ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.mesh[0]);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.cao);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.alpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.accuracy, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  Tcl_AppendResult(interp, "} {coulomb epsilon ", (char *) NULL);
  if (p3m.epsilon == P3M_EPSILON_METALLIC)
    Tcl_AppendResult(interp, " metallic ", (char *) NULL);
  else {
    Tcl_PrintDouble(interp, p3m.epsilon, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  }
  sprintf(buffer,"%d",p3m.inter);
  Tcl_AppendResult(interp, "n_interpol ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.mesh_off[0], buffer);
  Tcl_AppendResult(interp, "mesh_off ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.mesh_off[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.mesh_off[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
#endif

  return TCL_OK;
}

int tclprint_to_result_DipolarP3M(Tcl_Interp *interp)
{
#ifdef MAGNETOSTATICS
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, p3m.Dr_cut, buffer);
  Tcl_AppendResult(interp, "p3m ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.Dmesh[0]);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.Dcao);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.Dalpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.Daccuracy, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  Tcl_AppendResult(interp, "} {magnetic epsilon ", (char *) NULL);
  if (p3m.Depsilon == P3M_EPSILON_METALLIC)
    Tcl_AppendResult(interp, "metallic ", (char *) NULL);
  else {
    Tcl_PrintDouble(interp, p3m.Depsilon, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  }
  sprintf(buffer,"%d",p3m.Dinter);
  Tcl_AppendResult(interp, "n_interpol ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.Dmesh_off[0], buffer);
  Tcl_AppendResult(interp, "mesh_off ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.Dmesh_off[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.Dmesh_off[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
#endif 

  return TCL_OK;
}


double P3M_calc_kspace_forces(int force_flag, int energy_flag)
{
  double k_space_energy=0;
  
  #ifdef ELECTROSTATICS
    k_space_energy+=P3M_calc_kspace_forces_for_charges( force_flag,  energy_flag);  /* see p3m-charges.c */
  #endif
    
  #ifdef MAGNETOSTATICS
    k_space_energy+=P3M_calc_kspace_forces_for_dipoles( force_flag,  energy_flag); /* see p3m-dipoles.c */
  #endif
  
  return k_space_energy;
}

/************************************************/

//Note: this function P3m_exit is actually never called !!!
void   P3M_exit()
{
  int i;
  /* free memory */

#ifdef ELECTROSTATICS
  free(ca_frac);
  free(ca_fmp);
  free(send_grid);
  free(recv_grid);
  free(rs_mesh);
  free(ks_mesh); 
  for(i=0; i<p3m.cao; i++) free(int_caf[i]);
#endif
  
#ifdef MAGNETOSTATICS
  for (i=0;i<3;i++) free(Drs_mesh_dip[i]);
  free(Dca_frac);
  free(Dca_fmp);
  free(Dsend_grid);
  free(Drecv_grid);
  free(Drs_mesh);
  free(Dks_mesh); 
#endif
}


/************************************************
 * Debug functions printing p3m structures 
 ************************************************/

void p3m_print_p3m_struct(p3m_struct ps) {
#ifdef ELECTROSTATICS  
  fprintf(stderr,"%d: p3m_struct: \n",this_node);
  fprintf(stderr,"   alpha_L=%f, r_cut_iL=%f \n",
	  ps.alpha_L,ps.r_cut_iL);
  fprintf(stderr,"   mesh=(%d,%d,%d), mesh_off=(%.4f,%.4f,%.4f)\n",
	  ps.mesh[0],ps.mesh[1],ps.mesh[2],
	  ps.mesh_off[0],ps.mesh_off[1],ps.mesh_off[2]);
  fprintf(stderr,"   cao=%d, inter=%d, epsilon=%f\n",
	  ps.cao,ps.inter,ps.epsilon);
  fprintf(stderr,"   cao_cut=(%f,%f,%f)\n",
	  ps.cao_cut[0],ps.cao_cut[1],ps.cao_cut[2]);
  fprintf(stderr,"   a=(%f,%f,%f), ai=(%f,%f,%f)\n",
	  ps.a[0],ps.a[1],ps.a[2],ps.ai[0],ps.ai[1],ps.ai[2]);
#endif	  
	  
#ifdef MAGNETOSTATICS
  fprintf(stderr,"%d: dipolar p3m_struct: \n",this_node);
  fprintf(stderr,"   Dalpha_L=%f, Dr_cut_iL=%f \n",
	  ps.Dalpha_L,ps.Dr_cut_iL);
  fprintf(stderr,"   Dmesh=(%d,%d,%d), Dmesh_off=(%.4f,%.4f,%.4f)\n",
	  ps.Dmesh[0],ps.Dmesh[1],ps.Dmesh[2],
	  ps.Dmesh_off[0],ps.Dmesh_off[1],ps.Dmesh_off[2]);
  fprintf(stderr,"   Dcao=%d, Dinter=%d, Depsilon=%f\n",
	  ps.Dcao,ps.Dinter,ps.Depsilon);
  fprintf(stderr,"   Dcao_cut=(%f,%f,%f)\n",
	  ps.Dcao_cut[0],ps.Dcao_cut[1],ps.Dcao_cut[2]);
  fprintf(stderr,"   Da=(%f,%f,%f), Dai=(%f,%f,%f)\n",
	  ps.Da[0],ps.Da[1],ps.Da[2],ps.Dai[0],ps.Dai[1],ps.Dai[2]);
#endif	  
}

#undef P3M_C_CURRENT

#endif /* of ELP3M */

