// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file p3m.c  P3M algorithms, main file to handle all the long range interactions via P3M algorithms.
 *               Charge-Charge and magnetic dipole-dipole interactions are installed at this moment.
 *               Coming soon: charge-electrical dipole interaction.
 *     
 *  For more information about the different p3m algorithms,
 *  see \ref p3m.h "p3m.h"
 *  see \ref p3m-charges.c  "p3m-charges.c"
 *  see \ref p3m-charges.h  "p3m-charges.h"
 *  see \ref p3m-dipoles.c  "p3m-dipoles.c"
 *  see \ref p3m-dipoles.h  "p3m-dipoles.h"
 *  see \ref p3m-assignment.c  "p3m-assignment.c"
 
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

/* defined only within this header file, for checking that system (electro/magnetostatic)
   specific files are only included from here. */
#define P3M_C_CURRENT

/************************************************
 * DEFINES
 ************************************************/
   
/** increment size of charge assignment fields. */
#define CA_INCREMENT 32       
/** precision limit for the r_cut zero */
#define P3M_RCUT_PREC 1e-3
/** granularity of the time measurement */
#define P3M_TIME_GRAN 2

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
  0, P3M_N_INTERPOL, 0.0, P3M_EPSILON, 
  {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}, 0.0, 0.0, 0, 0, {0, 0, 0},
#endif
};

/** \name Private declarations of common grid related functions */
/************************************************************/
/*@{*/

/** print local mesh content. 
    \param l local mesh structure.
*/
void p3m_print_local_mesh(local_mesh l);

/** print send mesh content. 
 *  \param sm send mesh structure.
*/
void p3m_print_send_mesh(send_mesh sm);

/** Add values of a 3d-grid input block (size[3]) to values of 3d-grid
 *  ouput array with dimension dim[3] at start position start[3].  
 *
 *  \param in          Pointer to first element of input block data.
 *  \param out         Pointer to first element of output grid.
 *  \param start       Start position of block in output grid.
 *  \param size        Dimensions of the block
 *  \param dim         Dimensions of the output grid.
*/
void add_block(double *in, double *out, int start[3], int size[3], int dim[3]);

/** One of the aliasing sums used by \ref P3M_k_space_error. 
    (fortunately the one which is most important (because it converges
    most slowly, since it is not damped exponentially)) can be
    calculated analytically. The result (which depends on the order of
    the spline interpolation) can be written as an even trigonometric
    polynomial. The results are tabulated here (The employed formula
    is Eqn. 7.66 in the book of Hockney and Eastwood). */
double analytic_cotangent_sum(int n, double mesh_i, int cao);

/*@}*/

/*********************** miscelanea of functions *************************************/

int printP3MToResult(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

#ifdef ELECTROSTATICS
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
  
#ifdef MAGNETOSTATICS
  Tcl_PrintDouble(interp, p3m.Dr_cut, buffer);
  Tcl_AppendResult(interp, "dipolar p3m ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.Dmesh[0]);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.Dcao);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.Dalpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.Daccuracy, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  Tcl_AppendResult(interp, "} {dipolar epsilon ", (char *) NULL);
  if (p3m.Depsilon == P3M_EPSILON_METALLIC)
    Tcl_AppendResult(interp, " metallic ", (char *) NULL);
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

/************************************************/

void   P3M_init() {
#ifdef ELECTROSTATICS
  P3M_init_charges();
#endif
  
#ifdef  MAGNETOSTATICS  
  P3M_init_dipoles();
#endif    
}

/******************* common p3m assignment functions *****************************************/

#include "p3m-assignment.c"
/* include at this point the functions commmon
   to all p3m-algorithms related to the mapping to a lattice
*/

/****** specific functions for the different p3m algorithms installed ************************/

#ifdef ELECTROSTATICS
#include "p3m-charges.c"  /* for charge-charge interaction via p3m */
#endif

#ifdef MAGNETOSTATICS
#include "p3m-dipoles.c"   /* for  magnetic dipole-dipole interaction via p3m */
#endif 

/******** reciprocal part  (the real one is placed in the header files) *****************************************/

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

void p3m_print_local_mesh(local_mesh l) 
{
  fprintf(stderr,"%d: local_mesh: dim=(%d,%d,%d), size=%d\n",this_node,
	  l.dim[0],l.dim[1],l.dim[2],l.size);
  fprintf(stderr,"%d:    ld_ind=(%d,%d,%d), ld_pos=(%f,%f,%f)\n",this_node,
	  l.ld_ind[0],l.ld_ind[1],l.ld_ind[2],
	  l.ld_pos[0],l.ld_pos[1],l.ld_pos[2]);
  fprintf(stderr,"%d:    inner=(%d,%d,%d) [(%d,%d,%d)-(%d,%d,%d)]\n",this_node,
	  l.inner[0],l.inner[1],l.inner[2],
	  l.in_ld[0],l.in_ld[1],l.in_ld[2],
	  l.in_ur[0],l.in_ur[1],l.in_ur[2]);
  fprintf(stderr,"%d:    margin = (%d,%d, %d,%d, %d,%d)\n",this_node,
	  l.margin[0],l.margin[1],l.margin[2],l.margin[3],l.margin[4],l.margin[5]);
  fprintf(stderr,"%d:    r_margin=(%d,%d, %d,%d, %d,%d)\n",this_node,
	  l.r_margin[0],l.r_margin[1],l.r_margin[2],l.r_margin[3],l.r_margin[4],l.r_margin[5]);
}

/************************************************************/

void p3m_print_send_mesh(send_mesh sm) 
{
  int i;
  fprintf(stderr,"%d: send_mesh: max=%d\n",this_node,sm.max);
  for(i=0;i<6;i++) {
    fprintf(stderr,"%d:  dir=%d: s_dim (%d,%d,%d)  s_ld (%d,%d,%d) s_ur (%d,%d,%d) s_size=%d\n",this_node,i,sm.s_dim[i][0],sm.s_dim[i][1],sm.s_dim[i][2],sm.s_ld[i][0],sm.s_ld[i][1],sm.s_ld[i][2],sm.s_ur[i][0],sm.s_ur[i][1],sm.s_ur[i][2],sm.s_size[i]);
    fprintf(stderr,"%d:         r_dim (%d,%d,%d)  r_ld (%d,%d,%d) r_ur (%d,%d,%d) r_size=%d\n",this_node,sm.r_dim[i][0],sm.r_dim[i][1],sm.r_dim[i][2],sm.r_ld[i][0],sm.r_ld[i][1],sm.r_ld[i][2],sm.r_ur[i][0],sm.r_ur[i][1],sm.r_ur[i][2],sm.r_size[i]);
  }
}

/************************************************************/

void add_block(double *in, double *out, int start[3], int size[3], int dim[3])
{
  /* fast,mid and slow changing indices */
  int f,m,s;
  /* linear index of in grid, linear index of out grid */
  int li_in=0,li_out=0;
  /* offsets for indizes in output grid */
  int m_out_offset,s_out_offset;

  li_out = start[2] + ( dim[2]*( start[1] + (dim[1]*start[0]) ) );
  m_out_offset  = dim[2] - size[2];
  s_out_offset  = (dim[2] * (dim[1] - size[1]));

  for(s=0 ;s<size[0]; s++) {
    for(m=0; m<size[1]; m++) {
      for(f=0; f<size[2]; f++) {
	out[li_out++] += in[li_in++];
      }
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}

/************************************************************/

double analytic_cotangent_sum(int n, double mesh_i, int cao)
{
  double c, res=0.0;
  c = SQR(cos(PI*mesh_i*(double)n));

  switch (cao) {
  case 1 : { 
    res = 1; 
    break; }
  case 2 : { 
    res = (1.0+c*2.0)/3.0; 
    break; }
  case 3 : { 
    res = (2.0+c*(11.0+c*2.0))/15.0; 
    break; }
  case 4 : { 
    res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0; 
    break; }
  case 5 : { 
    res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0; 
    break; }
  case 6 : { 
    res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0; 
    break; }
  case 7 : { 
    res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0; 
    break; }
  default : {
    fprintf(stderr,"%d: INTERNAL_ERROR: The value %d for the interpolation order should not occur!\n",this_node, cao);
    errexit();
  }
  }
  
  return res;
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

