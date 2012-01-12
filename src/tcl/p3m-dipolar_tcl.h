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
#ifndef _P3M_MAGNETOSTATICS_TCL_H
#define _P3M_MAGNETOSTATICS_TCL_H
/** \file p3m-magnetostatics.h P3M algorithm for long range magnetic dipole-dipole interaction.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the dipolar Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm. The file p3m contains only the
 *  Particle-Mesh part.
 *
 *  Further reading: 
 *  <ul>
 *  <li> J.J. Cerda, P3M for dipolar interactions. J. Chem. Phys, 129, xxx ,(2008).
 *  </ul>
 *
 */
#include "p3m-common.h"
#include "interaction_data.h"
#include "p3m-dipolar.h"

#ifdef DP3M

//typedef struct {
//  p3m_parameter_struct params;
//
//  /** local mesh. */
//  p3m_local_mesh local_mesh;
//  /** real space mesh (local) for CA/FFT.*/
//  double *rs_mesh;
//  /** real space mesh (local) for CA/FFT of the dipolar field.*/
//  double *rs_mesh_dip[3];
//  /** k space mesh (local) for k space calculation and FFT.*/
//  double *ks_mesh;
//
//  /** number of dipolar particles (only on master node). */
//  int sum_dip_part; 
//  /** Sum of square of magnetic dipoles (only on master node). */
//  double sum_mu2;
//
//  /** interpolation of the charge assignment function. */
//  double *int_caf[7];
//
//  /** position shift for calc. of first assignment mesh point. */
//  double pos_shift;
//  /** help variable for calculation of aliasing sums */
//  double *meshift;
//
//  /** Spatial differential operator in k-space. We use an i*k differentiation. */
//  double *d_op;
//  /** Force optimised influence function (k-space) */
//  double *g_force;
//  /** Energy optimised influence function (k-space) */
//  double *g_energy;
//
//  /** number of charged particles on the node. */
//  int ca_num;
//
//  /** Charge fractions for mesh assignment. */
//  double *ca_frac;
//  /** index of first mesh point for charge assignment. */
//  int *ca_fmp;
//  /** number of permutations in k_space */
//  int ks_pnum;
//
//  /** send/recv mesh sizes */
//  p3m_send_mesh  sm;
//
//  /** Field to store grid points to send. */
//  double *send_grid; 
//  /** Field to store grid points to recv */
//  double *recv_grid;
//
//  /* Stores the value of the energy correction due to MS effects */
//  double  energy_correction;
//
//  /** Flag to know if we should calculate the constants for the energy 
//      (If you neither compute the energy, is a waste of time
//      spendig circa 3 or 4 min computing such constants) **/
//  int flag_constants_energy_dipolar;
//
//} dp3m_data_struct;
//
/** dipolar P3M parameters. */
extern dp3m_data_struct dp3m;

/** \name Exported Functions */
/************************************************************/
/*@{*/
/** dipolar p3m parser */
int tclcommand_inter_magnetic_parse_dp3m(Tcl_Interp * interp, int argc, char ** argv);

/** dipolar p3m parser, optional parameters */
int tclcommand_inter_magnetic_parse_dp3m_opt_params(Tcl_Interp * interp, int argc, char ** argv);

/** print the p3m parameters to the interpreters result */
int tclprint_to_result_dp3m(Tcl_Interp *interp);



#endif /* DP3M */
#endif /* _P3M_DIPOLES_H */
