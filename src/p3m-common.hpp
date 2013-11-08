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
#ifndef _P3M_COMMON_H 
#define _P3M_COMMON_H
/** \file p3m-common.hpp   common functions for dipolar and charge p3m.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm. The file p3m contains only the
 *  Particle-Mesh part.
 *
 *  Further reading: 
 *  <ul>
 *  <li> P.P. Ewald,
 *       <i>Die Berechnung optischer und elektrostatischer Gitterpotentiale</i>,
 *       Ann. Phys. (64) 253-287, 1921
 *  <li> R. W. Hockney and J. W. Eastwood, 
 *       <i>Computer Simulation Using Particles</i>,
 *       IOP, London, 1988
 *  <li> M. Deserno and C. Holm,
 *       <i>How to mesh up {E}wald sums. I. + II.</i>,
 *       J. Chem. Phys. (109) 7678, 1998; (109) 7694, 1998
 *  <li> M. Deserno, C. Holm and H. J. Limbach,
 *       <i>How to mesh up {E}wald sums. </i>,
 *       in Molecular Dynamics on Parallel Computers,
 *       Ed. R. Esser et al., World Scientific, Singapore, 2000
 *  <li> M. Deserno,
 *       <i>Counterion condensation for rigid linear polyelectrolytes</i>,
 *       PhdThesis, Universit{\"a}t Mainz, 2000
 *  <li> J.J. Cerda, P3M for dipolar interactions. J. Chem. Phys, 129, xxx ,(2008).
 *  </ul>
 *
 */
#include "config.hpp"
#include "utils.hpp"

#if defined(P3M) || defined(DP3M)

/** This value for p3m.epsilon indicates metallic boundary conditions. */
#define P3M_EPSILON_METALLIC 0.0

/** increment size of charge assignment fields. */
#define CA_INCREMENT 32       
/** precision limit for the r_cut zero */
#define P3M_RCUT_PREC 1e-3
/** granularity of the time measurement */
#define P3M_TIME_GRAN 2

/** whether the P3M charge assignment fraction is stored or not */
#define P3M_STORE_CA_FRAC

/************************************************
 * data types
 ************************************************/

/** Structure for local mesh parameters. */
typedef struct {
  /* local mesh characterization. */
  /** dimension (size) of local mesh. */
  int dim[3];
  /** number of local mesh points. */
  int size;
  /** index of lower left corner of the 
      local mesh in the global mesh. */
  int ld_ind[3];
  /** position of the first local mesh point. */
  double ld_pos[3];
  /** dimension of mesh inside node domain. */
  int inner[3];
  /** inner left down grid point */
  int in_ld[3];
  /** inner up right grid point + (1,1,1)*/
  int in_ur[3];
  /** number of margin mesh points. */
  int margin[6];
  /** number of margin mesh points from neighbour nodes */
  int r_margin[6];
  /** offset between mesh lines of the last dimension */
  int q_2_off;
  /** offset between mesh lines of the two last dimensions */
  int q_21_off;
} p3m_local_mesh;

/** Structure for send/recv meshs. */
typedef struct {
  /** dimension of sub meshs to send. */
  int s_dim[6][3];
  /** left down corners of sub meshs to send. */
  int s_ld[6][3];
  /** up right corners of sub meshs to send. */
  int s_ur[6][3];
  /** sizes for send buffers. */
  int s_size[6];
  /** dimensionof sub meshs to recv. */
  int r_dim[6][3];
  /** left down corners of sub meshs to recv. */
  int r_ld[6][3];
  /** up right corners of sub meshs to recv. */
  int r_ur[6][3];
  /** sizes for recv buffers. */
  int r_size[6];
  /** maximal size for send/recv buffers. */
  int max;
} p3m_send_mesh;

/** Structure to hold P3M parameters and some dependend variables. */
typedef struct {
  /** Ewald splitting parameter (0<alpha<1), rescaled to alpha_L = alpha * box_l. */
  double alpha_L;
  /** Cutoff radius for real space electrostatics (>0), rescaled to r_cut_iL = r_cut * box_l_i. */
  double r_cut_iL;
  /** number of mesh points per coordinate direction (>0). */
  int    mesh[3];
  /** offset of the first mesh point (lower left corner) from the
      coordinate origin ([0,1[). */
  double mesh_off[3];
  /** charge assignment order ([0,7]). */
  int    cao;
  /** number of interpolation points for charge assignment function */
  int    inter;
  /** Accuracy of the actual parameter set. */
  double accuracy;

  /** epsilon of the "surrounding dielectric". */
  double epsilon;
  /** Cutoff for charge assignment. */
  double cao_cut[3];
  /** mesh constant. */
  double a[3];
  /** inverse mesh constant. */
  double ai[3];
  /** unscaled \ref alpha_L for use with fast inline functions only */
  double alpha;
  /** unscaled \ref r_cut_iL for use with fast inline functions only */
  double r_cut;
  /** full size of the interpolated assignment function */
  int inter2;
  /** number of points unto which a single charge is interpolated, i.e. p3m.cao^3 */
  int cao3;
  /** additional points around the charge assignment mesh, for method like dielectric ELC
      creating virtual charges. */
  double additional_mesh[3];
} p3m_parameter_struct;

/** initialize the parameter struct */
void p3m_common_parameter_pre_init(p3m_parameter_struct *params);

/** print local mesh content. 
    \param l local mesh structure.
*/
void p3m_p3m_print_local_mesh(p3m_local_mesh l);

/** print send mesh content. 
 *  \param sm send mesh structure.
*/
void p3m_p3m_print_send_mesh(p3m_send_mesh sm);

/** Add values of a 3d-grid input block (size[3]) to values of 3d-grid
 *  ouput array with dimension dim[3] at start position start[3].  
 *
 *  \param in          Pointer to first element of input block data.
 *  \param out         Pointer to first element of output grid.
 *  \param start       Start position of block in output grid.
 *  \param size        Dimensions of the block
 *  \param dim         Dimensions of the output grid.
*/
void p3m_add_block(double *in, double *out, int start[3], int size[3], int dim[3]);

/** One of the aliasing sums used by \ref p3m_k_space_error. 
    (fortunately the one which is most important (because it converges
    most slowly, since it is not damped exponentially)) can be
    calculated analytically. The result (which depends on the order of
    the spline interpolation) can be written as an even trigonometric
    polynomial. The results are tabulated here (The employed formula
    is Eqn. 7.66 in the book of Hockney and Eastwood). */
double p3m_analytic_cotangent_sum(int n, double mesh_i, int cao);

/** Computes the  assignment function of for the \a i'th degree
    at value \a x. */
double p3m_caf(int i, double x,int cao_value);

#endif /* P3M || DP3M */

#endif /* _P3M_COMMON_H */
