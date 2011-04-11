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
#ifndef _P3M_H 
#define _P3M_H
/** \file p3m.h   main header-file for P3M algorithms intended to deal with long-range forces.
 *
 *  We use a P3M (Particle-Particle Particle-Mesh) method based on the
 *  Ewald summation. Details of the used method can be found in
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
 */

#include "utils.h"
#include "interaction_data.h"
#include "integrate.h"

#ifdef ELP3M

/* defined only within this header file, for checking that system (electro/magnetostatic)
   specific files are only included from here. */
#define P3M_H_CURRENT

/** This value for p3m.epsilon indicates metallic boundary conditions. */
#define P3M_EPSILON_METALLIC 0.0

/************************************************
 * data types
 ************************************************/

/** Structure to hold P3M parameters and some dependend variables. */
typedef struct {
#ifdef ELECTROSTATICS
  /** Ewald splitting parameter (0<alpha<1), rescaled to alpha_L = alpha * box_l. */
  double alpha_L;
  /** Cutoff radius for real space electrostatics (>0), rescaled to r_cut_iL = r_cut * box_l_i. */
  double r_cut_iL;
  /** number of mesh points per coordinate direction (>0). */
  int    mesh[3];
  /** offset of the first mesh point (lower left 
      corner) from the coordinate origin ([0,1[). */
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
#endif  
#ifdef MAGNETOSTATICS 
    /** Ewald splitting parameter (0<alpha<1), rescaled to alpha_L = alpha * box_l. */
  double Dalpha_L;
  /** Cutoff radius for real space electrostatics (>0), rescaled to r_cut_iL = r_cut * box_l_i. */
  double Dr_cut_iL;
  /** number of mesh points per coordinate direction (>0). */
  int    Dmesh[3];
  /** offset of the first mesh point (lower left 
      corner) from the coordinate origin ([0,1[). */
  double Dmesh_off[3];
  /** charge assignment order ([0,7]). */
  int    Dcao;
  /** number of interpolation points for charge assignment function */
  int    Dinter;
  /** Accuracy of the actual parameter set. */
  double Daccuracy;

  /** epsilon of the "surrounding dielectric". */
  double Depsilon;
  /** Cutoff for charge assignment. */
  double Dcao_cut[3];
  /** mesh constant. */
  double Da[3];
  /** inverse mesh constant. */
  double Dai[3];
  /** unscaled \ref alpha_L for use with fast inline functions only */
  double Dalpha;
  /** unscaled \ref r_cut_iL for use with fast inline functions only */
  double Dr_cut;
  /** full size of the interpolated assignment function */
  int Dinter2;
  /** number of points unto which a single charge is interpolated, i.e. p3m.Dcao^3 */
  int Dcao3;
  /** additional points around the charge assignment mesh, for method like dielectric ELC
      creating virtual charges. */
  double Dadditional_mesh[3];
#endif
} p3m_struct;


/** P3M parameters. */
extern p3m_struct p3m;

/** \name Exported Functions */
/************************************************************/
/*@{*/

/// print the p3m parameters to the interpreters result
int tclprint_to_result_ChargeP3M(Tcl_Interp *interp);
int tclprint_to_result_DipolarP3M(Tcl_Interp *interp);

/** Calculate the k-space contribution to the coulomb interaction
     forces for electro- and/or magnetostatic interactions, if
    compiled in. */ 
double P3M_calc_kspace_forces(int force_flag, int energy_flag);

/** Clean up P3M memory allocations. */
void P3M_exit();

/*@}*/


#endif /* of ifdef ELP3M */

#endif  /*of ifndef P3M_H */
