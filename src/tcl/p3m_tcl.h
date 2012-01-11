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
#ifndef _P3M_TCL_H 
#define _P3M_TCL_H
/** \file p3m.h P3M algorithm for long range coulomb interaction.
 *
 *  We use a P3M (Particle-Particle Particle-Mesh) method based on the
 *  Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm.
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

#include "p3m-common.h"
#include "interaction_data.h"

#ifdef P3M

/************************************************
 * data types
 ************************************************/

/** P3M parameters. */
extern p3m_data_struct p3m;


/** \name Exported Functions */
/************************************************************/
/*@{*/

/// parse the basic p3m parameters
int tclcommand_inter_coulomb_parse_p3m(Tcl_Interp * interp, int argc, char ** argv);

/// parse the optimization parameters of p3m and the tuner
int tclcommand_inter_coulomb_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv);

/** Tune P3M parameters to desired accuracy.

    Usage:
    \verbatim inter coulomb <bjerrum> p3m tune accuracy <value> [r_cut <value> mesh <value> cao <value>] \endverbatim

    The parameters are tuned to obtain the desired accuracy in best
    time, by running mpi_integrate(0) for several parameter sets.

    The function utilizes the analytic expression of the error estimate 
    for the P3M method in the book of Hockney and Eastwood (Eqn. 8.23) in 
    order to obtain the rms error in the force for a system of N randomly 
    distributed particles in a cubic box.
    For the real space error the estimate of Kolafa/Perram is used. 

    Parameter range if not given explicit values: For \ref p3m_struct::r_cut_iL
    the function uses the values (\ref min_local_box_l -\ref #skin) /
    (n * \ref box_l), n being an integer (this implies the assumption that \ref
    p3m_struct::r_cut_iL is the largest cutoff in the system!). For \ref
    p3m_struct::mesh the function uses the two values which matches best the
    equation: number of mesh point = number of charged particles. For
    \ref p3m_struct::cao the function considers all possible values.

    For each setting \ref p3m_struct::alpha_L is calculated assuming that the
    error contributions of real and reciprocal space should be equal.

    After checking if the total error fulfils the accuracy goal the
    time needed for one force calculation (including verlet list
    update) is measured via \ref mpi_integrate (0).

    The function returns a log of the performed tuning.

    The function is based on routines of the program HE_Q.c written by M. Deserno.
 */
int tclcommand_inter_coulomb_p3m_print_tune_parameters(Tcl_Interp *interp);

/// print the p3m parameters to the interpreters result
int tclprint_to_result_p3m(Tcl_Interp *interp);

/** assign the physical charges using the tabulated charge assignment function.
    If store_ca_frac is true, then the charge fractions are buffered in cur_ca_fmp and
    cur_ca_frac. */
// extern void p3m_charge_assign();

/** assign a single charge into the current charge grid. cp_cnt gives the a running index,
    which may be smaller than 0, in which case the charge is assumed to be virtual and is not
    stored in the ca_frac arrays. */
//extern void p3m_assign_charge(double q,
//		       double real_pos[3],
//		       int cp_cnt);
//
/** shrink wrap the charge grid */
//extern void p3m_shrink_wrap_charge_grid(int n_charges);
/** Clean up P3M memory allocations. */
//extern void p3m_free();

#endif /* of ifdef P3M */

#endif  /*of ifndef P3M_H */
