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
#ifndef _P3M_TCL_H 
#define _P3M_TCL_H

#include "parser.hpp"

#ifdef P3M

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

    Parameter range if not given explicit values: For \ref p3m_parameter_struct::r_cut_iL
    the function uses the values (\ref min_local_box_l -\ref #skin) /
    (n * \ref box_l), n being an integer (this implies the assumption that \ref
    p3m_parameter_struct::r_cut_iL is the largest cutoff in the system!). For \ref p3m_parameter_struct::mesh
	the function uses the two values which matches best the
    equation: number of mesh point = number of charged particles. For
    \ref p3m_parameter_struct::cao the function considers all possible values.

    For each setting \ref p3m_parameter_struct::alpha_L is calculated assuming that the
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

#endif /* of ifdef P3M */

#endif  /*of ifndef P3M_H */
