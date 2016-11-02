/*
 Copyright (C) 2010,2011 Florian Fahrenberger
 Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

/** \file maggs.hpp
 *  Maxwell Equations Molecular Dynamics (MEMD) method for electrostatic
 *  interactions.
 *
 *  We use a local update scheme to propagate artificial B-fields on a
 *  lattice in the system. In principal, the algorithm simulates full
 *  electrodynamics, but with a tunable speed of light.
 *
 *  The method is very usable for large particle numbers or highly
 *  parallel architectures, since it is local and scales linearly.
 *  It is not suited for high-precision calculation of forces, since
 *  the simple interpolation scheme produces errors in the order of
 *  10^-5 in the force.
 *
 *  The chosen mesh should roughly be of the size of the particles.
 *
 *  Further reading on the algorithm:
 *  <ul>
 *  <li> I. Pasichnyk and B. Dunweg, Coulomb interaction via local dynamics: a molecular-dynamics algorithm. J. Phys: Condens. Matter, 16 ,p. 1399-4020, (2004).
 *  </ul>
 *  
 */


/** protect header file: */
#ifndef _MAGGS_TCL_H
#define _MAGGS_TCL_H

#include "parser.hpp"

#ifdef ELECTROSTATICS

/*****************************/
/** \name External functions */
/*****************************/

/*@{*/

/** parce TCL command for dielectrics. Checks the parameters and calls the according
 functions.
 @return 0 for success, -1 otherwise
 @param interp  TCL Interpreter handle
 @param argc    number of TCL arguments after "inter coulomb $bjerrum maggs"
 @param argv    array of TCL arguments after "inter coulomb $bjerrum maggs"
 */
int tclcommand_localeps(Tcl_Interp* interp, int argc, char** argv);

/** parce TCL command for dielectrics. Checks the parameters and calls the according
 functions.
 @return 0 for success, -1 otherwise
 @param interp  TCL Interpreter handle
 @param argc    number of TCL arguments after "inter coulomb $bjerrum maggs"
 @param argv    array of TCL arguments after "inter coulomb $bjerrum maggs"
 */
int tclcommand_adaptive_eps(Tcl_Interp* interp, int argc, char** argv);

/** parse TCL command. The number of parameters is checked and
    maggs_set_parameters function is called.
    @return 0 for success, -1 otherwise
    @param interp  TCL Interpreter handle
    @param argc    number of TCL arguments after "inter coulomb $bjerrum maggs"
    @param argv    array of TCL arguments after "inter coulomb $bjerrum maggs"
*/
int tclcommand_inter_coulomb_parse_maggs(Tcl_Interp * interp, int argc, char ** argv);

/** Print out the results.
    @return 0 for success, -1 otherwise
    @param interp TCL Interpreter handle
*/
int tclprint_to_result_Maggs(Tcl_Interp *interp);

/*@}*/

#endif
#endif
