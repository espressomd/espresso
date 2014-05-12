/*
  Copyright (C) 2010,2011 Florian Fahrenberger
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

/** \file maggs.cpp
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


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "interaction_data.hpp"
#include "maggs.hpp"
#include "parser.hpp"

#ifdef ELECTROSTATICS

int tclcommand_localeps(Tcl_Interp* interp, int argc, char** argv)
{
    int mesh = maggs_get_mesh_1D();
    int node_x, node_y, node_z, direction;
    double relative_epsilon;
	
    /* number of arguments has to be 8 */
    if(argc != 9) {
        Tcl_AppendResult(interp, "Wrong number of paramters. Usage: \n", (char *) NULL);
        Tcl_AppendResult(interp, "inter coulomb <bjerrum> memd localeps node <x> <y> <z> dir <X/Y/Z> eps <epsilon>", (char *) NULL);
        return TCL_ERROR;
    }
	
    /* first argument should be "node" */
    if(! ARG_IS_S(1, "node")) return TCL_ERROR;
    
    /* arguments 2-4 should be integers */
    if(! ARG_IS_I(2, node_x)) {
        Tcl_AppendResult(interp, "integer expected", (char *) NULL);
        return TCL_ERROR; }
    if(! ARG_IS_I(3, node_y)) {
        Tcl_AppendResult(interp, "integer expected", (char *) NULL);
        return TCL_ERROR; }
    if(! ARG_IS_I(4, node_z)) {
        Tcl_AppendResult(interp, "integer expected", (char *) NULL);
        return TCL_ERROR; }
    /* check if mesh position is in range */
    if ( (node_x < 0) || (node_y < 0) || (node_z < 0) || (node_x > mesh) || (node_y > mesh) || (node_z > mesh) ) {
        char buffer[TCL_INTEGER_SPACE];
        sprintf(buffer, "%d", mesh);
        Tcl_AppendResult(interp, "epsilon position out of mesh range. Mesh in each dimension is ", buffer, ".", (char *) NULL);
        return TCL_ERROR;
    }
    
    /* parse fifth and sixth argument (e.g. dir X) */
    if(! ARG_IS_S(5, "dir")) return TCL_ERROR;
    if ( (! ARG_IS_S(6, "X")) && (! ARG_IS_S(6, "Y")) && (! ARG_IS_S(6, "Z")) ) {
        Tcl_AppendResult(interp, "Parameter dir should be 'X', 'Y' or 'Z'.", (char *) NULL);
        return TCL_ERROR; }
    if(ARG_IS_S(6, "X")) direction = 0;
    if(ARG_IS_S(6, "Y")) direction = 1;
    if(ARG_IS_S(6, "Z")) direction = 2;
    
    /* parse seventh and eight argument (e.g. eps 0.5) */
    if(! ARG_IS_S(7, "eps")) return TCL_ERROR;
    if ( (! ARG_IS_D(8, relative_epsilon)) || (relative_epsilon < 0.0) ) {
        Tcl_AppendResult(interp, "eps expects a positive double", (char *) NULL);
        return TCL_ERROR; }
    
    double eps_before = maggs_set_permittivity(node_x, node_y, node_z, direction, relative_epsilon);
    
    if (eps_before == 1.0) return TCL_OK;
    else return TCL_OK;
}

/** parse TCL command.
    number of parameters is checked and maggs_set_parameters function is called.
    @return zero if successful
    @param interp  TCL interpreter handle
    @param argc    number of arguments given
    @param argv    array of arguments given
*/
int tclcommand_inter_coulomb_parse_maggs(Tcl_Interp * interp, int argc, char ** argv)
{
    int mesh;
    double f_mass;
    double epsilon = 1.0;
    int finite_epsilon_flag = 1;
	
    /* if the command is localeps, call function */
    if ( (argc > 0) && (ARG_IS_S(0, "localeps")) )
        return tclcommand_localeps(interp, argc, argv);
    
    if(argc < 2) {
        Tcl_AppendResult(interp, "Not enough parameters: inter coulomb <bjerrum> memd <f_mass> <mesh>", (char *) NULL);
        return TCL_ERROR;
    }
	
    if(! ARG_IS_D(0, f_mass))
        return TCL_ERROR;
	
    if(! ARG_IS_I(1, mesh)) {
        Tcl_AppendResult(interp, "integer expected", (char *) NULL);
        return TCL_ERROR;
    }
	
    if(argc > 4) {
        Tcl_AppendResult(interp, "Too many parameters: inter coulomb memd <f_mass> <mesh> [epsilon <eps>]", (char *) NULL);
        return TCL_ERROR;
    }
    if(argc == 3) {
        Tcl_AppendResult(interp, "Usage: inter coulomb memd <f_mass> <mesh> [epsilon <eps>]", (char *) NULL);
        return TCL_ERROR;
    }
    if(argc == 4) {
        if (ARG_IS_S(2, "epsilon")) {
            if(! (ARG_IS_D(3, epsilon) && epsilon > 0.0)) {
                Tcl_AppendResult(interp, "epsilon expects a positive double",
                             (char *) NULL);
                return TCL_ERROR;
            }
        }
    } else finite_epsilon_flag=1;

  coulomb.method = COULOMB_MAGGS;
	
  int res = maggs_set_parameters(coulomb.bjerrum, f_mass, mesh,
                                 finite_epsilon_flag, epsilon);
  switch (res) {
  case -1:
    Tcl_AppendResult(interp, "mass of the field is negative", (char *)NULL);
    return TCL_ERROR;
  case -2:
    Tcl_AppendResult(interp, "mesh must be positive", (char *) NULL);
    return TCL_ERROR;
  case ES_OK:
    return TCL_OK;
  }
  Tcl_AppendResult(interp, "unknown error", (char *) NULL);
  return TCL_ERROR;
}


int tclprint_to_result_Maggs(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
	
  Tcl_PrintDouble(interp, maggs.f_mass, buffer);
  Tcl_AppendResult(interp, "maggs ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",maggs.mesh);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL); 
	
  return TCL_OK;
}

#endif // ELECTROSTATICS
