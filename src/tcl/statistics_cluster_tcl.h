/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
#ifndef STATISTICS_CLUSTER_TCL_H
#define STATISTICS_CLUSTER_TCL_H
#include "parser.h"

/** Parser for the necklace cluster algorithm

    \verbatim analyze necklace <pearl_treshold> <back_dist> <space_dist> <first> <length> \endverbatim 

    For more information see: Limbach H.J. and Holm C.  Single-Chain
    Properties of Polyelectrolytes in Poor Solvent J. Phys. Chem. B,
    107 (32), 8041 -8055, AUG 14 2003.  The first three parameters are
    tune parameters for the algorithm: pearl_treshold is the minimal
    number of monomers in a pearl. back_dist is the number of monomers
    along the chain backbone which are excluded from the space
    distance criterion to form clusters. space_dist is the distance
    between two monomers up to which they are considered to belong to
    the same clusters. The three parameters may be connected by
    scaling arguments. Make sure that your results are only weakly
    dependent on the exact choice of your parameters. For the
    algorithm the coordinates stored in \ref partCfg are used. The
    chain itself is defined by the identity first of its first monomer
    and the chain length length.
*/
int tclcommand_analyze_parse_necklace(Tcl_Interp *interp, int argc, char **argv);

/** Parser for Hole cluster algorithm
    \verbatim analyze holes <prob_part_type_number> <mesh_size> \endverbatim

    Identifies free space in the simulation box via a mesh based 
    cluster algorithm. Free space is defined via a probe particle 
    and its interactions with other particles which have to be 
    defined through LJ interactions with the other existing particle
    types via the inter command before calling this routine.
*/
int tclcommand_analyze_parse_holes(Tcl_Interp *interp, int argc, char **argv);

#endif
