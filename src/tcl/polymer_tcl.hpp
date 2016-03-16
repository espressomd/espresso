/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef POLYMER_TCL_H
#define POLYMER_TCL_H
#include "parser.hpp"

/************************************************************* 
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/

/** Implementation of the tcl-command <br>
    polymer \<N_P\> \<MPC\> \<bond_length\> [start \<part_id\>] [pos \<x\> \<y\> \<z\>] [mode { SAW | RW | PSAW } [\<shield\> [\<max_try\>]]] 
            [charge \<val_cM\>] [distance \<cM_dist\>] [types \<type_nM\> [\<type_cM\>]] [FENE \<type_FENE\>] [angle \<angle\>
	    [\<angle2\> [\<x2\> \<y2\> \<z2\>]]] <br>
    Creates some polymer chains within the simulation box,
    and returns how often the attempt to place a monomer failed in the worst case.
    Parameters: <br>
                 \<N_P\>         = how many polymers to create <br>
                 \<MPC\>         = monomers per chain <br>
                 \<bond_length\> = length of the bonds between two monomers <br>
                 \<part_id\>     = particle number of the start monomer (defaults to '0') <br>
	         \<pos\>         = sets the position of the start monomer of the first chain (defaults to a randomly chosen value) <br>
	         \<mode\>        = selects setup mode: (Pseudo) self avoiding walk ([P]SAW) or plain random walk (RW) (defaults to 'SAW') <br>
	         \<shield\>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0') <br>
	         \<max_try\>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000') <br>
	         \<val_cM\>      = valency of charged monomers (defaults to '0.0') <br>
	         \<cM_dist\>     = distance between two charged monomers' indices (defaults to '1') <br>
	         \<type_{n|c}P\> = type number of {neutral|charged} monomers to be used with "part" (default to '0' and '1') <br>
	         \<type_FENE\>   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0') <br>
		 \<angle\>       = freely rotating bond-angle to be fixed <br>
	         \<angle2\>      = second spherical bond-angle (for setting up helixes or planar polymers) <br>
		 \<x2,y2,z2\>    = sets the position of the 2nd monomer of the first chain <br>
		 \<constr\>      = shall constraints be respected when setting up polymer? (0=no, 1=yes, default: 0)
    <br>For more informations on the parameters see \ref polymerC. */
int tclcommand_polymer (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implementation of the tcl-command <br>
    counterions \<N_CI\> [start \<part_id\>] [mode { SAW | RW } [\<shield\> [\<max_try\>]]] [charge \<val_CI\>] [type \<type_CI\>] <br>
    Creates counterions of charge \<val_CI\> within the simulation box,
    and returns how often the attempt to place a counterion failed in the worst case.
    Parameters: <br>
                   \<N_CI\>        = number of counterions to create <br>
                   \<part_id\>     = particle number of the first counterion (defaults to 'n_total_particles') <br>
		   \<mode\>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW') <br>
		   \<shield\>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0') <br>
		   \<max_try\>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000') <br>
		   \<val_CI\>      = valency of the counterions (defaults to '-1.0') <br>
		   \<type_CI\>     = type number of the counterions to be used with "part" (default to '2') <br>
    <br>For more informations on the parameters see \ref counterionsC. */
int tclcommand_counterions (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implementation of the tcl-command <br>
    salt \<N_pS\> \<N_nS\> [start \<part_id\>] [mode { SAW | RW } [\<shield\> [\<max_try\>]]] [charges \<val_pS\> [\<val_nS\>]]
         [types \<type_pS\> [\<type_nS\>]] [rad \<rad\>]<br>
    Creates \<N_pS\> positively and \<N_nS\> negatively charged salt ions of charge \<val_pS\> and \<val_nS\> within the simulation box,
    and returns how often the attempt to place a salt ion failed in the worst case.
    Parameters: <br>
                   \<N_pS\>/\<N_nS\> = number of salt ions to create <br>
		   \<part_id\>     = particle number of the first salt ion (defaults to 'n_total_particles') <br>
		   \<mode\>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW') <br>
		   \<shield\>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0') <br>
		   \<max_try\>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000') <br>
		   \<val_{p|n}S\>  = valencies of the salt ions (default to '1' and '-1', respectively); if \<val_nS\> is not given, \<val_nS\> = -1*\<val_pS\> <br>
		   \<type_{p|n}S\> = type numbers to be used with "part" (default to '3' and '4'); if \<type_nS\> is not given, \<type_nS\> = \<type_pS\> is assumed. <br>
		   \<rad\>         = radius of the cell for the cell model. <br>
    <br>For more informations on the parameters see \ref saltC. */
int tclcommand_salt (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implementation of the tcl-command <br>
    velocities \<v_max\> [start \<part_id\>] [count \<N_T\>] <br>
    Sets the velocities of \<N_T\> particles to a random value [-vmax,vmax],
    and returns the averaged velocity assigned.
    Parameters: <br>
                   \<v_max\>       = maximum velocity to be used <br>
		   \<part_id\>     = particle number of the first of the \<N_T\> particles (defaults to '0') <br>
		   \<N_T\>         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id') <br>
    <br>For more informations on the parameters see \ref velocitiesC. */
int tclcommand_velocities (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implementation of the tcl-command <br>
    maxwell_velocities [start \<part_id\>] [count \<N_T\>] <br>
    Sets the velocities of \<N_T\> particles to a random value with maxwell distribution,
    and returns the averaged velocity assigned.
    Parameters: <br>
                   \<part_id\>     = particle number of the first of the \<N_T\> particles (defaults to '0') <br>
		   \<N_T\>         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id') <br>
    <br>For more informations on the parameters see \ref maxwell_velocitiesC. */
int tclcommand_maxwell_velocities (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implementation of the tcl-command <br>
    crosslink \<N_P\> \<MPC\> [start \<part_id\>] [catch \<r_catch\>] [distLink \<link_dist\>] [distChain \<chain_dist\>] [FENE \<type_FENE\>] [trials \<max_try\>] <br>
    Evaluates the current configuration and connects each chain's end to a random monomer of another chain at most \<r_catch\> away,
    if the next crosslink from there is at least \<link_dist\> monomers away,
    and returns how many ends are now successfully linked.
    Parameters: <br>
                   \<N_P\>         = number of polymer chains <br>
                   \<MPC\>         = monomers per chain <br>
                   \<part_id\>     = particle number of the start monomer (defaults to '0') <br>
		   \<r_catch\>     = maximum length of a crosslink (defaults to '1.9') <br>
		   \<link_dist\>   = minimum distance between the indices of two crosslinked monomers with different binding partners (defaults to '2') <br>
		   \<chain_dist\>  = same as \<link_dist\>, but for monomers of the same bond (defaults to \<MPC\> =\> no bonds to the same chain allowed) <br>
		   \<type_FENE\>   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0') <br>
		   \<max_try\>     = how often crosslinks should be removed if they are too close to other links (defaults to '30000') <br>
    <br>For more informations on the parameters see \ref crosslinkC. */
int tclcommand_crosslink (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implementation of the tcl-command <br>
    diamond \<a\> \<bond_length\> \<MPC\> [counterions \<N_CI\>] [charges \<val_nodes\> \<val_cM\> \<val_CI\>] [distance \<cM_dist\>] [nonet] */
int tclcommand_diamond (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implementation of the tcl-command <br>
    icosaeder \<a\> \<MPC\> [counterions \<N_CI\>] [charges \<val_cM\> \<val_CI\>] [distance \<cM_dist\>] */
int tclcommand_icosaeder (ClientData data, Tcl_Interp *interp, int argc, char **argv);

#endif

