// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef POLYMER_H
#define POLYMER_H
/** \file polymer.h

    This file contains everything needed to create a start-up
    configuration of (partially charged) polymer chains with
    counterions and salt molecules, assigning velocities to the
    particles and crosslinking the polymers if necessary.
 
    For more information on polymer, see \ref polymer.c "polymer.c"
  
    <b>Responsible:</b>
    <a href="mailto:mann@mpip-mainz.mpg.de">BAM</a>
*/

#include <tcl.h>
#include "particle_data.h"

/************************************************************* 
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/


/** C implementation of 'mindist \<posx\> \<posy\> \<posz\>',<br>
    which returns the minimum distance of all current particles
    to position (\<posx\>, \<posy\>, \<posz\>) as a double.<br>
    If it fails, return value equals -1. */
double mindist4(double pos[3]);

/** C implementation of 'mindist \<part_id\> \<r_catch\>',<br>
    which returns the size of an array \<ids\> of indices of particles which are 
    less than \<r_catch\> away from the position of the particle \<part_id\>. */
int mindist3(int part_id, double r_catch, int *ids);

/** Checks whether a particle at coordinates (\<posx\>, \<posy\>, \<posz\>) collides
    with any other particle due to a minimum image distance smaller than \<shield\>. 
    @return Returns '1' if there is a collision, '0' otherwise. */
int collision(double pos[3], double shield);

/** Implementation of the tcl-command <br>
    polymer \<N_P\> \<MPC\> \<bond_length\> [start \<part_id\>] [pos \<x\> \<y\> \<z\>] [mode { SAW | RW } [\<shield\> [\<max_try\>]]] 
            [charge \<val_cM\>] [distance \<cM_dist\>] [types \<type_nM\> [\<type_cM\>]] [FENE \<type_FENE\>] [angle \<angle\>
	    [\<angle2\>]] <br>
    Creates some polymer chains within the simulation box,
    and returns how often the attempt to place a monomer failed in the worst case.
    Parameters: <br>
                 \<N_P\>         = how many polymers to create <br>
                 \<MPC\>         = monomers per chain <br>
                 \<bond_length\> = length of the bonds between two monomers <br>
                 \<part_id\>     = particle number of the start monomer (defaults to '0') <br>
	         \<pos\>         = sets the position of the start monomer of the first chain (defaults to a randomly chosen value) <br>
	         \<mode\>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW') <br>
	         \<shield\>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0') <br>
	         \<max_try\>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000') <br>
	         \<val_cM\>      = valency of charged monomers (defaults to '0.0') <br>
	         \<cM_dist\>     = distance between two charged monomers' indices (defaults to '1') <br>
	         \<type_{n|c}P\> = type number of {neutral|charged} monomers to be used with "part" (default to '0' and '1') <br>
	         \<type_FENE\>   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0') <br>
		 \<angle\>       = desired bond-angle to be fixed <br>
	         \<angle2\>      = second spherical bond-angle <br>
    <br>For more informations on the parameters see \ref polymerC. */
int polymer (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'polymer \<N_P\> \<MPC\> \<bond_length\> [options]', which returns how often the attempt to place a monomer failed in the worst case.
    @param  N_P         = how many polymers to create <br>
    @param  MPC         = monomers per chain <br>
    @param  bond_length = length of the bonds between two monomers <br>
    @param  part_id     = particle number of the start monomer (defaults to '0') <br>
    @param  posed       = sets the position of the start monomer of the first chain (defaults to a randomly chosen value) <br>
    @param  mode        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW') <br>
    @param  shield      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0') <br>
    @param  max_try     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000') <br>
    @param  val_cM      = valency of charged monomers (defaults to '0.0') <br>
    @param  cM_dist     = distance between two charged monomers' indices (defaults to '1') <br>
    @param  type_nM     = type number of neutral monomers (defaults to '0') <br>
    @param  type_cM     = type number of charged monomers (default to '1') <br>
    @param  type_FENE   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0') <br>
    @param  angle       = desired bond-angle to be fixed <br>
    @param  angle2      = secon spherical bond-angle (relative to the plane s*(preceding bondvectore)+t*[0,0,1])<br>
    @return Returns how often the attempt to place a monomer failed in the worst case. <br>
    If val_cM \< 1e-10, the charge is assumed to be zero, and type_cM = type_nM.  */
int polymerC(int N_P, int MPC, double bond_length, int part_id, double *posed, int mode, double shield, int max_try, 
	     double val_cM, int cM_dist, int type_nM, int type_cM, int type_FENE, double angle, double angle2);

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
int counterions (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'counterions \<N_CI\> [options]'.
    @param  N_CI        = number of counterions to create
    @param  part_id     = particle number of the first counterion (defaults to 'n_total_particles')
    @param  mode        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
    @param  shield      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0')
    @param  max_try     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
    @param  val_CI      = valency of the counterions (defaults to '-1.0')
    @param  type_CI     = type number of the counterions to be used with "part" (default to '2')
    @return Returns how often the attempt to place a particle failed in the worst case. */
int counterionsC(int N_CI, int part_id, int mode, double shield, int max_try, double val_CI, int type_CI);


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
int salt (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'salt \<N_pS\> \<N_nS\> [options]',
    @param  N_pS        = number of positive salt ions to create
    @param  N_nS        = number of negative salt ions to create
    @param  part_id     = particle number of the first salt ion (defaults to 'n_total_particles')
    @param  mode        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
    @param  shield      = shield around each particle another particle's position may not enter if using SAW (defaults to '0')
    @param  max_try     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
    @param  val_pS      = valency of the positive salt ions (defaults to '1')
    @param  val_nS      = valencies of the negative salt ions (default to '-1')
    @param  type_pS     = type numbers to be used for positive salt ions (defaults to '3')
    @param  type_nS     = type numbers to be used for negative salt ions (defaults to '4')
    @param  rad         = radius of the sphere for the cell model. If value is >0 the ions are distributed in a sphere centered
                          in the simulation box. 
    @return Returns how often the attempt to place a particle failed in the worst case. */
int saltC(int N_pS, int N_nS, int part_id, int mode, double shield, int max_try, double val_pS, double val_nS, int type_pS, int type_nS, double rad);

/** Implementation of the tcl-command <br>
    velocities \<v_max\> [start \<part_id\>] [count \<N_T\>] <br>
    Sets the velocities of \<N_T\> particles to a random value [-vmax,vmax],
    and returns the averaged velocity assigned.
    Parameters: <br>
                   \<v_max\>       = maximum velocity to be used <br>
		   \<part_id\>     = particle number of the first of the \<N_T\> particles (defaults to '0') <br>
		   \<N_T\>         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id') <br>
    <br>For more informations on the parameters see \ref velocitiesC. */
int velocities (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'velocities \<v_max\> [options]',
    @param  v_max       = maximum velocity to be used
    @param  part_id     = particle number of the first of the \<N_T\> particles (defaults to '0') 
    @param  N_T         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id')
    @return Returns the averaged velocity when done. */
double velocitiesC(double v_max, int part_id, int N_T);

/** Implementation of the tcl-command <br>
    maxwell_velocities [start \<part_id\>] [count \<N_T\>] <br>
    Sets the velocities of \<N_T\> particles to a random value with maxwell distribution,
    and returns the averaged velocity assigned.
    Parameters: <br>
                   \<part_id\>     = particle number of the first of the \<N_T\> particles (defaults to '0') <br>
		   \<N_T\>         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id') <br>
    <br>For more informations on the parameters see \ref maxwell_velocitiesC. */
int maxwell_velocities (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'maxwell_velocities [options]',
    @param  part_id     = particle number of the first of the \<N_T\> particles (defaults to '0') 
    @param  N_T         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id')
    @return Returns the averaged velocity when done. */
double maxwell_velocitiesC(int part_id, int N_T);

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
int crosslink (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Collects the bonds leading to and from the ending monomers of the chains (mode == 1) or
    all the bonds leading to and from each monomer (mode == 2). 
    @return  Returns '0' upon success, '-2' otherwise. */
int collectBonds(int mode, int part_id, int N_P, int MPC, int type_bond, int **bond_out, int ***bonds_out);

/** C implementation of 'crosslink \<N_P\> \<MPC\> [options]',
    @param  N_P         = number of polymer chains
    @param  MPC         = monomers per chain
    @param  part_id     = particle number of the start monomer (defaults to '0')
    @param  r_catch     = maximum length of a crosslink (defaults to '1.9')
    @param  link_dist   = minimum distance between the indices of two crosslinked monomers with different binding partners (defaults to '2')
    @param  chain_dist  = same as \<link_dist\>, but for monomers of the same bond (defaults to \<MPC\> =\> no bonds to the same chain allowed)
    @param  type_FENE   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0') 
    @param  max_try     = how often crosslinks should be removed if they are too close to other links (defaults to '30000')
    @return Returns how many ends are now successfully linked. */
int crosslinkC(int N_P, int MPC, int part_id, double r_catch, int link_dist, int chain_dist, int type_FENE, int max_try);

/** Implementation of the tcl-command <br>
    diamond \<a\> \<bond_length\> \<MPC\> [counterions \<N_CI\>] [charges \<val_nodes\> \<val_cM\> \<val_CI\>] [distance \<cM_dist\>] [nonet] */
int diamond (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'diamond \<a\> \<bond_length\> \<MPC\> [options]' */
int diamondC(double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, int nonet);

/** Implementation of the tcl-command <br>
    icosaeder \<a\> \<MPC\> [counterions \<N_CI\>] [charges \<val_cM\> \<val_CI\>] [distance \<cM_dist\>] */
int icosaeder (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'icosaeder \<a\> \<bond_length\> \<MPC\> [options]' */
int icosaederC(double ico_a, int MPC, int N_CI, double val_cM, double val_CI, int cM_dist);


#endif

