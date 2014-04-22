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
#ifndef POLYMER_H
#define POLYMER_H
/** \file polymer.hpp

    This file contains everything needed to create a start-up
    configuration of (partially charged) polymer chains with
    counterions and salt molecules, assigning velocities to the
    particles and crosslinking the polymers if necessary.
 
    For more information on polymer, see \ref polymer.cpp "polymer.c"
*/

#include "particle_data.hpp"

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
    @param pos coordinates of particle to check
    @param shield minimum distance before it is defined as collision
    @param n_add number of additional coordinates to check
    @param add additional coordinates to check
    @return Returns '1' if there is a collision, '0' otherwise. */
int collision(double pos[3], double shield, int n_add, double *add);

/** Function used by polymerC to determine wether a constraint has been violated while setting up a polymer. Currently only "wall", "sphere" and "cylinder" constraints are respected.
    @param p1           = position of first particle given as double-array of lenght 3
    @param p2           = position of second particle given as double-array of length 3
    @return Returns 1 if p1 and p2 sit on opposite sites of any constraint currently defined in the system and 0 otherwise
 */
int constraint_collision(double *p1,double *p2);


/** C implementation of 'polymer \<N_P\> \<MPC\> \<bond_length\> [options]', which returns how often the attempt to place a monomer failed in the worst case.
    @param  N_P         = how many polymers to create <br>
    @param  MPC         = monomers per chain <br>
    @param  bond_length = length of the bonds between two monomers <br>
    @param  part_id     = particle number of the start monomer (defaults to '0') <br>
    @param  posed       = sets the position of the start monomer of the first chain (defaults to a randomly chosen value) <br>
    @param  mode        = selects setup mode: (Pseudo) self avoiding walk ([P]SAW) or plain random walk (RW) (defaults to 'SAW') <br>
    @param  shield      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0') <br>
    @param  max_try     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000') <br>
    @param  val_cM      = valency of charged monomers (defaults to '0.0') <br>
    @param  cM_dist     = distance between two charged monomers' indices (defaults to '1') <br>
    @param  type_nM     = type number of neutral monomers (defaults to '0') <br>
    @param  type_cM     = type number of charged monomers (default to '1') <br>
    @param  type_FENE   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0') <br>
    @param  angle       = desired bond-angle to be fixed <br>
    @param  angle2      = secon spherical bond-angle <br>
    @param  posed2      = sets the position of the 2nd monomer of the first chain <br>
    @param  constr      = shall constraints be respected when setting up polymer?  (0=no, 1=yes, default: 0)
    @return Returns how often the attempt to place a monomer failed in the worst case. <br>
    If val_cM \< 1e-10, the charge is assumed to be zero, and type_cM = type_nM.  */
int polymerC(int N_P, int MPC, double bond_length, int part_id, double *posed, int mode, double shield, int max_try, 
	     double val_cM, int cM_dist, int type_nM, int type_cM, int type_FENE, double angle, double angle2, double* posed2, int constr);

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

/** C implementation of 'velocities \<v_max\> [options]',
    @param  v_max       = maximum velocity to be used
    @param  part_id     = particle number of the first of the \<N_T\> particles (defaults to '0') 
    @param  N_T         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id')
    @return Returns the averaged velocity when done. */
double velocitiesC(double v_max, int part_id, int N_T);

/** C implementation of 'maxwell_velocities [options]',
    @param  part_id     = particle number of the first of the \<N_T\> particles (defaults to '0') 
    @param  N_T         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id')
    @return Returns the averaged velocity when done. */
double maxwell_velocitiesC(int part_id, int N_T);

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

/** C implementation of 'diamond \<a\> \<bond_length\> \<MPC\> [options]' */
int diamondC(double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, int nonet);

/** C implementation of 'icosaeder \<a\> \<bond_length\> \<MPC\> [options]' */
int icosaederC(double ico_a, int MPC, int N_CI, double val_cM, double val_CI, int cM_dist);

#endif

