#ifndef POLYMER_H
#define POLYMER_H
/** \file polymer.h
    This file contains everything needed to create a start-up configuration
    of (partially charged) polymer chains with counterions and salt molecules,
    assigning velocities to the particles and crosslinking the polymers if necessary.
 
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


/** C implementation of 'mindist <posx> <posy> <posz>',
    which returns the minimum distance of all current particles
    to position (<posx>, <posy>, <posz>) as a double.
    If it fails, return value equals -1. */
double mindist4(double pos[3]);

/** C implementation of 'mindist <part_id> <r_catch>',
    which returns the size of an array <ids> of indices of particles which are 
    less than <r_catch> away from the position of the particle <part_id>. */
int mindist3(int part_id, double r_catch, int *ids);

/** Checks whether a particle at coordinates (<posx>, <posy>, <posz>) collides
    with any other particle due to a minimum image distance smaller than <shield>. 
    @return Returns '1' if there is a collision, '0' otherwise. */
int collision(double pos[3], double shield);

/** Implementation of the tcl-command
    polymer <N_P> <MPC> <bond_length> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] 
                                      [charge <val_cM>] [distance <cM_dist>] [types <type_nM> [<type_cM>]] [FENE <type_FENE>]
    Creates some polymer chains within the simulation box.
    @param  <N_P>         = how many polymers to create
    @param  <MPC>         = monomers per chain
    @param  <bond_length> = length of the bonds between two monomers
    @param  <box_length>  = length of the simulation box
    @param  <part_id>     = particle number of the start monomer (defaults to '0')
    @param  <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
    @param  <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0')
    @param  <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
    @param  <val_cM>      = valency of charged monomers (defaults to '0.0')
    @param  <cM_dist>     = distance between two charged monomers' indices (defaults to '1')
    @param  <type_{n|c}P> = type number of {neutral|charged} monomers to be used with "part" (default to '0' and '1')
    @param  <type_FENE>   = type number of the FENE-typed bonded interaction bonds to be set between the monomers (defaults to '0')
    @return Returns how often the attempt to place a monomer failed in the worst case.
    If <val_cM> < 1e-10, the charge is assumed to be zero, and <type_cM> = <type_nM>.  */
int polymer (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'polymer <N_P> <MPC> <bond_length> [options]',
    which returns how often the attempt to place a monomer failed in the worst case. */
int polymerC(int N_P, int MPC, double bond_length, int part_id, int mode, double shield, int max_try, 
	     double val_cM, int cM_dist, int type_nM, int type_cM, int type_FENE);

/** Implementation of the tcl-command
    counterions <N_CI> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] [charge <val_CI>] [type <type_CI>]
    Creates counterions of charge <val_CI> within the simulation box.
    @param  <N_CI>        = number of counterions to create
    @param  <part_id>     = particle number of the first counterion (defaults to 'n_total_particles')
    @param  <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
    @param  <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0.0')
    @param  <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
    @param  <val_CI>      = valency of the counterions (defaults to '-1.0')
    @param  <type_CI>     = type number of the counterions to be used with "part" (default to '2')
    @return Returns how often the attempt to place a particle failed in the worst case. */
int counterions (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'counterions <N_CI> [options]',
    which returns how often the attempt to place a counterion failed in the worst case. */
int counterionsC(int N_CI, int part_id, int mode, double shield, int max_try, double val_CI, int type_CI);

/** Implementation of the tcl-command
    salt <N_pS> <N_nS> [start <part_id>] [mode { SAW | RW } [<shield> [<max_try>]]] [charges <val_pS> [<val_nS>]] [types <type_pS> [<type_nS>]]
    Creates <N_pS> positively and <N_nS> negatively charged salt ions of charge <val_pS> and <val_nS> within the simulation box.
    @param  <N_pS>/<N_nS> = number of salt ions to create
    @param  <box_length>  = length of the simulation box
    @param  <part_id>     = particle number of the first salt ion (defaults to '[setmd npart]')
    @param  <mode>        = selects setup mode: Self avoiding walk (SAW) or plain random walk (RW) (defaults to 'SAW')
    @param  <shield>      = shield around each particle another particle's position may not enter if using SAW (defaults to '0')
    @param  <max_try>     = how often a monomer should be reset if current position collides with a previous particle (defaults to '30000')
    @param  <val_{p|n}S>  = valencies of the salt ions (default to '1' and '-1', respectively); if <val_nS> is not given, <val_nS> = -1*<val_pS>
    @param  <type_{p|n}S> = type numbers to be used with "part" (default to '3' and '4'); if <type_nS> is not given, <type_nS> = <type_pS> is assumed. 
    @return Returns how often the attempt to place a particle failed in the worst case. */
int salt (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'salt <N_pS> <N_nS> [options]',
    which returns how often the attempt to place a salt ion failed in the worst case. */
int saltC(int N_pS, int N_nS, int part_id, int mode, double shield, int max_try, double val_pS, double val_nS, int type_pS, int type_nS);

/** Implementation of the tcl-command
    velocities <v_max> [start <part_id>] [count <N_T>]
    Sets the velocities of <N_T> particles to a random value [-vmax,vmax].
    @param  <v_max>       = maximum velocity to be used
    @param  <part_id>     = particle number of the first of the <N_T> particles (defaults to '0') 
    @param  <N_T>         = number of particles of which the velocities should be set (defaults to 'n_total_particles - part_id')
    @return Returns the averaged velocity when done. */
int velocities (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** C implementation of 'velocities <v_max> [options]',
    which returns the averaged velocity assigned. */
double velocitiesC(double v_max, int part_id, int N_T);

/** Implementation of the tcl-command
    crosslink <N_P> <MPC> [start <part_id>] [catch <r_catch>] [distLink <link_dist>] [distChain <chain_dist>] [FENE <type_FENE>] [trials <max_try>]
    Evaluates the current configuration and connects each chain's end to a random monomer of another chain at most <r_catch> away,
    if the next crosslink from there is at least <link_dist> monomers away.
    @param  <N_P>         = number of polymer chains
    @param  <MPC>         = monomers per chain
    @param  <part_id>     = particle number of the start monomer (defaults to '0')
    @param  <r_catch>     = maximum length of a crosslink (defaults to '1.9')
    @param  <link_dist>   = minimum distance between the indices of two crosslinked monomers with different binding partners (defaults to '2')
    @param  <chain_dist>  = same as <link_dist>, but for monomers of the same bond (defaults to <MPC> => no bonds to the same chain allowed)
    @param  <max_try>     = how often crosslinks should be removed if they are too close to other links (defaults to '30000')
    @return Returns how many ends are now successfully linked. */
int crosslink (ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Collects the bonds leading to and from the ending monomers of the chains (mode == 1) or
    all the bonds leading to and from each monomer (mode == 2). 
    @return  Returns '0' upon success, '-2' otherwise. */
int collectBonds(int mode, int part_id, int N_P, int MPC, int type_bond, int **bond_out, int ***bonds_out);

/** C implementation of 'crosslink <N_P> <MPC> [options]',
    which returns how many ends are now successfully linked.   */
int crosslinkC(int N_P, int MPC, int part_id, double r_catch, int link_dist, int chain_dist, int type_FENE, int max_try);


#endif

