#ifndef STATISTICS_H
#define STATISTICS_H
/** \file statistics.h
    This file contains the code for simply statistics on the data.

    <b>Responsible:</b>
    <a href="mailto:mann@mpip-mainz.mpg.de">BAM</a>

*/

#include <tcl.h>
#include "particle_data.h"

/** \name Data Types */
/************************************************************/
/*@{*/

typedef struct {
  /** Status flag for energy calculation. 0 re initialize energy
      struct, 1 every thing is fine, calculation can start. */
  int init_status;

  /** Array for energies on each node. */
  DoubleList node;
  /** Array for energies summed over all nodes. */
  DoubleList sum;

  /** number of energies. */
  int n;
  /** number of energies befor specific interaction energies. */
  int n_pre;
  /** number of energies for bonded interactions. */
  int n_bonded;
  /** number of energies for non bonded interactions. */
  int n_non_bonded;
  /** number of energies for coulomb interaction. */
  int n_coulomb;

  /** analyze specified energy. */
  int ana_num;
} Energy_stat;

/*@}*/

/** \name Exported Variables */
/************************************************************/
/*@{*/
extern Energy_stat energy;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Implements the Tcl command 'analyze <what> [<structure info>] [...]' for basic analysis.
    Possible arguments for <what> are:
    <li> 'analyze mindist' \\
         returns the minimal distance of two particles (needs no structure info).
    <li> 'analyze nbhood <part_id> <r_catch>' \\
         returns all particles within a given radius <r_catch> around the position of particle <part_id>.
    <li> 'analyze distto <posx> <posy> <posz>' \\
         returns the minimal distance of a particle to coordinates (<posx>, <posy>, <posz>).
    <li> 'analyze set <structure info>' \\
         defines the structure. The second argument defines the topology to set, i. e. chain at the moment.
	 <ul>
	 Possible values for <structure info> are
	 <li> 'analyze set chains [<chain_start> <n_chains> <chain_length>]'
              A set of equal-length chains. If no parameters are given, the ones currently stored are returned.
	 </ul>
    <li> 'analyze energy [interaction]' \\
         returns the energies of the system. output is blockfile format: \\
	 { energy <value> } { kinetic <value> } { interaction <value> } ... \\
	 if you specify an interaction, e.g. fene <type_num> or lj <type1> <type2> or coulomb or kinetic
	 it returns just that energy.
    All tasks below need the particles to be stored consecutively starting with identity 0.
    <li> 'analyze re [<chain_start> <n_chains> <chain_length>]' \\
         returns the quadratic end-to-end-distance averaged over all polymers (requires chain structure to be set).
    <li> 'analyze rg [<chain_start> <n_chains> <chain_length>]' \
         returns the radius of gyration averaged over all chains (requires chain structure to be set).
    <li> 'analyze rh [<chain_start> <n_chains> <chain_length>]' \\
         returns the hydrodynamic radius (requires chain structure to be set).
    <li> 'analyze g123 [[-init] <chain_start> <n_chains> <chain_length>]' \\
         returns the mean-square displacement g1(t) of a monomer,
                 the mean-square displacement g2(t) in the center of gravity of the chain itself, and
                 the motion of the center of mass g3(t)
	 as a tcl-list {g1(t) g2(t) g3(t)} (requires chain structure to be set). \\
	 If before the structure info you give '-init', the current configuration is stored as reference config.
*/
int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv);

void init_energies();
void calc_energy();
/*@}*/

#endif
