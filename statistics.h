#ifndef STATISTICS_H
#define STATISTICS_H
/** \file statistics.h
    This file contains the code for simply statistics on the data.

    <b>Responsible:</b>
    <a href="mailto:mann@mpip-mainz.mpg.de">BAM</a>

*/

#include <tcl.h>
#include "particle_data.h"
#include "utils.h"

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

/** Implements the Tcl command \ref tcl_analyze. This allows for basic system analysis,
    both online and offline.
*/
int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv);

void init_energies();

void calc_energy();

/** the minimal distance of two particles.
    @return the minimal distance of two particles */
double mindist();

/** returns all particles within a given radius r_catch around a position.
    @param posx x-coordinate of point
    @param posy y-coordinate of point
    @param posz z-coordinate of point
    @param r_catch the radius around the position
    @param il the list where to store the particles indices */
void nbhood(double posx, double posy, double posz, double r_catch, IntList *il);

/** minimal distance to point.
    @param posx x-coordinate of point
    @param posy y-coordinate of point
    @param posz z-coordinate of point
    @param pid  if a valid particle id, this particle is omitted from minimazation
                (this is a good idea if the posx, posy, posz is the position of a particle).
    @return the minimal distance of a particle to coordinates (<posx>, <posy>, <posz>). */
double distto(double posx, double posy, double posz, int pid);

/** calculate the end-to-end-distance. chain information \ref chain_start etc. must be set!
    @return the end-to-end-distance */
double calc_re();

/** calculate the end-to-end-distance averaged over all configurations stored in 'configs'. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged end-to-end-distance */
double calc_re_av();

/** calculate the radius of gyration. chain information \ref chain_start etc. must be set!
    @return the radius of gyration */
double calc_rg();

/** calculate the radius of gyration averaged over all configurations stored in 'configs'. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged radius of gyration */
double calc_rg_av();

/** calculate the hydrodynamic radius. chain information \ref chain_start etc. must be set!
    @return the hydrodynamic radius */
double calc_rh();

/** calculate the hydrodynamic radius averaged over all configurations stored in 'configs'. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged hydrodynamic radius */
double calc_rh_av();

/** calculate g123. chain information \ref chain_start etc. must be set!
    @param g1 contains g1
    @param g2 contains g2
    @param g3 contains g3
*/
void calc_g123(double *g1, double *g2, double *g3);

/** calculate <g1> averaged over all configurations stored in 'configs'. 
    Chain information \ref chain_start etc. must be set!
    @param g1 contains g1[0],...,g1[n_configs-1]
*/
void calc_g1_av(double **g1);

/** calculate <g2> averaged over all configurations stored in 'configs'. 
    Chain information \ref chain_start etc. must be set!
    @param g2 contains g2[0],...,g2[n_configs-1]
*/
void calc_g2_av(double **g2);

/** calculate <g3> averaged over all configurations stored in 'configs'. 
    Chain information \ref chain_start etc. must be set!
    @param g3 contains g3[0],...,g3[n_configs-1]
*/
void calc_g3_av(double **g3);

/** set the start configuration for g123.
    chain information \ref chain_start etc. must be set!
*/
void init_g123();

/** appends particles' positions in 'partCfg' to 'configs' */
void analyze_append();

/** appends the configuration stored in 'config[3*count]' to 'configs'
    @param config the configuration which should be added 
    @param count  how many particles in 'config' */
void analyze_configs(double *config, int count);

/** removes configs[0], pushes all entries forward, appends current 'partCfg' to last spot */
void analyze_push();

/** replaces configs[ind] with current 'partCfg'
    @param ind the entry in 'configs' to be replaced */
void analyze_replace(int ind);

/** removes configs[ind] and shrinks the array accordingly
    @param ind the entry in 'configs' to be removed */
void analyze_remove(int ind);

/*@}*/

#endif
