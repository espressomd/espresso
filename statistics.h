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
  /** Status flag for energy calculation. 0 re-initialize energy
      struct, 1 every thing is fine, calculation can start. */
  int init_status;

  /** Array for energies on each node. */
  DoubleList node;
  /** Array for energies summed over all nodes. */
  DoubleList sum;

  /** number of energies. */
  int n;
  /** number of energies before specific interaction energies. */
  int n_pre;
  /** number of energies for bonded interactions. */
  int n_bonded;
  /** number of energies for non-bonded interactions. */
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
extern Energy_stat virials;
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

/** calculate the end-to-end-distance averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged end-to-end-distance */
double calc_re_av();

/** calculate the radius of gyration. chain information \ref chain_start etc. must be set!
    @return the radius of gyration */
double calc_rg();

/** calculate the radius of gyration averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged radius of gyration */
double calc_rg_av();

/** calculate the hydrodynamic radius. chain information \ref chain_start etc. must be set!
    @return the hydrodynamic radius */
double calc_rh();

/** calculate the hydrodynamic radius averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged hydrodynamic radius */
double calc_rh_av();

/** calculates the internal distances within a chain. Chain information \ref chain_start etc. must be set!
    @param idf contains <tt>idf[0],...,idf[chain_length-1]</tt> */
void calc_internal_dist(double **idf);

/** calculates the internal distances within a chain averaged over all configurations stored in \ref #configs.
    Chain information \ref chain_start etc. must be set!
    @param idf contains <tt>idf[0],...,idf[chain_length-1]</tt> */
void calc_internal_dist_av(double **idf);

/** calculate g123. chain information \ref chain_start etc. must be set!
    @param g1 contains g1
    @param g2 contains g2
    @param g3 contains g3
*/
void calc_g123(double *g1, double *g2, double *g3);

/** calculate <g1> averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @param g1 contains <tt>g1[0],...,g1[n_configs-1]</tt>
*/
void calc_g1_av(double **g1);

/** calculate <g2> averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @param g2 contains <tt>g2[0],...,g2[n_configs-1]</tt>
*/
void calc_g2_av(double **g2);

/** calculate <g3> averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @param g3 contains <tt>g3[0],...,g3[n_configs-1]</tt>
*/
void calc_g3_av(double **g3);

/** set the start configuration for g123.
    chain information \ref chain_start etc. must be set!
*/
void init_g123();

/** appends particles' positions in 'partCfg' to \ref #configs */
void analyze_append();

/** appends the configuration stored in 'config[3*count]' to \ref #configs
    @param config the configuration which should be added 
    @param count  how many particles in 'config' */
void analyze_configs(double *config, int count);

/** removes configs[0], pushes all entries forward, appends current 'partCfg' to last spot */
void analyze_push();

/** replaces configs[ind] with current 'partCfg'
    @param ind the entry in \ref #configs to be replaced */
void analyze_replace(int ind);

/** removes configs[ind] and shrinks the array accordingly
    @param ind the entry in \ref #configs to be removed */
void analyze_remove(int ind);

/** Calculates the distribution of particles around others. 
    Calculates the distance distribution of particles with types given
    in the p1_types list around particles with types given in the
    p2_types list. The distances range from r_min to r_max, binned
    into r_bins bins which are either aequidistant (log_flag==0)or
    logarithmically aequidistant (log_flag==1). The result is stored
    in the array dist.
    @param p1_types list with types of particles to find the distribution for.
    @param n_p1     length of p1_types.
    @param p2_types list with types of particles the others are distributed around.
    @param n_p2     length of p2_types.
    @param r_min    Minimal distance for the distribution.
    @param r_max    Maximal distance for the distribution.
    @param r_bins   Number of bins.
    @param log_flag Wether the bins are (logarithmically) aequidistant.
    @param low      particles closer than r_min
    @param dist     Array to store the result (size: r_bins).
 */
void calc_part_distribution(int *p1_types, int n_p1, int *p2_types, int n_p2, 
			    double r_min, double r_max, int r_bins, int log_flag,
			    double *low, double *dist);
/** Calculates the radial distribution function.

    Calculates the radial distribution function of particles with
    types given in the p1_types list around particles with types given
    in the p2_types list. The range is given by r_min and r_max and
    the distribution function is binned into r_bin bins, which are
    equidistant. The result is stored in the array rdf.

    @param p1_types list with types of particles to find the distribution for.
    @param n_p1     length of p1_types.
    @param p2_types list with types of particles the others are distributed around.
    @param n_p2     length of p2_types.
    @param r_min    Minimal distance for the distribution.
    @param r_max    Maximal distance for the distribution.
    @param r_bins   Number of bins.
    @param rdf     Array to store the result (size: r_bins).
*/
void calc_rdf(int *p1_types, int n_p1, int *p2_types, int n_p2, 
	      double r_min, double r_max, int r_bins, double *rdf);

/** Initializes extern Energy_stat \ref #virials to be used by \ref calc_virials. */
void init_virials();

/** Calculates the virials of the system in parallel (hence it should be called by \ref mpi_gather_stats with job=2).<BR>
    Due to the nature of a virial being <tt>Sum(i=0..n_total_particles)(Sum(j=i+1..n_total_particles)(r_ij*F_ij))</tt>
    this function is based on a merge of \ref force_calc into \ref calc_energy. */
void calc_virials();

/** Calculates the pressure in the system from a virial expansion using the terms from \ref calc_virials.<BR>
    Output is stored in the \ref #virials array, in which (on the first node) each component carries the corresponding pressure,
    while <tt>virials.sum.e[0]</tt> contains the total pressure, <tt>virials.node.e[0]</tt> the sum of all squared pressure components,
    <tt>virials.sum.e[1]</tt> the pressure of the ideal gas, and <tt>virials.node.e[1]</tt> the kinetic energy.
*/
void calc_pressure(void);

/** Derives the spherically averaged formfactor S(q) = 1/chain_length * Sum(i,j=1..chain_length)[sin(q*r_ij)/q*r_ij] of a single chain,
    averaged over all \ref chain_n_chains currently allocated (-> chain information must be set!).
    @param qmin  smallest q-vector to look at (qmin > 0)
    @param qmax  biggest q-vector to look at (qmax > qmin)
    @param qbins decides how many S(q) are derived (note that the qbins+1 values will be logarithmically spaced)
    @param _ff   contains S(q) as an array of size qbins */
void analyze_formfactor(double qmin, double qmax, int qbins, double **_ff);

/** Derives the spherically averaged formfactor S(q) = 1/chain_length * Sum(i,j=1..chain_length)[sin(q*r_ij)/q*r_ij] of a single chain,
    averaged over all \ref chain_n_chains of all \ref n_configs stored configurations in \ref #configs.
    @param qmin  smallest q-vector to look at (qmin > 0)
    @param qmax  biggest q-vector to look at (qmax > qmin)
    @param qbins decides how many S(q) are derived (note that the qbins+1 values will be logarithmically spaced)
    @param _ff   contains S(q) as an array of size qbins */
void analyze_formfactor_av(double qmin, double qmax, int qbins, double **_ff);

/*@}*/

#endif
