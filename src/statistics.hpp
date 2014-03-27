/*
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
#ifndef _STATISTICS_H
#define _STATISTICS_H
/** \file statistics.hpp
    This file contains the code for statistics on the data.
*/

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "utils.hpp"
#include "topology.hpp"

/** \name Data Types */
/************************************************************/
/*@{*/

typedef struct {
  /** Status flag for observable calculation.  For 'analyze energy': 0
      re-initialize observable struct, else everything is fine,
      calculation can start.  For 'analyze pressure' and 'analyze
      p_inst': 0 or !(1+v_comp) re-initialize, else all OK. */
  int init_status;

  /** Array for observables on each node. */
  DoubleList data;

  /** number of coulomb interactions */
  int n_coulomb;
  /** number of coulomb interactions */
  int n_dipolar;
  /** number of non bonded interactions */
  int n_non_bonded;
  /** Number of virtual sites relative (rigid body) conributions */
  int n_vs_relative;

  /** start of bonded interactions. Right after the special ones */
  double *bonded;
  /** start of observables for non-bonded interactions. */
  double *non_bonded;
  /** start of observables for coulomb interaction. */
  double *coulomb;
  /** start of observables for coulomb interaction. */
  double *dipolar;
  /** Start of observables for virtual sites relative (rigid bodies) */
  double *vs_relative;

  /** number of doubles per data item */
  int chunk_size;
} Observable_stat;

/** Structure used only in the pressure and stress tensor calculation to distinguish 
    non-bonded intra- and inter- molecular contributions. */
typedef struct {
  /** Status flag for observable calculation.
      For 'analyze energy': 0 re-initialize observable struct, else every thing is fine, calculation can start.
      For 'analyze pressure' and 'analyze p_inst': 0 or !(1+v_comp) re-initialize, else all OK. */
  int init_status_nb;

  /** Array for observables on each node. */
  DoubleList data_nb;

  /** number of non bonded interactions */
  int n_nonbonded;

  /** start of observables for non-bonded intramolecular interactions. */
  double *non_bonded_intra;
  /** start of observables for non-bonded intermolecular interactions. */
  double *non_bonded_inter;

  /** number of doubles per data item */
  int chunk_size_nb;
} Observable_stat_non_bonded;

/*@}*/

/** \name Exported Variables
    Previous particle configurations (needed for offline analysis
    and correlation analysis in \ref tclcommand_analyze)
*/
/************************************************************/
/*@{*/
extern double **configs;
extern int n_configs;
extern int n_part_conf;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** the minimal distance of two particles with types in set1 rsp. set2.
    @param set1 types of particles
    @param set2 types of particles
    @return the minimal distance of two particles */
double mindist(IntList *set1, IntList *set2);

/** calculate the aggregate distribution for molecules.
    @param dist_criteria2 distance criteria squared
    @param min_contact minimum number of contacts 
    @param s_mol_id start molecule id
    @param f_mol_id finish molecule id
    @param head_list 
    @param link_list
    @param agg_id_list
    @param agg_num
    @param agg_size
    @param agg_max
    @param agg_min
    @param agg_avg
    @param agg_std
    @param charge_criteria
*/
int aggregation(double dist_criteria2, int min_contact, int s_mol_id, int f_mol_id, int *head_list, int *link_list,
		int *agg_id_list, int *agg_num, int *agg_size, int *agg_max, int *agg_min, int *agg_avg, int *agg_std, int charge_criteria);

/** returns all particles within a given radius r_catch around a position.
    @param pos position of sphere of point
    @param r_catch the radius around the position
    @param il the list where to store the particles indices
    @param planedims orientation of coordinate system
*/
void nbhood(double pos[3], double r_catch, IntList *il, int planedims[3]);

/** minimal distance to point.
    @param pos point
    @param pid  if a valid particle id, this particle is omitted from minimization
                (this is a good idea if the posx, posy, posz is the position of a particle).
    @return the minimal distance of a particle to coordinates (\<posx\>, \<posy\>, \<posz\>). */
double distto(double pos[3], int pid);

/** numerical solution for the integration constant \f$\gamma\f$ in the cell model, determined by 
    \f[\gamma\,\ln\frac{R}{r_0}=\arctan\frac{1}{\gamma}+\arctan\frac{\xi_M-1}{\gamma}\f]
    from which the second integration constant, the Manning radius \f$R_M\f$, follows to
    \f[R_M = R\cdot\exp\left(-\frac{1}{\gamma}\cdot\arctan\frac{1}{\gamma}\right)\f]
    Any value \f$\xi_M\geq 0\f$ is allowed, the function will automatically ensure the 
    analytical continuation required for \f$\xi_M<\ln(R/r_0)/(1+\ln(R/r_0))\f$, in which case 
    \f$\gamma\f$ becomes imaginary.
    @param xi_m   Manning parameter \f$\xi_M=\ell_B/a\f$ (with Bjerrum-length \f$\ell_B\f$ and charge distance \f$a\f$)
    @param Rc     outer radius \f$R_C\f$ of the cylindrical cell around each polyelectrolyte
    @param ro     inner radius \f$r_0\f$ of the cylindrical cell around each polyelectrolyte
    @param gacc   the accuracy up to which \f$\gamma\f$ should be determined
    @param maxtry maximum number of interations to find a solution 
    @param result pointer to double array containing \f$\gamma\f$ and \f$R_M\f$, 
                  and a third entry which is -1.0 if \f$\gamma\f$ is imaginary, +1.0 else. */
void calc_cell_gpb(double xi_m, double Rc, double ro, double gacc, int maxtry, double *result);

/** appends particles' positions in 'partCfg' to onfigs */
void analyze_append();

/** appends the configuration stored in 'config[3*count]' to configs
    @param config the configuration which should be added 
    @param count  how many particles in 'config' */
void analyze_configs(double *config, int count);

/** Docs missing!
\todo Docs missing
*/
void analyze_activate(int ind);

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


/** Calculates the radial distribution function averaged over last n_conf configurations.

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
    @param rdf      Array to store the result (size: r_bins).
    @param n_conf   Number of configurations from the last stored configuration.
*/
void calc_rdf_av(int *p1_types, int n_p1, int *p2_types, int n_p2,
	      double r_min, double r_max, int r_bins, double *rdf, int n_conf);

/** Calculates the intermolecular radial distribution function averaged over last n_conf configurations.

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
    @param rdf      Array to store the result (size: r_bins).
    @param n_conf   Number of configurations from the last stored configuration.
*/

void calc_rdf_intermol_av(int *p1_types, int n_p1, int *p2_types, int n_p2,
	      double r_min, double r_max, int r_bins, double *rdf, int n_conf);


/** Calculates the van Hove auto correlation function and as a side product the mean sqaure displacement (msd).

    Calculates the van Hove auto correlation function (acf)  G(r,t) which is the probability that a particle has moved
    a distance r after time t. In the case of a random walk G(r,t)/(4 pi r*r) is a gaussian. The mean square 
    displacement (msd) is connected to the van Hove acf via sqrt(msd(t)) = int G(r,t) dr. This is very useful for
    the investigation of diffusion processes.
    calc_vanhove does the calculation for one particle type ptype and stores the functions specified by rmin, rmax and
    rbins in the arrays msd and vanhove.

    @param ptype    particle type for which the analysis should be performed
    @param rmin     minimal distance for G(r,t)
    @param rmax     maximal distance for G(r,t)
    @param rbins    number of bins for the r distribution in G(r,t)
    @param tmax     max time, for which G(r,t) is computed, if omitted or set to zero, default tmax=n_configs-1 is used
    @param msd      array to store the mean square displacement (size n_configs-1)
    @param vanhove  array to store G(r,t) (size (n_configs-1)*(rbins))

*/
double calc_vanhove(int ptype, double rmin, double rmax, int rbins, int tmax, double *msd, double **vanhove);


/** Calculates the spherically averaged structure factor.

    Calculates the spherically averaged structure factor of particles of a
    given type. The possible wave vectors are given by q = 2PI/L sqrt(nx^2 + ny^2 + nz^2).
    The S(q) is calculated up to a given length measured in 2PI/L (the recommended order of
    the wave vector is less than 20).
    The data is stored starting with q=1, and contains alternatingly S(q-1) and the number
    of wave vectors l with l^2=q. Only if the second number is nonzero, the first is meaningful.
    This means the q=1 entries are sf[0]=S(1) and sf[1]=1. For q=7, there are no possible wave vectors,
    so sf[2*(7-1)]=sf[2*(7-1)+1]=0.
    
    @param type   the type of the particles to be analyzed
    @param order  the maximum wave vector length in 2PI/L
    @param sf     pointer to hold the base of the array containing the result (size: 2*order^2).
*/

void calc_structurefactor(int type, int order, double **sf);
	  

/** Calculates the density profile in dir direction */
void density_profile_av(int n_conf, int n_bin, double density, int dir, double *rho_ave, int type);

int calc_radial_density_map (int xbins,int ybins,int thetabins,double xrange,double yrange, double axis[3], double center[3], IntList *beadids, DoubleList *density_map, DoubleList *density_profile);

void calc_diffusion_profile(int dir, double xmin, double xmax, int nbins, int n_part, int n_conf, int time, int type, double *bins) ;  

/** returns the minimal squared distance between two positions in the perhaps periodic
    simulation box.
 *  \param pos1  Position one.
 *  \param pos2  Position two.
 */
double min_distance2(double pos1[3], double pos2[3]);

/** returns the minimal distance between two positions in the perhaps periodic
    simulation box.
 *  \param pos1  Position one.
 *  \param pos2  Position two.
 */
inline double min_distance(double pos1[3], double pos2[3]) {
  return sqrt(min_distance2(pos1, pos2));
}

/** calculate the center of mass of a special type of the current configuration
 *  \param type  type of the particle
 *  \param com   center of mass position
 */
void centermass(int type, double *com);

/** Docs missing
\todo Docs missing
*/
void centermass_vel(int type, double *com);

/** calculate the angular momentum of a special type of the current configuration
 *  \param type  type of the particle
 *  \param com   angular momentum vector
 */
void angularmomentum(int type, double *com);


/** calculate the center of mass of a special type of a saved configuration
 *  \param k       number of the saved configuration
 *  \param type_1  type of the particle, -1 for all
 *  \param com     center of mass position
 */
void centermass_conf(int k, int type_1, double *com);


void momentofinertiamatrix(int type, double *MofImatrix);
void calc_gyration_tensor(int type, double **gt);
void calculate_verlet_neighbors();

/** returns the momentum of the particles in the simulation box.
 * \param result Momentum of particles.
 */
void predict_momentum_particles(double *result);

/** Docs missing
\todo Docs missing
*/
void momentum_calc(double *momentum);

inline double *obsstat_bonded(Observable_stat *stat, int j)
{
  return stat->bonded + stat->chunk_size*j;
}

inline double *obsstat_nonbonded(Observable_stat *stat, int p1, int p2)
{
  int tmp;
  if (p1 > p2) {
    tmp = p2;
    p2 = p1;
    p1 = tmp;
  }
  return stat->non_bonded + stat->chunk_size*(((2 * n_particle_types - 1 - p1) * p1) / 2  +  p2);
}

inline double *obsstat_nonbonded_intra(Observable_stat_non_bonded *stat, int p1, int p2)
{
/*  return stat->non_bonded_intra + stat->chunk_size*1; */
  int tmp;
  if (p1 > p2) {
    tmp = p2;
    p2 = p1;
    p1 = tmp;
  }
  return stat->non_bonded_intra + stat->chunk_size_nb*(((2 * n_particle_types - 1 - p1) * p1) / 2  +  p2);
}

inline double *obsstat_nonbonded_inter(Observable_stat_non_bonded *stat, int p1, int p2)
{
/*  return stat->non_bonded_inter + stat->chunk_size*1; */
  int tmp;
  if (p1 > p2) {
    tmp = p2;
    p2 = p1;
    p1 = tmp;
  }
  return stat->non_bonded_inter + stat->chunk_size_nb*(((2 * n_particle_types - 1 - p1) * p1) / 2  +  p2);
}

void invalidate_obs();

/** Docs missing
\todo Docs missing
*/
void mark_neighbours(int type,int pa_nr,double dist,int *list);

void obsstat_realloc_and_clear(Observable_stat *stat, int n_pre, int n_bonded, int n_non_bonded,
			       int n_coulomb, int n_dipolar, int n_vsr, int chunk_size);

void obsstat_realloc_and_clear_non_bonded(Observable_stat_non_bonded *stat_nb, int n_nonbonded, int chunk_size_nb);

/*@}*/

#endif
