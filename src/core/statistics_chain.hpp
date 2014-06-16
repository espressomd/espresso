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
#ifndef STATISTICS_CHAIN_H
#define STATISTICS_CHAIN_H
/** \file statistics_chain.hpp 

    This file contains the code for statistics on the data using the
    molecule information set with analyse set chains.
*/

/** \name Exported Variables */
/************************************************************/
/** Particles' initial positions (needed for g1(t), g2(t), g3(t) in \ref tclcommand_analyze) */
/*@{*/
extern float *partCoord_g;
extern float *partCM_g;
extern int n_part_g;
extern int n_chains_g;
/*@}*/

/** data for a system consisting of chains. TBRS. */
/*@{*/
extern int chain_start;
extern int chain_n_chains;
extern int chain_length;
/*@}*/


/** \name Exported Functions */
/************************************************************/
/*@{*/

/** calculate the end-to-end-distance. chain information \ref chain_start etc. must be set!
    @return the end-to-end-distance */
void calc_re(double **re);

/** calculate the end-to-end-distance averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged end-to-end-distance */
void calc_re_av(double **re);

/** calculate the radius of gyration. chain information \ref chain_start etc. must be set!
    @return the radius of gyration */
void calc_rg(double **rg);

/** calculate the radius of gyration averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged radius of gyration */
void calc_rg_av(double **rg);

/** calculate the hydrodynamic radius (ref. Kirkwood-Zimm theory). chain information \ref chain_start etc. must be set!
    @return the hydrodynamic radius */
void calc_rh(double **rh);

/** calculate the hydrodynamic radius averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @return the averaged hydrodynamic radius */
void calc_rh_av(double **rh);

/** calculates the internal distances within a chain. Chain information \ref chain_start etc. must be set!
    @param idf contains <tt>idf[0],...,idf[chain_length-1]</tt> */
void calc_internal_dist(double **idf);

/** calculates the internal distances within a chain averaged over all configurations stored in \ref #configs.
    Chain information \ref chain_start etc. must be set!
    @param idf contains <tt>idf[0],...,idf[chain_length-1]</tt> */
void calc_internal_dist_av(double **idf);

/** calculates the bond length between two neighbouring monomers (i.e. idf[1] in \ref calc_internal_dist).
    Chain information \ref chain_start etc. must be set!
    @param bond_l returns the bond length */
void calc_bond_l(double **bond_l);

/** calculates the averaged bond length between two neighbouring monomers (i.e. idf[1] in \ref calc_internal_dist_av).
    Chain information \ref chain_start etc. must be set!
    @param bond_l returns the bond length */
void calc_bond_l_av(double **bond_l);

/** calculates the internal distances within a chain measured from monomer \<ind_n\>.
    Chain information \ref chain_start etc. must be set!
    @param bdf   contains <tt>bdf[0],...,bdf[(chain_length-1) - ind_n]</tt> 
    @param ind_n the index of the monomer from where all distances are taken */
void calc_bond_dist(double **bdf, int ind_n);

/** calculates the internal distances within a chain measured from monomer \<ind_n\> averaged over all configurations stored in \ref #configs.
    Chain information \ref chain_start etc. must be set!
    @param bdf contains <tt>bdf[0],...,bdf[(chain_length-1) - ind_n]</tt> 
    @param ind_n the index of the monomer from where all distances are taken */
void calc_bond_dist_av(double **bdf, int ind_n);

/** calculate g123. chain information \ref chain_start etc. must be set!
    @param g1 contains g1
    @param g2 contains g2
    @param g3 contains g3
*/
void calc_g123(double *g1, double *g2, double *g3);

/** calculate \<g1\> averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @param _g1     contains <tt>g1[0],...,g1[n_configs-1]</tt>
    @param window  if large than 0, the window size for a sliding window analysis
    @param weights weights for the different coordinates, basically to allow to calculate 2d g1
*/
void calc_g1_av(double **_g1, int window, double weights[3]);

/** calculate \<g2\> averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @param _g2 contains <tt>g2[0],...,g2[n_configs-1]</tt>
    @param window  if large than 0, the window size for a sliding window analysis
    @param weights weights for the different coordinates, basically to allow to calculate 2d g1
*/
void calc_g2_av(double **_g2, int window, double weights[3]);

/** calculate \<g3\> averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @param _g3 contains <tt>g3[0],...,g3[n_configs-1]</tt>
    @param window  if large than 0, the window size for a sliding window analysis
    @param weights weights for the different coordinates, basically to allow to calculate 2d g1
*/
void calc_g3_av(double **_g3, int window, double weights[3]);
//void calc_g3_av(double **g3);

/** set the start configuration for g123.
    chain information \ref chain_start etc. must be set!
*/
void init_g123();

/** Derives the spherically averaged formfactor S(q) = 1/chain_length * Sum(i,j=1..chain_length)[sin(q*r_ij)/q*r_ij] of a single chain,
    averaged over all \ref chain_n_chains currently allocated (-\> chain information must be set!).
    @param qmin  smallest q-vector to look at (qmin \> 0)
    @param qmax  biggest q-vector to look at (qmax \> qmin)
    @param qbins decides how many S(q) are derived (note that the qbins+1 values will be logarithmically spaced)
    @param _ff   contains S(q) as an array of size qbins */
void analyze_formfactor(double qmin, double qmax, int qbins, double **_ff);

/** Derives the spherically averaged formfactor S(q) = 1/chain_length * Sum(i,j=1..chain_length)[sin(q*r_ij)/q*r_ij] of a single chain,
    averaged over all \ref chain_n_chains of all \ref n_configs stored configurations in \ref #configs.
    @param qmin  smallest q-vector to look at (qmin \> 0)
    @param qmax  biggest q-vector to look at (qmax \> qmin)
    @param qbins decides how many S(q) are derived (note that the qbins+1 values will be logarithmically spaced)
    @param _ff   contains S(q) as an array of size qbins */
void analyze_formfactor_av(double qmin, double qmax, int qbins, double **_ff);

/** Calculates monomer-monomer distribution between monomers of different chains. 
    @param r_min   minimal distance for the distribution.
    @param r_max   maximal distance for the distribution.
    @param r_bins  the number of bins
    @param _rdf    contains the monomer-monomer distribution 
    @param _rdf_cm contains the distribution of centers of mass of the chains
    @param _rdf_d  contains the distribution of closest distances between the chains
    */
void analyze_rdfchain(double r_min, double r_max, int r_bins, double **_rdf, double **_rdf_cm, double **_rdf_d);

#ifdef ELECTROSTATICS
/** Calculates the (charge weighted) velocity auto-correlation function from the stored configurations.
 *  The charge weighted velocity auto-correlation function is used to determine
 *  the electrophoretic mobility of a chain using Green-Kubo relation.
 * 
 *  cwvac(tau) = < sum_i^N ( q_i*v_i(t0)*v_CM(t0+tau) ) >
  @param maxtau maximal tau
  @param interval step between t0, sampling frequency
  @param _avac contains the averaged(=over all chains in the system) velocity auto-correlation
  @param _evac contains the error associated with the averaged velocity auto-correlation function
  */
void analyze_cwvac(int maxtau, int interval, double **_avac, double **_evac); 
#endif

/** sets the particle mol_id according to the chain_structure info*/
void update_mol_ids_setchains();


#endif
