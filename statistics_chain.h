// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef STATISTICS_CHAIN_H
#define STATISTICS_CHAIN_H
/** \file statistics_chain.h 

    This file contains the code for statistics on the data using the
    molecule information set with analyse set chains.

    <b>Responsible:</b>
    <a href="mailto:mann@mpip-mainz.mpg.de">BAM</a>
*/

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

/** calculate the hydrodynamic radius. chain information \ref chain_start etc. must be set!
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
    @param g1 contains <tt>g1[0],...,g1[n_configs-1]</tt>
*/
void calc_g1_av(double **g1);

/** calculate \<g2\> averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @param g2 contains <tt>g2[0],...,g2[n_configs-1]</tt>
*/
void calc_g2_av(double **g2);

/** calculate \<g3\> averaged over all configurations stored in \ref #configs. 
    Chain information \ref chain_start etc. must be set!
    @param g3 contains <tt>g3[0],...,g3[n_configs-1]</tt>
*/
void calc_g3_av(double **g3);

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


///
int print_chain_structure_info(Tcl_Interp *interp);

/** this function scans the arguments for a description of the chain structure,
    i.e. start of chains, number of chains and chain length. Since this structure
    requires the particles to be sorted, this is performed, too. */
int parse_chain_structure_info(Tcl_Interp *interp, int argc, char **argv);

/** sets the particle mol_id according to the chain_structure info*/
void update_mol_ids_setchains();

/** same as \ref parse_chain_structure_info, but also allows for chain_structure info
    not to be given. You either get an error or chain topology is safe to be used. */
int check_and_parse_chain_structure_info(Tcl_Interp *interp, int argc, char **argv);

///
int parse_re(Tcl_Interp *interp, int average, int argc, char **argv);
///
int parse_rg(Tcl_Interp *interp, int average, int argc, char **argv);
///
int parse_rh(Tcl_Interp *interp, int average, int argc, char **argv);
///
int parse_intdist(Tcl_Interp *interp, int average, int argc, char **argv);
///
int parse_bond_l(Tcl_Interp *interp, int average, int argc, char **argv);
///
int parse_bond_dist(Tcl_Interp *interp, int average, int argc, char **argv);
///
int parse_g123(Tcl_Interp *interp, int average, int argc, char **argv);
///
int parse_g_av(Tcl_Interp *interp, int average, int argc, char **argv);
///
int parse_formfactor(Tcl_Interp *interp, int average, int argc, char **argv);

/*@}*/

#endif
