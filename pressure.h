// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef PRESSURE_H
#define PRESSURE_H
#include "statistics.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat virials;
extern Observable_stat p_tensor;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

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

/** implementation of 'analyze pressure' */
int parse_and_print_pressure(Tcl_Interp *interp, int argc, char **argv);

/** Does the binning for calc_p_tensor
    @param r_min minimum distance for binning
    @param r_max maximum distance for binning
    @param r_bins number of bins
    @param center 3 dim pointer to sphere origin 
*/
void calc_bins_sphere(int *_new_bin,int *_elements,double *_volumes,double r_min,double r_max,int r_bins, double *center);

/** Implementation of 'analyze bins' */
int parse_bins(Tcl_Interp *interp, int argc, char **argv);

/** Initializes extern Energy_stat \ref #p_tensor to be used by \ref calc_p_tensor. */
void init_p_tensor();

/** Calculates the pressure in the system from a virial expansion as a tensor.<BR>
    Output is stored in the \ref #p_tensor array, in which the <tt>p_tensor.sum</tt>-components contain the sum of the component-tensors
    stored in <tt>p_tensor.node</tt>. The function is executed on the master node only and uses sort of a N^2-loop to calculate the virials,
    so it is rather slow.
    @param volume the volume of the bin to be considered
    @param p_list contains the list of particles to look at
    @param flag   decides whether to derive the interactions of the particles in p_list to <i>all</i> other particles (=1) or not (=0).
*/
void calc_p_tensor(double volume, IntList *p_list, int flag);

/** implementation of 'analyze p_IK1' */
int parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
