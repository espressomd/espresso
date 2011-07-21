/*
  Copyright (C) 2010 The ESPResSo project
  
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
/* * @File Header file for the correlation class
 *
 * This module allow to compute correlations (and other two time averages) on
 * the fly an from files.
 *
 * The basic idea is that the user can write arbitrary function A and B that
 * can depend on e.g. particle coordinates or whatever state of the MD box.
 * These function can be vector valued, always indicated by dim_A and dim_B.
 *
 * The way they can be correlated is can be formulated as a (vector valued)
 * operation on A and B. One example would be to calculate the product
 * component by component, and if A and B are both the particle velocities,
 * then one would obtain 
 * { <v_1x v_1x(t)>  <v_1y v_1y(t)>  <v_1z v_1z(t)>  <v_2x v_2x(t)>  <v_2y v_2y(t)>  <v_2z v_2z(t)> ... }
 *
 * The idea of the algorithm is not to keep all As and Bs in memory, but
 * feed the algorithm (or the container object) successively by new As and Bs
 * and new correlation estimates are calculated on the fly. Those As and Bs
 * which have exceeded a certain age are (first compressed, see below) and then
 * discarded. 
 *
 * To save memory, increase statistics and make the calculation possible for
 * many orders of magnitude in time, the blocking algorithm in Frenkels book (p.91)
 * is applied. Thus not all As and Bs of the whole "past" are stored but 
 * some of them are blocked. In this implementation always a blocking based
 * on 2 is applied: All As and Bs not older than a certain tau_lin are stored
 * as they were, those which are older are not entire stored, but only
 * their compressed average value. All As and Bs older than 2*tau_lin
 * are compressed in blocks of four etc.
 *
 * This leads to a hiearchical "history": on level 0 the last tau_lin values
 * are stored. This is done in a cyclic array: The newest is always appended
 * at the end of the array, and if the array length is reached, values
 * are appended at the beginning, overwriting older values. We therefore
 * have to carry the index newest[0] which tells, what is the index of the
 * newest value of A and B. 
 *
 * As soon as the first row of As and Bs is full, two of them are 
 * compressed by calculating the arithmetic mean, and stored in the
 * first element on level 1. Then we can overwrite these two values
 * on level 0, because we have saved them. Always if necessary 
 * we produce space on level 0 by compressing to level 1. On 
 * level 1 also an array with tau_lin entries is available to store
 * the level-1-compresssed values. It is sucessivly filled
 * and also cyclic. When it is filled, the values are stored
 * on level 2 and so on.
 *
 * This allows to have a "history" over many orders of magnitude
 * in time, without the full memory effort. 
 *
 * Correlations are only calculated on each level: For
 * tau=1,2,..,tau_lin the values are taken from level 1.
 * For tau=tau_lin, tau_lin+2, .., 2*tau_lin we take the values
 * from level 2. On level 2 we halso have values for 0,..tau_lin-2, 
 * but these are discared as we have "better" estimates on level 1.
 *
 * The functions A and B can get extra arguments. This can e.g. be an
 * integer describing to the "type" of the particles to be considered,
 * or in case of files, it is a file_data_source pointer, which tells
 * the function from where to read. These arguments have to be a single
 * pointer to an object that carries all relevant information
 * to obtain A and B (from the MD Box or somewhere else).
 *
 * The correlation has to be initialized with all necessary information,
 * i.e. all function pointers, the dimensions of A and B and their dimensions etc.
 * When new A and B values should be processed, a single call of
 * double_correlation_get_data(..) with the correlation object as an 
 * argument is enough to update A and B and the correlation estimate.
 *
 * Eventually the correlation can be printed. 
 *
 *
 * There is a lot of stuff to do:
 * ==============================
 * -> Write input functions that take TCL arrays for A and B to
 *  make the method available for TCL-coded variables
 * -> Expand the file_data_source so that one can specify which
 *  columns of the file are to be processed
 * -> Write more correlation operations (scalar product)
 * -> Write more observable
 * -> calculate an estimate of average values. This might be 
 *  even necessary to calculate <(A-<A>)(B(tau)-<B>), which
 *  is often probably what people want
 * -> Use the A_args to calculate As and Bs only for particular 
 *  particle types (especially and example, so that other people can follow)
 * -> Use the A_args to calculate molecular stuff in combination with
 *  the topology concept
 * -> Tidy up the parsers (just a bit)
 * -> Write some nice samples 
 * -> Write a destructor
 * -> Finally: write the users guide
 *
 *
 */
#ifndef STATISTICS_CORRELATION_H
#define STATISTICS_CORRELATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tcl.h"
#include "parser.h"
#include "statistics_observable.h"

#define MAXLINELENGTH 2048

#define MIN(a,b) ((a)>(b)?(b):(a))

// IDs of different correlations
#define CORR_TYPE_GENERIC       0
#define CORR_TYPE_SF            1
#define CORR_TYPE_AVERAGED_SF   2
#define CORR_TYPE_VACF          3
#define CORR_TYPE_MSD           4


void autoupdate_correlations();

/** The struct that is used to calculate correlations *
 * 
 * Data organization:
 * We use a ring-like way to manage the data: at the beginning we have a linear array
 * which we fill from index 0 to tau_lin. The index newest[i] always indicates the latest
 * entry of the hierarchic "past" For every new entry in is incremented and if tau_lin is reached, 
 * it starts again from the beginning. 
 */
typedef struct {
  unsigned int autocorrelation;    // autocorrelation flag
  unsigned int correlation_type;   // to help keeping track what we are actually correlating
  unsigned int finalized;          // non-zero of correlation is finialized
  unsigned int hierarchy_depth;    // maximum level of data compression
  unsigned int dim_A;              // dimensionality of A
  unsigned int dim_B;
  unsigned int dim_corr;
  unsigned int *n_sweeps;          // number of correlation sweeps at a particular value of tau
  unsigned int *n_vals;            // number of data values already present at a particular value of tau
  unsigned int t;                  // global time in number of frames
  double dt;                       // time interval at which samples arrive
  int update_frequency;              // time distance between updates in MD timesteps 
  double tau_max;                  // maximum time, for which the correlation should be calculated
  unsigned int tau_lin;            // number of frames in the linear correlation
  unsigned int* newest;            // index of the newest entry in each hierarchy level
  unsigned int window_distance; 

  // Convenience pointers to our stored data
  // indices: A[level][tau_i][component]
  int* tau;                       // time differences
  double*** A;                     // input quantity 1
  double*** B;                     // input quantity 2
  double** result;                // output quantity
  unsigned int n_result;           // the total number of result values
  
  // The actual allocated storage space
  double* A_data;
  double* B_data;
  double* result_data;
  int* tau_data;                     // just for double-checking, store tau for all results

  double* A_accumulated_average;     // all A values are added up here
  double* B_accumulated_average;     // all B values are added up here
  double* A_accumulated_variance;    // all A**2 values are added up here
  double* B_accumulated_variance;    // all B**2 values are added up here
  unsigned int n_data;               // a counter to calculated averages and variances

  // compressing functions
  int (*compressA)( double* A1, double*A2, double* A_compressed, unsigned int dim_A );
  int (*compressB)( double* B1, double*B2, double* A_compressed, unsigned int dim_B );

  // correlation function
  int (*corr_operation)  ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr );

  // Functions producing observables A and B from the input data
  observable* A_obs;
  observable* B_obs;

  int is_from_file;
  int autoupdate;
  double last_update;

} double_correlation;

extern int correlations_autoupdate;
extern double_correlation* correlations;




/** This struct allow to use a file as input for the correlation.
 * It is initalized and then just passed as extra argument to the file_data_source_readline 
 * function that extracts floating point values from the particular file.
 *
 */
typedef struct {
  FILE* f;
  IntList requested_columns;
  int n_columns;
  char last_line[MAXLINELENGTH];
  int data_left;
} file_data_source;

/** The TCL command parser 
 */
int tclcommand_correlation(ClientData data, Tcl_Interp* interp, int argc, char** argv);
int tclcommand_correlation_parse_corr(Tcl_Interp* interp, int no, int argc, char** argv);

int correlation_print_usage(Tcl_Interp* interp);
// parsing calls to pre-defined correlations
int parse_structure_factor (Tcl_Interp* interp, int argc, char** argv, int*  change, void** A_args, int *tau_lin_p, double *tau_max_p, double* delta_t_p);
// TESTING
//void print_sf_params(sf_params *params);
// parsing generic correlation call
int tclcommand_correlation_parse_observable(Tcl_Interp* interp, int argc, char** argv, observable** obs);
int parse_corr_operation(Tcl_Interp* interp, int argc, char** argv, int* change, int (**corr_fun)( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ), int* dim_corr, int dim_A, int dim_B);

  

/**
 * The initialization procedure for the correlation object. All important parameters have to be speciefied
 * at the same time. They can not be change later, so every instance of the correlation class
 * has to be fed with correct data from the very beginning.
 *
 * @param self: The pointer to the correlation class instance that should be initialized
 * @param dt: The time interval between subsequent updates
 * @param @param tau_lin: The linear part of the correlation function. 
 * @param hierarchy_depth: The depth of the hierarchichal memory
 * @param window_distance: The distance in time domain between update of the correlation estimate
 * @param dim_A: The dimension of the A vector
 * @param dim_B: The dimension of the B vector
 * @param dim_corr: The dimension of the correlation to be calculated
 * @param A_fun: A pointer to the function that is used to calculate the value of A
 * @param A_args: The arguments passed to function A
 * @param B_fun: A pointer to the function that is used to calculate the value of B
 * @param B_args: The arguments passed to function B
 * @param compressA: A pointer to the function that indicates how the A values should be compressed (usually 
 *     the linear compression method)
 * @param compressB: A pointer to the function that indicates how the B values should be compressed (usually 
 *     the linear compression method)
 *
 */
int double_correlation_init(double_correlation* self, double dt, unsigned int tau_lin, unsigned int hierarchy_depth, 
                  unsigned int window_distance, unsigned int dim_A, unsigned int dim_B, unsigned int dim_corr, 
                  observable* A, observable* B, void* corr_operation, 
                  void* compressA, void* compressB,
		  int correlation_type, int autocorrelation);


/** The function to process a new datapoint of A and B
 *  
 * First the function finds out if it necessary to make some space for the new entries of A and B.
 * Then, if necessary, it compresses old Values of A and B to make for the new value. Finally
 * The new values of A and B are stored in A[newest[0]] and B[newest[0]], where the newest indices
 * have been increased before. Finally the correlation estimate is updated. TODO: Not all
 * the correlation estimates have to be updated.
 *
 */
int double_correlation_get_data(  double_correlation* self );

/** At the end of data collection, go through the whole hierarchy and correlate data left there
 *  
 * This works pretty much the same as get_data, but does not feed on new data, just uses what
 * is already available.
 *
 */
int double_correlation_get_data(double_correlation* self);

/** writes the correlation to the TCL console */
int double_correlation_print_correlation( double_correlation* self, Tcl_Interp* interp); 
int double_correlation_write_to_file( double_correlation* self, char* filename); 

/** Prints spherically averaged SF */
int double_correlation_print_spherically_averaged_sf(double_correlation* self, Tcl_Interp* interp);

int file_data_source_init(file_data_source* self, char* filename, IntList* columns);
int file_data_source_readline(void* xargs, double* A, int dim_A); 

int identity ( double* input, unsigned int n_input, double* A, unsigned int dim_A);

/* *************************
*
* Functions for compressing data
*
**************************/
/** The minimal version of compression function */
int compress_do_nothing( double* A1, double*A2, double* A_compressed, unsigned int dim_A );

/** Compress computing arithmetic mean: A_compressed=(A1+A2)/2 */
int compress_linear( double* A1, double*A2, double* A_compressed, unsigned int dim_A );

/** Compress discarding the 1st argument and return the 2nd */
int compress_discard1( double* A1, double*A2, double* A_compressed, unsigned int dim_A );

/** Compress discarding the 2nd argument and return the 1st */
int compress_discard2( double* A1, double*A2, double* A_compressed, unsigned int dim_A );

/* *************************
*
* Functions for computing observables
*
**************************/
/** The minimal version of observable computing function */

int scalar_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr );

int componentwise_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ); 

int complex_conjugate_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ); 

int square_distance_componentwise ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr );

int square_distance( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr );

int square_distance_cond_1d_int ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr );


/**************  Functions that calculate A and B from MD state ************/

/** Obtain the particle velocities.
 */ 
//
//typedef struct {
//  Tcl_Interp* interp;
//  int argc;
//  char** argv;
//} tcl_input_data;
//
//int tcl_input(void* data, double* A, unsigned int n_A);
#endif
