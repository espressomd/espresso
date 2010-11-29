
/** @File Header file for the correlation classe
 *
 * TODO: Add a return value to all voids to check if everything was working allright.
 * TODO: Add a destructor.
 *
 */

#ifndef STATISTICS_CORRALTION_H
#define STATISTICS_CORRALTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tcl.h"
#include "parser.h"

#define MAXLINELENGTH 2048

#define MIN(a,b) ((a)>(b)?(b):(a))

/** The TCL command parser 
 */
int parse_correlation(Tcl_Interp* interp, int argc, char** argv); 
int correlation_parse_corr(Tcl_Interp* interp, int no, int argc, char** argv);
int correlation_print_usage(Tcl_Interp* interp);
int parse_observable(Tcl_Interp* interp, int argc, char** argv, int* change, int (**A_fun)  ( void* A_args, double* A, unsigned int dim_A), int* dim_A, void** A_args);
int parse_corr_operation(Tcl_Interp* interp, int argc, char** argv, int* change, int (**corr_fun)( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ), int* dim_corr, int dim_A, int dim_B);

  

/** The struct that is used to calculate correlations *
 * 
 * Data organization:
 * We use a ring-like way to manage the data: at the beginning we have a linear array
 * which we fill from index 0 to tau_lin. The index newest[i] always indicates the latest
 * entry of the hierarchic "past" For every new entry in is incremented and if tau_lin is reached, 
 * it starts again from the beginning. 
 */
typedef struct {
  unsigned int hierarchy_depth;    // maximum level of data compression
  unsigned int dim_A;              // dimensionality of A
  unsigned int dim_B;
  unsigned int dim_corr;
  unsigned int *n_sweeps;          // number of correlation sweeps at a particular value of tau
  unsigned int *n_vals;            // number of data values already present at a particular value of tau
  unsigned int t;                  // global time in number of frames
  double dt;                       // time interval at which samples arrive
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

  // compressing functions
  int (*compressA)( double* A1, double*A2, double* A_compressed, unsigned int dim_A );
  int (*compressB)( double* B1, double*B2, double* A_compressed, unsigned int dim_B );

  // correlation function
  int (*corr_operation)  ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr );

  // Functions producing observables A and B from the input data
  int (*A_fun)  ( void* A_args, double* A, unsigned int dim_A);
  void* A_args;
  int(*B_fun)  ( void* B_args, double* B, unsigned int dim_B);
  void* B_args;

  int is_from_file;

} double_correlation;


typedef struct {
  FILE* f;
  IntList requested_columns;
  int n_columns;
  char last_line[MAXLINELENGTH];
  int data_left;
} file_data_source;


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
                  void* A_fun, void* A_args, void* B_fun, void* B_args, void* corr_operation, 
                  void* compressA, void* compressB);


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

/** writes the correlation to the TCL console */
int double_correlation_print_correlation( double_correlation* self, Tcl_Interp* interp); 

int file_data_source_init(file_data_source* self, char* filename, IntList* columns);
int file_data_source_readline(void* xargs, double* A, int dim_A); 



int identity ( double* input, unsigned int n_input, double* A, unsigned int dim_A);

int compress_linear( double* A1, double*A2, double* A_compressed, unsigned int dim_A );

int scalar_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr );

int componentwise_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ); 

int square_distance_componentwise ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr );

int particle_velocities(void* typelist, double* A, unsigned int n_A);
int particle_positions(void* typelist, double* A, unsigned int n_A);


#endif
