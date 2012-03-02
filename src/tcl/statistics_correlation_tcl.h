/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
  
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
#ifndef _STATISTICS_CORRELATION_TCL_H
#define _STATISTICS_CORRELATION_TCL_H

#include "config.h"
#include "tcl.h"
#include "parser.h"
#include "statistics_correlation.h"

#define MAXLINELENGTH 2048

#define MIN(a,b) ((a)>(b)?(b):(a))

/** This struct allow to use a file as input for the correlation.
 * It is initalized and then just passed as extra argument to the file_data_source_readline 
 * function that extracts floating point values from the particular file.
 *
 */
//typedef struct {
//  FILE* f;
//  IntList requested_columns;
//  int n_columns;
//  char last_line[MAXLINELENGTH];
//  int data_left;
//} file_data_source;

/** The TCL command parser 
 */
int tclcommand_correlation(ClientData data, Tcl_Interp* interp, int argc, char** argv);
int tclcommand_correlation_parse_corr(Tcl_Interp* interp, int no, int argc, char** argv);

int correlation_print_usage(Tcl_Interp* interp);
// parsing calls to pre-defined correlations
int parse_structure_factor (Tcl_Interp* interp, int argc, char** argv, int*  change, void** A_args, int *tau_lin_p, double *tau_max_p, double* delta_t_p);
// parsing generic correlation call
int tclcommand_correlation_parse_observable(Tcl_Interp* interp, int argc, char** argv, observable** obs);
//int parse_corr_operation(Tcl_Interp* interp, int argc, char** argv, int* change, int (**corr_fun)( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ), unsigned int* dim_corr, unsigned int dim_A, unsigned int dim_B);
int parse_corr_operation(Tcl_Interp* interp, int argc, char** argv, int* change, char **corr_operation_name, unsigned int* dim_corr, unsigned int dim_A, unsigned int dim_B, void **args);

  
/** writes the correlation to the TCL console */
int double_correlation_print_correlation( double_correlation* self, Tcl_Interp* interp); 
int double_correlation_write_to_file( double_correlation* self, char* filename); 

int file_data_source_init(file_data_source* self, char* filename, IntList* columns);
int file_data_source_readline(void* xargs, double* A, unsigned int dim_A); 


#endif
