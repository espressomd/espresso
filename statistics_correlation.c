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
//comentario siguiente
#include "statistics_correlation.h"
#include "particle_data.h"
#include "parser.h"
#include "integrate.h"

double_correlation* correlations=0;
unsigned int n_correlations = 0;

const char init_errors[][64] = {
  "",
  "No valid correlation given" ,
  "delta_t must be > 0",
  "tau_lin must be >= 2",
  "hierarchy_depth must be >=1",
  "window_distance must be >1",
  "dimension of A was not >1",
  "dimension of B was not >1",
  "dim_corr was not >1",
  "no proper function for first observable given",
  "no proper function for second observable given",
  "no proper function for correlation operation given",
  "no proper function for compression of first observable given",
  "no proper function for compression of second observable given"
};

const char file_data_source_init_errors[][64] = {
  "",
  "No valid filename given." ,
  "File could not be opened.",
  "No line found that was not commented out"
};

const char double_correlation_get_data_errors[][64] = {
  "",
  "Error calculating variable A" ,
  "Error calculating variable B" ,
  "Error calculating correlation\n", 
  "Error allocating temporary memory\n", 
  "Error in correlation operation: The vector sizes do not match\n"
};

int correlation_update(unsigned int no) {
  if (n_correlations > no)
    return double_correlation_get_data(&correlations[no]);
  else 
    return 1;
}

int correlation_update_from_file(unsigned int no) {
  if (!correlations[no].is_from_file)
    return 1;
  while ( ! double_correlation_get_data(&correlations[no]) ) {
  }
  return 0;
}


/* We can track several correlation functions at a time
*  identified by their ids.
*/
int tclcommand_analyze_parse_correlation(Tcl_Interp* interp, int argc, char** argv) {
  int no;
  if (argc < 1)
    return correlation_print_usage(interp);
  if (ARG0_IS_I(no)) {
    argc-=1;
    argv+=1;
    return correlation_parse_corr(interp, no, argc, argv);
  } else {
    return correlation_print_usage(interp);
  }
}

int correlation_print_usage(Tcl_Interp* interp) {
  Tcl_AppendResult(interp, "You don't know how to use correlation.", (char *)NULL);
  return TCL_ERROR;
}

int correlation_parse_corr(Tcl_Interp* interp, int no, int argc, char** argv) {
  int (*A_fun)  ( void* A_args, double* A, unsigned int dim_A) = 0;
  int(*compressA)  ( double* A1, double*A2, double* A_compressed, unsigned int dim_A ) = 0;
  void* A_args = 0;
  int dim_A;
  int(*B_fun)  ( void* B_args, double* B, unsigned int dim_B) = 0;
  int(*compressB)  ( double* B1, double*B2, double* B_compressed, unsigned int dim_B ) = 0;
  void* B_args = 0;
  int dim_B;
  int tau_lin = 16; // default values
  /* It would be perhaps better, if the user should specify 
  tau_max and hiererchy depth was computed internally 
  PK */
  int hierarchy_depth = 4; 
  double delta_t = 1;
  double tau_max=0.0;
  int(*corr_operation)  ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) = 0;
  int dim_corr;
  int change; // how many tcl argmuents are "consumed" by the parsing of arguments
  int error;
  tcl_input_data* tcl_input_p;
  char buffer[TCL_INTEGER_SPACE];
  static int *correlation_types=NULL; // list of already set up correlations
  int correlation_type; // correlation type of correlation which is currently being created

  // Check if ID is negative
  if ( no < 0 ) {
    Tcl_AppendResult(interp, "Correlation IDs must be positive", (char *)NULL);
    return TCL_ERROR;
  }
  // If correlation already exists then, the next argument must be print or update, because we
  // can not handle anything else yet.
  if ( no < n_correlations ) {
    if (argc > 0)
      if (ARG0_IS_S("print")) {
        return double_correlation_print_correlation(&correlations[no], interp);
      } else if (ARG0_IS_S("update") && correlations[no].A_fun==&tcl_input && correlations[no].B_fun==&tcl_input ) {
        correlations[no].A_args = tcl_input_p;
        tcl_input_p->interp=interp;
        tcl_input_p->argv=argv+1;
        tcl_input_p->argc=argc-1;
        correlations[no].B_args = tcl_input_p;
        tcl_input_p->interp=interp;
        tcl_input_p->argv=argv+2;
        tcl_input_p->argc=argc-2;
        error = double_correlation_get_data(&correlations[no]);
        if (error) {
          Tcl_AppendResult(interp, "Error reading tclinput", (char *)NULL);
        } else {
          return TCL_OK;
        }
      } else if (ARG0_IS_S("update")) {
        error = double_correlation_get_data(&correlations[no]);
        if (error) {
          Tcl_AppendResult(interp, double_correlation_get_data_errors[error], (char *)NULL);
        } else {
          return TCL_OK;
        }
      } else if (ARG0_IS_S("update_from_file")) {
        error = correlation_update_from_file(no);
        if (error) {
          Tcl_AppendResult(interp, "error in update_from_file", (char *)NULL);
        } else {
          return TCL_OK;
        }
      } else if (ARG0_IS_S("write_to_file")) {
        if (argc <2) {
          Tcl_AppendResult(interp, "You must pass a filename as argument of write_to_file", (char *)NULL);
          return TCL_ERROR;
        }
        error = double_correlation_write_to_file(&correlations[no], argv[1]);
        if (error) {
          Tcl_AppendResult(interp, "error in write to file", (char *)NULL);
        } else {
          return TCL_OK;
        }
      } else {
        Tcl_AppendResult(interp, "Usage for an already existing correlation:", (char *)NULL);
        Tcl_AppendResult(interp, "analyze correlation $ID [ print | update ]", (char *)NULL);
        return TCL_ERROR;
      }
    else
      Tcl_AppendResult(interp, "Usage for an already existing correlation:", (char *)NULL);
      Tcl_AppendResult(interp, "analyze correlation $ID [ print | update ]", (char *)NULL);
      return TCL_ERROR;

  } else if ( no == n_correlations ) {  
  
    // Else we must parse the other arguments and see if we can construct a fully
    // working correlation class instance from that.
    correlation_type=CORR_TYPE_GENERIC; // this is the default
    while (argc > 0) {
      if ( ARG_IS_S_EXACT(0,"vacf") ) {
        Tcl_AppendResult(interp, "Cannot use vacf yet. Use generic instead.", (char *)NULL);
	error=1;
        //error=parse_vacf(interp, argc, argv, &change, &A_fun, &dim_A, &A_args);
        argc -= change;
        argv += change;
        if (error)
          return TCL_ERROR; 
	else 
	  correlation_type=CORR_TYPE_VACF;
      } 
      else if ( ARG_IS_S_EXACT(0,"msd") ) {
        Tcl_AppendResult(interp, "Cannot use msd yet. Use generic instead.", (char *)NULL);
	error=1;
        //error=parse_msd(interp, argc, argv, &change, &A_fun, &dim_A, &A_args);
        argc -= change;
        argv += change;
        if (error)
          return TCL_ERROR;
	else 
	  correlation_type=CORR_TYPE_MSD;
      } 
      else if ( ARG_IS_S_EXACT(0,"structure_factor") ) {
        Tcl_AppendResult(interp, "Cannot use structure_factor yet. Use generic instead.", (char *)NULL);
	error=1;
        error=parse_structure_factor(interp, argc, argv, &change,  &A_args,  &tau_lin, &tau_max,  &delta_t);
        argc -= change;
        argv += change;
        if (error)
          return TCL_ERROR;
	else 
	  correlation_type=CORR_TYPE_SF;
      } 
      else if ( ARG0_IS_S("first_obs") ) {
        argc -= 1;
        argv += 1;
        error=parse_observable(interp, argc, argv, &change, &A_fun, &dim_A, &A_args);
	// by default, set second_observable to the same
        error = parse_observable(interp, argc, argv, &change, &B_fun, &dim_B, &B_args); 
        argc -= change;
        argv += change;
        if (error)
          return TCL_ERROR;
      } else if ( ARG0_IS_S("second_obs") ) {
        argc -= 1;
        argv += 1;
        error = parse_observable(interp, argc, argv, &change, &B_fun, &dim_B, &B_args); 
        argc -= change;
        argv += change;
        if (error)
          return TCL_ERROR;
      } else if ( ARG0_IS_S("corr_operation") ) {
        argc -= 1;
        argv += 1;
        if ( parse_corr_operation(interp, argc, argv, &change, &corr_operation, &dim_corr, dim_A, dim_B) ) 
          return TCL_ERROR;
        argc -= change;
        argv += change;
      } else if ( ARG0_IS_S("tau_lin") ) {
        if ( argc < 2 || !ARG1_IS_I(tau_lin))
          Tcl_AppendResult(interp, "Usage: analyze correlation ... tau_lin $tau_lin", (char *)NULL);
        else { 
          argc -= 2;
          argv += 2;
        } // tau_lin is already set
      } else if ( ARG0_IS_S("hierarchy_depth") ) {
        if ( argc < 2 || !ARG1_IS_I(hierarchy_depth))
          Tcl_AppendResult(interp, "Usage: analyze correlation ... hierarchy_depth $hierarchy_depth", (char *)NULL);
        else { 
          argc -= 2;
          argv += 2;
        } // hierarchy_depth is already set
      } else if ( ARG0_IS_S("delta_t") ) {
        if ( argc < 2 || !ARG1_IS_D(delta_t))
          Tcl_AppendResult(interp, "Usage: analyze correlation ... delta_t $delta_t", (char *)NULL);
        else { 
        argc -= 2;
        argv += 2;
        } // delta_t is already set
      } else if ( ARG_IS_S_EXACT(0,"compress1") ) {
        if ( ARG_IS_S_EXACT(1,"linear") )  compressA=compress_linear; 
        else if (ARG_IS_S_EXACT(1,"discard1")) compressA=compress_discard1;
        else if (ARG_IS_S_EXACT(1,"discard2")) compressA=compress_discard2;
	else {
	  Tcl_AppendResult(interp, "Compression function ", argv[1], (char *)NULL);
	  Tcl_AppendResult(interp, " is not implemented. ", (char *)NULL);
	  return TCL_ERROR;
	}
        argc -= 2;
        argv += 2; 
      } else if ( ARG_IS_S_EXACT(0,"compress2") ) {
        if ( ARG_IS_S_EXACT(1,"linear") )  compressB=compress_linear; 
        else if (ARG_IS_S_EXACT(1,"discard1")) compressB=compress_discard1;
        else if (ARG_IS_S_EXACT(1,"discard2")) compressB=compress_discard2;
	else {
	  Tcl_AppendResult(interp, "Compression function ", argv[1], (char *)NULL);
	  Tcl_AppendResult(interp, " is not implemented. ", (char *)NULL);
	  return TCL_ERROR;
	}
        argc -= 2;
        argv += 2; 
      } else if ( ARG0_IS_S("update") ) {
        sprintf(buffer,"%d",no);
        Tcl_AppendResult(interp, "Correlation error: cannot update correlation ", buffer, (char *)NULL);
        Tcl_AppendResult(interp, " It must be set up first", (char *)NULL);
        return TCL_ERROR;
      } else if ( ARG0_IS_S("print") ) {
        sprintf(buffer,"%d",no);
        Tcl_AppendResult(interp, "Correlation error: cannot print correlation ", buffer, (char *)NULL);
        Tcl_AppendResult(interp, " It must be set up first", (char *)NULL);
        return TCL_ERROR;
      } else {
        Tcl_AppendResult(interp, "unknown argument ", argv[0], (char *)NULL);
        return TCL_ERROR;
      }
    }
  } else {
      sprintf(buffer,"%d",no);
      Tcl_AppendResult(interp,"Correlation error: cannot set up correlation with id ",buffer,(char *)NULL);
      sprintf(buffer,"%d",n_correlations);
      Tcl_AppendResult(interp," Next correlation must have id ",buffer,(char *)NULL);
      return TCL_ERROR;
  }


  // Now we should find out, if this is enough to make a correlation.
  // Unfortunately still nothing happens here!
  if(hierarchy_depth>0) {
    if(compressA==NULL) {
        Tcl_AppendResult(interp, "You have not chosen the compression function! Correlation cannot be performed.", (char *)NULL);
        return TCL_ERROR; 
    }
    if(compressB==NULL) {
        Tcl_AppendResult(interp, "You have not chosen the compressB function. Taking compressB=compressA by default.", (char *)NULL);
	compressB=compressA;
    }
  }
  
  // Let us just realloc space. Does not matter even if we can not make a proper correlation out of that.
  correlations=(double_correlation*) realloc(correlations, (n_correlations+1)*sizeof(double_correlation)); 

  // Now initialize the new correlation
  /*
  error = double_correlation_init(&correlations[n_correlations], delta_t, tau_lin, hierarchy_depth, 1, 
      dim_A, dim_B, dim_corr, A_fun, A_args, B_fun, B_args,
      corr_operation, &compress_linear, &compress_linear);
      */
  error = double_correlation_init(&correlations[n_correlations], delta_t, tau_lin, hierarchy_depth, 1, 
      dim_A, dim_B, dim_corr, A_fun, A_args, B_fun, B_args,
      corr_operation, compressA, compressB);
  if ( error == 0 ) {
    n_correlations++;
    correlation_types=realloc(correlation_types,n_correlations*sizeof(int));
    correlation_types[n_correlations-1]=correlation_type;
    return TCL_OK;
  } else {
    printf("Error number %d\n", error);
    Tcl_AppendResult(interp, "correlation could not be corretly initialized: ", init_errors[error], "\n", (char *)NULL);
    return TCL_ERROR;
  }

}

int observable_usage(Tcl_Interp* interp) {
  Tcl_AppendResult(interp, "Usage: analyze correlation [ velocities | com_velocity | position | com_position ] [ first_obs | second_obs ] " , (char *)NULL);
  Tcl_AppendResult(interp, "   type [ int ... int | * ]  \n" , (char *)NULL);
  Tcl_AppendResult(interp, "or id [ int ... int ]  \n" , (char *)NULL);
  Tcl_AppendResult(interp, "or molecule int \n" , (char *)NULL);
  return TCL_ERROR;
}

int parse_structure_factor (Tcl_Interp* interp, int argc, char** argv, int* change, void** A_args, int *tau_lin_p, double *tau_max_p, double* delta_t_p) {
  sf_params* params;
  int order,order2,tau_lin;
  int i,j,k,l,n;
  double delta_t,tau_max;
  char buffer[TCL_DOUBLE_SPACE];
  if (ARG0_IS_I(order)) {
    if(order>1) {
      params->order=order;
      order2=order*order;
    } else {
      sprintf(buffer, "%d", order);
      Tcl_AppendResult(interp, "order must be > 1, got ", buffer);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, "problem reading order\n", "usage: structure_factor order delta_t tau_max  tau_lin\n");
    return TCL_ERROR; 
  }
  if (ARG1_IS_D(delta_t)) {
    if (delta_t > 0.0) *delta_t_p=delta_t;
    else {
      sprintf(buffer, "%lf", delta_t);
      Tcl_AppendResult(interp, "delta_t must be > 0.0, got ", buffer);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, "problem reading delta_t\n", "usage: structure_factor order delta_t tau_max  tau_lin\n");
    return TCL_ERROR; 
  }
  if (ARG_IS_D(2,tau_max)) {
    if (tau_max > 2.0*delta_t) *tau_max_p=tau_max;
    else {
      sprintf(buffer, "%lf", tau_max);
      Tcl_AppendResult(interp, "tau_max must be > 2.0*delta_t, got ", buffer);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, "problem reading tau_max\n", "usage: structure_factor order delta_t tau_max  tau_lin\n");
    return TCL_ERROR; 
  }
  if (ARG_IS_I(3,tau_lin)) {
    if (tau_lin > 2 && tau_lin < (tau_max/delta_t+1)) *tau_lin_p=tau_lin;
    else {
      sprintf(buffer, "%d", tau_lin);
      Tcl_AppendResult(interp, "tau_lin must be < tau_max/delta_t+1, got ", buffer);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, "problem reading tau_lin\n", "usage: structure_factor order delta_t tau_max  tau_lin\n");
    return TCL_ERROR; 
  }
  // compute the vectors
  l=0;
  for(i=-order; i<=order; i++) 
      for(j=-order; j<=order; j++) 
        for(k=-order; k<=order; k++) {
          n = i*i + j*j + k*k;
          if ((n<=order2) && (n>0)) 
            l += 2;
        }
  params->dim_sf=l;
  *A_args=(void*)params;
  *change=4; // if we reach this point, we have parsed 3 arguments, if not error is returned anyway
  return 0;
}

int parse_id_list(Tcl_Interp* interp, int argc, char** argv, int* change, IntList** ids ) {
  int i;
  char** temp_argv; int temp_argc;
  int temp;
  IntList* input=malloc(sizeof(IntList));
  init_intlist(input);
  alloc_intlist(input,1);


  if (ARG0_IS_S("id")) {
    if (!parse_int_list(interp, argv[1],input)) {
      Tcl_AppendResult(interp, "Error parsing id list\n", (char *)NULL);
      return TCL_ERROR;
    } 
    *ids=input;
    for (i=0; i<input->n; i++) {
      if (input->e[i] >= n_total_particles) {
        Tcl_AppendResult(interp, "Error parsing ID list. Given particle ID exceeds the number of existing particles\n", (char *)NULL);
        return TCL_ERROR;
      }
    }
    *change=2;
    return TCL_OK;

  } else if ( ARG0_IS_S("type") ) {
      Tcl_AppendResult(interp, "Error parsing ID list. Specifying particles by IDs is not possible yet.", (char *)NULL);
      return TCL_ERROR;
  }
  Tcl_AppendResult(interp, "unknown keyword given to observable: ", argv[0] , (char *)NULL);
  return TCL_ERROR;
}

int parse_observable(Tcl_Interp* interp, int argc, char** argv, int* change, int (**A_fun)  ( void* A_args, double* A, unsigned int dim_A), int* dim_A, void** A_args) {


  file_data_source* fds;
  IntList* types=0;
  int error=0;
  int temp;
  int order=0;
  int *order_p;
  IntList* ids;
//  if (ARG0_IS_S("particle_velocities") || ARG0_IS_S("particle_positions") || ARG0_IS_S("structure_factor"))  {
  if (ARG0_IS_S("particle_velocities")) {
    *A_fun = &particle_velocities;
    if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
      return TCL_ERROR;
    *A_args=(void*)ids;
    *dim_A=3*ids->n;
    *change=1+temp;
    return TCL_OK;
  } 
  if (ARG0_IS_S("com_velocity")) {
    *A_fun = &com_velocity;
    if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
      return TCL_ERROR;
    *A_args=(void*)ids;
    *dim_A=3;
    *change=1+temp;
    return TCL_OK;
  } 
  if (ARG0_IS_S("particle_positions")) {
    *A_fun = &particle_positions;
    *A_args=0;
    *dim_A=3*n_total_particles;
    *change=1;
    return TCL_OK;
  }
  if (ARG0_IS_S("structure_factor") ) {
    if (argc > 1 && ARG1_IS_I(order)) {
      *A_fun = &structure_factor;
      order_p=malloc(sizeof(int));
      *order_p=order;
      *A_args=(void*) order_p;
      int order2,i,j,k,l,n ; 
      order2=order*order ;
      l=0;
      // lets counter the number of entries for the DSF
      for(i=-order; i<=order; i++) 
        for(j=-order; j<=order; j++) 
          for(k=-order; k<=order; k++) {
	          n = i*i + j*j + k*k;
	          if ((n<=order2) && (n>=1)) 
              l=l+2;
	  }
    *dim_A=l;
    *change=2;
    return TCL_OK;
    } else { 
      Tcl_AppendResult(interp, "usage: structure_factor $order\n" , (char *)NULL);
      return TCL_ERROR; 
    }
  }
  if (ARG0_IS_S("textfile")) {
    // We still can only handle full files
    if ( argc>1 ) {
      fds= malloc(sizeof(file_data_source));
      error = file_data_source_init(fds, argv[1], 0);
      *change=2;
      if (!error) {
        *A_args=(void*) fds;
        *dim_A = fds->n_columns;
        *A_fun = (void*)&file_data_source_readline;
        return TCL_OK;
      } else {
        Tcl_AppendResult(interp, "Error reading file ", argv[1] ,"\n", (char *)NULL);
        Tcl_AppendResult(interp, file_data_source_init_errors[error] ,"\n", (char *)NULL);
        return TCL_ERROR;
      }
    } else {
      Tcl_AppendResult(interp, "Error in parse_observable textfile: no filename given" , (char *)NULL);
      return TCL_ERROR;
    }
    return TCL_OK;
  } 
  if (ARG0_IS_S("tclinput")) {
    if (argc>1 && ARG1_IS_I(temp)) {
      *dim_A = temp;
      *A_fun = &tcl_input;
      *A_args=malloc(sizeof(tcl_input_data));
      *change =2;
      return TCL_OK;
    } else {
      Tcl_AppendResult(interp, "\nError in parse_observable tclinfo. You must pass the dimension of the observable." , (char *)NULL);
      return TCL_ERROR;
    }
  }else {
    return observable_usage(interp);
  }
  return 0 ;
}


int parse_corr_operation(Tcl_Interp* interp, int argc, char** argv, int* change, int (**corr_fun)( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ), int* dim_corr, int dim_A, int dim_B) {
  if (ARG0_IS_S("componentwise_product")) {
    *corr_fun = &componentwise_product;
    *dim_corr = dim_A;
    *change=1;
    return TCL_OK;
  } else if (ARG0_IS_S("square_distance_componentwise")) {
    *corr_fun = &square_distance_componentwise;
    *dim_corr = dim_A;
    *change=1;
    return TCL_OK;
  } else if (ARG0_IS_S("scalar_product")) {
    *corr_fun = &scalar_product;
    *dim_corr = 1;
    *change=1;
    return TCL_OK;
  } else {
    Tcl_AppendResult(interp, "Unknown correlation operation: ", argv[0], "\n" , (char *)NULL);
    return TCL_ERROR;
  }
}

int double_correlation_init(double_correlation* self, double dt, unsigned int tau_lin, unsigned int hierarchy_depth, 
                  unsigned int window_distance, unsigned int dim_A, unsigned int dim_B, unsigned int dim_corr, 
                  void* A_fun, void* A_args, void* B_fun, void* B_args, void* corr_operation, 
                  void* compressA, void* compressB) {
  unsigned int i,j,k;
  self->dt = dt;
  self->tau_lin=tau_lin;
  self->hierarchy_depth = hierarchy_depth;
  self->dim_A = dim_A;
  self->dim_B = dim_B;
  self->dim_corr = dim_corr;
  self->A_fun = A_fun;
  self->A_args = A_args;
  self->B_fun = B_fun;
  self->B_args = B_args;
  self->corr_operation = corr_operation;
  self->compressA = compressA;
  self->compressB = compressB;
  self->window_distance = window_distance;
  self->t = 0;

  if (self==0)
    return 1;
  if (dt <= 0)
    return 2;
  if (tau_lin<2)
    return 3;
  if (hierarchy_depth<1)
    return 4;
  if (window_distance<1)
    return 5;
  if (dim_A<1)
    return 6;
  if (dim_B<1)
    return 7;
  if (dim_corr<1)
    return 8;
  if (A_fun == 0)
    return 9;
  if (B_fun == 0)
    return 10;
  if (corr_operation==0)
    return 11;
  if (compressA==0)
    return 12;
  if (compressB==0)
    return 13;

  if (A_fun == &file_data_source_readline && B_fun == &file_data_source_readline) {
    self->is_from_file = 1;
  }

  self->A_data = (double*)malloc((tau_lin+1)*hierarchy_depth*dim_A*sizeof(double));
  self->B_data = (double*)malloc((tau_lin+1)*hierarchy_depth*dim_B*sizeof(double));


  self->n_result=tau_lin+1 + (tau_lin+1)/2*(hierarchy_depth-1);
  self->tau = (int*)                malloc(self->n_result*sizeof(int));
  self->n_sweeps = (unsigned int*)  malloc(self->n_result*sizeof(int));
  self->result  = (double**)        malloc(self->n_result*sizeof(double*));
  self->result_data  = (double*)    malloc(self->n_result*dim_corr*sizeof(double));

  self->A = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->B = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->n_vals = (unsigned int*) malloc(hierarchy_depth*sizeof(unsigned int));

  for (i=0; i<self->hierarchy_depth; i++) {
    self->A[i] = (double**) malloc((self->tau_lin+1)*sizeof(double*));
    self->B[i] = (double**) malloc((self->tau_lin+1)*sizeof(double*));
  }
  for (i=0; i<self->hierarchy_depth; i++) {
    self->n_vals[i]=0;
    for (j=0; j<self->tau_lin+1; j++) {
      self->A[i][j] = &self->A_data[(i*(tau_lin+1))*dim_A+j*dim_A];
      self->B[i][j] = &self->B_data[(i*(tau_lin+1))*dim_B+j*dim_B];
      for (k=0; k<dim_A; k++) 
        self->A[i][j][k] = 0.;
      for (k=0; k<dim_B; k++) 
        self->B[i][j][k] = 0.;
    }
  }

  for (i=0; i<self->n_result; i++) {
    self->n_sweeps[i]=0;
    self->result[i]=&self->result_data[i*self->dim_corr];
    for (j=0; j<self->dim_corr; j++) 
      self->result[i][j]=0;
  }

  self->newest = (unsigned int *)malloc(hierarchy_depth*sizeof(unsigned int));
  for ( i = 0; i<self->hierarchy_depth; i++ ) {
    self->newest[i]= self->tau_lin;
  }
  for (i=0; i < tau_lin+1; i++) {
    self->tau[i] = i;
  }
  for (j=1; j < self->hierarchy_depth; j++)
    for (k=0; k < self->tau_lin/2; k++) {
      self->tau[self->tau_lin + 1 + (j-1)*tau_lin/2+k] = (k+(self->tau_lin/2)+1)*(1<<j); 
    }
  return 0;
}

int double_correlation_get_data( double_correlation* self ) {
  // We must now go through the hierarchy and make sure there is space for the new 
  // datapoint. For every hierarchy level we have to decide if it necessary to move 
  // something
  int i,j,k;
  int highest_level_to_compress;
  unsigned int index_new, index_old, index_res;
  int error;
  
  self->t++;

  highest_level_to_compress=-1;
  i=0;
  j=1;
  // Lets find out how far we have to go back in the hierarchy to make space for the new value
  while (1) {
    if ( ( (self->t - ((self->tau_lin + 1)*((1<<(i+1))-1) + 1) )% (1<<(i+1)) == 0) ) {
      if ( i < (self->hierarchy_depth - 1) && self->n_vals[i]> self->tau_lin) {

        highest_level_to_compress+=1;
        i++;
      } else break;
    } else break;
  }

  // Now we know we must make space on the levels 0..highest_level_to_compress
  // Now lets compress the data level by level.

  for ( i = highest_level_to_compress; i >= 0; i-- ) {
    // We increase the index indicating the newest on level i+1 by one (plus folding)
    self->newest[i+1] = (self->newest[i+1] + 1) % (self->tau_lin+1);
    self->n_vals[i+1]+=1;
//    printf("t %d compressing level %d no %d and %d to level %d no %d, nv %d\n",self->t, i, (self->newest[i]+1) % (self->tau_lin+1),
//(self->newest[i]+2) % (self->tau_lin+1), i+1, self->newest[i+1], self->n_vals[i]);
    (*self->compressA)(self->A[i][(self->newest[i]+1) % (self->tau_lin+1)],  
                       self->A[i][(self->newest[i]+2) % (self->tau_lin+1)], 
                       self->A[i+1][self->newest[i+1]],self->dim_A);
    (*self->compressB)(self->B[i][(self->newest[i]+1) % (self->tau_lin+1)],  
                       self->B[i][(self->newest[i]+2) % (self->tau_lin+1)], 
                       self->B[i+1][self->newest[i+1]],self->dim_B);
  }

  self->newest[0] = ( self->newest[0] + 1 ) % (self->tau_lin +1); 
  self->n_vals[0]++;

  if ( (*self->A_fun)(self->A_args, self->A[0][self->newest[0]], self->dim_A) != 0 )
    return 1;
  if ( (*self->B_fun)(self->B_args, self->B[0][self->newest[0]], self->dim_B) != 0 )
    return 2;

  double* temp = malloc(self->dim_corr*sizeof(double));
  if (!temp)
    return 4;
// Now update the lowest level correlation estimates
  for ( j = 0; j < MIN(self->tau_lin+1, self->n_vals[0]); j++) {
    index_new = self->newest[0];
    index_old =  (self->newest[0] - j + self->tau_lin + 1) % (self->tau_lin + 1);
//    printf("old %d new %d\n", index_old, index_new);
    error = (self->corr_operation)(self->A[0][index_old], self->dim_A, self->B[0][index_new], self->dim_B, temp, self->dim_corr);
    if ( error != 0)
      return error;
    self->n_sweeps[j]++;
    for (k = 0; k < self->dim_corr; k++) {
      self->result[j][k] += temp[k];
    }
  }
// Now for the higher ones
  for ( i = 1; i < highest_level_to_compress+2; i++) {
    for ( j = (self->tau_lin+1)/2+1; j < MIN(self->tau_lin+1, self->n_vals[i]); j++) {
      index_new = self->newest[i];
      index_old = (self->newest[i] - j + self->tau_lin + 1) % (self->tau_lin + 1);
      index_res = self->tau_lin + (i-1)*self->tau_lin/2 + (j - self->tau_lin/2+1) -1;
      error=(self->corr_operation)(self->A[i][index_old], self->dim_A, self->B[i][index_new], self->dim_B, temp, self->dim_corr);
      if ( error != 0)
        return error;
      self->n_sweeps[index_res]++;
      for (k = 0; k < self->dim_corr; k++) {
        self->result[index_res][k] += temp[k];
      }
    }
  }
  free(temp);
  return 0;
}

int double_correlation_print_correlation( double_correlation* self, Tcl_Interp* interp) {

  int j, k;
  double dt=self->dt;
  char buffer[TCL_DOUBLE_SPACE];

  for (j=0; j<self->n_result; j++) {
     Tcl_AppendResult(interp, " { ", (char *)NULL);
     Tcl_PrintDouble(interp, self->tau[j]*dt, buffer);
     Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
     for (k=0; k< self->dim_corr; k++) {
     if (self->n_sweeps[j] == 0 ) {
       Tcl_PrintDouble(interp, 0., buffer);
       Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
     }
     else {
       Tcl_PrintDouble(interp, self->result[j][k]/ (double) self->n_sweeps[j], buffer);
       Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
     }
     }
     Tcl_AppendResult(interp, " } \n", (char *)NULL);
  }
  return 0;
}

int double_correlation_write_to_file( double_correlation* self, char* filename) {
  FILE* file=0;
  int j, k;
  double dt=self->dt;
  file=fopen(filename, "w");
  if (!file) {
    return 1;
  }
  for (j=0; j<self->n_result; j++) {
    fprintf(file, "%f %d ", self->tau[j]*dt, self->n_sweeps[j]);
    for (k=0; k< self->dim_corr; k++) {
      if (self->n_sweeps[j] == 0 )
        fprintf(file, "%f ", 0.);
      else 
        fprintf(file, "%f ", self->result[j][k]/ (double) self->n_sweeps[j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  return 0;
}

int identity ( double* input, unsigned int n_input, double* A, unsigned int dim_A) {
  int i; 
  if ( n_input != dim_A ) {
    return 5;
  }
  for ( i = 0; i < dim_A; i++ ) {
    A[i] = input[i];
  }
  return 0;
}

int compress_linear( double* A1, double*A2, double* A_compressed, unsigned int dim_A ) {
  unsigned int i;
  for ( i = 0; i < dim_A; i++ )
    A_compressed[i] = 0.5*(A1[i]+A2[i]);
  return 0;
}

int compress_discard1( double* A1, double*A2, double* A_compressed, unsigned int dim_A ) {
  unsigned int i;
  for ( i = 0; i < dim_A; i++ )
    A_compressed[i] = A2[i];
  return 0;
}

int compress_discard2( double* A1, double*A2, double* A_compressed, unsigned int dim_A ) {
  unsigned int i;
  for ( i = 0; i < dim_A; i++ )
    A_compressed[i] = A1[i];
  return 0;
}

int scalar_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) {
  double temp = 0;
  unsigned int i;
  if (!(dim_A == dim_B && dim_corr == 1 )) {
    printf("Error in scalar product: The vector sizes do not match");
    return 5;
  }
  for ( i = 0; i < dim_A; i++ ) {
    temp += A[i]*B[i];
  }
  C[0] = temp; 
  return 0;
}

int componentwise_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) {
  unsigned int i;
  if (!(dim_A == dim_B )) {
    printf("Error in componentwise product: The vector sizes do not match");
    return 5;
  }
  for ( i = 0; i < dim_A; i++ ) {
    C[i] = A[i]*B[i];
  }
  return 0;
}

int square_distance_componentwise ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) {
  unsigned int i;
  if (!(dim_A == dim_B )) {
    printf("Error in componentwise product: The vector sizes do not match\n");
    return 5;
  }
  for ( i = 0; i < dim_A; i++ ) {
    C[i] = (A[i]-B[i])*(A[i]-B[i]);
  }
  return 0;
}

int particle_velocities(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_total_particles)
      return 1;
    A[3*i + 0] = partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}

int com_velocity(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  double v_com[3] = { 0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_total_particles)
      return 1;
    v_com[0] += partCfg[ids->e[i]].m.v[0]/time_step;
    v_com[1] += partCfg[ids->e[i]].m.v[1]/time_step;
    v_com[2] += partCfg[ids->e[i]].m.v[2]/time_step;
  }
  A[0]=v_com[0]/ids->n;
  A[1]=v_com[1]/ids->n;
  A[2]=v_com[2]/ids->n;
  return 0;
}

int particle_positions(void* typelist, double* A, unsigned int n_A) {
  unsigned int i;
  sortPartCfg();
  if (typelist == NULL) {
    if (n_total_particles*3 != n_A) {
      return 1;
    } else {
      // Then everything seems to be alright
      for ( i = 0; i<n_total_particles; i++ ) {
        A[3*partCfg[i].p.identity + 0] = partCfg[i].r.p[0];
        A[3*partCfg[i].p.identity + 1] = partCfg[i].r.p[1];
        A[3*partCfg[i].p.identity + 2] = partCfg[i].r.p[2];
      }
      return 0;
    }
  } else 
    return 1;
}



int structure_factor(void* params_p, double* A, unsigned int n_A) {
  int i,j,k,l,p;
  int order, order2, n;
  double twoPI_L, C_sum, S_sum, qr;
  sf_params params;
  params = *(sf_params*)params_p;
  order = params.order;
  order2=order*order;
  twoPI_L = 2*PI/box_l[0];
  
  sortPartCfg();

    for(p=0; p<n_A; p++) {
       A[p]   = 0.0;
    }

    l=0;
    for(i=-order; i<=order; i++) {
      for(j=-order; j<=order; j++) {
        for(k=-order; k<=order; k++) {
	  n = i*i + j*j + k*k;
	  if ((n<=order2) && (n>=1)) {
	    C_sum = S_sum = 0.0;
	    for(p=0; p<n_total_particles; p++) {
	      qr = twoPI_L * ( i*partCfg[p].r.p[0] + j*partCfg[p].r.p[1] + k*partCfg[p].r.p[2] );
	      C_sum+= partCfg[p].p.scattering_length * cos(qr);
	      S_sum-= partCfg[p].p.scattering_length * sin(qr);
	    }
            A[l]   =C_sum;
            A[l+1] =S_sum;
            l=l+2;
	  }
	}
      }
    }
    return 0;
}


int file_data_source_init(file_data_source* self, char* filename, IntList* columns) {
  int counter=1;
  char* token;
  if (filename==0)
    return 1;
  self->f = fopen(filename, "r");
  if (! self->f )
    return 2;
  fgets(self->last_line, MAXLINELENGTH, self->f);
  while (self->last_line && self->last_line[0] == 35) {
    fgets(self->last_line, MAXLINELENGTH, self->f);
  }
  if (!self->last_line)
    return 3;
// Now lets count the tokens in the first line
  token=strtok(self->last_line, " \t\n");
  while (token) {
//    printf("reading token **%s**\n", token);
    token=strtok(NULL, " \t\n");
    counter++;
  }
  self->n_columns = counter;
  rewind(self->f);
  self->data_left=1;
//  printf("I found out that your file has %d columns\n", self->n_columns);
  if (columns !=0)
    /// Here we would like to check if we can handle the desired columns, but this has to be implemented!
    return -1;
  return 0;
}

int file_data_source_readline(void* xargs, double* A, int dim_A) {
  file_data_source* self = xargs;
  int counter=0;
  char* token;
  char* temp;

  temp=fgets(self->last_line, MAXLINELENGTH, self->f);
  while (temp!= NULL && self->last_line && self->last_line[0] == 35) {
    temp=fgets(self->last_line, MAXLINELENGTH, self->f);
  }
  if (!self->last_line || temp==NULL) {
//    printf("nothing left\n");
    self->data_left=0;
    return 3;
  }
  token=strtok(self->last_line, " \t\n");
  while (token) {
//    printf("reading token: ");
    A[counter]=atof(token);
//    printf("%f ", A[counter]);
    token=strtok(NULL, " \t\n");
    counter++;
    if (counter >= dim_A) {
//      printf("urgs\n");
      return 4;
    }
  }
//  printf("\n");
  return 0;
}

int tcl_input(void* data, double* A, unsigned int n_A) {
  tcl_input_data* input_data = (tcl_input_data*) data;
  int i, tmp_argc, res = 1;
  char  **tmp_argv;
  Tcl_SplitList(input_data->interp, input_data->argv[0], &tmp_argc, &tmp_argv);
  if (tmp_argc < n_A) {
    Tcl_AppendResult(input_data->interp, "Not enough arguments passed to analyze correlation update", (char *)NULL);
    return 1;
  }
  for (i = 0; i < n_A; i++) {
    if (Tcl_GetDouble(input_data->interp, tmp_argv[i], &A[i]) != TCL_OK) {
      Tcl_AppendResult(input_data->interp, "error parsing argument ", input_data->argv[i],"\n", (char *)NULL);
      return 1;
    }
  }
  return 0;
}

