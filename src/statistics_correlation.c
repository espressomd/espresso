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
int correlations_autoupdate=0;

const char init_errors[][64] = {
  "",                            // 0
  "No valid correlation given" , // 1
  "delta_t must be > 0",         // 2
  "tau_lin must be >= 2",        // 3
  "hierarchy_depth must be >=1", // 4
  "window_distance must be >1",  // 5
  "dimension of A was not >1",   // 6
  "dimension of B was not >1",   // 7
  "dim_corr was not >1",         // 8
  "no proper function for first observable given",                // 9
  "no proper function for second observable given",               //10 
  "no proper function for correlation operation given",           //11
  "no proper function for compression of first observable given", //12
  "no proper function for compression of second observable given",//13
  "tau_lin must be divisible by 2" // 14
  "dt is smaller than the MD timestep",//15
  "dt is not a multiple of the MD timestep"//16
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

int correlation_print_average1(double_correlation* self, Tcl_Interp* interp) {
  int i;
  char buffer[TCL_DOUBLE_SPACE];
  if (self->n_data < 1) {
    Tcl_AppendResult(interp, buffer, "Error in print average: No input data available", (char *)NULL);
    return TCL_ERROR;
  }
  for (i=0; i< self->dim_A; i++) {
    Tcl_PrintDouble(interp, self->A_accumulated_average[i]/self->n_data, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  return TCL_OK;
}

int correlation_print_variance1(double_correlation* self, Tcl_Interp* interp) {
  int i;
  char buffer[TCL_DOUBLE_SPACE];
  if (self->n_data < 1) {
    Tcl_AppendResult(interp, buffer, "Error in print variance: No input data available", (char *)NULL);
    return TCL_ERROR;
  }
  for (i=0; i< self->dim_A; i++) {
    Tcl_PrintDouble(interp, 
        self->A_accumulated_variance[i]/self->n_data
        - (self->A_accumulated_average[i]/self->n_data)*(self->A_accumulated_average[i]/self->n_data), buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  return TCL_OK;
}

int correlation_get_correlation_time(double_correlation* self, Tcl_Interp* interp) {
  int j, k;
  double dt=self->dt;
  char buffer[TCL_DOUBLE_SPACE];
  double* correlation_time = (double*) malloc(self->dim_corr*sizeof(double));
  if (self->dim_A != self->dim_corr) {
    Tcl_AppendResult(interp, buffer, "Error in print correlation_time: Only makes sense when the dimensions of  \
        the observables and correlation match (Isn't it?) ", (char *)NULL);
    return TCL_ERROR;
  }

// We calculate the correlation time for each dim_corr by normalizing the correlation,
// integrating it and finding out where C(tau)=tau;
  double C_tau;
  int ok_flag;
  for (j=0; j<self->dim_corr; j++) {
    correlation_time[j] = 0.; 
  }

  // here we still have to fix the stuff a bit!
  if (!self->autocorrelation) {
    Tcl_AppendResult(interp, buffer, "Does this make sense for an non-autocorrelation function? ", (char *)NULL);
    return TCL_ERROR;
  } else {
    for (j=0; j<self->dim_corr; j++) {
      C_tau=1*self->dt;
      ok_flag=0;
      for (k=1; k<self->n_result-1; k++) {
        if (self->n_sweeps[k]==0)
          break;
        C_tau+=(self->result[k][j]/ (double) self->n_sweeps[k] - self->A_accumulated_average[j]*self->B_accumulated_average[j]/self->n_data/self->n_data)/(self->result[0][j]/self->n_sweeps[0])*self->dt*(self->tau[k]-self->tau[k-1]);
        printf("C_tau %f tau %f exp(-W/tau) + 2*sqrt(W/N) %f corr %f deltat %f \n", 
            C_tau, 
            self->tau[k]*self->dt, 
            exp(-self->tau[k]*self->dt/C_tau)+2*sqrt(self->tau[k]*self->dt/self->n_data),
            self->result[k][j]/ (double) self->n_sweeps[k],
            self->dt*(self->tau[k]-self->tau[k-1])
            );

//        if (C_tau < i*self->tau[k]*self->dt) {
          if (exp(-self->tau[k]*self->dt/C_tau)+2*sqrt(self->tau[k]*self->dt/self->n_data)
              >exp(-self->tau[k-1]*self->dt/C_tau)+2*sqrt(self->tau[k-1]*self->dt/self->n_data)) {
          correlation_time[j]=C_tau*(1+(2*(double)self->tau[k]+1)/(double)self->n_data);
          printf("stopped at tau=>%f\n", self->tau[k]*self->dt);
          ok_flag=1;
          break;
         }
      }
      if (!ok_flag) {
        correlation_time[j]=-1;
      }
    }
  }

  for (j=0; j<self->dim_corr; j++) {
     Tcl_PrintDouble(interp, correlation_time[j], buffer);
     Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
  }
  free(correlation_time);
  return 0;
}

/* We can track several correlation functions at a time
*  identified by their ids.
*/
int tclcommand_analyze_parse_correlation(Tcl_Interp* interp, int argc, char** argv) {
  int no;
  char buffer[TCL_INTEGER_SPACE];
  if (argc < 1)
    return correlation_print_usage(interp);
  if (argc==1 && ARG0_IS_S("n_corr")) {
    sprintf(buffer,"%d ",n_correlations);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return TCL_OK;
  }
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

int tclcommand_correlation_parse_autoupdate(Tcl_Interp* interp, int no, int argc, char** argv) {
  int i;
  if (argc > 0 ) {
    if (ARG0_IS_S("start")) {
      if(correlations[no].update_frequency > 0) {
        correlations_autoupdate = 1;
        correlations[no].autoupdate=1;
        correlations[no].last_update=sim_time;
        return TCL_OK;
      } else {
        Tcl_AppendResult(interp, "Could not start autoupdate: update frequency not set\n", (char *)NULL);
        return TCL_ERROR;
      }
    } else if (ARG0_IS_S("stop")) {
      correlations_autoupdate=0;
      for (i=0; i<n_correlations; i++) {
        if (correlations[i].autoupdate)
          correlations_autoupdate=1;
      }
      return TCL_OK;
    }
  }
  Tcl_AppendResult(interp, "Usage: analyze correlation $n autoupdate [ start | stop ]\n", (char *)NULL);
  return TCL_ERROR;
}


int tclcommand_correlation_parse_print(Tcl_Interp* interp, int no, int argc, char** argv) {
  char buffer[TCL_INTEGER_SPACE];
  if(argc==0) {
      return double_correlation_print_correlation(&correlations[no], interp);
  } 
  if (ARG0_IS_S("spherically_averaged_sf")) {
  	  if(correlations[no].correlation_type==CORR_TYPE_SF) {
  	    if(argc!=2) { 
          Tcl_AppendResult(interp, "usage: analyze <correlation_id> print spherically_averaged_sf\n", (char *)NULL);
  	      return TCL_ERROR;
  	    }
  	  } else {
  	    sprintf(buffer,"%d ",correlations[no].correlation_type);
        Tcl_AppendResult(interp, "canot print spherically averaged sf for correlation type ", buffer, (char *)NULL);
  	    return TCL_ERROR;
  	  }
            return double_correlation_print_spherically_averaged_sf(&correlations[no], interp);
  } else if (ARG0_IS_S("average1")) {
    return correlation_print_average1(&correlations[no], interp);
  } else if (ARG0_IS_S("variance1")) {
    return correlation_print_variance1(&correlations[no], interp);
  } else if (ARG0_IS_S("correlation_time")) {
    return correlation_get_correlation_time(&correlations[no], interp);
  }
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
  int hierarchy_depth=0; 
  double delta_t = 1;
  double tau_max = 0;
  int(*corr_operation)  ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) = 0;
  int dim_corr;
  int change; // how many tcl argmuents are "consumed" by the parsing of arguments
  int error;
  tcl_input_data* tcl_input_p;
  char buffer[TCL_INTEGER_SPACE+2];
  int correlation_type; // correlation type of correlation which is currently being created
  int autocorrelation=1; // by default, we are doing autocorrelation
  int n_bins; // number of bins for spherically averaged sf
  int i;

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
          return tclcommand_correlation_parse_print(interp, no, argc-1, argv+1); 
      } else if (ARG0_IS_S("finalize")) {
        if(argc==1) {
          return double_correlation_finalize(&correlations[no]);
        } else {
          Tcl_AppendResult(interp, "Usage: analyze <correlation_id> finalize", buffer, (char *)NULL);
    	    return TCL_ERROR;
        }
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
      } else if (ARG0_IS_S("autoupdate")) {
        return tclcommand_correlation_parse_autoupdate(interp, no, argc-1, argv+1);
      } else if (ARG0_IS_S("update_from_file")) {
        error = correlation_update_from_file(no);
        if (error) {
          Tcl_AppendResult(interp, "error in update_from_file", (char *)NULL);
          return TCL_ERROR;
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
    else {
      Tcl_AppendResult(interp, "Usage for an already existing correlation:", (char *)NULL);
      Tcl_AppendResult(interp, "analyze correlation $ID [ print | update ]", (char *)NULL);
      return TCL_ERROR;
    }

  } else if ( no == n_correlations ) {  
  
    // Else we must parse the other arguments and see if we can construct a fully
    // working correlation class instance from that.
    correlation_type=CORR_TYPE_GENERIC; // this is the default
    while (argc > 0) {
      if ( ARG0_IS_S("first_obs") ) {
        argc -= 1;
        argv += 1;
        error=parse_observable(interp, argc, argv, &change, &A_fun, &dim_A, &A_args);
        argc -= change;
        argv += change;
        if (error)
          return TCL_ERROR;
      } else if ( ARG0_IS_S("second_obs") ) {
        argc -= 1;
        argv += 1;
        error = parse_observable(interp, argc, argv, &change, &B_fun, &dim_B, &B_args); 
	      autocorrelation=0;
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
      } else if ( ARG0_IS_S("tau_max") ) {
        if ( argc < 2 || !ARG1_IS_D(tau_max)) {
          Tcl_AppendResult(interp, "Usage: analyze correlation ... tau_max $hierarchy_depth", (char *)NULL);
	  return TCL_ERROR;
	} else { 
          argc -= 2;
          argv += 2;
        }
      } else if ( ARG0_IS_S("hierarchy_depth") ) {
        // FIXME This will be removed later
	Tcl_AppendResult(interp, "*** CHANGE *** specify tau_max instead of hierarchy_depth *** ", (char *)NULL);
	  return TCL_ERROR;
      } else if ( ARG0_IS_S("delta_t") ) {
        if ( argc < 2 || !ARG1_IS_D(delta_t)) {
          Tcl_AppendResult(interp, "Usage: analyze correlation ... delta_t $delta_t ", (char *)NULL);
        } else { 
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
	autocorrelation=0;
        argc -= 2;
        argv += 2; 
      } else if ( ARG0_IS_S("update") ) {
        sprintf(buffer,"%d ",no);
        Tcl_AppendResult(interp, "Correlation error: cannot update correlation ", buffer, (char *)NULL);
        Tcl_AppendResult(interp, " It must be set up first", (char *)NULL);
        return TCL_ERROR;
      } else if ( ARG0_IS_S("print") ) {
        sprintf(buffer,"%d ",no);
        Tcl_AppendResult(interp, "Correlation error: cannot print correlation ", buffer, (char *)NULL);
        Tcl_AppendResult(interp, " It must be set up first", (char *)NULL);
        return TCL_ERROR;
      } else {
        Tcl_AppendResult(interp, "unknown argument ", argv[0], (char *)NULL);
        return TCL_ERROR;
      }
    }
  } else {
      sprintf(buffer,"%d ",no);
      Tcl_AppendResult(interp,"Correlation error: cannot set up correlation with id ",buffer,(char *)NULL);
      sprintf(buffer,"%d ",n_correlations);
      Tcl_AppendResult(interp," Next correlation must have id ",buffer,(char *)NULL);
      return TCL_ERROR;
  }

  // Now we should find out, if this is enough to make a correlation.
  if ( ! (delta_t > 0.0) ) { 
    Tcl_AppendResult(interp, "delta_t must be > 0.0. Correlation cannot be performed.\n", (char *)NULL);
    return TCL_ERROR; 
  }
  if (tau_max <= delta_t) { 
    Tcl_AppendResult(interp, "tau_max must be >= delta_t. Correlation cannot be performed.\n", (char *)NULL);
    return TCL_ERROR; 
  } else {
    //set hierarchy depth which can  accomodate at least tau_max
    hierarchy_depth=(int)ceil( 1 + log( (tau_max/delta_t)/(tau_lin-1) ) / log(2.0) );
  }
  if(compressA==NULL) {
    Tcl_AppendResult(interp, "You have not chosen the compression function! Correlation cannot be performed.\n", (char *)NULL);
    return TCL_ERROR; 
  }
  if(autocorrelation) {
    dim_B=dim_A;
    B_fun=&obs_nothing;
    compressB=&compress_do_nothing;
  } else {
    printf("auto %d\n", autocorrelation);
    if(B_fun==NULL && compressB!=NULL) {
      Tcl_AppendResult(interp, "You have chosen compressB but not a function for computing observable B.\n", (char *)NULL);
       return TCL_ERROR; 
    } 
    if(compressB==NULL) {
        Tcl_AppendResult(interp, "You are not performing autocorrelation but have not chosen the compressB function. Taking compressB=compressA by default.\n", (char *)NULL);
	      compressB=compressA;
    }
  }
  
  // Let us just realloc space. Does not matter even if we can not make a proper correlation out of that.
  correlations=(double_correlation*) realloc(correlations, (n_correlations+1)*sizeof(double_correlation)); 

  // Now initialize the new correlation
  error = double_correlation_init(&correlations[n_correlations], delta_t, tau_lin, hierarchy_depth, 1, 
      dim_A, dim_B, dim_corr, A_fun, A_args, B_fun, B_args,
      corr_operation, compressA, compressB, correlation_type, autocorrelation);
  if ( error == 0 ) {
    n_correlations++;
    return TCL_OK;
  } else {
    printf("Error number %d\n", error);
    Tcl_AppendResult(interp, "correlation could not be corretly initialized: ", init_errors[error], "\n", (char *)NULL);
    return TCL_ERROR;
  }

}

int observable_usage(Tcl_Interp* interp) {
  Tcl_AppendResult(interp, "Usage: analyze correlation [ particle_velocities | com_velocity | particle_positions ", (char *) NULL);
#ifdef ELECTROSTATICS
  Tcl_AppendResult(interp, "particle_currents | currents ] ", (char *) NULL);
#endif
  Tcl_AppendResult(interp, "] [ first_obs | second_obs ] " , (char *)NULL);
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
  char ibuffer[TCL_INTEGER_SPACE + 2];
  char dbuffer[TCL_DOUBLE_SPACE];
  int *vals;
  double *q_density;
  params=(sf_params*)malloc(sizeof(sf_params));
  
  if(argc!=5) { 
    sprintf(ibuffer, "%d ", argc);
    Tcl_AppendResult(interp, "structure_factor  needs 5 arguments, got ", ibuffer, (char*)NULL);
    sf_print_usage(interp);
    return TCL_ERROR;
  }
  if (ARG_IS_I(1,order)) {
    sprintf(ibuffer, "%d ", order);
    if(order>1) {
      params->order=order;
      order2=order*order;
    } else {
      Tcl_AppendResult(interp, "order must be > 1, got ", ibuffer, (char*)NULL);
      sf_print_usage(interp);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, "problem reading order",(char*)NULL);
    return TCL_ERROR; 
  }
  if (ARG_IS_D(2,delta_t)) {
    if (delta_t > 0.0) *delta_t_p=delta_t;
    else {
      Tcl_PrintDouble(interp,delta_t,dbuffer);
      Tcl_AppendResult(interp, "delta_t must be > 0.0, got ", dbuffer,(char*)NULL);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, "problem reading delta_t, got ",argv[2],(char*)NULL);
    return TCL_ERROR; 
  }
  if (ARG_IS_D(3,tau_max)) {
    if (tau_max > 2.0*delta_t) *tau_max_p=tau_max;
    else {
      Tcl_PrintDouble(interp,tau_max,dbuffer);
      Tcl_AppendResult(interp, "tau_max must be > 2.0*delta_t, got ", dbuffer,(char*)NULL);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, "problem reading tau_max, got",argv[3],(char*)NULL);
    return TCL_ERROR; 
  }
  if (ARG_IS_I(4,tau_lin)) {
    if (tau_lin > 2 && tau_lin < (tau_max/delta_t+1)) *tau_lin_p=tau_lin;
    else {
      sprintf(ibuffer, "%d", tau_lin);
      Tcl_AppendResult(interp, "tau_lin must be < tau_max/delta_t+1, got ", ibuffer,(char*)NULL);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, "problem reading tau_lin, got",argv[4],(char*)NULL);
    sf_print_usage(interp);
    return TCL_ERROR; 
  }
  // compute the number of vectors
  l=0;
  for(i=-order; i<=order; i++) 
      for(j=-order; j<=order; j++) 
        for(k=-order; k<=order; k++) {
          n = i*i + j*j + k*k;
          if ((n<=order2) && (n>0)) {
            l++;
	  }
        }
  params->dim_sf=l;
  params->q_vals=(int*)malloc(3*l*sizeof(double));
  q_density=(double*)malloc(order2*sizeof(double));
  for(i=0;i<order2;i++) q_density[i]=0.0;
  l=0;
  // Store their values and density
  for(i=-order; i<=order; i++) 
      for(j=-order; j<=order; j++) 
        for(k=-order; k<=order; k++) {
          n = i*i + j*j + k*k;
          if ((n<=order2) && (n>0)) {
	    params->q_vals[3*l  ]=i;
	    params->q_vals[3*l+1]=j;
	    params->q_vals[3*l+2]=k;
	    q_density[n-1]+=1.0;
            l++;
	  }
        }
  for(i=0;i<order2;i++) q_density[i]/=(double)l;
  params->q_density=q_density;
  *A_args=(void*)params;
  *change=5; // if we reach this point, we have parsed 5 arguments, if not, error is returned anyway
  return 0;
}

// just a test function, will be removed later
void print_sf_params(sf_params *params) {
  int i, imax;
  int *vals;
  printf("order: %d\n",params->order);
  printf("dim_sf: %d\n",params->dim_sf);
  //printf("n_bins: %d\n",params->n_bins);
  //printf("qmax: %g\n",params->qmax);
  //printf("q2max2: %g\n",params->q2max);
  printf("q_vals: \n");
  imax=params->dim_sf;
  vals=params->q_vals;
  for(i=0;i<imax;i++)
    printf("i:%d %d %d %d\n",i,vals[3*i],vals[3*i+1],vals[3*i+ 2]);
  printf("End of sf_params\n");
  return;
}


int sf_print_usage(Tcl_Interp* interp) {
  Tcl_AppendResult(interp, "\nusage: structure_factor order delta_t tau_max  tau_lin", (char *)NULL);
  return TCL_ERROR;
}


static int convert_types_to_ids(IntList * type_list, IntList * id_list){ 
      int i,j,n_ids=0,flag;
      sortPartCfg();
      for ( i = 0; i<n_total_particles; i++ ) {
         if(type_list==NULL) { 
		/* in this case we select all particles */
               flag=1 ;
         } else {  
                flag=0;
                for ( j = 0; j<type_list->n ; j++ ) {
                    if(partCfg[i].p.type == type_list->e[j])  flag=1;
	        }
         }
	 if(flag==1){
              realloc_intlist(id_list, id_list->n=n_ids+1);
	      id_list->e[n_ids] = i;
	      n_ids++;
	 }
      }
      return n_ids;
}

int parse_id_list(Tcl_Interp* interp, int argc, char** argv, int* change, IntList** ids ) {
  int i,ret;
  char** temp_argv; int temp_argc;
  int temp;
  IntList* input=malloc(sizeof(IntList));
  IntList* output=malloc(sizeof(IntList));
  init_intlist(input);
  alloc_intlist(input,1);
  init_intlist(output);
  alloc_intlist(output,1);


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

  } else if ( ARG0_IS_S("types") ) {
    if (!parse_int_list(interp, argv[1],input)) {
      Tcl_AppendResult(interp, "Error parsing types list\n", (char *)NULL);
      return TCL_ERROR;
    } 
    if( (ret=convert_types_to_ids(input,output))<=0){
        Tcl_AppendResult(interp, "Error parsing types list. No particle with given types found.\n", (char *)NULL);
        return TCL_ERROR;
    } else { 
      *ids=output;
    }
    *change=2;
    return TCL_OK;
  } else if ( ARG0_IS_S("all") ) {
    if( (ret=convert_types_to_ids(NULL,output))<=0){
        Tcl_AppendResult(interp, "Error parsing keyword all. No particle found.\n", (char *)NULL);
        return TCL_ERROR;
    } else { 
      *ids=output;
    }
    *change=1;
    return TCL_OK;
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
  double cutoff;
  int *order_p;
  IntList* ids;
  IntList* ids1;
  iw_params* iw_params_p;
  
  profile_data* pdata=0;
  radial_profile_data* rpdata=0;

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
    if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
       return TCL_ERROR;
    *A_args=(void*)ids;
    *dim_A=3*ids->n;
    *change=1+temp;
    return TCL_OK;
  }
  if (ARG0_IS_S("density_profile")) {
    *A_fun = &density_profile;
    if (! tclcommand_parse_profile(interp, argc-1, argv+1, &temp, &dim_A, &pdata) == TCL_OK ) 
       return TCL_ERROR;
    *A_args=(void*)pdata;
    *dim_A=pdata->nbins;
    *change=1+temp;
    return TCL_OK;
  }
  if (ARG0_IS_S("radial_density_profile")) {
    *A_fun = &radial_density_profile;
    if (! tclcommand_parse_radial_profile(interp, argc-1, argv+1, &temp, &dim_A, &rpdata) == TCL_OK ) 
       return TCL_ERROR;
    *A_args=(void*)rpdata;
    *dim_A=rpdata->nbins;
    *change=1+temp;
    return TCL_OK;
  }
#ifdef ELECTROSTATICS
  if (ARG0_IS_S("particle_currents")) {
    *A_fun = &particle_currents;
    if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
      return TCL_ERROR;
    *A_args=(void*)ids;
    *dim_A=3*ids->n;
    *change=1+temp;
    return TCL_OK;
  }

  if (ARG0_IS_S("currents")) {
    *A_fun = &currents;
    if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
      return TCL_ERROR;
    *A_args=(void*)ids;
    *dim_A=3;
    *change=1+temp;
    return TCL_OK;
  } 
#endif

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
      sf_print_usage(interp);
      return TCL_ERROR; 
    }
  }
  if (ARG0_IS_S("interacts_with") ) {
    *A_fun = &interacts_with;
    ids=(IntList*)malloc(2*sizeof(IntList));
    if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
      return TCL_ERROR;
    iw_params_p=(iw_params*)malloc(sizeof(iw_params));
    iw_params_p->ids1=ids;
    *change=1+temp;
    if (! parse_id_list(interp, argc-3, argv+3, &temp, &ids1) == TCL_OK ) {
      return TCL_ERROR;
    }
    *change+=temp;
    iw_params_p->ids2=ids1;
    if ( argc < 5 || !ARG_IS_D(5,cutoff)) {
          Tcl_AppendResult(interp, "Usage: analyze correlation ... interacts_with id_list1 id_list2 cutoff", (char *)NULL);
	  return TCL_ERROR;
    }
    *change+=1;
    iw_params_p->cutoff=cutoff;
    *A_args=(void*)iw_params_p;
    *dim_A=ids->n; // number of ids from the 1st argument
    return TCL_OK;
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
  if (ARG_IS_S_EXACT(0,"componentwise_product")) {
    *corr_fun = &componentwise_product;
    *dim_corr = dim_A;
    *change=1;
    return TCL_OK;
  } else if (ARG_IS_S_EXACT(0,"complex_conjugate_product")) {
    *corr_fun = &complex_conjugate_product;
    *dim_corr = dim_A;
    *change=1;
    return TCL_OK;
  } else if (ARG_IS_S_EXACT(0,"square_distance_componentwise")) {
    *corr_fun = &square_distance_componentwise;
    *dim_corr = dim_A;
    *change=1;
    return TCL_OK;
  } else if (ARG_IS_S_EXACT(0,"square_distance")) {
    *corr_fun = &square_distance;
    *dim_corr = 1;
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

int tclcommand_parse_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, profile_data** pdata_) {
  int temp;
  *change=0;
  profile_data* pdata=(profile_data*)malloc(sizeof(profile_data));
  *pdata_ = pdata;
  pdata->id_list=0;
  pdata->startz=1e100;
  pdata->stopz=1e100;
  pdata->nbins=0;
  if (ARG0_IS_S("id") || ARG0_IS_S("type") || ARG0_IS_S("all")) {
    if (!parse_id_list(interp, argc, argv, &temp, &pdata->id_list )==TCL_OK) {
      Tcl_AppendResult(interp, "Error reading profile: Error parsing particle id information\n" , (char *)NULL);
      return TCL_ERROR;
    } else {
      *change+=temp;
      argc-=temp;
      argv+=temp;
    }
  } 
  if ( ARG0_IS_S("startz")){
    if (argc>1 && ARG1_IS_D(pdata->startz)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in profile: could not read startz\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if ( ARG0_IS_S("stopz") ) {
    if (argc>1 && ARG1_IS_D(pdata->stopz)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in profile: could not read stopz\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if (ARG0_IS_S("nbins")) {
    if (argc>1 && ARG1_IS_I(pdata->nbins)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in profile: could not read nbins\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  }
  
  temp=0;
  if (pdata->id_list==0) {
    Tcl_AppendResult(interp, "Error in profile: particle ids/types not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->startz>1e90) {
    Tcl_AppendResult(interp, "Error in profile: startz not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->stopz>1e90) {
    Tcl_AppendResult(interp, "Error in profile: stopz not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->nbins<1) {
    Tcl_AppendResult(interp, "Error in profile: nbins not specified\n" , (char *)NULL);
    temp=1;
  }
  if (temp)
    return TCL_ERROR;
  else
    return TCL_OK;
}

int tclcommand_parse_radial_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, radial_profile_data** pdata_) {
  int temp;
  *change=0;
  radial_profile_data* pdata=(radial_profile_data*)malloc(sizeof(radial_profile_data));
  *pdata_ = pdata;
  pdata->id_list=0;
  pdata->stopr=1e100;
  pdata->center[0]=1e100;pdata->center[1]=1e100;pdata->center[2]=1e100;
  pdata->nbins=0;
  if (ARG0_IS_S("id") || ARG0_IS_S("type") || ARG0_IS_S("all")) {
    if (!parse_id_list(interp, argc, argv, &temp, &pdata->id_list )==TCL_OK) {
      Tcl_AppendResult(interp, "Error reading profile: Error parsing particle id information\n" , (char *)NULL);
      return TCL_ERROR;
    } else {
      *change+=temp;
      argc-=temp;
      argv+=temp;
    }
  } 
  if ( ARG0_IS_S("center")){
    if (argc>3 && ARG1_IS_D(pdata->center[0]) && ARG_IS_D(2,pdata->center[1]) && ARG_IS_D(3,pdata->center[2])) {
      argc-=4;
      argv+=4;
      *change+=4;
    } else {
      Tcl_AppendResult(interp, "Error in radial_profile: could not read center\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if ( ARG0_IS_S("stopr") ) {
    if (argc>1 && ARG1_IS_D(pdata->stopr)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in radial_profile: could not read stopr\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  } 
  if (ARG0_IS_S("nbins")) {
    if (argc>1 && ARG1_IS_I(pdata->nbins)) {
      argc-=2;
      argv+=2;
      *change+=2;
    } else {
      Tcl_AppendResult(interp, "Error in radial_profile: could not read nbins\n" , (char *)NULL);
      return TCL_ERROR;
    } 
  }
  
  temp=0;
  if (pdata->id_list==0) {
    Tcl_AppendResult(interp, "Error in radial_profile: particle ids/types not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->center[0]>1e90) {
    Tcl_AppendResult(interp, "Error in radial_profile: center not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->stopr>1e90) {
    Tcl_AppendResult(interp, "Error in radial_profile: stopr not specified\n" , (char *)NULL);
    temp=1;
  }
  if (pdata->nbins<1) {
    Tcl_AppendResult(interp, "Error in radial_profile: nbins not specified\n" , (char *)NULL);
    temp=1;
  }
  if (temp)
    return TCL_ERROR;
  else
    return TCL_OK;
}

int double_correlation_init(double_correlation* self, double dt, unsigned int tau_lin, unsigned int hierarchy_depth, 
                  unsigned int window_distance, unsigned int dim_A, unsigned int dim_B, unsigned int dim_corr, 
                  void* A_fun, void* A_args, void* B_fun, void* B_args, void* corr_operation, 
                  void* compressA, void* compressB,
		  int correlation_type, int autocorrelation) {
  unsigned int i,j,k;
  self->dt = dt;
  if ((dt-time_step)<-1e-6*time_step) {
    return 15;
  }
  // check if dt is a multiple of the md timestep
  if ( abs(dt/time_step - round(dt/time_step)>1e-6 ) ) 
    return 16;
  self->update_frequency = (int) floor(dt/time_step);
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
  self->autocorrelation=autocorrelation;
  self->correlation_type=correlation_type;

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
  if (tau_lin%2)
    return 14;

  if (A_fun == &file_data_source_readline && (B_fun == &file_data_source_readline|| autocorrelation)) {
    self->is_from_file = 1;
  }

  self->A_data = (double*)malloc((tau_lin+1)*hierarchy_depth*dim_A*sizeof(double));
  if (autocorrelation) 
    self->B_data = self->A_data;
  else 
    self->B_data = (double*)malloc((tau_lin+1)*hierarchy_depth*dim_B*sizeof(double));
  
  self->n_data=0;
  self->A_accumulated_average = (double*)malloc(dim_A*sizeof(double));
  self->A_accumulated_variance= (double*)malloc(dim_A*sizeof(double));
  for (k=0; k<dim_A; k++) {
    self->A_accumulated_average[k]=0;
    self->A_accumulated_variance[k]=0;
  }
  if (autocorrelation) {
    self->B_accumulated_average =  self->A_accumulated_average;
    self->B_accumulated_variance = self->B_accumulated_variance;
  } else {
    self->B_accumulated_average = (double*)malloc(dim_B*sizeof(double));
    self->B_accumulated_variance = (double*)malloc(dim_B*sizeof(double));
    for (k=0; k<dim_B; k++) {
      self->B_accumulated_average[k]=0;
      self->B_accumulated_variance[k]=0;
    }
  }


  self->n_result=tau_lin+1 + (tau_lin+1)/2*(hierarchy_depth-1);
  self->tau = (int*)                malloc(self->n_result*sizeof(int));
  self->n_sweeps = (unsigned int*)  malloc(self->n_result*sizeof(int));
  self->result  = (double**)        malloc(self->n_result*sizeof(double*));
  self->result_data  = (double*)    malloc(self->n_result*dim_corr*sizeof(double));

  self->A = (double***)malloc(hierarchy_depth*sizeof(double**));
  if(autocorrelation) self->B = self->A;
  else self->B = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->n_vals = (unsigned int*) malloc(hierarchy_depth*sizeof(unsigned int));

  for (i=0; i<self->hierarchy_depth; i++) {
    self->A[i] = (double**) malloc((self->tau_lin+1)*sizeof(double*));
    if(!autocorrelation) self->B[i] = (double**) malloc((self->tau_lin+1)*sizeof(double*));
  }
  for (i=0; i<self->hierarchy_depth; i++) {
    self->n_vals[i]=0;
    for (j=0; j<self->tau_lin+1; j++) {
      self->A[i][j] = &self->A_data[(i*(tau_lin+1))*dim_A+j*dim_A];
      for (k=0; k<dim_A; k++) 
        self->A[i][j][k] = 0.;
      if(!autocorrelation) {
        self->B[i][j] = &self->B_data[(i*(tau_lin+1))*dim_B+j*dim_B];
        for (k=0; k<dim_B; k++) 
          self->B[i][j][k] = 0.;
      }
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

  // Now we update the cumulated averages and variances of A and B
  self->n_data++;
  for (k=0; k<self->dim_A; k++) {
    self->A_accumulated_average[k]+=self->A[0][self->newest[0]][k];
    self->A_accumulated_variance[k]+=self->A[0][self->newest[0]][k]*self->A[0][self->newest[0]][k];
  }
  // Here we check if it is an autocorrelation
  if (!self->autocorrelation) {
    for (k=0; k<self->dim_B; k++) {
      self->B_accumulated_average[k]+=self->B[0][self->newest[0]][k];
      self->B_accumulated_variance[k]+=self->B[0][self->newest[0]][k]*self->B[0][self->newest[0]][k];
    }
  } 

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

int double_correlation_finalize( double_correlation* self ) {
  // We must now go through the hierarchy and make sure there is space for the new 
  // datapoint. For every hierarchy level we have to decide if it necessary to move 
  // something
  int i,j,k;
  int ll=0; // current lowest level
  int vals_ll=0; // number of values remaining in the lowest level
  int dt=2;
  int highest_level_to_compress;
  unsigned int index_new, index_old, index_res;
  int error;
  //int compress;
  int tau_lin=self->tau_lin;
  int hierarchy_depth=self->hierarchy_depth;
  
  double* temp = malloc(self->dim_corr*sizeof(double));
  if (!temp)
    return 4;
  //printf ("tau_lin:%d, hierarchy_depth: %d\n",tau_lin,hierarchy_depth); 
  //for(ll=0;ll<hierarchy_depth;ll++) printf("n_vals[l=%d]=%d\n",ll, self->n_vals[ll]);
  for(ll=0;ll<hierarchy_depth-1;ll++) {
    if(self->n_vals[ll] > tau_lin+1 ) vals_ll = tau_lin + self->n_vals[ll]%2;
    else vals_ll=self->n_vals[ll];
    //printf("\nfinalizing level %d with %d vals initially\n",ll,vals_ll);
    
    while(vals_ll) {
      // Check, if we will want to push the value from the lowest level
      if(vals_ll % 2)  {
        highest_level_to_compress=ll; 
      } else {
        highest_level_to_compress=-1;
      }
      //printf("have %d values in ll=%d ",vals_ll,ll);
      //if(highest_level_to_compress<0) printf("Do NOT ");
      //printf("Compress\n"); fflush(stdout);
  
      i=ll+1; // lowest level, for which we have to check for compression 
      j=1; 
      // Lets find out how far we have to go back in the hierarchy to make space for the new value 
      while (highest_level_to_compress>-1) { 
        //printf("test level %d for compression, n_vals=%d ... ",i,self->n_vals[i]);
        if ( self->n_vals[i]%2 ) { 
	  if ( i < (hierarchy_depth - 1) && self->n_vals[i]> tau_lin) { 
            //printf("YES\n");
  	    highest_level_to_compress+=1; 
  	    i++; 
  	  } else {
            //printf("NO\n");
	    break;
	  }
        } else {
          //printf("NO\n");
	  break; 
        }
      }
      vals_ll-=1; 
      //printf("ll: %d, highest_lvevel_to_compress:%d\n",ll, highest_level_to_compress); fflush(stdout);
    
      // Now we know we must make space on the levels 0..highest_level_to_compress 
      // Now lets compress the data level by level.  
    
      for ( i = highest_level_to_compress; i >=ll; i-- ) { 
        // We increase the index indicating the newest on level i+1 by one (plus folding) 
        self->newest[i+1] = (self->newest[i+1] + 1) % (tau_lin+1); 
        self->n_vals[i+1]+=1; 
        //printf("compressing level %d no %d and %d to level %d no %d, nv %d\n",i, (self->newest[i]+1) % (tau_lin+1), (self->newest[i]+2) % (tau_lin+1), i+1, self->newest[i+1], self->n_vals[i]); 
        (*self->compressA)(self->A[i][(self->newest[i]+1) % (tau_lin+1)],  
	                   self->A[i][(self->newest[i]+2) % (tau_lin+1)], 
                           self->A[i+1][self->newest[i+1]],self->dim_A);
        (*self->compressB)(self->B[i][(self->newest[i]+1) % (tau_lin+1)],  
                           self->B[i][(self->newest[i]+2) % (tau_lin+1)], 
                           self->B[i+1][self->newest[i+1]],self->dim_B);
      } 
      self->newest[ll] = (self->newest[ll] + 1) % (tau_lin+1); 

      // We only need to update correlation estimates for the higher levels
      for ( i = ll+1; i < highest_level_to_compress+2; i++) {
        for ( j = (tau_lin+1)/2+1; j < MIN(tau_lin+1, self->n_vals[i]); j++) {
          index_new = self->newest[i];
          index_old = (self->newest[i] - j + tau_lin + 1) % (tau_lin + 1);
          index_res = tau_lin + (i-1)*tau_lin/2 + (j - tau_lin/2+1) -1;
          error=(self->corr_operation)(self->A[i][index_old], self->dim_A, self->B[i][index_new], self->dim_B, temp, self->dim_corr);
          if ( error != 0)
            return error;
          self->n_sweeps[index_res]++;
          for (k = 0; k < self->dim_corr; k++) {
            self->result[index_res][k] += temp[k];
          }
        }
      }
      // lowest level exploited, go upwards
      if(!vals_ll) { 
        //printf("reached end of level %d, go up\n",ll); 
	//fflush(stdout);
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
  char ibuffer[TCL_INTEGER_SPACE+2];

  for (j=0; j<self->n_result; j++) {
     Tcl_AppendResult(interp, " { ", (char *)NULL);
     Tcl_PrintDouble(interp, self->tau[j]*dt, buffer);
     Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
     sprintf(buffer, "%d ", self->n_sweeps[j]);
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

// TODO print it out
// Look at how static SF is printed
int double_correlation_print_spherically_averaged_sf(double_correlation* self, Tcl_Interp* interp) {

  int i,j,k,n_samp;
  int qi,qj,qk,qn, dim_sf, order2;
  double dt=self->dt;
  sf_params* params=(sf_params*)self->A_args;
  char buffer[TCL_DOUBLE_SPACE];
  char ibuffer[ 3*(TCL_INTEGER_SPACE+1) + 1 ];
  int *q_vals;
  double *q_density;
  double *av_sf_Re;
  double *av_sf_Im;
  
  q_vals=params->q_vals;
  q_density=params->q_density;
  dim_sf=params->dim_sf;
  order2=params->order*params->order;
  if(dim_sf!=self->dim_corr) printf("problem: dim_sf=%d != dim_corr=%d\n",dim_sf,self->dim_corr); 
  else printf("OK: dim_sf=%d != dim_corr=%d\n",dim_sf,self->dim_corr); fflush(stdout);
  fflush(stdout);
  
  
  av_sf_Re=(double*)malloc(order2*sizeof(double));
  av_sf_Im=(double*)malloc(order2*sizeof(double));

  // compute spherically averaged sf

  for (j=0; j<self->n_result; j++) {
 //   Tcl_AppendResult(interp, " { \n", (char *)NULL);
    // compute the spherically averaged sf for current dt
    for(k=0;k<order2;k++) av_sf_Re[k]=av_sf_Im[k]=0.0;
    if(j==0) printf("q_vals (n^2, i, j, k:\n");
    for(k=0;k<dim_sf;k++) {
      qi=q_vals[3*k  ];
      qj=q_vals[3*k+1];
      qk=q_vals[3*k+2];
      qn= qi*qi + qj*qj + qk*qk;
      av_sf_Re[qn-1]+=self->result[j][2*k  ];
      av_sf_Im[qn-1]+=self->result[j][2*k+1];
      if(j==0) printf("%d: %d, %d, %d, %d\n",k, qn, qi, qj, qk);
    }
    if(j==0) printf("q_density:\n");
    for(k=0;k<order2;k++) { 
      if(q_density[k]>0.0) {
        av_sf_Re[k]/=q_density[k];
        av_sf_Im[k]/=q_density[k];
      }
      if(j==0) printf("%d %lf\n",k, q_density[k]);
      // note: if q_density[k]==0, we did not add anything to av_sf_Xx[k], so it is 0.0
    }
    if(j==0) printf("\n\n");
    fflush(stdout);
    // now print what we obtained
    for(k=0;k<order2;k++) { 
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      Tcl_PrintDouble(interp, self->tau[j]*dt, buffer);
      Tcl_AppendResult(interp, buffer, " } { ",(char *)NULL);
      Tcl_PrintDouble(interp, (double) self->n_sweeps[j], buffer);
      Tcl_AppendResult(interp, buffer, " } { ", (char *)NULL);
      Tcl_PrintDouble(interp, sqrt((double)k+1.0), buffer);
      Tcl_AppendResult(interp, buffer, " } { ", (char *)NULL);
      if (self->n_sweeps[j] == 0 ) {
        Tcl_PrintDouble(interp, 0., buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
        Tcl_PrintDouble(interp, 0., buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      }
      else {
        Tcl_PrintDouble(interp, av_sf_Re[k]/(double) self->n_sweeps[j], buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
        Tcl_PrintDouble(interp, av_sf_Im[k]/(double) self->n_sweeps[j], buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      }
      Tcl_AppendResult(interp, " } \n", (char *)NULL);
    }
//    Tcl_AppendResult(interp, " } \n", (char *)NULL);
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

int compress_do_nothing( double* A1, double*A2, double* A_compressed, unsigned int dim_A ) {
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

int obs_nothing (void* params, double* A, unsigned int n_A) {
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
  unsigned int i,j;
  if (!(dim_A == dim_B )) {
    printf("Error in componentwise product: The vector sizes do not match");
    return 5;
  } 
  if (!(dim_A == dim_corr )) {
    printf("Error in componentwise product: The vector sizes do not match");
    return 5;
  }
  for ( i = 0; i < dim_A; i++ )
    C[i] = A[i]*B[i];
  return 0;
}

int complex_conjugate_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) {
  unsigned int i,j;
  if (!(dim_A == dim_B )) {
    printf("Error in complex_conjugate product: The vector sizes do not match");
    return 5;
  }
  j=0;
  for ( i = 0; i < dim_A/2; i++ ) {
    C[j] = A[j]*B[j] + A[j+1]*B[j+1];
    C[j+1] = A[j+1]*B[j] - A[j]*B[j+1];
    j=j+2;
  }
  return 0;
}

int square_distance ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) {
  unsigned int i;
  double tmp=0.0;
  if (!(dim_A == dim_B )) {
    printf("Error in square distance: The vector sizes do not match\n");
    return 5;
  }
  for ( i = 0; i < dim_A; i++ ) {
    tmp += (A[i]-B[i])*(A[i]-B[i]);
  }
  C[0]=tmp;
  return 0;
}


int square_distance_componentwise ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) {
  unsigned int i;
  if (!(dim_A == dim_B )) {
    printf("Error in square distance componentwise: The vector sizes do not match\n");
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
    if (ids->e[i] >= n_total_particles)
      return 1;
    A[3*i + 0] = partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}

#ifdef ELECTROSTATICS
int particle_currents(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  double charge;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    A[3*i + 0] = charge * partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = charge * partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = charge * partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}
int currents(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_total_particles)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    j[0] += charge * partCfg[ids->e[i]].m.v[0]/time_step;
    j[1] += charge * partCfg[ids->e[i]].m.v[1]/time_step;
    j[2] += charge * partCfg[ids->e[i]].m.v[2]/time_step;
  }
  A[0]=j[0];
  A[1]=j[1];
  A[2]=j[2];
  return 0;
}
#endif

int com_velocity(void* idlist, double* A, unsigned int n_A) {
/* TODO: this does not work with MASS ... */
  unsigned int i;
  double v_com[3] = { 0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    v_com[0] += partCfg[ids->e[i]].m.v[0]/time_step;
    v_com[1] += partCfg[ids->e[i]].m.v[1]/time_step;
    v_com[2] += partCfg[ids->e[i]].m.v[2]/time_step;
  }
  A[0]=v_com[0]/ids->n;
  A[1]=v_com[1]/ids->n;
  A[2]=v_com[2]/ids->n;
  printf("v_com %f %f %f\n", A[0], A[1], A[2]);
  return 0;
}

int density_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i;
  int bin;
  double ppos[3];
  int img[3];
  IntList* ids;
  profile_data* pdata;
  sortPartCfg();
  pdata=(profile_data*) pdata_;
  ids=pdata->id_list;
    
  for ( i = 0; i<pdata->nbins; i++ ) {
    A[i]=0;
  }
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
/* We use folded coordinates here */
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    bin= (int) floor( pdata->nbins*  (ppos[2]-pdata->startz)/(pdata->stopz-pdata->startz));
/* uncomment this line for unfolded coordinates */
//    bin= (int) floor( pdata->nbins*  (partCfg[ids->e[i]].r.p[2]-pdata->startz)/(pdata->stopz-pdata->startz));
    if (bin>=0 && bin < pdata->nbins) {
      A[bin] += 1./box_l[0]/box_l[1]/(pdata->stopz - pdata->startz)*pdata->nbins;
    }
  }
  return 0;
}

int radial_density_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i;
  int bin;
  double ppos[3];
  double r;
  int img[3];
  double bin_volume;
  IntList* ids;
  radial_profile_data* pdata;
  sortPartCfg();
  pdata=(radial_profile_data*) pdata_;
  ids=pdata->id_list;
    
  for ( i = 0; i<pdata->nbins; i++ ) {
    A[i]=0;
  }
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
/* We use folded coordinates here */
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    r=sqrt( (ppos[0]-pdata->center[0])*(ppos[0]-pdata->center[0])+(ppos[1]-pdata->center[1])*(ppos[1]-pdata->center[1]));
    bin= (int) floor( pdata->nbins*  r/pdata->stopr);

    if (bin>=0 && bin < pdata->nbins) {
      bin_volume=PI*((bin+1)*(bin+1)-bin*bin)*pdata->stopr/pdata->nbins*box_l[2];
      A[bin] += 1./bin_volume;
    }
  }
  return 0;
}

int particle_positions(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
      A[3*i + 0] = partCfg[ids->e[i]].r.p[0];
      A[3*i + 1] = partCfg[ids->e[i]].r.p[1];
      A[3*i + 2] = partCfg[ids->e[i]].r.p[2];
  }
  return 0;
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
    //printf("n_A: %d, dim_sf: %d\n",n_A, params.dim_sf); fflush(stdout);
    for(i=-order; i<=order; i++) {
      for(j=-order; j<=order; j++) {
        for(k=-order; k<=order; k++) {
	  n = i*i + j*j + k*k;
	  if ((n<=order2) && (n>=1)) {
	    C_sum = S_sum = 0.0;
            //printf("l: %d, n: %d %d %d\n",l,i,j,k); fflush(stdout);
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
    //printf("finished calculating sf\n"); fflush(stdout);
    return 0;
}

int interacts_with (void* params_p, double* A, unsigned int n_A) {
  iw_params *params=(iw_params*)params_p;
  IntList* ids1;
  IntList* ids2;
  int i,j;
  double dx,dy,dz,dist2;
  double cutoff2=params->cutoff*params->cutoff;
  double pos1[3], pos2[3], dist[3];
  ids1=params->ids1;
  ids2=params->ids2;
  sortPartCfg();
  for ( i = 0; i<ids1->n; i++ ) {
    if (ids1->e[i] >= n_total_particles)
      return 1;
    pos1[0]=partCfg[ids1->e[i]].r.p[0];
    pos1[1]=partCfg[ids1->e[i]].r.p[1];
    pos1[2]=partCfg[ids1->e[i]].r.p[2];
    for ( j = 0; j<ids2->n; j++ ) {
      if (ids2->e[j] >= n_total_particles)
        return 1;
      A[i] = 0;
      pos2[0]=partCfg[ids2->e[j]].r.p[0];
      pos2[1]=partCfg[ids2->e[j]].r.p[1];
      pos2[2]=partCfg[ids2->e[j]].r.p[2];
      get_mi_vector(dist,pos1,pos2);
      dist2= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if(dist2<cutoff2) {
        A[i] = 1;
	break;
	// interaction found for i, go for next
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
  const char  **tmp_argv;
  Tcl_SplitList(input_data->interp, input_data->argv[0], &tmp_argc, &tmp_argv);
  // function prototype from man page:
  // int Tcl_SplitList(interp, list, argcPtr, argvPtr)
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

void autoupdate_correlations() {
  int i;
  for (i=0; i<n_correlations; i++) {
//    printf("checking correlation %d autoupdate is %d \n", i, correlations[i].autoupdate);
    if (correlations[i].autoupdate && sim_time-correlations[i].last_update>correlations[i].dt*0.99999) {
//      printf("updating %d\n", i);
      correlations[i].last_update=sim_time;
      double_correlation_get_data(&correlations[i]);
    }
  }
}
