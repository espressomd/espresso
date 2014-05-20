 /*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  
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
#include "statistics_correlation.hpp"
#include "particle_data.hpp"
#include "integrate.hpp"
#include <cstring>

/* global variables */
double_correlation* correlations=0;
unsigned int n_correlations = 0;

/* forward declarations */
int correlation_get_correlation_time(double_correlation* self, double* correlation_time);
int double_correlation_finalize( double_correlation* self );
int read_until_end_of_line(FILE *f);


/* Error codes */
const char init_errors[][64] = {
  "",                            // 0
  "No valid correlation given" , // 1
  "delta_t must be specified and > 0",         // 2
  "tau_lin must be >= 2",        // 3
  "tau_max must be >= delta_t", //4
  "window_distance must be >1",  // 5
  "dimension of A was not >1",   // 6
  "dimension of B was not >1",   // 7
  "dimension of B must match dimension of A ",   // 8
  "no proper function for first observable given",                // 9
  "no proper function for second observable given",               //10 
  "no proper function for correlation operation given",           //11
  "no proper function for compression of first observable given", //12
  "no proper function for compression of second observable given",//13
  "tau_lin must be divisible by 2", // 14
  "dt is smaller than the MD timestep",//15
  "dt is not a multiple of the MD timestep", //16
  "cannot set compress2 for autocorrelation", //17
  "dim_A must be divisible by 3 for fcs_acf", //18
  "fcs_acf requires 3 additional parameters"  //19
};

const char file_data_source_init_errors[][64] = {
  "",
  "No valid filename given." ,
  "File could not be opened.",
  "No line found that was not commented out"
};

const char double_correlation_get_data_errors[][64] = {
  "",                                                                // 0
  "Error calculating variable A\n" ,                                   // 2
  "Error calculating variable B\n" ,                                   // 3
  "Error calculating correlation\n",                                 // 4
  "Error allocating temporary memory\n",                             // 4
  "Error in corr_operation: observable dimensions do not match\n",   // 5
  "Error: dim_corr and dim_A do not match for fcs_acf\n"             // 6
};

int correlations_autoupdate=0;


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



int correlation_get_correlation_time(double_correlation* self, double* correlation_time) {
  // We calculate the correlation time for each dim_corr by normalizing the correlation,
  // integrating it and finding out where C(tau)=tau;
  double C_tau;
  int ok_flag;
  for (unsigned j=0; j<self->dim_corr; j++) {
    correlation_time[j] = 0.; 
  }

  // here we still have to fix the stuff a bit!
  for (unsigned j=0; j<self->dim_corr; j++) {
    C_tau=1*self->dt;
    ok_flag=0;
    for (unsigned k=1; k<self->n_result-1; k++) {
      if (self->n_sweeps[k]==0)
        break;
      C_tau+=(self->result[k][j]/ (double) self->n_sweeps[k] - self->A_accumulated_average[j]*self->B_accumulated_average[j]/self->n_data/self->n_data)/(self->result[0][j]/self->n_sweeps[0])*self->dt*(self->tau[k]-self->tau[k-1]);
//        printf("C_tau %f tau %f exp(-W/tau) + 2*sqrt(W/N) %f corr %f deltat %f \n", 
//            C_tau, 
//            self->tau[k]*self->dt, 
//            exp(-self->tau[k]*self->dt/C_tau)+2*sqrt(self->tau[k]*self->dt/self->n_data),
//            self->result[k][j]/ (double) self->n_sweeps[k],
//            self->dt*(self->tau[k]-self->tau[k-1])
//            );

//        if (C_tau < i*self->tau[k]*self->dt) {
        if (exp(-self->tau[k]*self->dt/C_tau)+2*sqrt(self->tau[k]*self->dt/self->n_data)
            >exp(-self->tau[k-1]*self->dt/C_tau)+2*sqrt(self->tau[k-1]*self->dt/self->n_data)) {
        correlation_time[j]=C_tau*(1+(2*(double)self->tau[k]+1)/(double)self->n_data);
//          printf("stopped at tau=>%f\n", self->tau[k]*self->dt);
        ok_flag=1;
        break;
       }
    }
    if (!ok_flag) {
      correlation_time[j]=-1;
    }
  }

  return 0;
}


int double_correlation_init(double_correlation* self, double dt, unsigned int tau_lin, double tau_max,
                  unsigned int window_distance, unsigned int dim_A, unsigned int dim_B, unsigned int dim_corr, 
                  observable* A, observable* B, char* corr_operation_name, 
                  char* compressA_name, char* compressB_name, void *args) {
  unsigned int i,j,k;
  unsigned int hierarchy_depth=0;

  if (self==0)
    return 1;
  // first input-independent values
  self->t = 0;
  self->finalized=0;
  self->autoupdate=0;
  self->autocorrelation=1; // the default may change later if dim_B != 0
  
  // then input-dependent ones
  if (dt <= 0)
    return 2;
  if ((dt-time_step)<-1e-6*time_step) {
    return 15;
  }
  // check if dt is a multiple of the md timestep
  if ( abs(dt/time_step - round(dt/time_step)) >1e-6 ) 
    return 16;
  self->dt = dt;
  self->update_frequency = (int) floor(dt/time_step);
  
  if ( tau_lin == 1 ) { // use the default
    tau_lin=(int)ceil(tau_max/dt);
    printf("tau_lin: %d\n", tau_lin);
    if (tau_lin%2) tau_lin+=1;
  }
  if (tau_lin<2)
    return 3;
  if (tau_lin%2)
    return 14;
  self->tau_lin=tau_lin;
  
  if (tau_max <= dt) { 
    return 4;
  } else { //set hierarchy depth which can  accomodate at least tau_max
    if ( (tau_max/dt) < tau_lin ) 
      hierarchy_depth = 1;
    else 
      hierarchy_depth=(unsigned int)ceil( 1 + log( (tau_max/dt)/(tau_lin-1) ) / log(2.0) );
  }
  self->tau_max=tau_max;
  self->hierarchy_depth = hierarchy_depth;
  
  if (window_distance<1)
    return 5;
  self->window_distance = window_distance;
  
  if (dim_A<1)
    return 6;
  self->dim_A = dim_A;
  
  if (dim_B==0)
    dim_B = dim_A;
  else if (dim_B>0)
    self->autocorrelation=0;
  else if (dim_A != dim_B) 
    return 8;
    // currently there is no correlation function able to handel observables of different dimensionalities
  else 
    return 7;
  self->dim_B = dim_B;

  if (A == 0)
    return 9;
  self->A_obs = A;
  
  if (B == 0 && !self->autocorrelation)
    return 10;
  self->B_obs = B;
  

  // choose the correlation operation 
  if (corr_operation_name==0) { 
    return 11; // there is no reasonable default
  } else if ( strcmp(corr_operation_name,"componentwise_product") == 0 ) {
    dim_corr = dim_A;
    self->corr_operation = &componentwise_product;
    self->args = NULL;
  } else if ( strcmp(corr_operation_name,"complex_conjugate_product") == 0 ) {
    dim_corr = dim_A;
    self->corr_operation = &complex_conjugate_product;
    self->args = NULL;
  } else if ( strcmp(corr_operation_name,"square_distance_componentwise") == 0 ) {
    dim_corr = dim_A;
    self->corr_operation = &square_distance_componentwise;
    self->args = NULL;
  } else if ( strcmp(corr_operation_name,"fcs_acf") == 0 ) {
    if (dim_A %3 )
      return 18;
    dim_corr = dim_A/3;
    self->corr_operation = &fcs_acf;
    self->args = args;
// square_distance will be removed -- will be replaced by strides and blocks
//  } else if ( strcmp(corr_operation_name,"square_distance") == 0 ) {
//    self->corr_operation = &square_distance;
  } else if ( strcmp(corr_operation_name,"scalar_product") == 0 ) {
    dim_corr=1;
    self->corr_operation = &scalar_product;
    self->args = NULL;
  } else {
    return 11; 
  }
  self->dim_corr = dim_corr;
  self->corr_operation_name = corr_operation_name;
  
  // Choose the compression function
  if (compressA_name==0) { // this is the default
    compressA_name=strdup("discard2");
    self->compressA=&compress_discard2;
  } else if ( strcmp(compressA_name,"discard2") == 0 ) {
    self->compressA=&compress_discard2;
  } else if ( strcmp(compressA_name,"discard1") == 0 ) {
    self->compressA=&compress_discard1;
  } else if ( strcmp(compressA_name,"linear") == 0 ) {
    self->compressA=&compress_linear;
  } else {
    return 12;
  }
  self->compressA_name=compressA_name; 
  
  if (compressB_name==0) { 
    if(self->autocorrelation) { // the default for autocorrelation
      compressB_name=strdup("none"); 
      self->compressB=&compress_do_nothing;
    } else { // the default for corsscorrelation
      compressB_name=self->compressA_name;
      self->compressB=self->compressA;
    } 
  } else if ( self->autocorrelation ) {
    return 17;
  } else if ( strcmp(compressB_name,"discard2") == 0 ) {
    self->compressB=&compress_discard2;
  } else if ( strcmp(compressB_name,"discard1") == 0 ) {
    self->compressB=&compress_discard1;
  } else if ( strcmp(compressB_name,"linear") == 0 ) {
    self->compressB=&compress_linear;
  } else {
    return 13;
  }
  self->compressB_name=compressB_name; 
 
//  if (A_fun == &file_data_source_readline && (B_fun == &file_data_source_readline|| autocorrelation)) {
//    self->is_from_file = 1;
//  } else {
//    self->is_from_file = 0;
//  }

  self->A_data = (double*)malloc((tau_lin+1)*hierarchy_depth*dim_A*sizeof(double));
  if (self->autocorrelation) 
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
  if (self->autocorrelation) {
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
  self->result_data  = (double*)    malloc(self->n_result*dim_corr*sizeof(double));
  self->n_vals = (unsigned int*) malloc(hierarchy_depth*sizeof(unsigned int));

  // allocate space for convenience pointer to A and B buffers
  self->A = (double***)malloc(hierarchy_depth*sizeof(double**));
  if(self->autocorrelation) self->B = self->A;
  else self->B = (double***)malloc(hierarchy_depth*sizeof(double**));

  for (i=0; i<self->hierarchy_depth; i++) {
    self->A[i] = (double**) malloc((self->tau_lin+1)*sizeof(double*));
    if(!self->autocorrelation) self->B[i] = (double**) malloc((self->tau_lin+1)*sizeof(double*));
  }

  // and initialize the values
  for (i=0; i<self->hierarchy_depth; i++) {
    self->n_vals[i]=0;
    for (j=0; j<self->tau_lin+1; j++) {
      self->A[i][j] = &self->A_data[(i*(tau_lin+1))*dim_A+j*dim_A];
      for (k=0; k<dim_A; k++) 
        self->A[i][j][k] = 0.;
      if(!self->autocorrelation) {
        self->B[i][j] = &self->B_data[(i*(tau_lin+1))*dim_B+j*dim_B];
        for (k=0; k<dim_B; k++) 
          self->B[i][j][k] = 0.;
      }
    }
  }

  // allocate space for convenience pointer to the result
  self->result  = (double**)        malloc(self->n_result*sizeof(double*));
  for (i=0; i<self->n_result; i++) {
    self->n_sweeps[i]=0;
    self->result[i]=&self->result_data[i*self->dim_corr];
    for (j=0; j<self->dim_corr; j++) 
      // and initialize the values
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
  int i,j;
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
      if ( i < (int(self->hierarchy_depth) - 1) && self->n_vals[i]> self->tau_lin) {

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
    if (!self->autocorrelation)
      (*self->compressB)(self->B[i][(self->newest[i]+1) % (self->tau_lin+1)],  
                       self->B[i][(self->newest[i]+2) % (self->tau_lin+1)], 
                       self->B[i+1][self->newest[i+1]],self->dim_B);
  }

  self->newest[0] = ( self->newest[0] + 1 ) % (self->tau_lin +1); 
  self->n_vals[0]++;

  if ( observable_calculate(self->A_obs) != 0 )
    return 1;
  // copy the result:
  memcpy(self->A[0][self->newest[0]], self->A_obs->last_value, self->dim_A*sizeof(double));

  if (!self->autocorrelation) {
    if ( observable_calculate(self->B_obs) != 0 )
      return 2;
    memcpy(self->B[0][self->newest[0]], self->B_obs->last_value, self->dim_B*sizeof(double));
  }

  // Now we update the cumulated averages and variances of A and B
  self->n_data++;
  for (unsigned k=0; k<self->dim_A; k++) {
    self->A_accumulated_average[k]+=self->A[0][self->newest[0]][k];
    self->A_accumulated_variance[k]+=self->A[0][self->newest[0]][k]*self->A[0][self->newest[0]][k];
  }
  // Here we check if it is an autocorrelation
  if (!self->autocorrelation) {
    for (unsigned k=0; k<self->dim_B; k++) {
      self->B_accumulated_average[k]+=self->B[0][self->newest[0]][k];
      self->B_accumulated_variance[k]+=self->B[0][self->newest[0]][k]*self->B[0][self->newest[0]][k];
    }
  } 

  double* temp = (double*)malloc(self->dim_corr*sizeof(double));
  if (!temp)
    return 4;
// Now update the lowest level correlation estimates
  for ( j = 0; j < int(MIN(self->tau_lin+1, self->n_vals[0]) ); j++) {
    index_new = self->newest[0];
    index_old =  (self->newest[0] - j + self->tau_lin + 1) % (self->tau_lin + 1);
//    printf("old %d new %d\n", index_old, index_new);
    error = (self->corr_operation)(self->A[0][index_old], self->dim_A, self->B[0][index_new], self->dim_B, temp, self->dim_corr, self->args);
    if ( error != 0)
      return error;
    self->n_sweeps[j]++;
    for (unsigned k = 0; k < self->dim_corr; k++) {
      self->result[j][k] += temp[k];
    }
  }
// Now for the higher ones
  for ( int i = 1; i < highest_level_to_compress+2; i++) {
    for ( unsigned j = (self->tau_lin+1)/2+1; j < MIN(self->tau_lin+1, self->n_vals[i]); j++) {
      index_new = self->newest[i];
      index_old = (self->newest[i] - j + self->tau_lin + 1) % (self->tau_lin + 1);
      index_res = self->tau_lin + (i-1)*self->tau_lin/2 + (j - self->tau_lin/2+1) -1;
      error=(self->corr_operation)(self->A[i][index_old], self->dim_A, self->B[i][index_new], self->dim_B, temp, self->dim_corr, self->args);
      if ( error != 0)
        return error;
      self->n_sweeps[index_res]++;
      for (unsigned k = 0; k < self->dim_corr; k++) {
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
  int i,j;
  int ll=0; // current lowest level
  int vals_ll=0; // number of values remaining in the lowest level
  int highest_level_to_compress;
  unsigned int index_new, index_old, index_res;
  int error;
  //int compress;
  unsigned tau_lin=self->tau_lin;
  int hierarchy_depth=self->hierarchy_depth;

  double* temp = (double*)malloc(self->dim_corr*sizeof(double));

  // make a flag that the correlation is finalized
  self->finalized=1;

  if (!temp)
    return 4;
  //printf ("tau_lin:%d, hierarchy_depth: %d\n",tau_lin,hierarchy_depth); 
  //for(ll=0;ll<hierarchy_depth;ll++) printf("n_vals[l=%d]=%d\n",ll, self->n_vals[ll]);
  for(ll=0;ll<hierarchy_depth-1;ll++) {
    if (self->n_vals[ll] > tau_lin+1 ) vals_ll = tau_lin + self->n_vals[ll]%2;
    else vals_ll=self->n_vals[ll];
    //printf("\nfinalizing level %d with %d vals initially\n",ll,vals_ll);
    
    while (vals_ll) {
      // Check, if we will want to push the value from the lowest level
      if (vals_ll % 2)  {
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
	  if ( i < (hierarchy_depth-1) && self->n_vals[i]> tau_lin) { 
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
        for ( j = (tau_lin+1)/2+1; j < int(MIN(tau_lin+1, self->n_vals[i])); j++) {
          index_new = self->newest[i];
          index_old = (self->newest[i] - j + tau_lin + 1) % (tau_lin + 1);
          index_res = tau_lin + (i-1)*tau_lin/2 + (j - tau_lin/2+1) -1;
          error=(self->corr_operation)(self->A[i][index_old], self->dim_A, self->B[i][index_new], self->dim_B, temp, self->dim_corr, self->args);
          if ( error != 0)
            return error;
          self->n_sweeps[index_res]++;
          for (unsigned k = 0; k < self->dim_corr; k++) {
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


int identity ( double* input, unsigned int n_input, double* A, unsigned int dim_A) {
  if ( n_input != dim_A ) {
    return 5;
  }
  for (unsigned i = 0; i < dim_A; i++ ) {
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

int scalar_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, void *args ) {
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

int componentwise_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, void *args ) {
  unsigned int i;
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

int complex_conjugate_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, void *args ) {
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


int square_distance_componentwise ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, void *args ) {
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


int fcs_acf ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, void *args ) {
  DoubleList *wsquare = (DoubleList*)args;
  if (args == NULL )
    return 1;
  if (!(dim_A == dim_B )) {
    return 5;
  }
  if ( dim_A / dim_corr != 3) {
    return 6; 
  }
  for (unsigned i = 0; i < dim_corr; i++ ) 
    C[i] = 0;
  for (unsigned i = 0; i < dim_A; i++ ) {
    C [i/3] -= ( (A[i]-B[i])*(A[i]-B[i]) ) / wsquare->e[i%3];
  }
  for (unsigned i = 0; i < dim_corr; i++ ) 
    C[i] = exp(C[i]);
  return 0;
}


void autoupdate_correlations() {
  for (unsigned i=0; i<n_correlations; i++) {
    if (correlations[i].autoupdate && sim_time-correlations[i].last_update>correlations[i].dt*0.99999) {
      correlations[i].last_update=sim_time;
      double_correlation_get_data(&correlations[i]);
    }
  }
}


