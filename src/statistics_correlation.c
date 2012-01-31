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

/* global variables */
double_correlation* correlations=0;
unsigned int n_correlations = 0;

/* forward declarations */
int correlation_get_correlation_time(double_correlation* self, double* correlation_time);
int double_correlation_finalize( double_correlation* self );
int read_until_end_of_line(FILE *f);


/* Error codes */
const char init_from_checkpoint_errors[][64] = {
  "                                                ",                 //0
  "unknown error",                                                    //1
  "cannot open the file",                                             //2
  "failed reading autocorrelation",                                   //3
  "failed reading finalized",                                         //4
  "failed reading hierarchy_depth",                                   //5
  "failed reading tau_lin",                                           //6
  "failed reading dim_A",                                             //7
  "dim_A from checkpoint != dim_A from observable",                   //8
  "failed reading dim_B",                                             //9
  "failed reading dim_corr",                                          //10
  "failed reading t",                                                 //11
  "failed reading dt",                                                //12
  "failed reading tau_max",                                           //13
  "failed reading update_frequency",                                  //14
  "failed reading window_distance",                                   //15
  "failed reading n_result",                                          //16
  "failed reading n_data",                                            //17
  "failed reading compressA_name",                                    //18
  "failed reading compressB_name",                                    //19
  "failed reading corr_operation_name",                               //20
  "failed reading is_from_file",                                      //21
  "failed reading autoupdate",                                        //22
  "failed reading last_update",                                       //23
  "failed reading tau",                                               //24
  "failed reading n_sweeps",                                          //25
  "failed reading result_data",                                       //26
  "failed reading n_vals",                                            //27
  "failed reading newest",                                            //28
  "failed reading A_accumulated_average",                             //29
  "failed reading A_accumulated_variance",                            //30
  "failed reading A_data",                                            //31
  "failed reading B_accumulated_average",                             //32
  "failed reading B_accumulated_variance",                            //33
  "failed reading B_data",                                            //34
  "dim_A does not match dim_corr",                                    //35
  "unknown corr_operation",                                           //36
  "unknown compression operation for obs1",                           //37
  "unknown compressB has to be none for autocorrelation",             //38
  "unknown compression operation for obs2",                           //39
  "you need to specify obs2",                                         //40
  "dim_B from checkpoint != dim_B from observable",                   //41
  "no error"                                                          //end of the eror codes
};

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
  "dim_A and dim_corr do not match for conditional correlation", //18
  "unknown error, dim_corr should not be 0 at this point" //19
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
  int j, k;

// We calculate the correlation time for each dim_corr by normalizing the correlation,
// integrating it and finding out where C(tau)=tau;
  double C_tau;
  int ok_flag;
  for (j=0; j<self->dim_corr; j++) {
    correlation_time[j] = 0.; 
  }

  // here we still have to fix the stuff a bit!
  for (j=0; j<self->dim_corr; j++) {
    C_tau=1*self->dt;
    ok_flag=0;
    for (k=1; k<self->n_result-1; k++) {
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
                  char* compressA_name, char* compressB_name) {
  unsigned int i,j,k;
  unsigned int hierarchy_depth=0;

  // FIXME After shuffling consistency checks around the source file, some error messages may be still messed up!
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
  if ( abs(dt/time_step - round(dt/time_step)>1e-6 ) ) 
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
    hierarchy_depth=(int)ceil( 1 + log( (tau_max/dt)/(tau_lin-1) ) / log(2.0) );
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
  } else if ( strcmp(corr_operation_name,"complex_conjugate_product") == 0 ) {
    dim_corr = dim_A;
    self->corr_operation = &complex_conjugate_product;
  } else if ( strcmp(corr_operation_name,"square_distance_componentwise") == 0 ) {
    dim_corr = dim_A;
    self->corr_operation = &square_distance_componentwise;
// square_distance will be removed -- will be replaced by strides and blocks
//  } else if ( strcmp(corr_operation_name,"square_distance") == 0 ) {
//    self->corr_operation = &square_distance;
  } else if ( strcmp(corr_operation_name,"square_distance_cond_chain") == 0 ) {
    if ( (dim_A-1) % (dim_corr + 1) ) 
      return 18;
    dim_corr = (dim_A-1) / (dim_corr + 1) * dim_corr;
    self->corr_operation = &square_distance_cond_chain;
  } else if ( strcmp(corr_operation_name,"square_distance_cond") == 0 ) {
    if (dim_A % (dim_corr + 1) ) 
      return 18;
    dim_corr = dim_A / (dim_corr + 1) * dim_corr;
    self->corr_operation = &square_distance_cond;
  } else if ( strcmp(corr_operation_name,"scalar_product") == 0 ) {
    dim_corr=1;
    self->corr_operation = &scalar_product;
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
    if (!self->autocorrelation)
      (*self->compressB)(self->B[i][(self->newest[i]+1) % (self->tau_lin+1)],  
                       self->B[i][(self->newest[i]+2) % (self->tau_lin+1)], 
                       self->B[i+1][self->newest[i+1]],self->dim_B);
  }

  self->newest[0] = ( self->newest[0] + 1 ) % (self->tau_lin +1); 
  self->n_vals[0]++;

  if ( (*self->A_obs->fun)(self->A_obs->args, self->A[0][self->newest[0]], self->dim_A) != 0 )
    return 1;
  if (!self->autocorrelation)
    if ( (*self->B_obs->fun)(self->B_obs->args, self->B[0][self->newest[0]], self->dim_B) != 0 )
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
  int highest_level_to_compress;
  unsigned int index_new, index_old, index_res;
  int error;
  //int compress;
  int tau_lin=self->tau_lin;
  int hierarchy_depth=self->hierarchy_depth;

  double* temp = malloc(self->dim_corr*sizeof(double));

  // make a flag that the correlation is finalized
  self->finalized=1;

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

/* to be removed
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
*/


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

/** 
Compute the correlation if the condition is fulfiled, or zero otherwise

Assume that the observable array is composed of the following sub-units
arranged linearly one after another:
     condition; 
     position[ dim_A / (dim_A - dim_corr) ];

     correlate only values within the same block and if condition > 0
*/
int square_distance_cond ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) {
  const double tiny=0.00001; // to avoid roundoff errors
  unsigned int i, j, stepA, stepC, n_vals;
  double dist;
  n_vals = dim_A - dim_corr; // A contains the values in d dimensions plus condition per each value, corr contains just the values
  stepA=dim_A/n_vals;
  stepC=dim_corr/n_vals;
  for ( i = 0; i < n_vals; i++ ) { 
    // if both conditions are positive and the same and positive
    if ( A[i*stepA] > 0 &&  fabs( A[i*stepA] - B[i*stepA] ) < tiny  ) {
      for (j=0; j<stepC; j++) { 
        dist = A[i*stepA+j+1] - B[i*stepA+j+1];
        C[i*stepC+j] = dist*dist; 
      }
    } else { 
      for (j=0; j<stepC; j++) {
        C[i*stepC+j] = 0.0;
      }
    }
  }
  return 0;
}

/** 
Compute the correlation if the condition is fulfiled, or zero otherwise

Assume that the observable array is composed of the following sub-units
arranged linearly one after another:
     condition; 
     position;

     correlate only values within the same block and if condition > 0
     the last entry in the observable array is the chain length.
     All computed distances are folded to be <= chain_length/2;
     Warning: this function is very problem-specific.
     
*/
int square_distance_cond_chain ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) {
  const double tiny=0.00001; // to void roundoff errors in double->int conversion
  unsigned int i, imax;
  double dist, distmax, halfmax;
  distmax=(double)(A[dim_A-1]);
  halfmax=0.5*distmax;
  imax=(dim_A-1)/2;
  for ( i = 0; i < imax; i++ ) { 
    // if both conditions are positive and the same and positive
    if ( A[2*i] > 0 &&  fabs(A[2*i] - B[2*i]) < tiny  ) {
      dist= fabs(A[2*i+1] - B[2*i+1]);
      if (dist > halfmax) 
	dist -= distmax;
      C[i] = dist < tiny ? 0.0 : dist*dist;
    } else 
      C[i] = 0.0;
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

int double_correlation_write_checkpoint( double_correlation* self, char* filename) {
  FILE* f=0;
  int i, imax; // indexing variables
  f=fopen(filename, "w");
  if (!f) {
    return 1;
  }
  fprintf(f,"%d\t#autocorrelation \n", self->autocorrelation); // autocorrelation flag
  fprintf(f,"%d\t#finalized \n", self->finalized);       // non-zero if correlation is finialized
  fprintf(f,"%d\t#hierarchy_depth \n", self->hierarchy_depth); // maximum level of data compression
  fprintf(f,"%d\t#tau_lin \n", self->tau_lin);                 // number of frames in the linear correlation
  fprintf(f,"%d\t#dim_A \n", self->dim_A);                     // dimensionality of A
  fprintf(f,"%d\t#dim_B \n", self->dim_B);
  fprintf(f,"%d\t#dim_corr \n", self->dim_corr);
  fprintf(f,"%d\t#time (No of frames)\n", self->t);            // global time in number of frames
  fprintf(f,"%.6g\t#dt \n", self->dt);                         // time interval at which samples arrive
  fprintf(f,"%.6g\t#tau_max \n", self->tau_max);               // maximum time, for which the correlation should be calculated
  fprintf(f,"%d\t#update_frequency \n", self->update_frequency);// time distance between updates in MD timesteps 
  fprintf(f,"%d\t#window_distance \n", self->window_distance); 
  
  fprintf(f,"%d\t#n_result \n", self->n_result);                     // the total number of result values
  fprintf(f,"%d\t#n_data\n", self->n_data);               // a counter to calculated averages and variances

  fprintf(f,"%s\t#compressA_name\n", self->compressA_name);
  fprintf(f,"%s\t#compressB_name\n", self->compressB_name);
  fprintf(f,"%s\t#corr_operation_name\n", self->corr_operation_name);
  fprintf(f,"%d\t#is_from_file \n", self->is_from_file);
  fprintf(f,"%d\t#autoupate \n", self->autoupdate);
  fprintf(f,"%10g\t#last_update\n", self->last_update); // this requires high precision

  // print the buffers
  fprintf(f,"#tau values: \n"); // time differences                  
  for (i=0;i<self->n_result;i++)
    fprintf(f,"%d ", self->tau[i]); 
  fprintf(f,"\n#n_sweeps values: \n");                       
  for (i=0;i<self->n_result;i++)
    fprintf(f,"%d ", self->n_sweeps[i]); 
  fprintf(f,"\n#result_data values: \n");
  imax=self->n_result*self->dim_corr;
  for (i=0;i<imax;i++)
    fprintf(f,"%.6g ", self->result_data[i]); 
  fprintf(f,"\n#result_data values: \n");
  for (i=0;i<self->hierarchy_depth;i++)
    fprintf(f,"%d ", self->n_vals[i]); 
  fprintf(f,"\n#newest values: \n");
  for (i=0;i<self->hierarchy_depth;i++)
    fprintf(f,"%d ", self->newest[i]); 

//  fprintf(f,"%\t#comment \n", double*** A;                     // input quantity 1
//  fprintf(f,"%\t#comment \n", double*** B;                     // input quantity 2
//  fprintf(f,"%\t#comment \n", double** result;                // output quantity

  fprintf(f,"\n#A_accumulated_average:\n");
  for (i=0;i<self->dim_A;i++)  
      fprintf(f,"%.6g ", self->A_accumulated_average[i]);     // all A values are added up here
  fprintf(f,"\n#A_accumulated_variance:\n");                       
  for (i=0;i<self->dim_A;i++)  
    fprintf(f,"%.6g ", self->A_accumulated_variance[i]);    // all A**2 values are added up here
  fprintf(f,"\n#A_data:\n");
  imax=(self->tau_lin+1)*self->hierarchy_depth*self->dim_A;
  for (i=0;i<imax;i++)  
    fprintf(f,"%.6g ", self->A_data[i]);

  if (!self->autocorrelation) { 
      fprintf(f,"\n#B_accumulated_average:\n");
      for (i=0;i<self->dim_B;i++)  
          fprintf(f,"%.6g ", self->B_accumulated_average[i]);     // all B values are added up here 
      fprintf(f,"\n#B_accumulated_variance:\n");                       
      for (i=0;i<self->dim_B;i++)  
          fprintf(f,"%.6g ", self->B_accumulated_variance[i]);    // all B**2 values are added up here 
      fprintf(f,"\n#B_data:\n"); 
      imax=(self->tau_lin+1)*self->hierarchy_depth*self->dim_B; 
      for (i=0;i<imax;i++)  
          fprintf(f,"%.6g ", self->B_data[i]); 
  }
  fclose(f);
  return TCL_OK;
}

int double_correlation_init_from_checkpoint(double_correlation* self, char* filename, int dim_A, int dim_B, observable *A, observable *B) {
  FILE* f=0;
  char *buffer; // make sure that no names are longer than this!
  int i, imax; // indexing variables
  int j, jmax; // indexing variables
  //int tmp;
  f=fopen(filename, "r");
  if (!f) 
    return 2;
  if (fscanf(f,"%d", &self->autocorrelation) < 1) // autocorrelation flag
    return 3;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->finalized))       // non-zero if correlation is finialized
    return 4;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->hierarchy_depth)) // maximum level of data compression
    return 5;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->tau_lin) )                 // number of frames in the linear correlation
    return 6;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->dim_A) )                     // dimensionality of A
    return 7;
  //fprintf(stderr, "\n***\ndim_A: %d, self->dim_A: %d\n***n\n",dim_A, self->dim_A);
  if (dim_A != self->dim_A) 
    return 8;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->dim_B) )
    return 9;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->dim_corr) )
    return 10;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->t) )            // global time in number of frames
    return 11;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%lf", &self->dt) )                         // time interval at which samples arrive
    return 12;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%lf", &self->tau_max) )               // maximum time, for which the correlation should be calculated
    return 13;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->update_frequency) )// time distance between updates in MD timesteps 
    return 14;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->window_distance) ) 
    return 15;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->n_result) )                     // the total number of result values
    return 16;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->n_data) )               // a counter to calculated averages and variances
    return 17;
  read_until_end_of_line(f);
  buffer=(char*)malloc(256*sizeof(char));
  if (1 > fscanf(f,"%s ", buffer) )
    return 18;
  self->compressA_name=strdup(buffer);
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%s", buffer) )
    return 19;
  self->compressB_name=strdup(buffer);
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%s", buffer) )
    return 20;
  self->corr_operation_name=strdup(buffer);
  free((void*)buffer);
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->is_from_file) )
    return 21;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%d", &self->autoupdate) )
    return 22;
  read_until_end_of_line(f);
  if (1 > fscanf(f,"%lf", &self->last_update) ) // this requires high precision
      return 23;
  read_until_end_of_line(f);
  
  // read the buffers
  imax=self->n_result;
  self->tau=(int*)malloc(imax*sizeof(int));
  for (i=0;i<imax;i++) { 
    if (1 > fscanf(f,"%d", &self->tau[i]) ) 
      return 24;
  }
  read_until_end_of_line(f);
  imax=self->n_result;
  self->n_sweeps=(unsigned int*)malloc(imax*sizeof(unsigned int));
  for (i=0;i<imax;i++) { 
    if (1 > fscanf(f,"%u", &self->n_sweeps[i]) )
      return 25;
  }
  read_until_end_of_line(f);
  imax=self->n_result*self->dim_corr;
  self->result_data=(double*)malloc(imax*sizeof(double));
  for (i=0;i<imax;i++) {
    if (1 > fscanf(f,"%lf", &self->result_data[i]) )
      return 26;
  }
  read_until_end_of_line(f);
  imax=self->hierarchy_depth;
  self->n_vals=(unsigned int*)malloc(imax*sizeof(unsigned int));
  for (i=0;i<imax;i++) { 
    if (1 > fscanf(f,"%u", &self->n_vals[i]) )
      return 27;
  }
  read_until_end_of_line(f);
  self->newest=(unsigned int*)malloc(self->hierarchy_depth*sizeof(unsigned int));
  for (i=0;i<self->hierarchy_depth;i++) {
    if (1 > fscanf(f,"%u ", &self->newest[i]) )
      return 28;
  }
  read_until_end_of_line(f);
  imax=self->dim_A;
  self->A_accumulated_average=(double*)malloc(imax*sizeof(double));
  for (i=0;i<imax;i++)  {
      if (1 > fscanf(f,"%lf", &self->A_accumulated_average[i]) )     // all A values are added up here
        return 29;
  }
  read_until_end_of_line(f);
  imax=self->dim_A;
  self->A_accumulated_variance=(double*)malloc(imax*sizeof(double));
  for (i=0;i<imax;i++)  {
      if (1 > fscanf(f,"%lf", &self->A_accumulated_variance[i]) )     // all A values are added up here
        return 30;
  }
  read_until_end_of_line(f);
  imax=(self->tau_lin+1)*self->hierarchy_depth*self->dim_A;
  self->A_data=(double*)malloc(imax*sizeof(double));
  for (i=0;i<imax;i++) { 
    if (1 > fscanf(f,"%lf", &self->A_data[i]) )
      return 31;
  }
  if (self->autocorrelation) { 
    self->B_accumulated_average=self->A_accumulated_average;
    self->B_accumulated_variance=self->A_accumulated_variance;
    self->B_data=self->A_data;
  } else {
    read_until_end_of_line(f);
    imax=self->dim_B;
    self->B_accumulated_average=(double*)malloc(imax*sizeof(double));
    for (i=0;i<imax;i++)  {
        if (1 > fscanf(f,"%lf", &self->B_accumulated_average[i]) )     // all B values are added up here
          return 32;
    }
    read_until_end_of_line(f);
    imax=self->dim_B;
    self->B_accumulated_variance=(double*)malloc(imax*sizeof(double));
    for (i=0;i<imax;i++)  {
        if (1 > fscanf(f,"%lf", &self->B_accumulated_variance[i]) )     
          return 33;
    }
    read_until_end_of_line(f);
    imax=(self->tau_lin+1)*self->hierarchy_depth*self->dim_B;
    self->B_data=(double*)malloc(imax*sizeof(double));
    for (i=0;i<imax;i++) { 
      if (1 > fscanf(f,"%lf", &self->B_data[i]) )
        return 34;
    }
  }
 
  // finally, get the function pointers
  // choose the correlation operation 
  if ( strcmp(self->corr_operation_name,"componentwise_product") == 0 ) {
    self->corr_operation = &componentwise_product;
  } else if ( strcmp(self->corr_operation_name,"complex_conjugate_product") == 0 ) {
    self->corr_operation = &complex_conjugate_product;
  } else if ( strcmp(self->corr_operation_name,"square_distance_componentwise") == 0 ) {
    self->corr_operation = &square_distance_componentwise;
// square_distance will be removed -- will be replaced by strides and blocks
//  } else if ( strcmp(corr_operation_name,"square_distance") == 0 ) {
//    self->corr_operation = &square_distance;
  } else if ( strcmp(self->corr_operation_name,"square_distance_cond_chain") == 0 ) {
    if ( (dim_A-1) % (self->dim_corr + 1) ) 
      return 35;
    self->corr_operation = &square_distance_cond_chain;
  } else if ( strcmp(self->corr_operation_name,"square_distance_cond") == 0 ) {
    self->corr_operation = &square_distance_cond;
  } else if ( strcmp(self->corr_operation_name,"scalar_product") == 0 ) {
    self->corr_operation = &scalar_product;
  } else {
    return 36; 
  }
  // Choose the compression function
  if ( strcmp(self->compressA_name,"discard2") == 0 ) {
    self->compressA=&compress_discard2;
  } else if ( strcmp(self->compressA_name,"discard1") == 0 ) {
    self->compressA=&compress_discard1;
  } else if ( strcmp(self->compressA_name,"linear") == 0 ) {
    self->compressA=&compress_linear;
  } else {
    return 37;
  }
  
  if ( strcmp(self->compressB_name,"none") == 0 ) {
      self->compressB=&compress_do_nothing;
  } else if ( self->autocorrelation ) {
    return 38;
  } else if ( strcmp(self->compressB_name,"discard2") == 0 ) {
    self->compressB=&compress_discard2;
  } else if ( strcmp(self->compressB_name,"discard1") == 0 ) {
    self->compressB=&compress_discard1;
  } else if ( strcmp(self->compressB_name,"linear") == 0 ) {
    self->compressB=&compress_linear;
  } else {
    return 39;
  }
  
  // set the observable pointer
  self->A_obs = A;
  // set the observable pointer
  if (B == 0 && !self->autocorrelation)
      return 40;
  self->B_obs = B;

  // finally restore the convenience pointers
  // to obs1 buffer
  imax=self->hierarchy_depth;
  jmax=self->tau_lin + 1;
  self->A = (double***)malloc(imax*sizeof(double**));
  for (i=0;i<imax;i++) 
      self->A[i] = (double**)malloc(jmax*sizeof(double*));
  for (i=0;i<imax;i++) 
      for (j=0;j<imax;j++) 
        self->A[i][j] = &self->A_data[(i*(self->tau_lin+1))*self->dim_A+j*self->dim_A];
  // and to the result
  imax=self->n_result;
  self->result  = (double**) malloc(imax*sizeof(double*));
  for (i=0; i<imax; i++)
    self->result[i]=&self->result_data[i*self->dim_corr];
    
  if ( self->autocorrelation ) {
    self->B = self->A; 
  } else {
    if (dim_B != self->dim_B) 
      return 41;
    // convenience pointer to obs2 buffer
    imax=self->hierarchy_depth;
    jmax=self->tau_lin + 1;
    self->B = (double***)malloc(imax*sizeof(double**));
    for (i=0;i<imax;i++) 
      self->B[i] = (double**)malloc(jmax*sizeof(double*));
    for (i=0;i<imax;i++) 
      for (j=0;j<imax;j++) 
        self->B[i][j] = &self->B_data[(i*(self->tau_lin+1))*self->dim_B+j*self->dim_B];
   
  }
  
  fclose(f);
  return TCL_OK;
}

int read_until_end_of_line(FILE *f) {
    char c;
    c=fgetc(f);
    while (c!='\n' && c!=EOF)  {
        c=fgetc(f); 
    }
    if (c=='\n') {
        c=fgetc(f);
        if (c=='#')
            return read_until_end_of_line(f);
        else {
            ungetc(c,f); 
            return 0;
        }
    } 
    else
        return 1;
}
