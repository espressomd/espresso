 /*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  
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
#include "Correlator.hpp"
#include "particle_data.hpp"
#include <cstring>
#include "integrate.hpp"
#include "utils.hpp" 

#include <limits>

namespace Correlators {

/* global variables */


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
//    return correlations.at(no)->get_data();
}

int correlation_update_from_file(unsigned int no) {
//  if (!correlations.at(no)->is_from_file)
//    return 1;
//  while ( ! correlations.at(no)->get_data()) {
//  }
//  return 0;
}



int Correlator::get_correlation_time(double* correlation_time) {
  // We calculate the correlation time for each dim_corr by normalizing the correlation,
  // integrating it and finding out where C(tau)=tau;
  double C_tau;
  int ok_flag;
  for (unsigned j=0; j<dim_corr; j++) {
    correlation_time[j] = 0.; 
  }

  // here we still have to fix the stuff a bit!
  for (unsigned j=0; j<dim_corr; j++) {
    C_tau=1*dt;
    ok_flag=0;
    for (unsigned k=1; k<n_result-1; k++) {
      if (n_sweeps[k]==0)
        break;
      C_tau+=(result[k][j]/ (double) n_sweeps[k] - A_accumulated_average[j]*B_accumulated_average[j]/n_data/n_data)/(result[0][j]/n_sweeps[0])*dt*(tau[k]-tau[k-1]);
//      C_tau+=(self->result[k][j]/ (double) self->n_sweeps[k] - self->A_accumulated_average[j]*self->B_accumulated_average[j]/self->n_data/self->n_data)/(self->result[0][j]/self->n_sweeps[0])*self->dt*(self->tau[k]-self->tau[k-1]);

//        if (C_tau < i*tau[k]*dt) {
        if (exp(-tau[k]*dt/C_tau)+2*sqrt(tau[k]*dt/n_data)
            >exp(-tau[k-1]*dt/C_tau)+2*sqrt(tau[k-1]*dt/n_data)) {
        correlation_time[j]=C_tau*(1+(2*(double)tau[k]+1)/(double)n_data);
//          printf("stopped at tau=>%f\n", tau[k]*dt);
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


Correlator::Correlator() :
			    t(0), finalized(0), autoupdate(0),autocorrelation(1),initialized(0),correlation_args{}
           {}

void Correlator::initialize() {
  unsigned int i,j,k;
  hierarchy_depth=0;
  // Class members are assigned via the initializer list

  
  // Input validation
  if (dt <= 0) {
    throw std::runtime_error( init_errors[2]);
  }

  if ((dt-time_step)<-1e-6*time_step) {
    throw std::runtime_error( init_errors[15]);
  }

  // check if dt is a multiple of the md timestep
  if ( std::abs(dt/time_step - round(dt/time_step)) > 1e-6 ) {
    throw std::runtime_error( init_errors[16]);
  }

  // Time steps and intervals
  update_frequency = std::floor(dt/time_step + std::numeric_limits<double>::round_error());

  if ( tau_lin == 1 ) { // use the default
    tau_lin=(int)ceil(tau_max/dt);
    if (tau_lin%2) tau_lin+=1;
    printf("tau_lin: %d\n", tau_lin);
  }

  if (tau_lin<2) {
    throw std::runtime_error( init_errors[3]);
  }

  if (tau_lin%2) {
    throw std::runtime_error( init_errors[14]);
  }

  
  if (tau_max <= dt) { 
    throw std::runtime_error( init_errors[4]);

  } else { //set hierarchy depth which can  accommodate at least tau_max
    if ( (tau_max/dt) < tau_lin ) {
      hierarchy_depth = 1;
    } else {
      hierarchy_depth=(unsigned int)ceil( 1 + log( (tau_max/dt)/(tau_lin-1) ) / log(2.0) );
    }
  }

  

  
  dim_A=0;
  dim_B=0;
  if (A_obs) dim_A=A_obs->n_values();
  if (B_obs.get()) dim_B=B_obs->n_values();

  if (dim_A<1) {
    throw std::runtime_error( init_errors[6]);
  }

  autocorrelation=1;
  if (dim_B==0) {
    dim_B = dim_A; 
  } else if (dim_B>0) {
    autocorrelation=0;
  } else {
    throw std::runtime_error( init_errors[7]);
  }

  

  // choose the correlation operation 
  if (corr_operation_name=="") { 
    throw std::runtime_error( init_errors[11]); // there is no reasonable default
  } else if ( corr_operation_name=="componentwise_product")  {
    dim_corr = dim_A;
    corr_operation = &componentwise_product;
    correlation_args = Vector3d{0,0,0};
  } else if ( corr_operation_name=="complex_conjugate_product")  {
    dim_corr = dim_A;
    corr_operation = &complex_conjugate_product;
    correlation_args = Vector3d{0,0,0};
  } else if ( corr_operation_name=="tensor_product") {
    dim_corr = dim_A*dim_B;
    corr_operation = &tensor_product;
    correlation_args = Vector3d{0,0,0};
  } else if ( corr_operation_name=="square_distance_componentwise") {
    dim_corr = dim_A;
    corr_operation = &square_distance_componentwise;
    correlation_args = Vector3d{0,0,0};
  } else if ( corr_operation_name=="fcs_acf")  {
	// note: user provides w=(wx,wy,wz) but we want to use wsquare=(wx^2,wy^2,wz^2)
    if(correlation_args[0] <= 0 || correlation_args[1] <= 0 || correlation_args[2] <= 0) {
		throw std::runtime_error("missing parameter for fcs_acf: w_x w_y w_z");
	}
	correlation_args[0]=correlation_args[0]*correlation_args[0];
    correlation_args[1]=correlation_args[1]*correlation_args[1];
    correlation_args[2]=correlation_args[2]*correlation_args[2];
	fprintf(stderr,"args2: %f %f %f\n",correlation_args[0],correlation_args[1],correlation_args[2]);
	if (dim_A %3 )
       throw std::runtime_error( init_errors[18]);
    dim_corr = dim_A/3;
    corr_operation = &fcs_acf;
// square_distance will be removed -- will be replaced by strides and blocks
//  } else if ( strcmp(corr_operation_name,"square_distance") == 0 ) {
//    corr_operation = &square_distance;
  } else if ( corr_operation_name=="scalar_product")  {
    dim_corr=1;
    corr_operation = &scalar_product;
    correlation_args = Vector3d{0,0,0};
  } else {
    throw std::runtime_error( init_errors[11]); 
  }
  
  
  // Choose the compression function
  if (compressA_name=="") { // this is the default
    compressA_name=strdup("discard2");
    compressA=&compress_discard2;
  } else if ( compressA_name=="discard2")  {
    compressA=&compress_discard2;
  } else if ( compressA_name=="discard1") {
    compressA=&compress_discard1;
  } else if ( compressA_name=="linear")  {
    compressA=&compress_linear;
  } else {
    throw std::runtime_error( init_errors[12]);
  }
  
  if (compressB_name=="") { 
    if(autocorrelation) { // the default for autocorrelation
      compressB_name=strdup("none"); 
      compressB=&compress_do_nothing;
    } else { // the default for corsscorrelation
      compressB_name=compressA_name;
      compressB=compressA;
    } 
  } else if ( autocorrelation ) {
    throw std::runtime_error( init_errors[17]);
  } else if ( compressB_name=="discard2")  {
    compressB=&compress_discard2;
  } else if ( compressB_name=="discard1")  {
    compressB=&compress_discard1;
  } else if ( compressB_name=="linear") {
    compressB=&compress_linear;
  } else {
    throw std::runtime_error( init_errors[13]);
  }
 
//  if (A_fun == &file_data_source_readline && (B_fun == &file_data_source_readline|| autocorrelation)) {
//    is_from_file = 1;
//  } else {
//    is_from_file = 0;
//  }

  
  // Memmory allocation
  A_data = (double*)Utils::malloc((tau_lin+1)*hierarchy_depth*dim_A*sizeof(double));
  if (autocorrelation) 
    B_data = A_data;
  else 
    B_data = (double*)Utils::malloc((tau_lin+1)*hierarchy_depth*dim_B*sizeof(double));
  
  n_data=0;
  A_accumulated_average = (double*)Utils::malloc(dim_A*sizeof(double));
  A_accumulated_variance= (double*)Utils::malloc(dim_A*sizeof(double));
  for (k=0; k<dim_A; k++) {
    A_accumulated_average[k]=0;
    A_accumulated_variance[k]=0;
  }
  if (autocorrelation) {
    B_accumulated_average =  A_accumulated_average;
    B_accumulated_variance = B_accumulated_variance;
  } else {
    B_accumulated_average = (double*)Utils::malloc(dim_B*sizeof(double));
    B_accumulated_variance = (double*)Utils::malloc(dim_B*sizeof(double));
    for (k=0; k<dim_B; k++) {
      B_accumulated_average[k]=0;
      B_accumulated_variance[k]=0;
    }
  }


  n_result=tau_lin+1 + (tau_lin+1)/2*(hierarchy_depth-1);
  tau = (int*)                Utils::malloc(n_result*sizeof(int));
  n_sweeps = (unsigned int*)  Utils::malloc(n_result*sizeof(int));
  result_data  = (double*)    Utils::malloc(n_result*dim_corr*sizeof(double));
  n_vals = (unsigned int*) Utils::malloc(hierarchy_depth*sizeof(unsigned int));

  // allocate space for convenience pointer to A and B buffers
  A = (double***)Utils::malloc(hierarchy_depth*sizeof(double**));
  if(autocorrelation) B = A;
  else B = (double***)Utils::malloc(hierarchy_depth*sizeof(double**));

  for (i=0; i<hierarchy_depth; i++) {
    A[i] = (double**) Utils::malloc((tau_lin+1)*sizeof(double*));
    if(!autocorrelation) B[i] = (double**) Utils::malloc((tau_lin+1)*sizeof(double*));
  }

  // and initialize the values
  for (i=0; i<hierarchy_depth; i++) {
    n_vals[i]=0;
    for (j=0; j<tau_lin+1; j++) {
      A[i][j] = &A_data[(i*(tau_lin+1))*dim_A+j*dim_A];
      for (k=0; k<dim_A; k++) 
        A[i][j][k] = 0.;
      if(!autocorrelation) {
        B[i][j] = &B_data[(i*(tau_lin+1))*dim_B+j*dim_B];
        for (k=0; k<dim_B; k++) 
          B[i][j][k] = 0.;
      }
    }
  }
  
  if (A == 0) {
    throw std::runtime_error( init_errors[9]);
  }

  if (B == 0 && !autocorrelation) {
    throw std::runtime_error( init_errors[10]);
  }


  // allocate space for convenience pointer to the result
  result  = (double**)        Utils::malloc(n_result*sizeof(double*));
  for (i=0; i<n_result; i++) {
    n_sweeps[i]=0;
    result[i]=&result_data[i*dim_corr];
    for (j=0; j<dim_corr; j++) 
      // and initialize the values
      result[i][j]=0;
  }

  newest = (unsigned int *)Utils::malloc(hierarchy_depth*sizeof(unsigned int));
  for ( i = 0; i<hierarchy_depth; i++ ) {
    newest[i]= tau_lin;
  }
  for (i=0; i < tau_lin+1; i++) {
    tau[i] = i;
  }
  for (j=1; j < hierarchy_depth; j++)
    for (k=0; k < tau_lin/2; k++) {
      tau[tau_lin + 1 + (j-1)*tau_lin/2+k] = (k+(tau_lin/2)+1)*(1<<j); 
    }


   initialized=1;
}

int Correlator::get_data() {
  // We must now go through the hierarchy and make sure there is space for the new 
  // datapoint. For every hierarchy level we have to decide if it necessary to move 
  // something
  int i,j;
  int highest_level_to_compress;
  unsigned int index_new, index_old, index_res;
  int error;
  
  t++;

  highest_level_to_compress=-1;
  i=0;
  j=1;
  // Lets find out how far we have to go back in the hierarchy to make space for the new value
  while (1) {
    if ( ( (t - ((tau_lin + 1)*((1<<(i+1))-1) + 1) )% (1<<(i+1)) == 0) ) {
      if ( i < (int(hierarchy_depth) - 1) && n_vals[i]> tau_lin) {

        highest_level_to_compress+=1;
        i++;
      } else break;
    } else break;
  }

  // Now we know we must make space on the levels 0..highest_level_to_compress
  // Now lets compress the data level by level.

  for ( i = highest_level_to_compress; i >= 0; i-- ) {
    // We increase the index indicating the newest on level i+1 by one (plus folding)
    newest[i+1] = (newest[i+1] + 1) % (tau_lin+1);
    n_vals[i+1]+=1;
//    printf("t %d compressing level %d no %d and %d to level %d no %d, nv %d\n",t, i, (newest[i]+1) % (tau_lin+1),
//(newest[i]+2) % (tau_lin+1), i+1, newest[i+1], n_vals[i]);
    (*compressA)(A[i][(newest[i]+1) % (tau_lin+1)],  
                       A[i][(newest[i]+2) % (tau_lin+1)], 
                       A[i+1][newest[i+1]],dim_A);
    if (!autocorrelation)
      (*compressB)(B[i][(newest[i]+1) % (tau_lin+1)],  
                       B[i][(newest[i]+2) % (tau_lin+1)], 
                       B[i+1][newest[i+1]],dim_B);
  }

  newest[0] = ( newest[0] + 1 ) % (tau_lin +1); 
  n_vals[0]++;

  if ( A_obs->calculate() != 0 )
    return 1;
  // copy the result:
  memmove(A[0][newest[0]], &(A_obs->last_value[0]), dim_A*sizeof(double));

  if (!autocorrelation) {
    if ( B_obs->calculate() != 0 )
      return 2;
    memmove(B[0][newest[0]], &(B_obs->last_value[0]), dim_B*sizeof(double));
  }

  // Now we update the cumulated averages and variances of A and B
  n_data++;
  for (unsigned k=0; k<dim_A; k++) {
    A_accumulated_average[k]+=A[0][newest[0]][k];
    A_accumulated_variance[k]+=A[0][newest[0]][k]*A[0][newest[0]][k];
  }
  // Here we check if it is an autocorrelation
  if (!autocorrelation) {
    for (unsigned k=0; k<dim_B; k++) {
      B_accumulated_average[k]+=B[0][newest[0]][k];
      B_accumulated_variance[k]+=B[0][newest[0]][k]*B[0][newest[0]][k];
    }
  } 

  double* temp = (double*)Utils::malloc(dim_corr*sizeof(double));
  if (!temp)
    return 4;
// Now update the lowest level correlation estimates
  for ( j = 0; j < int(MIN(tau_lin+1, n_vals[0]) ); j++) {
    index_new = newest[0];
    index_old =  (newest[0] - j + tau_lin + 1) % (tau_lin + 1);
//    printf("old %d new %d\n", index_old, index_new);
    error = (corr_operation)(A[0][index_old], dim_A, B[0][index_new], dim_B, temp, dim_corr, correlation_args);
    if ( error != 0)
      return error;
    n_sweeps[j]++;
    for (unsigned k = 0; k < dim_corr; k++) {
      result[j][k] += temp[k];
    }
  }
// Now for the higher ones
  for ( int i = 1; i < highest_level_to_compress+2; i++) {
    for ( unsigned j = (tau_lin+1)/2+1; j < MIN(tau_lin+1, n_vals[i]); j++) {
      index_new = newest[i];
      index_old = (newest[i] - j + tau_lin + 1) % (tau_lin + 1);
      index_res = tau_lin + (i-1)*tau_lin/2 + (j - tau_lin/2+1) -1;
      error=(corr_operation)(A[i][index_old], dim_A, B[i][index_new], dim_B, temp, dim_corr, correlation_args);
      if ( error != 0)
        return error;
      n_sweeps[index_res]++;
      for (unsigned k = 0; k < dim_corr; k++) {
        result[index_res][k] += temp[k];
      }
    }
  }
  free(temp);
  last_update=sim_time;
  return 0;
}


void write_double(FILE * fp, const double * data, unsigned int n, bool binary) {
  if (binary){
    fwrite(data,sizeof(double),n,fp);
  } else {
    for(unsigned int i=0;i<n-1;i++){
      fprintf(fp,"%e\t",data[i]);
    }
    fprintf(fp,"%e\n",data[n-1]);
  }
}
void write_uint(FILE * fp, const unsigned int * data, unsigned int n, bool binary){
  if (binary){
    fwrite(data,sizeof(unsigned int),n,fp);
  } else {
    for(unsigned int i=0;i<n-1;i++){
      fprintf(fp,"%u\t",data[i]);
    }
    fprintf(fp,"%u\n",data[n-1]);
  }
}
int Correlator::write_data_to_file(const char * filename, bool binary) const {
  FILE* fp=0;
  fp=fopen(filename, "w");
  if (!fp) {
    return 1;
  }
  for (unsigned int i=0; i<hierarchy_depth; i++) {
    for (unsigned int j=0; j<tau_lin+1; j++) {
      write_double(fp,A[i][j],dim_A,binary);
    }
  }
  if (!autocorrelation){
    for (unsigned int i=0; i<hierarchy_depth; i++) {
      for (unsigned int j=0; j<tau_lin+1; j++) {
	write_double(fp,B[i][j],dim_B,binary);
      }
    } 
  }
  for (unsigned int i=0; i<n_result; i++) {
    write_double(fp,result[i],dim_corr,binary);
  }
  write_uint(fp,n_sweeps,n_result       ,binary);
  write_uint(fp,n_vals  ,hierarchy_depth,binary);
  write_uint(fp,newest  ,hierarchy_depth,binary);
  
  write_double(fp,A_accumulated_average ,dim_A,binary);
  write_double(fp,A_accumulated_variance,dim_A,binary);
  if (!autocorrelation){
    write_double(fp,B_accumulated_average ,dim_B,binary);
    write_double(fp,B_accumulated_variance,dim_B,binary);
  }
  write_uint(fp,&(n_data),1,binary);
  write_uint(fp,&(t     ),1,binary);
  write_double(fp,&(last_update),1,binary);
  fclose(fp);
  return 0;
}

int read_double(FILE * fp, double * data, unsigned int n, bool binary){
  if (binary){
    size_t tmp = fread(data,sizeof(double),n,fp);
    if (tmp < n){
      return 1;
    }
  } else {
    for(unsigned int i=0;i<n;i++){
      if (fscanf(fp,"%le",data+i) == 0 ) return 1;
    }
  }
  return 0;
}
int read_uint(FILE * fp, unsigned int * data, unsigned int n, bool binary){
  if (binary){
    size_t tmp = fread(data,sizeof(unsigned int),n,fp);
    if (tmp < n){
      return 1;
    }
  } else {
    for(unsigned int i=0;i<n;i++){
      if (fscanf(fp,"%u",data+i) == 0)
	return 1;
    }
  }
  return 0;
}
int Correlator::read_data_from_file(const char * filename, bool binary){
  FILE* fp=0;
  fp=fopen(filename, "r");
  if (!fp) {
    return 2;
  }
  for (unsigned int i=0; i<hierarchy_depth; i++) {
    for (unsigned int j=0; j<tau_lin+1; j++) {
      if (read_double(fp,A[i][j],dim_A,binary)) return 1;
    }
  }
  if (!autocorrelation){
    for (unsigned int i=0; i<hierarchy_depth; i++) {
      for (unsigned int j=0; j<tau_lin+1; j++) {
	if (read_double(fp,B[i][j],dim_B,binary)) return 1;
      }
    } 
  }
  for (unsigned int i=0; i<n_result; i++) {
    if (read_double(fp,result[i],dim_corr,binary))return 1;
  }
  if (read_uint(fp,n_sweeps,n_result       ,binary))return 1;
  if (read_uint(fp,n_vals  ,hierarchy_depth,binary))return 1;
  if (read_uint(fp,newest  ,hierarchy_depth,binary))return 1;
  
  if (read_double(fp,A_accumulated_average ,dim_A,binary))return 1;
  if (read_double(fp,A_accumulated_variance,dim_A,binary))return 1;
  if (!autocorrelation){
    if (read_double(fp,B_accumulated_average ,dim_B,binary))return 1;
    if (read_double(fp,B_accumulated_variance,dim_B,binary))return 1;
  }
  if (read_uint(fp,&(n_data),1,binary))return 1;
  if (read_uint(fp,&(t     ),1,binary))return 1;
  if (read_double(fp,&(last_update),1,binary))return 1;
  fclose(fp);
  return 0;
}


int Correlator::finalize() {
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
  unsigned tau_lin=tau_lin;
  int hierarchy_depth=hierarchy_depth;

  double* temp = (double*)Utils::malloc(dim_corr*sizeof(double));

  // make a flag that the correlation is finalized
  finalized=1;

  if (!temp)
    return 4;
  //printf ("tau_lin:%d, hierarchy_depth: %d\n",tau_lin,hierarchy_depth); 
  //for(ll=0;ll<hierarchy_depth;ll++) printf("n_vals[l=%d]=%d\n",ll, n_vals[ll]);
  for(ll=0;ll<hierarchy_depth-1;ll++) {
    if (n_vals[ll] > tau_lin+1 ) vals_ll = tau_lin + n_vals[ll]%2;
    else vals_ll=n_vals[ll];
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
        //printf("test level %d for compression, n_vals=%d ... ",i,n_vals[i]);
        if ( n_vals[i]%2 ) { 
	  if ( i < (hierarchy_depth-1) && n_vals[i]> tau_lin) { 
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
        newest[i+1] = (newest[i+1] + 1) % (tau_lin+1); 
        n_vals[i+1]+=1; 
        //printf("compressing level %d no %d and %d to level %d no %d, nv %d\n",i, (newest[i]+1) % (tau_lin+1), (newest[i]+2) % (tau_lin+1), i+1, newest[i+1], n_vals[i]); 
        (*compressA)(A[i][(newest[i]+1) % (tau_lin+1)],  
	                   A[i][(newest[i]+2) % (tau_lin+1)], 
                           A[i+1][newest[i+1]],dim_A);
        (*compressB)(B[i][(newest[i]+1) % (tau_lin+1)],  
                           B[i][(newest[i]+2) % (tau_lin+1)], 
                           B[i+1][newest[i+1]],dim_B);
      } 
      newest[ll] = (newest[ll] + 1) % (tau_lin+1); 

      // We only need to update correlation estimates for the higher levels
      for ( i = ll+1; i < highest_level_to_compress+2; i++) {
        for ( j = (tau_lin+1)/2+1; j < int(MIN(tau_lin+1, n_vals[i])); j++) {
          index_new = newest[i];
          index_old = (newest[i] - j + tau_lin + 1) % (tau_lin + 1);
          index_res = tau_lin + (i-1)*tau_lin/2 + (j - tau_lin/2+1) -1;
          error=(corr_operation)(A[i][index_old], dim_A, B[i][index_new], dim_B, temp, dim_corr, correlation_args);
          if ( error != 0)
            return error;
          n_sweeps[index_res]++;
          for (unsigned k = 0; k < dim_corr; k++) {
            result[index_res][k] += temp[k];
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


void Correlator::start_auto_update() {
      if(update_frequency > 0) {
        correlations_autoupdate = 1;
        autoupdate=1;
        last_update=sim_time;
      } else {
        throw std::runtime_error("Could not start autoupdate: update frequency not set");
      }
}

void Correlator::stop_auto_update() {
  autoupdate=0;
  // Todo
  // Insert logic to determine if global correlations_auto_update can be set to 0
}


std::vector<double> Correlator::get_correlation() {
       std::vector<double> res;
       
       // time + n_sweeps + corr_1...corr_n
       int cols=2+dim_corr;
       res.resize(n_result*cols);

       for (int i=0;i<n_result; i++) {
         res[cols*i +0]=tau[i]*dt;
         res[cols*i +1]=n_sweeps[i];
         for (int k=0;k<dim_corr;k++) {
           if (n_sweeps[i]>0) {
             res[cols*i+2+k] =result[i][k]/n_sweeps[i];
           } 
           else { 
             res[cols*i+2+k] =0;
           }
         }
       }
       return res;
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

int scalar_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, Vector3d ) {
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

int componentwise_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, Vector3d ) {
  unsigned int i;
  if (!(dim_A == dim_B )) {
    printf("Error in componentwise product: The vector sizes do not match");
    return 5;
  } 
  if (!(dim_A == dim_corr )) {
    printf("Error in componentwise product: The vector sizes do not match");
    return 5;
  }
  for ( i = 0; i < dim_A; i++ ) {
    C[i] = A[i]*B[i];
  }
  return 0;
}

int complex_conjugate_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, Vector3d ) {
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

int tensor_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, Vector3d ) {
  unsigned int i,j;
  for ( i = 0; i < dim_A; i++ )
  {
    for ( j = 0; j < dim_B; j++ )
    {
      C[i*dim_B + j] = A[i]*B[j];
    }
  }
  return 0;
}

int square_distance_componentwise ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, Vector3d ) {
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

// note: the argument name wsquare denotes that it value is w^2 while the user sets w
int fcs_acf ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr, Vector3d wsquare) {
  if (!(dim_A == dim_B )) {
    return 5;
  }
  if ( dim_A / dim_corr != 3) {
    return 6; 
  }
  for (unsigned i = 0; i < dim_corr; i++ ) 
    C[i] = 0;
  for (unsigned i = 0; i < dim_A; i++ ) {
    C [i/3] -= ( (A[i]-B[i])*(A[i]-B[i]) ) / wsquare[i%3];
  }
  for (unsigned i = 0; i < dim_corr; i++ ) 
    C[i] = exp(C[i]);
  return 0;
}



} // Namespace Correlators



