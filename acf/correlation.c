// FIXME counting n_sweps is currently wrong
// FIXME last hierarchy level may be incomplete and in most cases it will be, unless we have exactly as many data points as to match our blocking scheme
// FIXME tau_lin may be confusig as it is the number of entries, and not their difference in time units; similar for tau_max
// Double_correlation = correlation of Real values of double precision

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef struct {
  unsigned int hierarchy_depth; // maximum level of data compression
  unsigned int dim_A; // dimensionality of A
  unsigned int dim_B;
  unsigned int dim_C;
  unsigned int *n_sweeps; // number of correlation sweeps at a particular value of tau
  unsigned int *n_vals; // number of data values already present at a particular value of tau
  unsigned int t; // global time in number of frames
  double dt; // time interval at which samples arrive
  double tau_max; // maximum time, for which the correlation should be calculated
  unsigned int tau_lin; // number of frames in the linear correlation
  unsigned int* tau_last;
  unsigned int window_distance; 

  // Convenience pointers to our stored data
  // indices: A[level][tau_i][component]
  int** tau; // time differences
  double*** A; // input quantity 1
  double*** B; // input quantity 2
  double*** result; // otput quantity
  
  // The actual allocated storage space
  double* A_data;
  double* B_data;
  double* result_data;
  int* tau_data; // just for double-checking, store tau for all results

  // compressing functions
  void (*compressA)( double* A1, double*A2, double* A_compressed, unsigned int dim_A );
  void (*compressB)( double* B1, double*B2, double* A_compressed, unsigned int dim_B );

  // correlation function
  void (*Corr_operation)  ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_C );

  // Functions producing observables A and B from the input data
  void (*A_fun)  ( double* input, unsigned int n_input, double* A, unsigned int dim_A);
  void (*B_fun)  ( double* input, unsigned int n_input, double* B, unsigned int dim_B);

} double_correlation;

struct {
  char* filename;
//  FILE* inputfile;
//  bool has_data;
  unsigned int n_columns;
  // more stuff for file handling
} datasource_double_file;


// Maybe we should think of using variable name "level" when referring to a particular level in the hierarchy?
void double_correlation_init(double_correlation* self, double dt, unsigned int tau_lin, unsigned int hierarchy_depth, unsigned int window_distance, unsigned int dim_A, unsigned int dim_B, unsigned int dim_C, void* A_fun, void* B_fun, void* Corr_operation, void* compressA, void* compressB) {
  unsigned int i,j;
  self->dt = dt;
  self->tau_lin=tau_lin;
  self->hierarchy_depth = hierarchy_depth;
  self->dim_A = dim_A;
  self->dim_B = dim_B;
  self->dim_C = dim_C;
  self->A_fun = A_fun;
  self->B_fun = B_fun;
  self->Corr_operation = Corr_operation;
  self->compressA = compressA;
  self->compressB = compressB;
  self->window_distance = window_distance;

  self->A_data = (double*)malloc(tau_lin*hierarchy_depth*dim_A*sizeof(double));
  self->B_data = (double*)malloc(tau_lin*hierarchy_depth*dim_B*sizeof(double));
  self->result_data = (double*)malloc(tau_lin*hierarchy_depth*dim_C*sizeof(double));
  self->tau_data = (int*)malloc(tau_lin*hierarchy_depth*sizeof(int));

  self->A = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->B = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->result = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->tau = (int**) malloc(hierarchy_depth*sizeof(double*));
  self->n_sweeps = (unsigned int*) malloc(hierarchy_depth*tau_lin*sizeof(unsigned int));
  self->n_vals = (unsigned int*) malloc(hierarchy_depth*sizeof(unsigned int));

  for (i=0; i<self->hierarchy_depth; i++) {
    self->A[i] = (double**) malloc(self->tau_lin*sizeof(double*));
    self->B[i] = (double**) malloc(self->tau_lin*sizeof(double*));
    self->result[i] = (double**) malloc(self->tau_lin*sizeof(double*));
    self->tau[i] = &self->tau_data[i*tau_lin];
  }
  for (i=0; i<self->hierarchy_depth; i++) {
    self->n_vals[i]=0;
    for (j=0; j<self->tau_lin; j++) {
      self->n_sweeps[i]=0;
      self->A[i][j] = &self->A_data[(i*hierarchy_depth+j)*tau_lin];
      self->B[i][j] = &self->B_data[(i*hierarchy_depth+j)*tau_lin];
      self->result[i][j] = &self->result_data[(i*hierarchy_depth+j)*tau_lin];
    }
  }

  /* This is wrong:
  printf("%d %d %ld\n", tau_lin, hierarchy_depth, sizeof(double));
  self->tau_last = (unsigned int *)malloc(hierarchy_depth*sizeof(unsigned int));
  self->t=0;
  for ( i = 0; i<self->hierarchy_depth; i++ ) {
    self->tau_last[i] = 0;
    for ( j = 0; j<self->tau_lin; j++ ) {
      printf("%f\n", (j+1)*(pow(tau_lin, i)));
      self->tau[i][j]=((1<<i)*tau_lin+j)*self->dt;
    }
  } */
}

// Data organization:
// We use a ring-like way to manage the data: at the beginning we have a linear array
// which we fill from index 0 to tau_lin. The index tau_last[i] always indicates the latest
// entry of the hierarchic "past" For every new entry in is incremented and if tau_lin is reached, 
// it starts again from the beginning. 

void double_correlation_get_data(  double_correlation* self, double* input, unsigned int n_input ) {
  // We must now go through the hierarchy and make sure there is space for the new 
  // datapoint. For every hierarchy level we have to decide if it necessary to move 
  // something
  double* C_temp = (double*)malloc(self->dim_C*sizeof(double));
  unsigned int i,j,k;

  self->t++;
  for ( i = self->hierarchy_depth-1; i > 0; i-- ) {
    if ( self->t % (1<<i) == self->tau_lin - 1 ) {
      self->tau_last[i] = (self->tau_last[i] + 1) % self->tau_lin;
      (*self->compressA)(self->A[i-1][self->tau_last[i-1]+1],  
                         self->A[i-1][self->tau_last[i-1]+2], 
                         self->A[i][self->tau_last[i]],self->dim_A);
      (*self->compressB)(self->B[i-1][self->tau_last[i-1]+1],  
                         self->B[i-1][self->tau_last[i-1]+2], 
                         self->B[i][self->tau_last[i]],self->dim_B);
    }
  }

  self->tau_last[0] = ( self->tau_last[0] + 1 ) % self->tau_lin; 

  (*self->A_fun)(input, n_input, self->A[0][self->tau_last[0]], self->dim_A);
  (*self->B_fun)(input, n_input, self->B[0][self->tau_last[0]], self->dim_A);

  // Now update the corellation estimate

  if ((self->t + 1) % self->window_distance == 0) {
    for ( i = 0; i < self->hierarchy_depth; i-- ) {
      for ( j = 0; j < self->tau_lin; j++ ) {
        (*self->Corr_operation)
             (self->A[0][self->tau_last[0]], self->dim_A, 
              self->B[i][(self->tau_last[i]+j)%self->tau_lin ], self->dim_B, 
              C_temp, self->dim_C);
        for ( k = 0; k < self->dim_C; k++ )
          self->result[i][j][k] += C_temp[k];
      }
    }
    self->n_sweeps[i*self->tau_lin+j]++;
  }
}

void double_correlation_write_corellation( double_correlation* self, char* outputfilename) {
  int j, k;
  int level; // current hierarchy level
  FILE* outputfile = fopen(outputfilename, "w");
  double dt=self->dt;
  for ( level = 0; level < self->hierarchy_depth; level++ ) 
    for ( j = 0; j < self->tau_lin; j++ ) {
      fprintf(outputfile, "%f ", self->tau[level][j]*dt);
      for ( k = 0; k < self->dim_C; k++ ) {
        fprintf(outputfile, "%f ", self->result[level][j][k]/self->n_sweeps[level]); // n_sweeps differs for each hierarchy level
      }
      fprintf(outputfile, "\n");
    }
  
}


void identity ( double* input, unsigned int n_input, double* A, unsigned int dim_A) {
  int i; 
  if ( n_input != dim_A ) {
    printf("Error in identiy: The dimensions do not match. \n");
    return;
  }
  for ( i = 0; i < dim_A; i++ ) {
    A[i] = input[i];
  }
}

void compress_linear( double* A1, double*A2, double* A_compressed, unsigned int dim_A ) {
  unsigned int i;
  for ( i = 0; i < dim_A; i++ )
    A_compressed[i] = 0.5*(A1[i]-A2[i]);
}

void scalar ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_C ) {
  double temp;
  unsigned int i;
  if (!(dim_A == dim_B && dim_C == 1 )) {
    printf("Error in scalar product: The vector sizes do not match");
    return;
  }
  for ( i = 0; i < dim_A; i++ ) {
    temp += A[i]*B[i];
  }
  C[0] = temp; 
}

void componentwise_product ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_C ) {
  unsigned int i;
  if (!(dim_A == dim_B )) {
    printf("Error in componentwise product: The vector sizes do not match");
    return;
  }
  for ( i = 0; i < dim_A; i++ ) {
    C[i] = A[i]*B[i];
  }
}

void square_distance ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_C ) {
  unsigned int i;
  if (!(dim_A == dim_B )) {
    printf("Error in componentwise product: The vector sizes do not match");
    return;
  }
  for ( i = 0; i < dim_A; i++ ) {
    C[i] += (A[i]-B[i])*(A[i]-B[i]);
  }
}

int main(int argc, char** argv) {
  double* r, * v, * f;
  unsigned int i,j;
  unsigned int n = 1;
  unsigned int dim = 3;
  double gamma = 0.;

  // These are the input parameters
  double dt=1.0; // time difference between consecutive samples
  unsigned int tau_lin=4; // number of entries in one linear correlation block
  unsigned int hierarchy_depth=3; // how deep in hierarchy to go
  unsigned int window_distance=1; // what is this?

  r = (double*)malloc(dim*n*sizeof(double));
  v = (double*)malloc(dim*n*sizeof(double));
  f = (double*)malloc(dim*n*sizeof(double));
  for ( i = 0; i < n*dim; i++) {
    r[i] = v[i] = f[i] = 0;
  }
  
  double_correlation cor;
  
  //double_correlation_init(&cor, 0.1, 10, 5, 1, n*dim, n*dim, n*dim, &identity, &identity, &componentwise_product, &compress_linear, &compress_linear);
  double_correlation_init(&cor, dt, tau_lin, hierarchy_depth, window_distance, n*dim, n*dim, n*dim, &identity, &identity, &componentwise_product, &compress_linear, &compress_linear);


  for (i=0; i<10000; i++) {
    for (j=0; j<n*dim; j++) {
      f[j] = -gamma*v[j]+gamma*((float) rand())/RAND_MAX-0.5;
      v[j] += f[j];
      r[j] += v[j];
      double_correlation_get_data(&cor, v, n*dim);
    }
  }

  double_correlation_write_corellation(&cor, "test_cor.dat");
  free(r);
  free(v);
  free(f);
  return 0;
}
