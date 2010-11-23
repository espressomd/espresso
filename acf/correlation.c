#include <stdio.h>
#include <stdlib.h>
typedef struct {
  unsigned int tau_lin;
  unsigned int tau_max;
  unsigned int hierarchy_depth;
  unsigned int n_A;
  unsigned int n_B;
  unsigned int n_C;
  unsigned int n_sweeps;
  unsigned int t;
  double tau_unit;
  unsigned int* tau_last;
  unsigned int window_distance;

  double** tau;
  double*** A;
  double*** B;
  double*** C;
  
  double* A_data;
  double* B_data;
  double* C_data;
  double* tau_data;

  void (*compressA)( double* A1, double*A2, double* A_compressed, unsigned int n_A );
  void (*compressB)( double* B1, double*B2, double* A_compressed, unsigned int n_B );

  void (*C_fun)  ( double* A, unsigned int n_A, double* B, unsigned int n_B, double* C, unsigned int n_C );

  void (*A_fun)  ( double* input, unsigned int n_input, double* A, unsigned int n_A);
  void (*B_fun)  ( double* input, unsigned int n_input, double* B, unsigned int n_B);

} double_correlation;

struct {
  char* filename;
//  FILE* inputfile;
//  bool has_data;
  unsigned int n_columns;
  // more stuff for file handling
} datasource_double_file;




void double_correlation_init(double_correlation* self, double tau_unit, unsigned int tau_lin, unsigned int hierarchy_depth, unsigned int window_distance, unsigned int n_A, unsigned int n_B, unsigned int n_C, void* A_fun, void* B_fun, void* C_fun, void* compressA, void* compressB) {
  unsigned int i,j;
  self->tau_unit = tau_unit;
  self->tau_lin=tau_lin;
  self->hierarchy_depth = hierarchy_depth;
  self->n_A = n_A;
  self->n_B = n_B;
  self->n_C = n_C;
  self->A_fun = A_fun;
  self->B_fun = B_fun;
  self->C_fun = C_fun;
  self->compressA = compressA;
  self->compressB = compressB;
  self->window_distance = window_distance;

  self->A_data = (double*)malloc(tau_lin*hierarchy_depth*n_A*sizeof(double));
  self->B_data = (double*)malloc(tau_lin*hierarchy_depth*n_B*sizeof(double));
  self->C_data = (double*)malloc(tau_lin*hierarchy_depth*n_C*sizeof(double));
  self->tau_data = (double*)malloc(tau_lin*hierarchy_depth*sizeof(double));

  self->A = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->B = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->C = (double***)malloc(hierarchy_depth*sizeof(double**));
  self->tau = (double**) malloc(hierarchy_depth*sizeof(double*));

  for (i=0; i<self->hierarchy_depth; i++) {
    self->A[i] = (double**) malloc(self->tau_lin*sizeof(double*));
    self->B[i] = (double**) malloc(self->tau_lin*sizeof(double*));
    self->C[i] = (double**) malloc(self->tau_lin*sizeof(double*));
    self->tau[i] = &self->tau_data[i*tau_lin];
  }
  for (i=0; i<self->hierarchy_depth; i++) {
    for (j=0; j<self->tau_lin; j++) {
      self->A[i][j] = &self->A_data[(i*hierarchy_depth+j)*tau_lin];
      self->B[i][j] = &self->B_data[(i*hierarchy_depth+j)*tau_lin];
      self->C[i][j] = &self->C_data[(i*hierarchy_depth+j)*tau_lin];
    }
  }

  printf("%d %d %d\n", tau_lin, hierarchy_depth, sizeof(double));
  self->tau_last = (unsigned int *)malloc(hierarchy_depth*sizeof(unsigned int));
  self->t=0;
  for ( i = 0; i<self->hierarchy_depth; i++ ) {
    self->tau_last[i] = 0;
    for ( j = 0; j<self->tau_lin; j++ ) {
      printf("%f\n", (j+1)*(pow(tau_lin, i)));
      self->tau[i][j]=((1<<i)*tau_lin+j)*self->tau_unit;
    }
  }
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
  double* C_temp = (double*)malloc(self->n_C*sizeof(double));
  unsigned int i,j,k;

  self->t++;
  for ( i = self->hierarchy_depth-1; i > 0; i-- ) {
    if ( self->t % (1<<i) == self->tau_lin - 1 ) {
      self->tau_last[i] = (self->tau_last[i] + 1) % self->tau_lin;
      (*self->compressA)(self->A[i-1][self->tau_last[i-1]+1],  
                         self->A[i-1][self->tau_last[i-1]+2], 
                         self->A[i][self->tau_last[i]],self->n_A);
      (*self->compressB)(self->B[i-1][self->tau_last[i-1]+1],  
                         self->B[i-1][self->tau_last[i-1]+2], 
                         self->B[i][self->tau_last[i]],self->n_B);
    }
  }

  self->tau_last[0] = ( self->tau_last[0] + 1 ) % self->tau_lin; 

  (*self->A_fun)(input, n_input, self->A[0][self->tau_last[0]], self->n_A);
  (*self->B_fun)(input, n_input, self->B[0][self->tau_last[0]], self->n_A);

  // Now update the corellation estimate

  if ((self->t + 1) % self->window_distance == 0) {
    for ( i = 0; i < self->hierarchy_depth; i-- ) {
      for ( j = 0; j < self->tau_lin; j++ ) {
        (*self->C_fun)
             (self->A[0][self->tau_last[0]], self->n_A, 
              self->B[i][(self->tau_last[i]+j)%self->tau_lin ], self->n_B, 
              C_temp, self->n_C);
        for ( k = 0; k < self->n_C; k++ )
          self->C[i][j][k] += C_temp[k];
      }
    }
    self->n_sweeps++;
  }
}

void double_correlation_write_corellation( double_correlation* self, char* outputfilename) {
  int i, j, k;
  FILE* outputfile = fopen(outputfilename, "w");
  for ( i = 0; i < self->hierarchy_depth; i++ ) 
    for ( j = 0; j < self->tau_lin; j++ ) {
      fprintf(outputfile, "%f ", self->tau[i][j]);
      for ( k = 0; k < self->n_C; k++ ) {
        fprintf(outputfile, "%f ", self->C[i][j][k]/self->n_sweeps);
      }
      fprintf(outputfile, "\n");
    }
  
}



void identity ( double* input, unsigned int n_input, double* A, unsigned int n_A) {
  int i; 
  if ( n_input != n_A ) {
    printf("Error in identiy: The dimensions do not match. \n");
    return;
  }
  for ( i = 0; i < n_A; i++ ) {
    A[i] = input[i];
  }
}

void compress_linear( double* A1, double*A2, double* A_compressed, unsigned int n_A ) {
  unsigned int i;
  for ( i = 0; i < n_A; i++ )
    A_compressed[i] = 0.5*(A1[i]-A2[i]);
}

void scalar ( double* A, unsigned int n_A, double* B, unsigned int n_B, double* C, unsigned int n_C ) {
  double temp;
  unsigned int i;
  if (!(n_A == n_B && n_C == 1 )) {
    printf("Error in scalar product: The vector sizes do not match");
    return;
  }
  for ( i = 0; i < n_A; i++ ) {
    temp += A[i]*B[i];
  }
  C[0] = temp;
}

void componentwise_product ( double* A, unsigned int n_A, double* B, unsigned int n_B, double* C, unsigned int n_C ) {
  unsigned int i;
  if (!(n_A == n_B )) {
    printf("Error in componentwise product: The vector sizes do not match");
    return;
  }
  for ( i = 0; i < n_A; i++ ) {
    C[i] = A[i]*B[i];
  }
}

int main(int argc, char** argv) {
  double* r, * v, * f;
  unsigned int i,j;
  unsigned int n = 1;
  unsigned int dim = 3;
  double gamma = 0.1;
  r = (double*)malloc(dim*n*sizeof(double));
  v = (double*)malloc(dim*n*sizeof(double));
  f = (double*)malloc(dim*n*sizeof(double));
  for ( i = 0; i < n*dim; i++) {
    r[i] = v[i] = f[i] = 0;
  }
  
  double_correlation cor;
  
  double_correlation_init(&cor, 0.1, 10, 5, 1, n*dim, n*dim, n*dim, &identity, &identity, &componentwise_product, &compress_linear, &compress_linear);


  for (i<0; i<10000; i++) {
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
