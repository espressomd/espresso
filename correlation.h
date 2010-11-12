struct {
  unsigned int tau_lin;
  unsigned int tau_max;
  unsigned int hierarchy_depth;
  unsigned int n_A;
  unsigned int n_B;
  unsigned int n_C;

  double*** A;
  double*** B;
  double*** C;

  void (*compressA)( double* A1, double*A2, double* A_compressed, unsigned int n_A );
  void (*compressB)( double* B1, double*B2, double* A_compressed, unsigned int n_B );

  void (*C_Value)  ( double* A, unsigned int n_A, double* B, unsigned int n_B, double* C, unsigned int n_C );

  void (*A_Value)  ( double* input, unsigned int n_input, double* A, unsigned int n_A);
  void (*B_Value)  ( double* input, unsigned int n_input, double* B, unsigned int n_B);

} double_correllation;

struct {
  char* filename;
  FILE* inputfile;
  bool has_data;
  unsigned int n_columns;
  // more stuff for file handling
} datasource_double_file;

void double_correllation_get_data(  double_correllation correlation, double* input, unsigned int n_input ) {
  double* A, B;
  malloc(A, n_input*sizeof(double));
  malloc(B, n_input*sizeof(double));
  free(A);
  free(B);
}


void identity ( double* input, unsigned int n_input, double* A, unsigned int n_A) {
  int i; 
  if ( n_input && n_A ) {
    printf("Error in identiy: The dimensions do not match.\n");
    return;
  }
  for ( i = 0; i < n_A; i++ ) {
    A[i] = input[i];
  }
}



void compress_linear( double* A1, double*A2, double* A_compressed, unsigned int n_A ) {
  unsigned int i;
  for ( i = 0; i < n_A; i++ )
    A_compressed[i] = 0.5*A1[i]-A2[i];
}

void scalar ( double* A, unsigned int n_A, double* B, unsigned int n_B, double* C, unsigned int n_C ) {
  double temp;
  if (!(n_A == n_B && n_C == 1 )) {
    printf("Error in scalar product: The vector sizes do not match");
    return;
  }
  for ( i = 0; i < n_A; i++ ) {
    temp += A[i]*B[i];
  }
  n_C[0] = temp;
}

void main(int argc, char** argv) {
  double* data;
  double_correllation* cor = malloc(sizeof(double_correllation));
  datasource_double_file* df = malloc(sizeof(datasource_double_file));
  datasource_double_file_init(df, "test.dat");
  double_correllation_init(10, 5, df.n_columns, &identiy, df.n_columns, &identiy, &scalar);
  while(df.has_data) {
    data = malloc(df.n_columns*sizeof(double));
    datasource_double_file_get_data(df);
    double_correllation_push_data(data);
  }
  double_correllation_write_corellation(df, "test_cor.dat");
  free(df);
}
