int ObservableStructureFactor::actual_calculate(PartCfg & partCfg) {
  double* A = last_value;
  // FIXME Currently scattering length is hardcoded as 1.0
  int l;
  int order, order2, n;
  double twoPI_L, C_sum, S_sum, qr; 
  //  DoubleList *scattering_length;
  observable_sf_params * params;
  params = (observable_sf_params*) container;
  //  scattering_length = params->scattering_length;
  const double scattering_length=1.0;
  order = params->order;
  order2=order*order;
  twoPI_L = 2*PI/box_l[0];
  
  if (!sortPartCfg()) {
    runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }

  
  l=0;
  float partCache[n_part*3];
  for(int p=0; p<n_part; p++) {
    for (int i=0;i<3;i++){
      partCache[3*p+i]=partCfg[p].r.p[i];
    }
  }
  //printf("n: %d, dim_sf: %d\n",n_A, params.dim_sf); fflush(stdout);
  for(int i=-order; i<=order; i++) {
    for(int j=-order; j<=order; j++) {
      for(int k=-order; k<=order; k++) {
	n = i*i + j*j + k*k;
	if ((n<=order2) && (n>=1)) {
	  C_sum = S_sum = 0.0;
	  //printf("l: %d, n: %d %d %d\n",l,i,j,k); fflush(stdout);
	  for(int p=0; p<n_part; p++) {
	    //qr = twoPI_L * ( i*partCfg[p].r.p[0] + j*partCfg[p].r.p[1] + k*partCfg[p].r.p[2] );
	    qr = twoPI_L * ( i*partCache[3*p+0] + j*partCache[3*p+1] + k*partCache[3*p+2] );
	    C_sum+= scattering_length * cos(qr);
	    S_sum-= scattering_length * sin(qr);
	  }
	  A[l]   =C_sum;
	  A[l+1] =S_sum;
	  l=l+2;
	}
      }
    }
  }
  l = 0;
  for(int k=0;k<n;k++) {
    //devide by the sqrt(number_of_particle) due to complex product and no k-vector averaging so far
    A[k] /= sqrt(n_part);
  }
  //printf("finished calculating sf\n"); fflush(stdout);
  return 0;
}
