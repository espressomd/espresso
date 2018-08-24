int ObservableStructureFactorFast::actual_calculate(PartCfg & partCfg) {
  //printf("calculating\n");
  double* A = last_value;
  // FIXME Currently scattering length is hardcoded as 1.0
  observable_sf_params * params = (observable_sf_params*) container;
  const int k_max = params->order * params->k_density;
  const double scattering_length=1.0;
  const double twoPI_L = 2*PI/box_l[0];
  
  if (!sortPartCfg()) {
    runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  
  for(int p=0; p<n; p++) {
    A[p]   = 0.0;
  }
  
  float partCache[n_part*3];
  for(int p=0; p<n_part; p++) {
    for (int i=0;i<3;i++){
      partCache[3*p+i]=partCfg[p].r.p[i];
    }
  }
  int k_density = params->k_density;
  int l=0;
  for(int k=0; k<k_max; k++) {
    int order=k/k_density+1;
    switch (k % k_density){
    case 0: // length sqrt(1)
      for (int dir=0;dir<3;dir++){
	double C_sum = 0;
	double S_sum = 0;
	for(int p=0; p<n_part; p++) {
	//double qr = twoPI_L * k  * ( ix*partCache[3*p+0] + iy*partCache[3*p+1] + iz*partCache[3*p+2] );
	  double qr = twoPI_L * order * ( partCache[3*p+dir]);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 1: // length sqrt(2)
      for (int dir=0;dir<6;dir++){
	int fac1,fac2,off1,off2;
	switch (dir){
	case 0: fac1= 1; off1=0; fac2= 1; off2=1; break;
	case 1: fac1= 1; off1=0; fac2= 1; off2=2; break;
	case 2: fac1= 1; off1=1; fac2= 1; off2=2; break;
	case 3: fac1=-1; off1=0; fac2= 1; off2=1; break;
	case 4: fac1=-1; off1=0; fac2= 1; off2=2; break;
	case 5: fac1=-1; off1=1; fac2= 1; off2=2; break;
	}
	double C_sum = 0;
	double S_sum = 0;
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+off1]*fac1+ partCache[3*p+off2]*fac2);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 2: // length sqrt(3)
      for (int dir=0;dir<4;dir++){
	double C_sum = 0;
	double S_sum = 0;
	int fac1=(1-2*(dir%2));
	int fac2=(1-2*(dir/2));
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+0]*fac1 + partCache[3*p+1]*fac2 + partCache[3*p+2]);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 3: // length sqrt(5) 
      for (int dir=0;dir<6;dir++){
	double C_sum = 0;
	double S_sum = 0;
	// pick 6 random vectors
	int fac1=(1-2*(dir/3));
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+(dir%3)]*fac1 + partCache[3*p+((dir+1)%3)]*2);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 4: // length sqrt(6) 
      for (int dir=0;dir<6;dir++){
	double C_sum = 0;
	double S_sum = 0;
	// pick 6 random vectors
	int fac1=(1-2*(dir/3))*2;
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+(dir%3)] + partCache[3*p+((dir+1)%3)]*fac1 + partCache[3*p+((dir+2)%3)]);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 5: // length sqrt(9)
      for (int dir=0;dir<6;dir++){
	double C_sum = 0;
	double S_sum = 0;
	// pick 6 random vectors
	int fac1=(1-2*(dir/3))*2;
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+(dir%3)]*fac1 + partCache[3*p+((dir+1)%3)]*2 + partCache[3*p+((dir+2)%3)]);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    default:
      runtimeErrorMsg() <<"so many samples per order not yet implemented";
      return -1;
    }
  }
  for(int l=0;l<n;l++) {
    //devide by the sqrt of number_of_particle, average later
    A[l] /= sqrt(n_part);
  }
  return 0;
}

