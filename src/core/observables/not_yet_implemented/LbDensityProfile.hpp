#ifdef LB
int ObservableLbVelocityProfile::actual_calculate(PartCfg & partCfg) {
  double* A= last_value;
  void* pdata_ = container;
  unsigned int n_A = n;
  unsigned int maxi, maxj, maxk;
  double xoffset, yoffset, zoffset;
  double x_incr, y_incr, z_incr;
  double p[3], v[3];
  int linear_index;
  profile_data* pdata;
  pdata=(profile_data*) container;


#ifdef LB_GPU
  if (lattice_switch & LATTICE_LB_GPU)
    return statistics_observable_lbgpu_velocity_profile((profile_data*) pdata_, A, n_A);
#endif
  if (lattice_switch & LATTICE_LB) {
    for ( int i = 0; i<n; i++ ) {
      A[i]=0;
    }
    double normalization_factor = 1.;
    if ( pdata->xbins == 1 ) {
      maxi = (int) floor(box_l[0]/lbpar.agrid);
      normalization_factor/=maxi;
      xoffset=0;
      x_incr=lbpar.agrid;
    } else {
      maxi = pdata->xbins;
      xoffset=pdata->minx;
      x_incr=(pdata->maxx-pdata->minx)/(pdata->xbins-1);
    }
    if ( pdata->ybins == 1 ) {
      maxj = (int) floor(box_l[1]/lbpar.agrid);
      normalization_factor/=maxj;
      yoffset=0;
      y_incr=lbpar.agrid;
    } else {
      maxj = pdata->ybins;
      yoffset=pdata->miny;
      y_incr=(pdata->maxy-pdata->miny)/(pdata->ybins-1);
    }
    if ( pdata->zbins == 1 ) {
      maxk = (int) floor(box_l[2]/lbpar.agrid);
      normalization_factor/=maxk;
      zoffset=0;
      z_incr=lbpar.agrid;
    } else {
      maxk = pdata->zbins;
      zoffset=pdata->minz;
      z_incr=(pdata->maxz-pdata->minz)/(pdata->zbins-1);
    }
    unsigned int i, j, k;
    for ( i = 0; i < maxi; i++ ) {
      for ( j = 0; j < maxj; j++ ) {
	for ( k = 0; k < maxk; k++ ) {
	  p[0]=xoffset + i*x_incr;
	  p[1]=yoffset + j*y_incr;
	  p[2]=zoffset + k*z_incr;
	  if (lb_lbfluid_get_interpolated_velocity(p, v)!=0)
	    return 1;
	  linear_index = 0;
	  if (pdata->xbins > 1)
	    linear_index += i*pdata->ybins*pdata->zbins;
	  if (pdata->ybins > 1)
	    linear_index += j*pdata->zbins;
	  if (pdata->zbins > 1)
	    linear_index +=k;

	  A[3*linear_index+0]+=v[0];
	  A[3*linear_index+1]+=v[1];
	  A[3*linear_index+2]+=v[2];
	}
      }
    }
  
    for ( int i = 0; i<n; i++ ) {
      A[i]*=normalization_factor;
    }

  
  }
  return 0;
}
#endif
