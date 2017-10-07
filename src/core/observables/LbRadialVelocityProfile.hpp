#ifndef OBSERVABLES_LBRADIALVELOCITYPROFILE
#define OBSERVABLES_LBRADIALVELOCITYPROFILE

#include "config.hpp"
#include "Observable.hpp"

namespace Observables {

#ifdef LB
//int ObservableLbRadialVelocityProfile::actual_calculate(PartCfg & partCfg) {
//  double* A = last_value;
//  void* pdata = container;
//  unsigned int n_A = n;
//
//#ifdef LB_GPU
//  if (lattice_switch & LATTICE_LB_GPU)
//    return statistics_observable_lbgpu_radial_velocity_profile((radial_profile_data*) pdata, A, n_A);
//#endif
//  
//  if (!(lattice_switch & LATTICE_LB))
//    return ES_ERROR;
//
//  if (n_nodes==1) {
//    mpi_observable_lb_radial_velocity_profile_parallel(pdata, A, n_A);
//    return ES_OK;
//  } else {
//    mpi_observable_lb_radial_velocity_profile();
//    MPI_Bcast(pdata, sizeof(radial_profile_data), MPI_BYTE, 0, comm_cart);
//    double* data = (double*) Utils::malloc(n_A*sizeof(double));
//    mpi_observable_lb_radial_velocity_profile_parallel(pdata, data, n_A);
//    MPI_Reduce(data, A, n_A, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
//    free(data);
//    return ES_OK;
//  }
//}

int mpi_observable_lb_radial_velocity_profile_parallel(void* pdata_, double* A, unsigned int n_A);

void mpi_observable_lb_radial_velocity_profile_slave_implementation() {
  radial_profile_data pdata;
  MPI_Bcast(&pdata, sizeof(radial_profile_data), MPI_BYTE, 0, comm_cart);
  unsigned int n_A=3*pdata.rbins*pdata.phibins*pdata.zbins;
  double* data = (double*) Utils::malloc(n_A*sizeof(double));
  mpi_observable_lb_radial_velocity_profile_parallel(&pdata, data, n_A);
  MPI_Reduce(data, 0, n_A, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
  free(data);

}

int mpi_observable_lb_radial_velocity_profile_parallel(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i, j, k;
  unsigned int maxi, maxj, maxk;
  double roffset, phioffset, zoffset;
  double r, phi, z;
  double r_incr, phi_incr, z_incr;
  double p[3], v[3];
  double v_r, v_phi, v_z;
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) pdata_;
  int linear_index;

    
  for ( i = 0; i<n_A; i++ ) {
    A[i]=0;
  }
  double normalization_factor = 1.;
  if ( pdata->rbins == 1 ) {
    return 1;
  } else {
    maxi = pdata->rbins;
    roffset=pdata->minr;
    r_incr=(pdata->maxr-pdata->minr)/(pdata->rbins-1);
  }
  if ( pdata->phibins == 1 ) {
    maxj = (int)floor( 2*3.1415*pdata->maxr/lbpar.agrid ) ; 
    normalization_factor/=maxj;
    phioffset=0;
    phi_incr=2*3.1415/maxj;
  } else {
    maxj = pdata->phibins;
    phioffset=pdata->minphi;
    phi_incr=(pdata->maxphi-pdata->minphi)/(pdata->phibins-1);
  }
  if ( pdata->zbins == 1 ) {
    maxk = (int) floor(box_l[2]/lbpar.agrid);
    normalization_factor/=maxk;
    zoffset=-pdata->center[2];
    z_incr=lbpar.agrid;
  } else {
    maxk = pdata->zbins;
    zoffset=pdata->minz;
    z_incr=(pdata->maxz-pdata->minz)/(pdata->zbins-1);
  }

  for ( i = 0; i < maxi; i++ ) {
    for ( j = 0; j < maxj; j++ ) {
      for ( k = 0; k < maxk; k++ ) {
        r = roffset + i*r_incr;
        phi = phioffset + j*phi_incr;
        z = zoffset + k*z_incr;
        p[0]=r*cos(phi)+pdata->center[0];
        p[1]=r*sin(phi)+pdata->center[1];
        p[2]=z+pdata->center[2];
        if (   p[0] < my_left[0] || p[0]>my_right[0] 
            || p[1] < my_left[1] || p[1]>my_right[1] 
            || p[2] < my_left[2] || p[2]>my_right[2] )
          continue;
          return 1;
        linear_index = 0;
        if (pdata->rbins > 1)
          linear_index += i*pdata->phibins*pdata->zbins;
        if (pdata->phibins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;
        if (r>0) {
          v_r = 1/r*((p[0]-pdata->center[0])*v[0] + (p[1]-pdata->center[1])*v[1]); 
          v_phi = 1/r/r*((p[0]-pdata->center[0])*v[1]-(p[1]-pdata->center[1])*v[0]);
        } else {
          v_r = 0;
          v_phi = 0;
        }
        v_z = v[2];

        A[3*linear_index+0]+=v_r;
        A[3*linear_index+1]+=v_phi;
        A[3*linear_index+2]+=v_z;
      }
    }
  }
  
  for ( i = 0; i<n_A; i++ ) {
    A[i]*=normalization_factor;
  }

  
  return 0;
}
#endif

void transform_to_cylinder_coordinates(double x, double y, double z_, double* r, double* phi, double* z) {
  *z =  z_;
  *r =  sqrt(x*x+y*y);
  *phi = atan2(y,x);
}
}
#endif

