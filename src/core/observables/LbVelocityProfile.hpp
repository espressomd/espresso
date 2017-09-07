#ifndef OBSERVABLES_LBVELOCITYPROFILE_HPP
#define OBSERVABLES_LBVELOCITYPROFILE_HPP

#include "ProfileObservable.hpp"
#include "particle_data.hpp" 
#include <vector>


namespace Observables {


class LBVelocityProfile : public ProfileObservable {
public:
    virtual int n_values() const override { return 3 * xbins *ybins*zbins; }
    virtual int actual_calculate(PartCfg & partCfg) override {
#ifdef LB
  unsigned int maxi, maxj, maxk;
  double xoffset, yoffset, zoffset;
  double x_incr, y_incr, z_incr;
  double p[3], v[3];
  int linear_index;


#ifdef LB_GPU
  if (lattice_switch & LATTICE_LB_GPU) {
//    return statistics_observable_lbgpu_velocity_profile((profile_data*) pdata_, A, n_A);
    throw std::runtime_error("The Lb Velocity profile observable is currently not available for lb glu");
    /* notes:
      * The proifle_data struct is no longer used. Instead, the values are
        stored in the Observable class. The profile code, however, needs the
        original struct on the gpu.
      * Instance a profile_Data struct from Observable.hpp here
      * fill it with info from ProfileObservable-class member variables
      * pass it and a &(last_value[0]) to the lbgpu profile function
    */
}
#endif
  if (lattice_switch & LATTICE_LB) {
    double normalization_factor = 1.;
    if ( xbins == 1 ) {
      maxi = (int) floor(box_l[0]/lbpar.agrid);
      normalization_factor/=maxi;
      xoffset=0;
      x_incr=lbpar.agrid;
    } else {
      maxi = xbins;
      xoffset=minx;
      x_incr=(maxx-minx)/(xbins-1);
    }
    if ( ybins == 1 ) {
      maxj = (int) floor(box_l[1]/lbpar.agrid);
      normalization_factor/=maxj;
      yoffset=0;
      y_incr=lbpar.agrid;
    } else {
      maxj = ybins;
      yoffset=miny;
      y_incr=(maxy-miny)/(ybins-1);
    }
    if ( zbins == 1 ) {
      maxk = (int) floor(box_l[2]/lbpar.agrid);
      normalization_factor/=maxk;
      zoffset=0;
      z_incr=lbpar.agrid;
    } else {
      maxk = zbins;
      zoffset=minz;
      z_incr=(maxz-minz)/(zbins-1);
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
	  if (xbins > 1)
	    linear_index += i*ybins*zbins;
	  if (ybins > 1)
	    linear_index += j*zbins;
	  if (zbins > 1)
	    linear_index +=k;

	  last_value[3*linear_index+0]+=v[0];
	  last_value[3*linear_index+1]+=v[1];
	  last_value[3*linear_index+2]+=v[2];
	}
}
}
for ( int i = 0; i<last_value.size(); i++ ) {
      last_value[i]*=normalization_factor;
}
}
#endif  
  return 0;
}
};

} // Namespace Observables

#endif

