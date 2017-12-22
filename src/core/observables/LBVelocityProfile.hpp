#ifndef OBSERVABLES_LBVELOCITYPROFILE_HPP
#define OBSERVABLES_LBVELOCITYPROFILE_HPP

#include "ProfileObservable.hpp"
#include "particle_data.hpp"
#include "lb.hpp"

#include <vector>

namespace Observables {

class LBVelocityProfile : public ProfileObservable {
public:
  virtual int n_values() const override {
    return 3 * n_x_bins * n_y_bins * n_z_bins;
  }
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
#ifdef LB
    unsigned int maxi, maxj, maxk;
    double xoffset, yoffset, zoffset;
    double x_incr, y_incr, z_incr;
    double p[3], v[3];
    int linear_index;

#ifdef LB_GPU
    if (lattice_switch & LATTICE_LB_GPU) {
      throw std::runtime_error("The Lb Velocity profile observable is "
                               "currently not available for LBGPU");
    }
#endif
    if (lattice_switch & LATTICE_LB) {
      double normalization_factor = 1.;
      if (n_x_bins == 1) {
        maxi = (int)floor(box_l[0] / lbpar.agrid);
        normalization_factor /= maxi;
        xoffset = 0;
        x_incr = lbpar.agrid;
      } else {
        maxi = n_x_bins;
        xoffset = min_x;
        x_incr = (max_x - min_x) / (n_x_bins - 1);
      }
      if (n_y_bins == 1) {
        maxj = (int)floor(box_l[1] / lbpar.agrid);
        normalization_factor /= maxj;
        yoffset = 0;
        y_incr = lbpar.agrid;
      } else {
        maxj = n_y_bins;
        yoffset = min_y;
        y_incr = (max_y - min_y) / (n_y_bins - 1);
      }
      if (n_z_bins == 1) {
        maxk = (int)floor(box_l[2] / lbpar.agrid);
        normalization_factor /= maxk;
        zoffset = 0;
        z_incr = lbpar.agrid;
      } else {
        maxk = n_z_bins;
        zoffset = min_z;
        z_incr = (max_z - min_z) / (n_z_bins - 1);
      }
      unsigned int i, j, k;
      for (i = 0; i < maxi; i++) {
        for (j = 0; j < maxj; j++) {
          for (k = 0; k < maxk; k++) {
            p[0] = xoffset + i * x_incr;
            p[1] = yoffset + j * y_incr;
            p[2] = zoffset + k * z_incr;
            if (lb_lbfluid_get_interpolated_velocity(p, v) != 0)
              throw std::runtime_error("LB velocity interpolation failed.");
            linear_index = 0;
            if (n_x_bins > 1)
              linear_index += i * n_y_bins * n_z_bins;
            if (n_y_bins > 1)
              linear_index += j * n_z_bins;
            if (n_z_bins > 1)
              linear_index += k;

            res[3 * linear_index + 0] += v[0];
            res[3 * linear_index + 1] += v[1];
            res[3 * linear_index + 2] += v[2];
          }
        }
      }
      for (int i = 0; i < res.size(); i++) {
        res[i] *= normalization_factor;
      }
    }
#endif
    return res;
  }
};

} // Namespace Observables

#endif
