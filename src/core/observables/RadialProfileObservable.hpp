#ifndef OBSERVABLES_RADIALPROFILEOBSERVABLE_HPP
#define OBSERVABLES_RADIALPROFILEOBSERVABLE_HPP

#include <cmath>

#include "PidObservable.hpp"
#include "particle_data.hpp"
#include "Vector.hpp"
#include "integrate.hpp"

namespace Observables {

class RadialProfileObservable : public PidObservable {
public:
  ::Vector<3, double> center;
  double min_r, max_r;
  double min_phi, max_phi;
  double min_z, max_z;
  // Number of bins for each coordinate.
  int n_r_bins, n_phi_bins, n_z_bins;
  double r_bin_size() {
      return (max_r - min_r) / n_r_bins;
  }
  double phi_bin_size() {
      return (max_phi - min_phi) / n_phi_bins;
  }
  double z_bin_size() {
      return (max_z - min_z) / n_z_bins;
  }
  virtual int n_values() const override { return n_r_bins * n_phi_bins * n_z_bins; };
  ::Vector<3, double> transform_to_cylinder_coordinates(::Vector<3, double> const &pos) const {
    double r = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    double phi = std::atan2(pos[1], pos[0]);
    return ::Vector<3, double> {r, phi, pos[2]};
  }
};

} // Namespace Observables
#endif
