#ifndef OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP

#include <cmath>

#include "PidObservable.hpp"
#include "Vector.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"

namespace Observables {

class CylindricalProfileObservable : public PidObservable {
public:
  ::Vector<3, double> center;
  std::string axis;
  double min_r, max_r;
  double min_phi, max_phi;
  double min_z, max_z;
  // Number of bins for each coordinate.
  int n_r_bins, n_phi_bins, n_z_bins;
  double r_bin_size() const { return (max_r - min_r) / n_r_bins; }
  double phi_bin_size() const { return (max_phi - min_phi) / n_phi_bins; }
  double z_bin_size() const { return (max_z - min_z) / n_z_bins; }
  virtual int n_values() const override {
    return n_r_bins * n_phi_bins * n_z_bins;
  };
};

} // Namespace Observables
#endif
