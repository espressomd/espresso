#ifndef OBSERVABLES_LBVELOCITYPROFILE_HPP
#define OBSERVABLES_LBVELOCITYPROFILE_HPP

#include "LBProfileObservable.hpp"
#include "lb.hpp"
#include "particle_data.hpp"

#include <vector>

namespace Observables {

class LBVelocityProfile : public LBProfileObservable {
public:
  virtual int n_values() const override {
    return 3 * n_x_bins * n_y_bins * n_z_bins;
  }
  virtual std::vector<double> operator()(PartCfg &partCfg) const override;
};

} // Namespace Observables

#endif
