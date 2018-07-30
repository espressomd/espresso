#ifndef OBSERVABLES_PROFILEOBSERVABLE_HPP
#define OBSERVABLES_PROFILEOBSERVABLE_HPP

#include "Observable.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include <vector>

namespace Observables {

// Observable which acts on a given list of particle ids
class ProfileObservable : virtual public Observable {
public:
  double min_x, max_x;
  double min_y, max_y;
  double min_z, max_z;
  int n_x_bins, n_y_bins, n_z_bins;
  virtual int n_values() const override {
    return n_x_bins * n_y_bins * n_z_bins;
  };
};

} // Namespace Observables
#endif
