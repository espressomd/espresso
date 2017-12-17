#ifndef OBSERVABLES_PROFILEOBSERVABLE_HPP
#define OBSERVABLES_PROFILEOBSERVABLE_HPP

#include "Observable.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include <vector>

namespace Observables {

// Observable which acts on a given list of particle ids
class ProfileObservable : public Observable {
public:
  std::vector<int> ids;
  double min_x;
  double max_x;
  double min_y;
  double max_y;
  double min_z;
  double max_z;
  int n_x_bins;
  int n_y_bins;
  int n_z_bins;
  std::vector<double> container;
  virtual int n_values() const override { return n_x_bins * n_y_bins * n_z_bins; };
};

} // Namespace Observables
#endif
