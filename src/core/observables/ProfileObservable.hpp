#ifndef OBSERVABLES_PROFILEOBSERVABLE_HPP
#define OBSERVABLES_PROFILEOBSERVABLE_HPP

#include "Observable.hpp"
#include "particle_data.hpp"
#include <vector>
#include "integrate.hpp"

namespace Observables {

// Observable which acts on a given list of particle ids
class ProfileObservable : public Observable {
public:
  std::vector<int> ids;
  double minx;
  double maxx;
  double miny;
  double maxy;
  double minz;
  double maxz;
  int xbins;
  int ybins;
  int zbins;
  std::vector<double> container;
  virtual int n_values() const override { return xbins * ybins * zbins; };
};

} // Namespace Observables
#endif
