#ifndef OBSERVABLES_PIDBLOCKEDOBSERVABLE_HPP
#define OBSERVABLES_PIDBLOCKEDOBSERVABLE_HPP

#include "Observable.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include <vector>

namespace Observables {

// Observable which acts on a given list of particle ids
class PidBlockedObservable : public Observable {
public:
  std::vector<int> ids;
  int blocks;
  virtual int n_values() const override { return 3 * blocks; };
};

} // Namespace Observables
#endif
