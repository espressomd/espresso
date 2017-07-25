#ifndef OBSERVABLES_PIDOBSERVABLE_HPP
#define OBSERVABLES_PIDOBSERVABLE_HPP

#include "Observable.hpp"
#include <vector>

namespace Observables {

// Observable which acts on a given list of particle ids
class PidObservable : public Observable {
public:
    std::vector<int> ids;
    virtual int n_values() const override { return 3*ids.size();};
};

} // Namespace Observables
#endif
  

