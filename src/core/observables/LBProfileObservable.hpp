#ifndef OBSERVABLES_LBPROFILEOBSERVABLE_HPP
#define OBSERVABLES_LBPROFILEOBSERVABLE_HPP

#include "Observable.hpp"
#include "ProfileObservable.hpp"

namespace Observables {

class LBProfileObservable : public LBObservable, public ProfileObservable {};

}

#endif
