#ifndef OBSERVABLES_LBPROFILEOBSERVABLE_HPP
#define OBSERVABLES_LBPROFILEOBSERVABLE_HPP

#include "LBObservable.hpp"
#include "ProfileObservable.hpp"

namespace Observables {

class LBProfileObservable : public LBObservable, public ProfileObservable {};

}

#endif
