#ifndef CORE_OBSERVABLES_LBPROFILEOBSERVABLE_HPP
#define CORE_OBSERVABLES_LBPROFILEOBSERVABLE_HPP

#include "LBObservable.hpp"
#include "ProfileObservable.hpp"

namespace Observables {

class LBProfileObservable : public LBObservable, public ProfileObservable {};

}
#endif
