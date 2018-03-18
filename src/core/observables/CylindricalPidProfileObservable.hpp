#ifndef OBSERVABLES_CYLINDRICALPIDPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALPIDPROFILEOBSERVABLE_HPP

#include <cmath>

#include "CylindricalProfileObservable.hpp"
#include "PidObservable.hpp"

namespace Observables {

class CylindricalPidProfileObservable : public PidObservable,
                                        public CylindricalProfileObservable {};

} // Namespace Observables
#endif
