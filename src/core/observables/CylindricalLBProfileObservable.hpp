#ifndef OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP

#include <cmath>

#include "CylindricalProfileObservable.hpp"
#include "LBObservable.hpp"
#include "Vector.hpp"

namespace Observables {

class CylindricalLBProfileObservable : public CylindricalProfileObservable, public LBObservable {};

} // Namespace Observables
#endif
