#ifndef OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP

#include <cmath>

#include "CylindricalProfileObservable.hpp"
#include "LBObservable.hpp"
#include "Vector.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "utils/coordinate_transformation.hpp"

namespace Observables {

class CylindricalLBProfileObservable : public CylindricalProfileObservable, public LBObservable {
};

} // Namespace Observables
#endif
