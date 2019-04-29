/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP

#include <cmath>

#include "CylindricalProfileObservable.hpp"
#include "LBObservable.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include <utils/Vector.hpp>
#include <utils/coordinate_transformation.hpp>

namespace Observables {

class CylindricalLBProfileObservable : public CylindricalProfileObservable,
                                       public LBObservable {};

} // Namespace Observables
#endif
