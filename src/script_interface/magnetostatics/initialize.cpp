/*
 * Copyright (C) 2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "initialize.hpp"

#include "config/config.hpp"

#ifdef DIPOLES

#include "Actor_impl.hpp"

#include "DipolarBarnesHutGpu.hpp"
#include "DipolarDirectSum.hpp"
#include "DipolarDirectSumGpu.hpp"
#include "DipolarLayerCorrection.hpp"
#include "DipolarP3M.hpp"
#include "DipolarScafacos.hpp"

#include "core/magnetostatics/dipoles.hpp"
#include "core/magnetostatics/registration.hpp"

#include "script_interface/auto_parameters/AutoParameter.hpp"

#endif // DIPOLES

#include <utils/Factory.hpp>

namespace ScriptInterface {
namespace Dipoles {

void initialize(Utils::Factory<ObjectHandle> *om) {
#ifdef DIPOLES
  om->register_new<DipolarDirectSum>("Dipoles::DipolarDirectSumCpu");
#ifdef DIPOLAR_DIRECT_SUM
  om->register_new<DipolarDirectSumGpu>("Dipoles::DipolarDirectSumGpu");
#endif
#ifdef DIPOLAR_BARNES_HUT
  om->register_new<DipolarBarnesHutGpu>("Dipoles::DipolarBarnesHutGpu");
#endif
#ifdef DP3M
  om->register_new<DipolarP3M>("Dipoles::DipolarP3M");
#endif
#ifdef SCAFACOS_DIPOLES
  om->register_new<DipolarScafacos>("Dipoles::DipolarScafacos");
#endif
  om->register_new<DipolarLayerCorrection>("Dipoles::DipolarLayerCorrection");
#endif // DIPOLES
}

} // namespace Dipoles
} // namespace ScriptInterface
