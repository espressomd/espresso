/*
  Copyright (C) 2015,2016 The ESPResSo project

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


#include "config.hpp"
#include "initialize.hpp"
#include "constraints/initialize.hpp"
#include "shapes/initialize.hpp"
#ifdef H5MD
#include "h5md/initialize.hpp"
#endif
#include "observables/initialize.hpp" 
#include "correlators/initialize.hpp" 
#include "lbboundaries/initialize.hpp"

#include "ComFixed.hpp"

#include "ParallelScriptInterface.hpp"
#include "VariantTester.hpp"

#include "core/communication.hpp"

namespace ScriptInterface {

void initialize() {
  ParallelScriptInterface::initialize(Communication::mpiCallbacks());

  Shapes::initialize();
  Constraints::initialize();
#ifdef H5MD
  Writer::initialize();
#endif
  Observables::initialize();
  Correlators::initialize();
  LBBoundaries::initialize();

  ScriptInterface::register_new<Testing::VariantTester>("Testing::VariantTester");
  ScriptInterface::register_new<ComFixed>();
}

} /* namespace ScriptInterface */
