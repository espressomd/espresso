/*
 * Copyright (C) 2015-2019 The ESPResSo project
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
#include "cluster_analysis/initialize.hpp"
#include "config.hpp"
#include "constraints/initialize.hpp"
#include "pair_criteria/initialize.hpp"
#include "shapes/initialize.hpp"
#ifdef H5MD
#include "h5md/initialize.hpp"
#endif
#include "accumulators/initialize.hpp"
#include "collision_detection/initialize.hpp"
#include "generic_dd/initialize.hpp"
#include "lbboundaries/initialize.hpp"
#include "mpiio/initialize.hpp"
#include "observables/initialize.hpp"

#include "ComFixed.hpp"

#include "ParallelScriptInterface.hpp"

#include "core/communication.hpp"
#include "virtual_sites/initialize.hpp"

namespace ScriptInterface {

void initialize() {
  ParallelScriptInterface::initialize(Communication::mpiCallbacks());

  Shapes::initialize();
  Constraints::initialize();
#ifdef H5MD
  Writer::initialize();
#endif
  Accumulators::initialize();
  Observables::initialize();
  ClusterAnalysis::initialize();
  LBBoundaries::initialize();
  PairCriteria::initialize();
  VirtualSites::initialize();
  MPIIO::initialize();
  CollisionDetection::initialize();
  GenericDD::initialize();

  ScriptInterface::register_new<ComFixed>("ComFixed");
}

} /* namespace ScriptInterface */
