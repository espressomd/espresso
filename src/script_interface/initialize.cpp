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
#include "ScriptInterface.hpp"
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
#include "lbboundaries/initialize.hpp"
#include "mpiio/initialize.hpp"
#include "observables/initialize.hpp"

#include "ComFixed.hpp"

#include "core/communication.hpp"
#include "virtual_sites/initialize.hpp"

#include "ObjectManager.hpp"

namespace ScriptInterface {
namespace {
std::unique_ptr<ObjectManager> m_om;
}

void initialize(Communication::MpiCallbacks &cb) {
  m_om = std::make_unique<ObjectManager>(&cb);
  auto om = m_om.get();

  ObjectHandle::initialize(om);

  Shapes::initialize(om);
  Constraints::initialize(om);
#ifdef H5MD
  Writer::initialize(om);
#endif
  Accumulators::initialize(om);
  Observables::initialize(om);
  ClusterAnalysis::initialize(om);
  LBBoundaries::initialize(om);
  PairCriteria::initialize(om);
  VirtualSites::initialize(om);
  MPIIO::initialize(om);
  CollisionDetection::initialize(om);

  om->register_new<ComFixed>("ComFixed");
}

} /* namespace ScriptInterface */
