/*
 * Copyright (C) 2015-2022 The ESPResSo project
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

#include "accumulators/initialize.hpp"
#include "analysis/initialize.hpp"
#include "bond_breakage/initialize.hpp"
#include "cell_system/initialize.hpp"
#include "cluster_analysis/initialize.hpp"
#include "code_info/initialize.hpp"
#include "collision_detection/initialize.hpp"
#include "constraints/initialize.hpp"
#include "electrostatics/initialize.hpp"
#include "galilei/initialize.hpp"
#include "h5md/initialize.hpp"
#include "integrators/initialize.hpp"
#include "interactions/initialize.hpp"
#include "lees_edwards/initialize.hpp"
#include "magnetostatics/initialize.hpp"
#include "math/initialize.hpp"
#include "mpiio/initialize.hpp"
#include "observables/initialize.hpp"
#include "pair_criteria/initialize.hpp"
#include "particle_data/initialize.hpp"
#include "reaction_methods/initialize.hpp"
#include "shapes/initialize.hpp"
#include "system/initialize.hpp"
#include "virtual_sites/initialize.hpp"
#include "walberla/initialize.hpp"

namespace ScriptInterface {
void initialize(Utils::Factory<ObjectHandle> *f) {
  Accumulators::initialize(f);
  Analysis::initialize(f);
  BondBreakage::initialize(f);
  CellSystem::initialize(f);
  ClusterAnalysis::initialize(f);
  CodeInfo::initialize(f);
  CollisionDetection::initialize(f);
  Constraints::initialize(f);
  Coulomb::initialize(f);
  Dipoles::initialize(f);
  Galilei::initialize(f);
  Integrators::initialize(f);
  Interactions::initialize(f);
  LeesEdwards::initialize(f);
  Math::initialize(f);
  MPIIO::initialize(f);
  Observables::initialize(f);
  PairCriteria::initialize(f);
  Particles::initialize(f);
  Shapes::initialize(f);
  System::initialize(f);
  VirtualSites::initialize(f);
  ReactionMethods::initialize(f);
#ifdef H5MD
  Writer::initialize(f);
#endif
#ifdef WALBERLA
  walberla::initialize(f);
#endif
}

} // namespace ScriptInterface
