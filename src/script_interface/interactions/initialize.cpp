/*
 * Copyright (C) 2015-2021 The ESPResSo project
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

#include "BondedInteraction.hpp"
#include "BondedInteractions.hpp"

#include <utils/Factory.hpp>

namespace ScriptInterface {
namespace Interactions {
void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<BondedInteractions>("Interactions::BondedInteractions");
  om->register_new<FeneBond>("Interactions::FeneBond");
  om->register_new<HarmonicBond>("Interactions::HarmonicBond");
  om->register_new<QuarticBond>("Interactions::QuarticBond");
  om->register_new<BondedCoulomb>("Interactions::BondedCoulomb");
  om->register_new<BondedCoulombSR>("Interactions::BondedCoulombSR");
  om->register_new<AngleHarmonicBond>("Interactions::AngleHarmonicBond");
  om->register_new<AngleCosineBond>("Interactions::AngleCosineBond");
  om->register_new<AngleCossquareBond>("Interactions::AngleCossquareBond");
  om->register_new<DihedralBond>("Interactions::DihedralBond");
  om->register_new<TabulatedDistanceBond>(
      "Interactions::TabulatedDistanceBond");
  om->register_new<TabulatedAngleBond>("Interactions::TabulatedAngleBond");
  om->register_new<TabulatedDihedralBond>(
      "Interactions::TabulatedDihedralBond");
  om->register_new<ThermalizedBond>("Interactions::ThermalizedBond");
  om->register_new<RigidBond>("Interactions::RigidBond");
  om->register_new<IBMTriel>("Interactions::IBMTriel");
  om->register_new<IBMVolCons>("Interactions::IBMVolCons");
  om->register_new<IBMTribend>("Interactions::IBMTribend");
  om->register_new<OifGlobalForcesBond>("Interactions::OifGlobalForcesBond");
  om->register_new<OifLocalForcesBond>("Interactions::OifLocalForcesBond");
  om->register_new<VirtualBond>("Interactions::VirtualBond");
}
} // namespace Interactions
} // namespace ScriptInterface
