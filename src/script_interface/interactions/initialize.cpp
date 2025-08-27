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

#include "BondedInteraction.hpp"
#include "BondedInteractions.hpp"
#include "NonBondedInteraction.hpp"
#include "NonBondedInteractions.hpp"

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

  om->register_new<NonBondedInteractions>(
      "Interactions::NonBondedInteractions");
  om->register_new<NonBondedInteractionHandle>(
      "Interactions::NonBondedInteractionHandle");
#ifdef ESPRESSO_LENNARD_JONES
  om->register_new<InteractionLJ>("Interactions::InteractionLJ");
#endif
#ifdef ESPRESSO_LENNARD_JONES_GENERIC
  om->register_new<InteractionLJGen>("Interactions::InteractionLJGen");
#endif
#ifdef ESPRESSO_LJCOS
  om->register_new<InteractionLJcos>("Interactions::InteractionLJcos");
#endif
#ifdef ESPRESSO_LJCOS2
  om->register_new<InteractionLJcos2>("Interactions::InteractionLJcos2");
#endif
#ifdef ESPRESSO_WCA
  om->register_new<InteractionWCA>("Interactions::InteractionWCA");
#endif
#ifdef ESPRESSO_HERTZIAN
  om->register_new<InteractionHertzian>("Interactions::InteractionHertzian");
#endif
#ifdef ESPRESSO_GAUSSIAN
  om->register_new<InteractionGaussian>("Interactions::InteractionGaussian");
#endif
#ifdef ESPRESSO_BMHTF_NACL
  om->register_new<InteractionBMHTF>("Interactions::InteractionBMHTF");
#endif
#ifdef ESPRESSO_MORSE
  om->register_new<InteractionMorse>("Interactions::InteractionMorse");
#endif
#ifdef ESPRESSO_BUCKINGHAM
  om->register_new<InteractionBuckingham>(
      "Interactions::InteractionBuckingham");
#endif
#ifdef ESPRESSO_SOFT_SPHERE
  om->register_new<InteractionSoftSphere>(
      "Interactions::InteractionSoftSphere");
#endif
#ifdef ESPRESSO_HAT
  om->register_new<InteractionHat>("Interactions::InteractionHat");
#endif
#ifdef ESPRESSO_GAY_BERNE
  om->register_new<InteractionGayBerne>("Interactions::InteractionGayBerne");
#endif
#ifdef ESPRESSO_TABULATED
  om->register_new<InteractionTabulated>("Interactions::InteractionTabulated");
#endif
#ifdef ESPRESSO_DPD
  om->register_new<InteractionDPD>("Interactions::InteractionDPD");
#endif
#ifdef ESPRESSO_THOLE
  om->register_new<InteractionThole>("Interactions::InteractionThole");
#endif
#ifdef ESPRESSO_SMOOTH_STEP
  om->register_new<InteractionSmoothStep>(
      "Interactions::InteractionSmoothStep");
#endif
}
} // namespace Interactions
} // namespace ScriptInterface
