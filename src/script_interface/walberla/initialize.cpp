/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef WALBERLA

#include "FluidNodeWalberla.hpp"
#include "FluidWalberla.hpp"

#include "EKContainer.hpp"
#include "EKFFT.hpp"
#include "EKNone.hpp"
#include "EKSpecies.hpp"
#include "EKSpeciesNode.hpp"

#include "EKReactant.hpp"
#include "EKReaction.hpp"
#include "EKReactions.hpp"

#include "VTKHandle.hpp"

#include "script_interface/ObjectHandle.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>

#include <utils/Factory.hpp>

#ifdef WALBERLA
#ifdef WALBERLA_STATIC_ASSERT
#error "waLberla headers should not be visible to the ESPResSo script interface"
#endif
#endif

namespace ScriptInterface::walberla {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<LatticeWalberla>("walberla::LatticeWalberla");

  om->register_new<FluidWalberla>("walberla::FluidWalberla");
  om->register_new<FluidNodeWalberla>("walberla::FluidNodeWalberla");
  om->register_new<VTKHandle>("walberla::VTKHandle");

  om->register_new<EKContainer>("walberla::EKContainer");
  om->register_new<EKSpecies>("walberla::EKSpecies");
  om->register_new<EKSpeciesNode>("walberla::EKSpeciesNode");
#ifdef WALBERLA_FFT
  om->register_new<EKFFT>("walberla::EKFFT");
#endif // WALBERLA_FFT
  om->register_new<EKNone>("walberla::None");
  om->register_new<EKVTKHandle>("walberla::EKVTKHandle");

  om->register_new<EKReactant>("walberla::EKReactant");
  om->register_new<EKBulkReaction>("walberla::EKBulkReaction");
  om->register_new<EKIndexedReaction>("walberla::EKIndexedReaction");
  om->register_new<EKReactions>("walberla::EKReactions");
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
