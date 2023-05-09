/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include "LatticeWalberla.hpp"

#include "LBFluid.hpp"
#include "LBFluidNode.hpp"
#include "LBFluidSlice.hpp"

#include "EKContainer.hpp"
#include "EKFFT.hpp"
#include "EKNone.hpp"

#include "EKSpecies.hpp"
#include "EKSpeciesNode.hpp"
#include "EKSpeciesSlice.hpp"

#include "EKReactant.hpp"
#include "EKReaction.hpp"
#include "EKReactions.hpp"

#include <script_interface/ObjectHandle.hpp>

#include <utils/Factory.hpp>

#ifdef WALBERLA_STATIC_ASSERT
#error "waLberla headers should not be visible to the ESPResSo script interface"
#endif

namespace ScriptInterface::walberla {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<LatticeWalberla>("walberla::LatticeWalberla");

  om->register_new<LBFluid>("walberla::LBFluid");
  om->register_new<LBFluidNode>("walberla::LBFluidNode");
  om->register_new<LBFluidSlice>("walberla::LBFluidSlice");
  om->register_new<LBVTKHandle>("walberla::LBVTKHandle");

  om->register_new<EKContainer>("walberla::EKContainer");
  om->register_new<EKSpecies>("walberla::EKSpecies");
  om->register_new<EKSpeciesNode>("walberla::EKSpeciesNode");
  om->register_new<EKSpeciesSlice>("walberla::EKSpeciesSlice");
#ifdef WALBERLA_FFT
  om->register_new<EKFFT>("walberla::EKFFT");
#endif // WALBERLA_FFT
  om->register_new<EKNone>("walberla::EKNone");
  om->register_new<EKVTKHandle>("walberla::EKVTKHandle");

  om->register_new<EKReactant>("walberla::EKReactant");
  om->register_new<EKBulkReaction>("walberla::EKBulkReaction");
  om->register_new<EKIndexedReaction>("walberla::EKIndexedReaction");
  om->register_new<EKReactions>("walberla::EKReactions");
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
