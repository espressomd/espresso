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

#include "config.hpp"

#ifdef LB_WALBERLA

#include "FluidNodeWalberla.hpp"
#include "FluidWalberla.hpp"
#include "LatticeWalberla.hpp"
#include "VTKHandle.hpp"

#include <script_interface/ObjectHandle.hpp>
#include <utils/Factory.hpp>

namespace ScriptInterface::walberla {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<LatticeWalberla>("walberla::LatticeWalberla");
  om->register_new<FluidWalberla>("walberla::FluidWalberla");
  om->register_new<VTKHandle>("walberla::VTKHandle");
  om->register_new<FluidNodeWalberla>("walberla::FluidNodeWalberla");
}

} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
