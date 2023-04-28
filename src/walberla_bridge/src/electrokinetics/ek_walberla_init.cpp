/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "EKinWalberlaImpl.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <memory>

std::shared_ptr<EKinWalberlaBase>
new_ek_walberla(std::shared_ptr<LatticeWalberla> const &lattice,
                double diffusion, double kT, double valency,
                Utils::Vector3d ext_efield, double density, bool advection,
                bool friction_coupling, bool single_precision) {
  if (single_precision) {
    return std::make_shared<walberla::EKinWalberlaImpl<13, float>>(
        lattice, diffusion, kT, valency, ext_efield, density, advection,
        friction_coupling);
  }

  return std::make_shared<walberla::EKinWalberlaImpl<13, double>>(
      lattice, diffusion, kT, valency, ext_efield, density, advection,
      friction_coupling);
}
