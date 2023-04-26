/*
 * Copyright (C) 2023 The ESPResSo project
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

#pragma once

#include <utils/Vector.hpp>

/** Checkpoint data for a EK node. */
struct EKWalberlaNodeState {
  double density;
  bool is_boundary_density;
  double density_boundary;
  bool is_boundary_flux;
  Utils::Vector3d flux_boundary;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &density &is_boundary_density &density_boundary &is_boundary_flux
        &flux_boundary;
  }
};
