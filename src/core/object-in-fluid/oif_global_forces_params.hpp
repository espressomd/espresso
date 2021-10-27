/*
 * Copyright (C) 2012-2021 The ESPResSo project
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
#ifndef OIF_GLOBAL_FORCES_PARAMS_HPP
#define OIF_GLOBAL_FORCES_PARAMS_HPP

#include <utils/Vector.hpp>

/** Parameters for OIF global forces
 *
 *  Characterize the distribution of the force of the global mesh deformation
 *  onto individual vertices of the mesh.
 */
struct OifGlobalForcesBond {
  /** Relaxed area of the mesh */
  double A0_g;
  /** Area coefficient */
  double ka_g;
  /** Relaxed volume of the mesh */
  double V0;
  /** Volume coefficient */
  double kv;

  double cutoff() const { return 0.; }

  static constexpr int num = 2;

  OifGlobalForcesBond(double A0_g, double ka_g, double V0, double kv) {
    this->ka_g = ka_g;
    this->A0_g = A0_g;
    this->V0 = V0;
    this->kv = kv;
  }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &A0_g;
    ar &ka_g;
    ar &V0;
    ar &kv;
  }
};

#endif
