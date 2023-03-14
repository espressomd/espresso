// kernel generated with pystencils v1.1.1, lbmpy v1.1.1,
// lbmpy_walberla/pystencils_walberla from commit
// e1fe2ad1dcbe8f31ea79d95e8a5a5cc0ee3691f3

/*
 * Copyright (C) 2021-2023 The ESPResSo project
 * Copyright (C) 2020 The waLBerla project
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

/**
 * @file
 * Lattice field accessors.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/a16141524c58ab88386e2a0f8fdd7c63c5edd704/python/lbmpy_walberla/templates/LatticeModel.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>
#include <core/math/Matrix3.h>
#include <core/math/Vector3.h>

#include <field/GhostLayerField.h>
#include <stencil/D3Q19.h>

#include <array>
#include <tuple>

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif

namespace walberla {
namespace lbm {
namespace accessor {

namespace Population {
inline std::array<double, 19u>
get(GhostLayerField<double, uint_t{19u}> const *pdf_field, Cell const &cell) {
  double const &xyz0 = pdf_field->get(cell, uint_t{0u});
  std::array<double, 19u> pop;
  pop[0] = pdf_field->getF(&xyz0, 0);
  pop[1] = pdf_field->getF(&xyz0, 1);
  pop[2] = pdf_field->getF(&xyz0, 2);
  pop[3] = pdf_field->getF(&xyz0, 3);
  pop[4] = pdf_field->getF(&xyz0, 4);
  pop[5] = pdf_field->getF(&xyz0, 5);
  pop[6] = pdf_field->getF(&xyz0, 6);
  pop[7] = pdf_field->getF(&xyz0, 7);
  pop[8] = pdf_field->getF(&xyz0, 8);
  pop[9] = pdf_field->getF(&xyz0, 9);
  pop[10] = pdf_field->getF(&xyz0, 10);
  pop[11] = pdf_field->getF(&xyz0, 11);
  pop[12] = pdf_field->getF(&xyz0, 12);
  pop[13] = pdf_field->getF(&xyz0, 13);
  pop[14] = pdf_field->getF(&xyz0, 14);
  pop[15] = pdf_field->getF(&xyz0, 15);
  pop[16] = pdf_field->getF(&xyz0, 16);
  pop[17] = pdf_field->getF(&xyz0, 17);
  pop[18] = pdf_field->getF(&xyz0, 18);
  return pop;
}

inline void set(GhostLayerField<double, uint_t{19u}> *pdf_field,
                std::array<double, 19u> const &pop, Cell const &cell) {
  double &xyz0 = pdf_field->get(cell, uint_t{0u});
  pdf_field->getF(&xyz0, 0) = pop[0];
  pdf_field->getF(&xyz0, 1) = pop[1];
  pdf_field->getF(&xyz0, 2) = pop[2];
  pdf_field->getF(&xyz0, 3) = pop[3];
  pdf_field->getF(&xyz0, 4) = pop[4];
  pdf_field->getF(&xyz0, 5) = pop[5];
  pdf_field->getF(&xyz0, 6) = pop[6];
  pdf_field->getF(&xyz0, 7) = pop[7];
  pdf_field->getF(&xyz0, 8) = pop[8];
  pdf_field->getF(&xyz0, 9) = pop[9];
  pdf_field->getF(&xyz0, 10) = pop[10];
  pdf_field->getF(&xyz0, 11) = pop[11];
  pdf_field->getF(&xyz0, 12) = pop[12];
  pdf_field->getF(&xyz0, 13) = pop[13];
  pdf_field->getF(&xyz0, 14) = pop[14];
  pdf_field->getF(&xyz0, 15) = pop[15];
  pdf_field->getF(&xyz0, 16) = pop[16];
  pdf_field->getF(&xyz0, 17) = pop[17];
  pdf_field->getF(&xyz0, 18) = pop[18];
}
} // namespace Population

namespace Vector {
inline Vector3<double> get(GhostLayerField<double, uint_t{3u}> const *vec_field,
                           Cell const &cell) {
  const double &xyz0 = vec_field->get(cell, 0);
  Vector3<double> vec;
  vec[0] = vec_field->getF(&xyz0, 0);
  vec[1] = vec_field->getF(&xyz0, 1);
  vec[2] = vec_field->getF(&xyz0, 2);
  return vec;
}

inline void set(GhostLayerField<double, uint_t{3u}> *vec_field,
                Vector3<double> const &vec, Cell const &cell) {
  double &xyz0 = vec_field->get(cell, 0);
  vec_field->getF(&xyz0, 0) = vec[0];
  vec_field->getF(&xyz0, 1) = vec[1];
  vec_field->getF(&xyz0, 2) = vec[2];
}

inline void add(GhostLayerField<double, uint_t{3u}> *vec_field,
                Vector3<double> const &vec, Cell const &cell) {
  double &xyz0 = vec_field->get(cell, 0);
  vec_field->getF(&xyz0, 0) += vec[0];
  vec_field->getF(&xyz0, 1) += vec[1];
  vec_field->getF(&xyz0, 2) += vec[2];
}

inline void broadcast(GhostLayerField<double, uint_t{3u}> *vec_field,
                      Vector3<double> const &vec) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
    double &xyz0 = vec_field->get(x, y, z, 0);
    vec_field->getF(&xyz0, 0) = vec[0];
    vec_field->getF(&xyz0, 1) = vec[1];
    vec_field->getF(&xyz0, 2) = vec[2];
  });
}

inline void add_to_all(GhostLayerField<double, uint_t{3u}> *vec_field,
                       Vector3<double> const &vec) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
    double &xyz0 = vec_field->get(x, y, z, 0);
    vec_field->getF(&xyz0, 0) += vec[0];
    vec_field->getF(&xyz0, 1) += vec[1];
    vec_field->getF(&xyz0, 2) += vec[2];
  });
}
} // namespace Vector

namespace EquilibriumDistribution {
inline double get(stencil::Direction const direction,
                  Vector3<double> const &u = Vector3<double>(double(0.0)),
                  double rho = double(1.0)) {

  using namespace stencil;
  switch (direction) {
  case C:
    return rho * -0.33333333333333331 * (u[0] * u[0]) +
           rho * -0.33333333333333331 * (u[1] * u[1]) +
           rho * -0.33333333333333331 * (u[2] * u[2]) +
           rho * 0.33333333333333331;
  case N:
    return rho * -0.16666666666666666 * (u[0] * u[0]) +
           rho * -0.16666666666666666 * (u[2] * u[2]) +
           rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[1] +
           rho * 0.16666666666666666 * (u[1] * u[1]);
  case S:
    return rho * -0.16666666666666666 * u[1] +
           rho * -0.16666666666666666 * (u[0] * u[0]) +
           rho * -0.16666666666666666 * (u[2] * u[2]) +
           rho * 0.055555555555555552 +
           rho * 0.16666666666666666 * (u[1] * u[1]);
  case W:
    return rho * -0.16666666666666666 * u[0] +
           rho * -0.16666666666666666 * (u[1] * u[1]) +
           rho * -0.16666666666666666 * (u[2] * u[2]) +
           rho * 0.055555555555555552 +
           rho * 0.16666666666666666 * (u[0] * u[0]);
  case E:
    return rho * -0.16666666666666666 * (u[1] * u[1]) +
           rho * -0.16666666666666666 * (u[2] * u[2]) +
           rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[0] +
           rho * 0.16666666666666666 * (u[0] * u[0]);
  case T:
    return rho * -0.16666666666666666 * (u[0] * u[0]) +
           rho * -0.16666666666666666 * (u[1] * u[1]) +
           rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[2] +
           rho * 0.16666666666666666 * (u[2] * u[2]);
  case B:
    return rho * -0.16666666666666666 * u[2] +
           rho * -0.16666666666666666 * (u[0] * u[0]) +
           rho * -0.16666666666666666 * (u[1] * u[1]) +
           rho * 0.055555555555555552 +
           rho * 0.16666666666666666 * (u[2] * u[2]);
  case NW:
    return rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[1] +
           rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[1] * u[1]);
  case NE:
    return rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
           rho * 0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.25 * u[0] * u[1];
  case SW:
    return rho * -0.083333333333333329 * u[0] +
           rho * -0.083333333333333329 * u[1] + rho * 0.027777777777777776 +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.25 * u[0] * u[1];
  case SE:
    return rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[0] * u[1] +
           rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[1] * u[1]);
  case TN:
    return rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.083333333333333329 * (u[2] * u[2]) +
           rho * 0.25 * u[1] * u[2];
  case TS:
    return rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[1] * u[2] +
           rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.083333333333333329 * (u[2] * u[2]);
  case TW:
    return rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[2] +
           rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[2] * u[2]);
  case TE:
    return rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
           rho * 0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[2] * u[2]) +
           rho * 0.25 * u[0] * u[2];
  case BN:
    return rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[1] * u[2] +
           rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.083333333333333329 * (u[2] * u[2]);
  case BS:
    return rho * -0.083333333333333329 * u[1] +
           rho * -0.083333333333333329 * u[2] + rho * 0.027777777777777776 +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.083333333333333329 * (u[2] * u[2]) +
           rho * 0.25 * u[1] * u[2];
  case BW:
    return rho * -0.083333333333333329 * u[0] +
           rho * -0.083333333333333329 * u[2] + rho * 0.027777777777777776 +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[2] * u[2]) +
           rho * 0.25 * u[0] * u[2];
  case BE:
    return rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[0] * u[2] +
           rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[2] * u[2]);
  default:
    WALBERLA_ABORT("Invalid Direction")
  }
}
} // namespace EquilibriumDistribution

namespace Equilibrium {
inline void set(GhostLayerField<double, uint_t{19u}> *pdf_field,
                Vector3<double> const &u, double const rho, Cell const &cell) {

  double &xyz0 = pdf_field->get(cell, 0);
  pdf_field->getF(&xyz0, 0) = rho * -0.33333333333333331 * (u[0] * u[0]) +
                              rho * -0.33333333333333331 * (u[1] * u[1]) +
                              rho * -0.33333333333333331 * (u[2] * u[2]) +
                              rho * 0.33333333333333331;
  pdf_field->getF(&xyz0, 1) = rho * -0.16666666666666666 * (u[0] * u[0]) +
                              rho * -0.16666666666666666 * (u[2] * u[2]) +
                              rho * 0.055555555555555552 +
                              rho * 0.16666666666666666 * u[1] +
                              rho * 0.16666666666666666 * (u[1] * u[1]);
  pdf_field->getF(&xyz0, 2) = rho * -0.16666666666666666 * u[1] +
                              rho * -0.16666666666666666 * (u[0] * u[0]) +
                              rho * -0.16666666666666666 * (u[2] * u[2]) +
                              rho * 0.055555555555555552 +
                              rho * 0.16666666666666666 * (u[1] * u[1]);
  pdf_field->getF(&xyz0, 3) = rho * -0.16666666666666666 * u[0] +
                              rho * -0.16666666666666666 * (u[1] * u[1]) +
                              rho * -0.16666666666666666 * (u[2] * u[2]) +
                              rho * 0.055555555555555552 +
                              rho * 0.16666666666666666 * (u[0] * u[0]);
  pdf_field->getF(&xyz0, 4) = rho * -0.16666666666666666 * (u[1] * u[1]) +
                              rho * -0.16666666666666666 * (u[2] * u[2]) +
                              rho * 0.055555555555555552 +
                              rho * 0.16666666666666666 * u[0] +
                              rho * 0.16666666666666666 * (u[0] * u[0]);
  pdf_field->getF(&xyz0, 5) = rho * -0.16666666666666666 * (u[0] * u[0]) +
                              rho * -0.16666666666666666 * (u[1] * u[1]) +
                              rho * 0.055555555555555552 +
                              rho * 0.16666666666666666 * u[2] +
                              rho * 0.16666666666666666 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 6) = rho * -0.16666666666666666 * u[2] +
                              rho * -0.16666666666666666 * (u[0] * u[0]) +
                              rho * -0.16666666666666666 * (u[1] * u[1]) +
                              rho * 0.055555555555555552 +
                              rho * 0.16666666666666666 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 7) =
      rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[1] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[1] * u[1]);
  pdf_field->getF(&xyz0, 8) =
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
      rho * 0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.25 * u[0] * u[1];
  pdf_field->getF(&xyz0, 9) =
      rho * -0.083333333333333329 * u[0] + rho * -0.083333333333333329 * u[1] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.25 * u[0] * u[1];
  pdf_field->getF(&xyz0, 10) =
      rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[0] * u[1] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[1] * u[1]);
  pdf_field->getF(&xyz0, 11) =
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[1] * u[1]) +
      rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[1] * u[2];
  pdf_field->getF(&xyz0, 12) =
      rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[1] * u[2] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[1] * u[1]) +
      rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 13) =
      rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[2] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 14) =
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
      rho * 0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[0] * u[2];
  pdf_field->getF(&xyz0, 15) =
      rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[1] * u[2] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * (u[1] * u[1]) +
      rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, 16) =
      rho * -0.083333333333333329 * u[1] + rho * -0.083333333333333329 * u[2] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[1] * u[1]) +
      rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[1] * u[2];
  pdf_field->getF(&xyz0, 17) =
      rho * -0.083333333333333329 * u[0] + rho * -0.083333333333333329 * u[2] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[0] * u[2];
  pdf_field->getF(&xyz0, 18) =
      rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[0] * u[2] +
      rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[2] * u[2]);
}
} // namespace Equilibrium

namespace Density {
inline double get(GhostLayerField<double, uint_t{19u}> const *pdf_field,
                  Cell const &cell) {
  const double &xyz0 = pdf_field->get(cell, 0);
  const double f_0 = pdf_field->getF(&xyz0, 0);
  const double f_1 = pdf_field->getF(&xyz0, 1);
  const double f_2 = pdf_field->getF(&xyz0, 2);
  const double f_3 = pdf_field->getF(&xyz0, 3);
  const double f_4 = pdf_field->getF(&xyz0, 4);
  const double f_5 = pdf_field->getF(&xyz0, 5);
  const double f_6 = pdf_field->getF(&xyz0, 6);
  const double f_7 = pdf_field->getF(&xyz0, 7);
  const double f_8 = pdf_field->getF(&xyz0, 8);
  const double f_9 = pdf_field->getF(&xyz0, 9);
  const double f_10 = pdf_field->getF(&xyz0, 10);
  const double f_11 = pdf_field->getF(&xyz0, 11);
  const double f_12 = pdf_field->getF(&xyz0, 12);
  const double f_13 = pdf_field->getF(&xyz0, 13);
  const double f_14 = pdf_field->getF(&xyz0, 14);
  const double f_15 = pdf_field->getF(&xyz0, 15);
  const double f_16 = pdf_field->getF(&xyz0, 16);
  const double f_17 = pdf_field->getF(&xyz0, 17);
  const double f_18 = pdf_field->getF(&xyz0, 18);
  const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const double vel1Term = f_1 + f_11 + f_15 + f_7;
  const double vel2Term = f_12 + f_13 + f_5;
  const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                     vel1Term + vel2Term;
  return rho;
}
} // namespace Density

namespace DensityAndVelocity {
inline std::tuple<double, Vector3<double>>
get(GhostLayerField<double, uint_t{19u}> const *pdf_field,
    GhostLayerField<double, uint_t{3u}> const *force_field, Cell const &cell) {
  const auto x = cell.x();
  const auto y = cell.y();
  const auto z = cell.z();
  const double &xyz0 = pdf_field->get(cell, 0);
  const double f_0 = pdf_field->getF(&xyz0, 0);
  const double f_1 = pdf_field->getF(&xyz0, 1);
  const double f_2 = pdf_field->getF(&xyz0, 2);
  const double f_3 = pdf_field->getF(&xyz0, 3);
  const double f_4 = pdf_field->getF(&xyz0, 4);
  const double f_5 = pdf_field->getF(&xyz0, 5);
  const double f_6 = pdf_field->getF(&xyz0, 6);
  const double f_7 = pdf_field->getF(&xyz0, 7);
  const double f_8 = pdf_field->getF(&xyz0, 8);
  const double f_9 = pdf_field->getF(&xyz0, 9);
  const double f_10 = pdf_field->getF(&xyz0, 10);
  const double f_11 = pdf_field->getF(&xyz0, 11);
  const double f_12 = pdf_field->getF(&xyz0, 12);
  const double f_13 = pdf_field->getF(&xyz0, 13);
  const double f_14 = pdf_field->getF(&xyz0, 14);
  const double f_15 = pdf_field->getF(&xyz0, 15);
  const double f_16 = pdf_field->getF(&xyz0, 16);
  const double f_17 = pdf_field->getF(&xyz0, 17);
  const double f_18 = pdf_field->getF(&xyz0, 18);
  const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
  const double vel1Term = f_1 + f_11 + f_15 + f_7;
  const double momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
  const double vel2Term = f_12 + f_13 + f_5;
  const double momdensity_2 =
      f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
  const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                     vel1Term + vel2Term;
  const double md_0 =
      force_field->get(x, y, z, 0) * 0.50000000000000000 + momdensity_0;
  const double md_1 =
      force_field->get(x, y, z, 1) * 0.50000000000000000 + momdensity_1;
  const double md_2 =
      force_field->get(x, y, z, 2) * 0.50000000000000000 + momdensity_2;

  const double conversion = double(1) / rho;
  Vector3<double> velocity;
  velocity[0] = md_0 * conversion;
  velocity[1] = md_1 * conversion;
  velocity[2] = md_2 * conversion;

  return {rho, velocity};
}

inline void set(GhostLayerField<double, uint_t{19u}> *pdf_field,
                GhostLayerField<double, uint_t{3u}> const *force_field,
                Vector3<double> const &u, double const rho_in,
                Cell const &cell) {
  const auto x = cell.x();
  const auto y = cell.y();
  const auto z = cell.z();
  const double rho = rho_in;
  const double delta_rho = rho - 1;
  const double u_0 =
      -force_field->get(x, y, z, 0) * 0.50000000000000000 / rho + u[0];
  const double u_1 =
      -force_field->get(x, y, z, 1) * 0.50000000000000000 / rho + u[1];
  const double u_2 =
      -force_field->get(x, y, z, 2) * 0.50000000000000000 / rho + u[2];

  Equilibrium::set(pdf_field, Vector3<double>(u_0, u_1, u_2), rho, cell);
}
} // namespace DensityAndVelocity

namespace DensityAndMomentumDensity {
inline std::tuple<double, Vector3<double>>
get(GhostLayerField<double, uint_t{19u}> const *pdf_field,
    GhostLayerField<double, uint_t{3u}> const *force_field, Cell const &cell) {
  const auto x = cell.x();
  const auto y = cell.y();
  const auto z = cell.z();
  const double &xyz0 = pdf_field->get(cell, 0);
  const double f_0 = pdf_field->getF(&xyz0, 0);
  const double f_1 = pdf_field->getF(&xyz0, 1);
  const double f_2 = pdf_field->getF(&xyz0, 2);
  const double f_3 = pdf_field->getF(&xyz0, 3);
  const double f_4 = pdf_field->getF(&xyz0, 4);
  const double f_5 = pdf_field->getF(&xyz0, 5);
  const double f_6 = pdf_field->getF(&xyz0, 6);
  const double f_7 = pdf_field->getF(&xyz0, 7);
  const double f_8 = pdf_field->getF(&xyz0, 8);
  const double f_9 = pdf_field->getF(&xyz0, 9);
  const double f_10 = pdf_field->getF(&xyz0, 10);
  const double f_11 = pdf_field->getF(&xyz0, 11);
  const double f_12 = pdf_field->getF(&xyz0, 12);
  const double f_13 = pdf_field->getF(&xyz0, 13);
  const double f_14 = pdf_field->getF(&xyz0, 14);
  const double f_15 = pdf_field->getF(&xyz0, 15);
  const double f_16 = pdf_field->getF(&xyz0, 16);
  const double f_17 = pdf_field->getF(&xyz0, 17);
  const double f_18 = pdf_field->getF(&xyz0, 18);
  const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
  const double vel1Term = f_1 + f_11 + f_15 + f_7;
  const double momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
  const double vel2Term = f_12 + f_13 + f_5;
  const double momdensity_2 =
      f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
  const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                     vel1Term + vel2Term;
  const double md_0 =
      force_field->get(x, y, z, 0) * 0.50000000000000000 + momdensity_0;
  const double md_1 =
      force_field->get(x, y, z, 1) * 0.50000000000000000 + momdensity_1;
  const double md_2 =
      force_field->get(x, y, z, 2) * 0.50000000000000000 + momdensity_2;

  Vector3<double> momentumDensity;
  momentumDensity[0] = md_0;
  momentumDensity[1] = md_1;
  momentumDensity[2] = md_2;

  return {rho, momentumDensity};
}
} // namespace DensityAndMomentumDensity

namespace MomentumDensity {
inline Vector3<double>
reduce(GhostLayerField<double, uint_t{19u}> const *pdf_field,
       GhostLayerField<double, uint_t{3u}> const *force_field) {
  Vector3<double> momentumDensity(double{0});
  WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
    const double &xyz0 = pdf_field->get(x, y, z, 0);
    const double f_0 = pdf_field->getF(&xyz0, 0);
    const double f_1 = pdf_field->getF(&xyz0, 1);
    const double f_2 = pdf_field->getF(&xyz0, 2);
    const double f_3 = pdf_field->getF(&xyz0, 3);
    const double f_4 = pdf_field->getF(&xyz0, 4);
    const double f_5 = pdf_field->getF(&xyz0, 5);
    const double f_6 = pdf_field->getF(&xyz0, 6);
    const double f_7 = pdf_field->getF(&xyz0, 7);
    const double f_8 = pdf_field->getF(&xyz0, 8);
    const double f_9 = pdf_field->getF(&xyz0, 9);
    const double f_10 = pdf_field->getF(&xyz0, 10);
    const double f_11 = pdf_field->getF(&xyz0, 11);
    const double f_12 = pdf_field->getF(&xyz0, 12);
    const double f_13 = pdf_field->getF(&xyz0, 13);
    const double f_14 = pdf_field->getF(&xyz0, 14);
    const double f_15 = pdf_field->getF(&xyz0, 15);
    const double f_16 = pdf_field->getF(&xyz0, 16);
    const double f_17 = pdf_field->getF(&xyz0, 17);
    const double f_18 = pdf_field->getF(&xyz0, 18);
    const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
    const double vel1Term = f_1 + f_11 + f_15 + f_7;
    const double momdensity_1 =
        -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
    const double vel2Term = f_12 + f_13 + f_5;
    const double momdensity_2 =
        f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
    const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                       vel1Term + vel2Term;
    const double md_0 =
        force_field->get(x, y, z, 0) * 0.50000000000000000 + momdensity_0;
    const double md_1 =
        force_field->get(x, y, z, 1) * 0.50000000000000000 + momdensity_1;
    const double md_2 =
        force_field->get(x, y, z, 2) * 0.50000000000000000 + momdensity_2;

    momentumDensity[0] += md_0;
    momentumDensity[1] += md_1;
    momentumDensity[2] += md_2;
  });
  return momentumDensity;
}
} // namespace MomentumDensity

namespace PressureTensor {
inline Matrix3<double>
get(GhostLayerField<double, uint_t{19u}> const *pdf_field, Cell const &cell) {
  const double &xyz0 = pdf_field->get(cell, 0);
  const double f_0 = pdf_field->getF(&xyz0, 0);
  const double f_1 = pdf_field->getF(&xyz0, 1);
  const double f_2 = pdf_field->getF(&xyz0, 2);
  const double f_3 = pdf_field->getF(&xyz0, 3);
  const double f_4 = pdf_field->getF(&xyz0, 4);
  const double f_5 = pdf_field->getF(&xyz0, 5);
  const double f_6 = pdf_field->getF(&xyz0, 6);
  const double f_7 = pdf_field->getF(&xyz0, 7);
  const double f_8 = pdf_field->getF(&xyz0, 8);
  const double f_9 = pdf_field->getF(&xyz0, 9);
  const double f_10 = pdf_field->getF(&xyz0, 10);
  const double f_11 = pdf_field->getF(&xyz0, 11);
  const double f_12 = pdf_field->getF(&xyz0, 12);
  const double f_13 = pdf_field->getF(&xyz0, 13);
  const double f_14 = pdf_field->getF(&xyz0, 14);
  const double f_15 = pdf_field->getF(&xyz0, 15);
  const double f_16 = pdf_field->getF(&xyz0, 16);
  const double f_17 = pdf_field->getF(&xyz0, 17);
  const double f_18 = pdf_field->getF(&xyz0, 18);
  const double p_0 =
      f_10 + f_13 + f_14 + f_17 + f_18 + f_3 + f_4 + f_7 + f_8 + f_9;
  const double p_1 = -f_10 - f_7 + f_8 + f_9;
  const double p_2 = -f_13 + f_14 + f_17 - f_18;
  const double p_3 = -f_10 - f_7 + f_8 + f_9;
  const double p_4 =
      f_1 + f_10 + f_11 + f_12 + f_15 + f_16 + f_2 + f_7 + f_8 + f_9;
  const double p_5 = f_11 - f_12 - f_15 + f_16;
  const double p_6 = -f_13 + f_14 + f_17 - f_18;
  const double p_7 = f_11 - f_12 - f_15 + f_16;
  const double p_8 =
      f_11 + f_12 + f_13 + f_14 + f_15 + f_16 + f_17 + f_18 + f_5 + f_6;

  Matrix3<double> pressureTensor;
  pressureTensor[0] = p_0;
  pressureTensor[1] = p_1;
  pressureTensor[2] = p_2;

  pressureTensor[3] = p_3;
  pressureTensor[4] = p_4;
  pressureTensor[5] = p_5;

  pressureTensor[6] = p_6;
  pressureTensor[7] = p_7;
  pressureTensor[8] = p_8;

  return pressureTensor;
}
} // namespace PressureTensor

} // namespace accessor
} // namespace lbm
} // namespace walberla

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif