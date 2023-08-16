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

// kernel generated with pystencils v1.3.1+2.g60e24c4, lbmpy
// v1.3.1+6.gcd1bc2f.dirty, lbmpy_walberla/pystencils_walberla from waLBerla
// commit 065ce5f311850371a97ac4766f47dbb5ca8424ba

/*
 * Lattice field accessors.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/a16141524c58ab88386e2a0f8fdd7c63c5edd704/python/lbmpy_walberla/templates/LatticeModel.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Matrix3.h>
#include <core/math/Vector3.h>

#include <field/GhostLayerField.h>
#include <stencil/D3Q19.h>

#include <array>
#include <cassert>
#include <tuple>
#include <vector>

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
  pop[0u] = pdf_field->getF(&xyz0, uint_t{0u});
  pop[1u] = pdf_field->getF(&xyz0, uint_t{1u});
  pop[2u] = pdf_field->getF(&xyz0, uint_t{2u});
  pop[3u] = pdf_field->getF(&xyz0, uint_t{3u});
  pop[4u] = pdf_field->getF(&xyz0, uint_t{4u});
  pop[5u] = pdf_field->getF(&xyz0, uint_t{5u});
  pop[6u] = pdf_field->getF(&xyz0, uint_t{6u});
  pop[7u] = pdf_field->getF(&xyz0, uint_t{7u});
  pop[8u] = pdf_field->getF(&xyz0, uint_t{8u});
  pop[9u] = pdf_field->getF(&xyz0, uint_t{9u});
  pop[10u] = pdf_field->getF(&xyz0, uint_t{10u});
  pop[11u] = pdf_field->getF(&xyz0, uint_t{11u});
  pop[12u] = pdf_field->getF(&xyz0, uint_t{12u});
  pop[13u] = pdf_field->getF(&xyz0, uint_t{13u});
  pop[14u] = pdf_field->getF(&xyz0, uint_t{14u});
  pop[15u] = pdf_field->getF(&xyz0, uint_t{15u});
  pop[16u] = pdf_field->getF(&xyz0, uint_t{16u});
  pop[17u] = pdf_field->getF(&xyz0, uint_t{17u});
  pop[18u] = pdf_field->getF(&xyz0, uint_t{18u});
  return pop;
}

inline void set(GhostLayerField<double, uint_t{19u}> *pdf_field,
                std::array<double, 19u> const &pop, Cell const &cell) {
  double &xyz0 = pdf_field->get(cell, uint_t{0u});
  pdf_field->getF(&xyz0, uint_t{0u}) = pop[0u];
  pdf_field->getF(&xyz0, uint_t{1u}) = pop[1u];
  pdf_field->getF(&xyz0, uint_t{2u}) = pop[2u];
  pdf_field->getF(&xyz0, uint_t{3u}) = pop[3u];
  pdf_field->getF(&xyz0, uint_t{4u}) = pop[4u];
  pdf_field->getF(&xyz0, uint_t{5u}) = pop[5u];
  pdf_field->getF(&xyz0, uint_t{6u}) = pop[6u];
  pdf_field->getF(&xyz0, uint_t{7u}) = pop[7u];
  pdf_field->getF(&xyz0, uint_t{8u}) = pop[8u];
  pdf_field->getF(&xyz0, uint_t{9u}) = pop[9u];
  pdf_field->getF(&xyz0, uint_t{10u}) = pop[10u];
  pdf_field->getF(&xyz0, uint_t{11u}) = pop[11u];
  pdf_field->getF(&xyz0, uint_t{12u}) = pop[12u];
  pdf_field->getF(&xyz0, uint_t{13u}) = pop[13u];
  pdf_field->getF(&xyz0, uint_t{14u}) = pop[14u];
  pdf_field->getF(&xyz0, uint_t{15u}) = pop[15u];
  pdf_field->getF(&xyz0, uint_t{16u}) = pop[16u];
  pdf_field->getF(&xyz0, uint_t{17u}) = pop[17u];
  pdf_field->getF(&xyz0, uint_t{18u}) = pop[18u];
}

inline void broadcast(GhostLayerField<double, uint_t{19u}> *pdf_field,
                      std::array<double, 19u> const &pop) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(pdf_field, {
    double &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
    pdf_field->getF(&xyz0, uint_t{0u}) = pop[0u];
    pdf_field->getF(&xyz0, uint_t{1u}) = pop[1u];
    pdf_field->getF(&xyz0, uint_t{2u}) = pop[2u];
    pdf_field->getF(&xyz0, uint_t{3u}) = pop[3u];
    pdf_field->getF(&xyz0, uint_t{4u}) = pop[4u];
    pdf_field->getF(&xyz0, uint_t{5u}) = pop[5u];
    pdf_field->getF(&xyz0, uint_t{6u}) = pop[6u];
    pdf_field->getF(&xyz0, uint_t{7u}) = pop[7u];
    pdf_field->getF(&xyz0, uint_t{8u}) = pop[8u];
    pdf_field->getF(&xyz0, uint_t{9u}) = pop[9u];
    pdf_field->getF(&xyz0, uint_t{10u}) = pop[10u];
    pdf_field->getF(&xyz0, uint_t{11u}) = pop[11u];
    pdf_field->getF(&xyz0, uint_t{12u}) = pop[12u];
    pdf_field->getF(&xyz0, uint_t{13u}) = pop[13u];
    pdf_field->getF(&xyz0, uint_t{14u}) = pop[14u];
    pdf_field->getF(&xyz0, uint_t{15u}) = pop[15u];
    pdf_field->getF(&xyz0, uint_t{16u}) = pop[16u];
    pdf_field->getF(&xyz0, uint_t{17u}) = pop[17u];
    pdf_field->getF(&xyz0, uint_t{18u}) = pop[18u];
  });
}

inline std::vector<double>
get(GhostLayerField<double, uint_t{19u}> const *pdf_field,
    CellInterval const &ci) {
  std::vector<double> out;
  out.reserve(ci.numCells() * uint_t(19u));
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        double const &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{0u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{1u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{2u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{3u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{4u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{5u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{6u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{7u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{8u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{9u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{10u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{11u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{12u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{13u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{14u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{15u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{16u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{17u}));
        out.emplace_back(pdf_field->getF(&xyz0, uint_t{18u}));
      }
    }
  }
  return out;
}

inline void set(GhostLayerField<double, uint_t{19u}> *pdf_field,
                std::vector<double> const &values, CellInterval const &ci) {
  assert(uint_c(values.size()) == ci.numCells() * uint_t(19u));
  auto values_ptr = values.data();
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        double &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
        pdf_field->getF(&xyz0, uint_t{0u}) = values_ptr[0u];
        pdf_field->getF(&xyz0, uint_t{1u}) = values_ptr[1u];
        pdf_field->getF(&xyz0, uint_t{2u}) = values_ptr[2u];
        pdf_field->getF(&xyz0, uint_t{3u}) = values_ptr[3u];
        pdf_field->getF(&xyz0, uint_t{4u}) = values_ptr[4u];
        pdf_field->getF(&xyz0, uint_t{5u}) = values_ptr[5u];
        pdf_field->getF(&xyz0, uint_t{6u}) = values_ptr[6u];
        pdf_field->getF(&xyz0, uint_t{7u}) = values_ptr[7u];
        pdf_field->getF(&xyz0, uint_t{8u}) = values_ptr[8u];
        pdf_field->getF(&xyz0, uint_t{9u}) = values_ptr[9u];
        pdf_field->getF(&xyz0, uint_t{10u}) = values_ptr[10u];
        pdf_field->getF(&xyz0, uint_t{11u}) = values_ptr[11u];
        pdf_field->getF(&xyz0, uint_t{12u}) = values_ptr[12u];
        pdf_field->getF(&xyz0, uint_t{13u}) = values_ptr[13u];
        pdf_field->getF(&xyz0, uint_t{14u}) = values_ptr[14u];
        pdf_field->getF(&xyz0, uint_t{15u}) = values_ptr[15u];
        pdf_field->getF(&xyz0, uint_t{16u}) = values_ptr[16u];
        pdf_field->getF(&xyz0, uint_t{17u}) = values_ptr[17u];
        pdf_field->getF(&xyz0, uint_t{18u}) = values_ptr[18u];
        values_ptr += 19u;
      }
    }
  }
}
} // namespace Population

namespace Vector {
inline Vector3<double> get(GhostLayerField<double, uint_t{3u}> const *vec_field,
                           Cell const &cell) {
  const double &xyz0 = vec_field->get(cell, uint_t{0u});
  Vector3<double> vec;
  vec[0] = vec_field->getF(&xyz0, uint_t{0u});
  vec[1] = vec_field->getF(&xyz0, uint_t{1u});
  vec[2] = vec_field->getF(&xyz0, uint_t{2u});
  return vec;
}

inline void set(GhostLayerField<double, uint_t{3u}> *vec_field,
                Vector3<double> const &vec, Cell const &cell) {
  double &xyz0 = vec_field->get(cell, uint_t{0u});
  vec_field->getF(&xyz0, uint_t{0u}) = vec[0u];
  vec_field->getF(&xyz0, uint_t{1u}) = vec[1u];
  vec_field->getF(&xyz0, uint_t{2u}) = vec[2u];
}

inline void add(GhostLayerField<double, uint_t{3u}> *vec_field,
                Vector3<double> const &vec, Cell const &cell) {
  double &xyz0 = vec_field->get(cell, uint_t{0u});
  vec_field->getF(&xyz0, uint_t{0u}) += vec[0u];
  vec_field->getF(&xyz0, uint_t{1u}) += vec[1u];
  vec_field->getF(&xyz0, uint_t{2u}) += vec[2u];
}

inline void broadcast(GhostLayerField<double, uint_t{3u}> *vec_field,
                      Vector3<double> const &vec) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
    double &xyz0 = vec_field->get(x, y, z, uint_t{0u});
    vec_field->getF(&xyz0, uint_t{0u}) = vec[0u];
    vec_field->getF(&xyz0, uint_t{1u}) = vec[1u];
    vec_field->getF(&xyz0, uint_t{2u}) = vec[2u];
  });
}

inline void add_to_all(GhostLayerField<double, uint_t{3u}> *vec_field,
                       Vector3<double> const &vec) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
    double &xyz0 = vec_field->get(x, y, z, uint_t{0u});
    vec_field->getF(&xyz0, uint_t{0u}) += vec[0u];
    vec_field->getF(&xyz0, uint_t{1u}) += vec[1u];
    vec_field->getF(&xyz0, uint_t{2u}) += vec[2u];
  });
}

inline std::vector<double>
get(GhostLayerField<double, uint_t{3u}> const *vec_field,
    CellInterval const &ci) {
  std::vector<double> out;
  out.reserve(ci.numCells() * uint_t(3u));
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        const double &xyz0 = vec_field->get(x, y, z, uint_t{0u});
        out.emplace_back(vec_field->getF(&xyz0, uint_t{0u}));
        out.emplace_back(vec_field->getF(&xyz0, uint_t{1u}));
        out.emplace_back(vec_field->getF(&xyz0, uint_t{2u}));
      }
    }
  }
  return out;
}

inline void set(GhostLayerField<double, uint_t{3u}> *vec_field,
                std::vector<double> const &values, CellInterval const &ci) {
  assert(uint_c(values.size()) == ci.numCells() * uint_t(3u));
  auto values_ptr = values.data();
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        double &xyz0 = vec_field->get(x, y, z, uint_t{0u});
        vec_field->getF(&xyz0, uint_t{0u}) = values_ptr[0u];
        vec_field->getF(&xyz0, uint_t{1u}) = values_ptr[1u];
        vec_field->getF(&xyz0, uint_t{2u}) = values_ptr[2u];
        values_ptr += 3u;
      }
    }
  }
}
} // namespace Vector

namespace EquilibriumDistribution {
inline double get(stencil::Direction const direction,
                  Vector3<double> const &u = Vector3<double>(double(0.0)),
                  double rho = double(0.0)) {
  double delta_rho = rho - double(1.0);

  using namespace stencil;
  switch (direction) {
  case C:
    return delta_rho * 0.33333333333333331 +
           rho * -0.33333333333333331 * (u[0] * u[0]) +
           rho * -0.33333333333333331 * (u[1] * u[1]) +
           rho * -0.33333333333333331 * (u[2] * u[2]);
  case N:
    return delta_rho * 0.055555555555555552 +
           rho * -0.16666666666666666 * (u[0] * u[0]) +
           rho * -0.16666666666666666 * (u[2] * u[2]) +
           rho * 0.16666666666666666 * u[1] +
           rho * 0.16666666666666666 * (u[1] * u[1]);
  case S:
    return delta_rho * 0.055555555555555552 +
           rho * -0.16666666666666666 * u[1] +
           rho * -0.16666666666666666 * (u[0] * u[0]) +
           rho * -0.16666666666666666 * (u[2] * u[2]) +
           rho * 0.16666666666666666 * (u[1] * u[1]);
  case W:
    return delta_rho * 0.055555555555555552 +
           rho * -0.16666666666666666 * u[0] +
           rho * -0.16666666666666666 * (u[1] * u[1]) +
           rho * -0.16666666666666666 * (u[2] * u[2]) +
           rho * 0.16666666666666666 * (u[0] * u[0]);
  case E:
    return delta_rho * 0.055555555555555552 +
           rho * -0.16666666666666666 * (u[1] * u[1]) +
           rho * -0.16666666666666666 * (u[2] * u[2]) +
           rho * 0.16666666666666666 * u[0] +
           rho * 0.16666666666666666 * (u[0] * u[0]);
  case T:
    return delta_rho * 0.055555555555555552 +
           rho * -0.16666666666666666 * (u[0] * u[0]) +
           rho * -0.16666666666666666 * (u[1] * u[1]) +
           rho * 0.16666666666666666 * u[2] +
           rho * 0.16666666666666666 * (u[2] * u[2]);
  case B:
    return delta_rho * 0.055555555555555552 +
           rho * -0.16666666666666666 * u[2] +
           rho * -0.16666666666666666 * (u[0] * u[0]) +
           rho * -0.16666666666666666 * (u[1] * u[1]) +
           rho * 0.16666666666666666 * (u[2] * u[2]);
  case NW:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[1] +
           rho * 0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[1] * u[1]);
  case NE:
    return delta_rho * 0.027777777777777776 +
           rho * 0.083333333333333329 * u[0] +
           rho * 0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.25 * u[0] * u[1];
  case SW:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[0] +
           rho * -0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.25 * u[0] * u[1];
  case SE:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[0] * u[1] +
           rho * 0.083333333333333329 * u[0] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[1] * u[1]);
  case TN:
    return delta_rho * 0.027777777777777776 +
           rho * 0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.083333333333333329 * (u[2] * u[2]) +
           rho * 0.25 * u[1] * u[2];
  case TS:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[1] * u[2] +
           rho * 0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.083333333333333329 * (u[2] * u[2]);
  case TW:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[2] +
           rho * 0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[2] * u[2]);
  case TE:
    return delta_rho * 0.027777777777777776 +
           rho * 0.083333333333333329 * u[0] +
           rho * 0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[2] * u[2]) +
           rho * 0.25 * u[0] * u[2];
  case BN:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[1] * u[2] +
           rho * 0.083333333333333329 * u[1] +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.083333333333333329 * (u[2] * u[2]);
  case BS:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[1] +
           rho * -0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[1] * u[1]) +
           rho * 0.083333333333333329 * (u[2] * u[2]) +
           rho * 0.25 * u[1] * u[2];
  case BW:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[0] +
           rho * -0.083333333333333329 * u[2] +
           rho * 0.083333333333333329 * (u[0] * u[0]) +
           rho * 0.083333333333333329 * (u[2] * u[2]) +
           rho * 0.25 * u[0] * u[2];
  case BE:
    return delta_rho * 0.027777777777777776 +
           rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[0] * u[2] +
           rho * 0.083333333333333329 * u[0] +
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
  double delta_rho = rho - double(1.0);

  double &xyz0 = pdf_field->get(cell, uint_t{0u});
  pdf_field->getF(&xyz0, uint_t{0u}) =
      delta_rho * 0.33333333333333331 +
      rho * -0.33333333333333331 * (u[0] * u[0]) +
      rho * -0.33333333333333331 * (u[1] * u[1]) +
      rho * -0.33333333333333331 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{1u}) =
      delta_rho * 0.055555555555555552 +
      rho * -0.16666666666666666 * (u[0] * u[0]) +
      rho * -0.16666666666666666 * (u[2] * u[2]) +
      rho * 0.16666666666666666 * u[1] +
      rho * 0.16666666666666666 * (u[1] * u[1]);
  pdf_field->getF(&xyz0, uint_t{2u}) =
      delta_rho * 0.055555555555555552 + rho * -0.16666666666666666 * u[1] +
      rho * -0.16666666666666666 * (u[0] * u[0]) +
      rho * -0.16666666666666666 * (u[2] * u[2]) +
      rho * 0.16666666666666666 * (u[1] * u[1]);
  pdf_field->getF(&xyz0, uint_t{3u}) =
      delta_rho * 0.055555555555555552 + rho * -0.16666666666666666 * u[0] +
      rho * -0.16666666666666666 * (u[1] * u[1]) +
      rho * -0.16666666666666666 * (u[2] * u[2]) +
      rho * 0.16666666666666666 * (u[0] * u[0]);
  pdf_field->getF(&xyz0, uint_t{4u}) =
      delta_rho * 0.055555555555555552 +
      rho * -0.16666666666666666 * (u[1] * u[1]) +
      rho * -0.16666666666666666 * (u[2] * u[2]) +
      rho * 0.16666666666666666 * u[0] +
      rho * 0.16666666666666666 * (u[0] * u[0]);
  pdf_field->getF(&xyz0, uint_t{5u}) =
      delta_rho * 0.055555555555555552 +
      rho * -0.16666666666666666 * (u[0] * u[0]) +
      rho * -0.16666666666666666 * (u[1] * u[1]) +
      rho * 0.16666666666666666 * u[2] +
      rho * 0.16666666666666666 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{6u}) =
      delta_rho * 0.055555555555555552 + rho * -0.16666666666666666 * u[2] +
      rho * -0.16666666666666666 * (u[0] * u[0]) +
      rho * -0.16666666666666666 * (u[1] * u[1]) +
      rho * 0.16666666666666666 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{7u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[0] +
      rho * -0.25 * u[0] * u[1] + rho * 0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[1] * u[1]);
  pdf_field->getF(&xyz0, uint_t{8u}) =
      delta_rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
      rho * 0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.25 * u[0] * u[1];
  pdf_field->getF(&xyz0, uint_t{9u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[0] +
      rho * -0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.25 * u[0] * u[1];
  pdf_field->getF(&xyz0, uint_t{10u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[1] +
      rho * -0.25 * u[0] * u[1] + rho * 0.083333333333333329 * u[0] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[1] * u[1]);
  pdf_field->getF(&xyz0, uint_t{11u}) =
      delta_rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[1] * u[1]) +
      rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[1] * u[2];
  pdf_field->getF(&xyz0, uint_t{12u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[1] +
      rho * -0.25 * u[1] * u[2] + rho * 0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[1] * u[1]) +
      rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{13u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[0] +
      rho * -0.25 * u[0] * u[2] + rho * 0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{14u}) =
      delta_rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] +
      rho * 0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[0] * u[2];
  pdf_field->getF(&xyz0, uint_t{15u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[2] +
      rho * -0.25 * u[1] * u[2] + rho * 0.083333333333333329 * u[1] +
      rho * 0.083333333333333329 * (u[1] * u[1]) +
      rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{16u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[1] +
      rho * -0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[1] * u[1]) +
      rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[1] * u[2];
  pdf_field->getF(&xyz0, uint_t{17u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[0] +
      rho * -0.083333333333333329 * u[2] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[0] * u[2];
  pdf_field->getF(&xyz0, uint_t{18u}) =
      delta_rho * 0.027777777777777776 + rho * -0.083333333333333329 * u[2] +
      rho * -0.25 * u[0] * u[2] + rho * 0.083333333333333329 * u[0] +
      rho * 0.083333333333333329 * (u[0] * u[0]) +
      rho * 0.083333333333333329 * (u[2] * u[2]);
}
} // namespace Equilibrium

namespace Density {
inline double get(GhostLayerField<double, uint_t{19u}> const *pdf_field,
                  Cell const &cell) {
  const double &xyz0 = pdf_field->get(cell, uint_t{0u});
  const double f_0 = pdf_field->getF(&xyz0, uint_t{0u});
  const double f_1 = pdf_field->getF(&xyz0, uint_t{1u});
  const double f_2 = pdf_field->getF(&xyz0, uint_t{2u});
  const double f_3 = pdf_field->getF(&xyz0, uint_t{3u});
  const double f_4 = pdf_field->getF(&xyz0, uint_t{4u});
  const double f_5 = pdf_field->getF(&xyz0, uint_t{5u});
  const double f_6 = pdf_field->getF(&xyz0, uint_t{6u});
  const double f_7 = pdf_field->getF(&xyz0, uint_t{7u});
  const double f_8 = pdf_field->getF(&xyz0, uint_t{8u});
  const double f_9 = pdf_field->getF(&xyz0, uint_t{9u});
  const double f_10 = pdf_field->getF(&xyz0, uint_t{10u});
  const double f_11 = pdf_field->getF(&xyz0, uint_t{11u});
  const double f_12 = pdf_field->getF(&xyz0, uint_t{12u});
  const double f_13 = pdf_field->getF(&xyz0, uint_t{13u});
  const double f_14 = pdf_field->getF(&xyz0, uint_t{14u});
  const double f_15 = pdf_field->getF(&xyz0, uint_t{15u});
  const double f_16 = pdf_field->getF(&xyz0, uint_t{16u});
  const double f_17 = pdf_field->getF(&xyz0, uint_t{17u});
  const double f_18 = pdf_field->getF(&xyz0, uint_t{18u});
  const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const double vel1Term = f_1 + f_11 + f_15 + f_7;
  const double vel2Term = f_12 + f_13 + f_5;
  const double delta_rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 +
                           vel0Term + vel1Term + vel2Term;
  const double rho = delta_rho + 1;
  return rho;
}

inline void set(GhostLayerField<double, uint_t{19u}> *pdf_field,
                double const rho_in, Cell const &cell) {
  double const delta_rho_in = rho_in - double(1.0);
  const double &xyz0 = pdf_field->get(cell, uint_t{0u});
  const double f_0 = pdf_field->getF(&xyz0, uint_t{0u});
  const double f_1 = pdf_field->getF(&xyz0, uint_t{1u});
  const double f_2 = pdf_field->getF(&xyz0, uint_t{2u});
  const double f_3 = pdf_field->getF(&xyz0, uint_t{3u});
  const double f_4 = pdf_field->getF(&xyz0, uint_t{4u});
  const double f_5 = pdf_field->getF(&xyz0, uint_t{5u});
  const double f_6 = pdf_field->getF(&xyz0, uint_t{6u});
  const double f_7 = pdf_field->getF(&xyz0, uint_t{7u});
  const double f_8 = pdf_field->getF(&xyz0, uint_t{8u});
  const double f_9 = pdf_field->getF(&xyz0, uint_t{9u});
  const double f_10 = pdf_field->getF(&xyz0, uint_t{10u});
  const double f_11 = pdf_field->getF(&xyz0, uint_t{11u});
  const double f_12 = pdf_field->getF(&xyz0, uint_t{12u});
  const double f_13 = pdf_field->getF(&xyz0, uint_t{13u});
  const double f_14 = pdf_field->getF(&xyz0, uint_t{14u});
  const double f_15 = pdf_field->getF(&xyz0, uint_t{15u});
  const double f_16 = pdf_field->getF(&xyz0, uint_t{16u});
  const double f_17 = pdf_field->getF(&xyz0, uint_t{17u});
  const double f_18 = pdf_field->getF(&xyz0, uint_t{18u});
  const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
  const double vel1Term = f_1 + f_11 + f_15 + f_7;
  const double momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
  const double vel2Term = f_12 + f_13 + f_5;
  const double delta_rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 +
                           vel0Term + vel1Term + vel2Term;
  const double momdensity_2 =
      f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
  const double rho = delta_rho + 1;

  // calculate current velocity (before density change)
  const double conversion = double(1) / rho;
  Vector3<double> velocity;
  velocity[0u] = momdensity_0 * conversion;
  velocity[1u] = momdensity_1 * conversion;
  velocity[2u] = momdensity_2 * conversion;

  Equilibrium::set(pdf_field, velocity, rho_in, cell);
}

inline std::vector<double>
get(GhostLayerField<double, uint_t{19u}> const *pdf_field,
    CellInterval const &ci) {
  std::vector<double> out;
  out.reserve(ci.numCells());
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        const double &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
        const double f_0 = pdf_field->getF(&xyz0, uint_t{0u});
        const double f_1 = pdf_field->getF(&xyz0, uint_t{1u});
        const double f_2 = pdf_field->getF(&xyz0, uint_t{2u});
        const double f_3 = pdf_field->getF(&xyz0, uint_t{3u});
        const double f_4 = pdf_field->getF(&xyz0, uint_t{4u});
        const double f_5 = pdf_field->getF(&xyz0, uint_t{5u});
        const double f_6 = pdf_field->getF(&xyz0, uint_t{6u});
        const double f_7 = pdf_field->getF(&xyz0, uint_t{7u});
        const double f_8 = pdf_field->getF(&xyz0, uint_t{8u});
        const double f_9 = pdf_field->getF(&xyz0, uint_t{9u});
        const double f_10 = pdf_field->getF(&xyz0, uint_t{10u});
        const double f_11 = pdf_field->getF(&xyz0, uint_t{11u});
        const double f_12 = pdf_field->getF(&xyz0, uint_t{12u});
        const double f_13 = pdf_field->getF(&xyz0, uint_t{13u});
        const double f_14 = pdf_field->getF(&xyz0, uint_t{14u});
        const double f_15 = pdf_field->getF(&xyz0, uint_t{15u});
        const double f_16 = pdf_field->getF(&xyz0, uint_t{16u});
        const double f_17 = pdf_field->getF(&xyz0, uint_t{17u});
        const double f_18 = pdf_field->getF(&xyz0, uint_t{18u});
        const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const double vel1Term = f_1 + f_11 + f_15 + f_7;
        const double vel2Term = f_12 + f_13 + f_5;
        const double delta_rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 +
                                 vel0Term + vel1Term + vel2Term;
        const double rho = delta_rho + 1;
        out.emplace_back(rho);
      }
    }
  }
  return out;
}

inline void set(GhostLayerField<double, uint_t{19u}> *pdf_field,
                std::vector<double> const &values, CellInterval const &ci) {
  assert(uint_c(values.size()) == ci.numCells());
  auto values_it = values.begin();
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        const double &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
        const double f_0 = pdf_field->getF(&xyz0, uint_t{0u});
        const double f_1 = pdf_field->getF(&xyz0, uint_t{1u});
        const double f_2 = pdf_field->getF(&xyz0, uint_t{2u});
        const double f_3 = pdf_field->getF(&xyz0, uint_t{3u});
        const double f_4 = pdf_field->getF(&xyz0, uint_t{4u});
        const double f_5 = pdf_field->getF(&xyz0, uint_t{5u});
        const double f_6 = pdf_field->getF(&xyz0, uint_t{6u});
        const double f_7 = pdf_field->getF(&xyz0, uint_t{7u});
        const double f_8 = pdf_field->getF(&xyz0, uint_t{8u});
        const double f_9 = pdf_field->getF(&xyz0, uint_t{9u});
        const double f_10 = pdf_field->getF(&xyz0, uint_t{10u});
        const double f_11 = pdf_field->getF(&xyz0, uint_t{11u});
        const double f_12 = pdf_field->getF(&xyz0, uint_t{12u});
        const double f_13 = pdf_field->getF(&xyz0, uint_t{13u});
        const double f_14 = pdf_field->getF(&xyz0, uint_t{14u});
        const double f_15 = pdf_field->getF(&xyz0, uint_t{15u});
        const double f_16 = pdf_field->getF(&xyz0, uint_t{16u});
        const double f_17 = pdf_field->getF(&xyz0, uint_t{17u});
        const double f_18 = pdf_field->getF(&xyz0, uint_t{18u});
        const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
        const double vel1Term = f_1 + f_11 + f_15 + f_7;
        const double momdensity_1 =
            -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
        const double vel2Term = f_12 + f_13 + f_5;
        const double delta_rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 +
                                 vel0Term + vel1Term + vel2Term;
        const double momdensity_2 =
            f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
        const double rho = delta_rho + 1;

        // calculate current velocity (before density change)
        const double conversion = double(1) / rho;
        Vector3<double> velocity;
        velocity[0u] = momdensity_0 * conversion;
        velocity[1u] = momdensity_1 * conversion;
        velocity[2u] = momdensity_2 * conversion;

        Equilibrium::set(pdf_field, velocity, *values_it, Cell{x, y, z});
        ++values_it;
      }
    }
  }
}
} // namespace Density

namespace Velocity {
inline void set(GhostLayerField<double, uint_t{19u}> *pdf_field,
                GhostLayerField<double, uint_t{3u}> const *force_field,
                Vector3<double> const &u, Cell const &cell) {
  const double &xyz0 = pdf_field->get(cell, uint_t{0u});
  const double f_0 = pdf_field->getF(&xyz0, uint_t{0u});
  const double f_1 = pdf_field->getF(&xyz0, uint_t{1u});
  const double f_2 = pdf_field->getF(&xyz0, uint_t{2u});
  const double f_3 = pdf_field->getF(&xyz0, uint_t{3u});
  const double f_4 = pdf_field->getF(&xyz0, uint_t{4u});
  const double f_5 = pdf_field->getF(&xyz0, uint_t{5u});
  const double f_6 = pdf_field->getF(&xyz0, uint_t{6u});
  const double f_7 = pdf_field->getF(&xyz0, uint_t{7u});
  const double f_8 = pdf_field->getF(&xyz0, uint_t{8u});
  const double f_9 = pdf_field->getF(&xyz0, uint_t{9u});
  const double f_10 = pdf_field->getF(&xyz0, uint_t{10u});
  const double f_11 = pdf_field->getF(&xyz0, uint_t{11u});
  const double f_12 = pdf_field->getF(&xyz0, uint_t{12u});
  const double f_13 = pdf_field->getF(&xyz0, uint_t{13u});
  const double f_14 = pdf_field->getF(&xyz0, uint_t{14u});
  const double f_15 = pdf_field->getF(&xyz0, uint_t{15u});
  const double f_16 = pdf_field->getF(&xyz0, uint_t{16u});
  const double f_17 = pdf_field->getF(&xyz0, uint_t{17u});
  const double f_18 = pdf_field->getF(&xyz0, uint_t{18u});
  const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const double vel1Term = f_1 + f_11 + f_15 + f_7;
  const double vel2Term = f_12 + f_13 + f_5;
  const double delta_rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 +
                           vel0Term + vel1Term + vel2Term;
  const double rho = delta_rho + 1;

  const auto x = cell.x();
  const auto y = cell.y();
  const auto z = cell.z();
  const double u_0 =
      -force_field->get(x, y, z, 0) * 0.50000000000000000 / rho + u[0];
  const double u_1 =
      -force_field->get(x, y, z, 1) * 0.50000000000000000 / rho + u[1];
  const double u_2 =
      -force_field->get(x, y, z, 2) * 0.50000000000000000 / rho + u[2];

  Equilibrium::set(pdf_field, Vector3<double>(u_0, u_1, u_2), rho, cell);
}
} // namespace Velocity

namespace MomentumDensity {
inline Vector3<double>
reduce(GhostLayerField<double, uint_t{19u}> const *pdf_field,
       GhostLayerField<double, uint_t{3u}> const *force_field) {
  Vector3<double> momentumDensity(double{0});
  WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
    const double &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
    const double f_0 = pdf_field->getF(&xyz0, uint_t{0u});
    const double f_1 = pdf_field->getF(&xyz0, uint_t{1u});
    const double f_2 = pdf_field->getF(&xyz0, uint_t{2u});
    const double f_3 = pdf_field->getF(&xyz0, uint_t{3u});
    const double f_4 = pdf_field->getF(&xyz0, uint_t{4u});
    const double f_5 = pdf_field->getF(&xyz0, uint_t{5u});
    const double f_6 = pdf_field->getF(&xyz0, uint_t{6u});
    const double f_7 = pdf_field->getF(&xyz0, uint_t{7u});
    const double f_8 = pdf_field->getF(&xyz0, uint_t{8u});
    const double f_9 = pdf_field->getF(&xyz0, uint_t{9u});
    const double f_10 = pdf_field->getF(&xyz0, uint_t{10u});
    const double f_11 = pdf_field->getF(&xyz0, uint_t{11u});
    const double f_12 = pdf_field->getF(&xyz0, uint_t{12u});
    const double f_13 = pdf_field->getF(&xyz0, uint_t{13u});
    const double f_14 = pdf_field->getF(&xyz0, uint_t{14u});
    const double f_15 = pdf_field->getF(&xyz0, uint_t{15u});
    const double f_16 = pdf_field->getF(&xyz0, uint_t{16u});
    const double f_17 = pdf_field->getF(&xyz0, uint_t{17u});
    const double f_18 = pdf_field->getF(&xyz0, uint_t{18u});
    const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
    const double vel1Term = f_1 + f_11 + f_15 + f_7;
    const double momdensity_1 =
        -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
    const double vel2Term = f_12 + f_13 + f_5;
    const double delta_rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 +
                             vel0Term + vel1Term + vel2Term;
    const double momdensity_2 =
        f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
    const double rho = delta_rho + 1;
    const double md_0 =
        force_field->get(x, y, z, 0) * 0.50000000000000000 + momdensity_0;
    const double md_1 =
        force_field->get(x, y, z, 1) * 0.50000000000000000 + momdensity_1;
    const double md_2 =
        force_field->get(x, y, z, 2) * 0.50000000000000000 + momdensity_2;

    momentumDensity[0u] += md_0;
    momentumDensity[1u] += md_1;
    momentumDensity[2u] += md_2;
  });
  return momentumDensity;
}
} // namespace MomentumDensity

namespace PressureTensor {
inline Matrix3<double>
get(GhostLayerField<double, uint_t{19u}> const *pdf_field, Cell const &cell) {
  const double &xyz0 = pdf_field->get(cell, uint_t{0u});
  const double f_0 = pdf_field->getF(&xyz0, uint_t{0u});
  const double f_1 = pdf_field->getF(&xyz0, uint_t{1u});
  const double f_2 = pdf_field->getF(&xyz0, uint_t{2u});
  const double f_3 = pdf_field->getF(&xyz0, uint_t{3u});
  const double f_4 = pdf_field->getF(&xyz0, uint_t{4u});
  const double f_5 = pdf_field->getF(&xyz0, uint_t{5u});
  const double f_6 = pdf_field->getF(&xyz0, uint_t{6u});
  const double f_7 = pdf_field->getF(&xyz0, uint_t{7u});
  const double f_8 = pdf_field->getF(&xyz0, uint_t{8u});
  const double f_9 = pdf_field->getF(&xyz0, uint_t{9u});
  const double f_10 = pdf_field->getF(&xyz0, uint_t{10u});
  const double f_11 = pdf_field->getF(&xyz0, uint_t{11u});
  const double f_12 = pdf_field->getF(&xyz0, uint_t{12u});
  const double f_13 = pdf_field->getF(&xyz0, uint_t{13u});
  const double f_14 = pdf_field->getF(&xyz0, uint_t{14u});
  const double f_15 = pdf_field->getF(&xyz0, uint_t{15u});
  const double f_16 = pdf_field->getF(&xyz0, uint_t{16u});
  const double f_17 = pdf_field->getF(&xyz0, uint_t{17u});
  const double f_18 = pdf_field->getF(&xyz0, uint_t{18u});
  const double p_0 = f_10 + f_13 + f_14 + f_17 + f_18 + f_3 + f_4 + f_7 + f_8 +
                     f_9 + 0.33333333333333333;
  const double p_1 = -f_10 - f_7 + f_8 + f_9;
  const double p_2 = -f_13 + f_14 + f_17 - f_18;
  const double p_3 = -f_10 - f_7 + f_8 + f_9;
  const double p_4 = f_1 + f_10 + f_11 + f_12 + f_15 + f_16 + f_2 + f_7 + f_8 +
                     f_9 + 0.33333333333333333;
  const double p_5 = f_11 - f_12 - f_15 + f_16;
  const double p_6 = -f_13 + f_14 + f_17 - f_18;
  const double p_7 = f_11 - f_12 - f_15 + f_16;
  const double p_8 = f_11 + f_12 + f_13 + f_14 + f_15 + f_16 + f_17 + f_18 +
                     f_5 + f_6 + 0.33333333333333333;

  Matrix3<double> pressureTensor;
  pressureTensor[0u] = p_0;
  pressureTensor[1u] = p_1;
  pressureTensor[2u] = p_2;

  pressureTensor[3u] = p_3;
  pressureTensor[4u] = p_4;
  pressureTensor[5u] = p_5;

  pressureTensor[6u] = p_6;
  pressureTensor[7u] = p_7;
  pressureTensor[8u] = p_8;

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