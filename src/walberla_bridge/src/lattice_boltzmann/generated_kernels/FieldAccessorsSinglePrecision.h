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

// kernel generated with pystencils v1.2, lbmpy v1.2,
// lbmpy_walberla/pystencils_walberla from waLBerla commit
// 4d10e7f2358fc4a4f7e99195d0f67f0b759ecb6f

/**
 * @file
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
inline std::array<float, 19u>
get(GhostLayerField<float, uint_t{19u}> const *pdf_field, Cell const &cell) {
  float const &xyz0 = pdf_field->get(cell, uint_t{0u});
  std::array<float, 19u> pop;
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

inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                std::array<float, 19u> const &pop, Cell const &cell) {
  float &xyz0 = pdf_field->get(cell, uint_t{0u});
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

inline void broadcast(GhostLayerField<float, uint_t{19u}> *pdf_field,
                      std::array<float, 19u> const &pop) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(pdf_field, {
    float &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
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

inline std::vector<float>
get(GhostLayerField<float, uint_t{19u}> const *pdf_field,
    CellInterval const &ci) {
  std::vector<float> out;
  out.reserve(ci.numCells() * uint_t(19u));
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        float const &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
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

inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                std::vector<float> const &values, CellInterval const &ci) {
  assert(uint_c(values.size()) == ci.numCells() * uint_t(19u));
  auto values_ptr = values.data();
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        float &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
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
inline Vector3<float> get(GhostLayerField<float, uint_t{3u}> const *vec_field,
                          Cell const &cell) {
  const float &xyz0 = vec_field->get(cell, uint_t{0u});
  Vector3<float> vec;
  vec[0] = vec_field->getF(&xyz0, uint_t{0u});
  vec[1] = vec_field->getF(&xyz0, uint_t{1u});
  vec[2] = vec_field->getF(&xyz0, uint_t{2u});
  return vec;
}

inline void set(GhostLayerField<float, uint_t{3u}> *vec_field,
                Vector3<float> const &vec, Cell const &cell) {
  float &xyz0 = vec_field->get(cell, uint_t{0u});
  vec_field->getF(&xyz0, uint_t{0u}) = vec[0u];
  vec_field->getF(&xyz0, uint_t{1u}) = vec[1u];
  vec_field->getF(&xyz0, uint_t{2u}) = vec[2u];
}

inline void add(GhostLayerField<float, uint_t{3u}> *vec_field,
                Vector3<float> const &vec, Cell const &cell) {
  float &xyz0 = vec_field->get(cell, uint_t{0u});
  vec_field->getF(&xyz0, uint_t{0u}) += vec[0u];
  vec_field->getF(&xyz0, uint_t{1u}) += vec[1u];
  vec_field->getF(&xyz0, uint_t{2u}) += vec[2u];
}

inline void broadcast(GhostLayerField<float, uint_t{3u}> *vec_field,
                      Vector3<float> const &vec) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
    float &xyz0 = vec_field->get(x, y, z, uint_t{0u});
    vec_field->getF(&xyz0, uint_t{0u}) = vec[0u];
    vec_field->getF(&xyz0, uint_t{1u}) = vec[1u];
    vec_field->getF(&xyz0, uint_t{2u}) = vec[2u];
  });
}

inline void add_to_all(GhostLayerField<float, uint_t{3u}> *vec_field,
                       Vector3<float> const &vec) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
    float &xyz0 = vec_field->get(x, y, z, uint_t{0u});
    vec_field->getF(&xyz0, uint_t{0u}) += vec[0u];
    vec_field->getF(&xyz0, uint_t{1u}) += vec[1u];
    vec_field->getF(&xyz0, uint_t{2u}) += vec[2u];
  });
}

inline std::vector<float>
get(GhostLayerField<float, uint_t{3u}> const *vec_field,
    CellInterval const &ci) {
  std::vector<float> out;
  out.reserve(ci.numCells() * uint_t(3u));
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        const float &xyz0 = vec_field->get(x, y, z, uint_t{0u});
        out.emplace_back(vec_field->getF(&xyz0, uint_t{0u}));
        out.emplace_back(vec_field->getF(&xyz0, uint_t{1u}));
        out.emplace_back(vec_field->getF(&xyz0, uint_t{2u}));
      }
    }
  }
  return out;
}

inline void set(GhostLayerField<float, uint_t{3u}> *vec_field,
                std::vector<float> const &values, CellInterval const &ci) {
  assert(uint_c(values.size()) == ci.numCells() * uint_t(3u));
  auto values_ptr = values.data();
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        float &xyz0 = vec_field->get(x, y, z, uint_t{0u});
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
inline float get(stencil::Direction const direction,
                 Vector3<float> const &u = Vector3<float>(float(0.0)),
                 float rho = float(1.0)) {

  using namespace stencil;
  switch (direction) {
  case C:
    return rho * -0.33333333333333331f * (u[0] * u[0]) +
           rho * -0.33333333333333331f * (u[1] * u[1]) +
           rho * -0.33333333333333331f * (u[2] * u[2]) +
           rho * 0.33333333333333331f;
  case N:
    return rho * -0.16666666666666666f * (u[0] * u[0]) +
           rho * -0.16666666666666666f * (u[2] * u[2]) +
           rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[1] +
           rho * 0.16666666666666666f * (u[1] * u[1]);
  case S:
    return rho * -0.16666666666666666f * u[1] +
           rho * -0.16666666666666666f * (u[0] * u[0]) +
           rho * -0.16666666666666666f * (u[2] * u[2]) +
           rho * 0.055555555555555552f +
           rho * 0.16666666666666666f * (u[1] * u[1]);
  case W:
    return rho * -0.16666666666666666f * u[0] +
           rho * -0.16666666666666666f * (u[1] * u[1]) +
           rho * -0.16666666666666666f * (u[2] * u[2]) +
           rho * 0.055555555555555552f +
           rho * 0.16666666666666666f * (u[0] * u[0]);
  case E:
    return rho * -0.16666666666666666f * (u[1] * u[1]) +
           rho * -0.16666666666666666f * (u[2] * u[2]) +
           rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[0] +
           rho * 0.16666666666666666f * (u[0] * u[0]);
  case T:
    return rho * -0.16666666666666666f * (u[0] * u[0]) +
           rho * -0.16666666666666666f * (u[1] * u[1]) +
           rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[2] +
           rho * 0.16666666666666666f * (u[2] * u[2]);
  case B:
    return rho * -0.16666666666666666f * u[2] +
           rho * -0.16666666666666666f * (u[0] * u[0]) +
           rho * -0.16666666666666666f * (u[1] * u[1]) +
           rho * 0.055555555555555552f +
           rho * 0.16666666666666666f * (u[2] * u[2]);
  case NW:
    return rho * -0.083333333333333329f * u[0] + rho * -0.25f * u[0] * u[1] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[1] * u[1]);
  case NE:
    return rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
           rho * 0.083333333333333329f * u[1] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.25f * u[0] * u[1];
  case SW:
    return rho * -0.083333333333333329f * u[0] +
           rho * -0.083333333333333329f * u[1] + rho * 0.027777777777777776f +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.25f * u[0] * u[1];
  case SE:
    return rho * -0.083333333333333329f * u[1] + rho * -0.25f * u[0] * u[1] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[1] * u[1]);
  case TN:
    return rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
           rho * 0.083333333333333329f * u[2] +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.083333333333333329f * (u[2] * u[2]) +
           rho * 0.25f * u[1] * u[2];
  case TS:
    return rho * -0.083333333333333329f * u[1] + rho * -0.25f * u[1] * u[2] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[2] +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.083333333333333329f * (u[2] * u[2]);
  case TW:
    return rho * -0.083333333333333329f * u[0] + rho * -0.25f * u[0] * u[2] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[2] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[2] * u[2]);
  case TE:
    return rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
           rho * 0.083333333333333329f * u[2] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[2] * u[2]) +
           rho * 0.25f * u[0] * u[2];
  case BN:
    return rho * -0.083333333333333329f * u[2] + rho * -0.25f * u[1] * u[2] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.083333333333333329f * (u[2] * u[2]);
  case BS:
    return rho * -0.083333333333333329f * u[1] +
           rho * -0.083333333333333329f * u[2] + rho * 0.027777777777777776f +
           rho * 0.083333333333333329f * (u[1] * u[1]) +
           rho * 0.083333333333333329f * (u[2] * u[2]) +
           rho * 0.25f * u[1] * u[2];
  case BW:
    return rho * -0.083333333333333329f * u[0] +
           rho * -0.083333333333333329f * u[2] + rho * 0.027777777777777776f +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[2] * u[2]) +
           rho * 0.25f * u[0] * u[2];
  case BE:
    return rho * -0.083333333333333329f * u[2] + rho * -0.25f * u[0] * u[2] +
           rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
           rho * 0.083333333333333329f * (u[0] * u[0]) +
           rho * 0.083333333333333329f * (u[2] * u[2]);
  default:
    WALBERLA_ABORT("Invalid Direction")
  }
}
} // namespace EquilibriumDistribution

namespace Equilibrium {
inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                Vector3<float> const &u, float const rho, Cell const &cell) {

  float &xyz0 = pdf_field->get(cell, uint_t{0u});
  pdf_field->getF(&xyz0, uint_t{0u}) =
      rho * -0.33333333333333331f * (u[0] * u[0]) +
      rho * -0.33333333333333331f * (u[1] * u[1]) +
      rho * -0.33333333333333331f * (u[2] * u[2]) + rho * 0.33333333333333331f;
  pdf_field->getF(&xyz0, uint_t{1u}) =
      rho * -0.16666666666666666f * (u[0] * u[0]) +
      rho * -0.16666666666666666f * (u[2] * u[2]) +
      rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[1] +
      rho * 0.16666666666666666f * (u[1] * u[1]);
  pdf_field->getF(&xyz0, uint_t{2u}) =
      rho * -0.16666666666666666f * u[1] +
      rho * -0.16666666666666666f * (u[0] * u[0]) +
      rho * -0.16666666666666666f * (u[2] * u[2]) +
      rho * 0.055555555555555552f + rho * 0.16666666666666666f * (u[1] * u[1]);
  pdf_field->getF(&xyz0, uint_t{3u}) =
      rho * -0.16666666666666666f * u[0] +
      rho * -0.16666666666666666f * (u[1] * u[1]) +
      rho * -0.16666666666666666f * (u[2] * u[2]) +
      rho * 0.055555555555555552f + rho * 0.16666666666666666f * (u[0] * u[0]);
  pdf_field->getF(&xyz0, uint_t{4u}) =
      rho * -0.16666666666666666f * (u[1] * u[1]) +
      rho * -0.16666666666666666f * (u[2] * u[2]) +
      rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[0] +
      rho * 0.16666666666666666f * (u[0] * u[0]);
  pdf_field->getF(&xyz0, uint_t{5u}) =
      rho * -0.16666666666666666f * (u[0] * u[0]) +
      rho * -0.16666666666666666f * (u[1] * u[1]) +
      rho * 0.055555555555555552f + rho * 0.16666666666666666f * u[2] +
      rho * 0.16666666666666666f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{6u}) =
      rho * -0.16666666666666666f * u[2] +
      rho * -0.16666666666666666f * (u[0] * u[0]) +
      rho * -0.16666666666666666f * (u[1] * u[1]) +
      rho * 0.055555555555555552f + rho * 0.16666666666666666f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{7u}) =
      rho * -0.083333333333333329f * u[0] + rho * -0.25f * u[0] * u[1] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[1] * u[1]);
  pdf_field->getF(&xyz0, uint_t{8u}) =
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
      rho * 0.083333333333333329f * u[1] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[1] * u[1]) + rho * 0.25f * u[0] * u[1];
  pdf_field->getF(&xyz0, uint_t{9u}) =
      rho * -0.083333333333333329f * u[0] +
      rho * -0.083333333333333329f * u[1] + rho * 0.027777777777777776f +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[1] * u[1]) + rho * 0.25f * u[0] * u[1];
  pdf_field->getF(&xyz0, uint_t{10u}) =
      rho * -0.083333333333333329f * u[1] + rho * -0.25f * u[0] * u[1] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[1] * u[1]);
  pdf_field->getF(&xyz0, uint_t{11u}) =
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
      rho * 0.083333333333333329f * u[2] +
      rho * 0.083333333333333329f * (u[1] * u[1]) +
      rho * 0.083333333333333329f * (u[2] * u[2]) + rho * 0.25f * u[1] * u[2];
  pdf_field->getF(&xyz0, uint_t{12u}) =
      rho * -0.083333333333333329f * u[1] + rho * -0.25f * u[1] * u[2] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[2] +
      rho * 0.083333333333333329f * (u[1] * u[1]) +
      rho * 0.083333333333333329f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{13u}) =
      rho * -0.083333333333333329f * u[0] + rho * -0.25f * u[0] * u[2] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[2] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{14u}) =
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
      rho * 0.083333333333333329f * u[2] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[2] * u[2]) + rho * 0.25f * u[0] * u[2];
  pdf_field->getF(&xyz0, uint_t{15u}) =
      rho * -0.083333333333333329f * u[2] + rho * -0.25f * u[1] * u[2] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[1] +
      rho * 0.083333333333333329f * (u[1] * u[1]) +
      rho * 0.083333333333333329f * (u[2] * u[2]);
  pdf_field->getF(&xyz0, uint_t{16u}) =
      rho * -0.083333333333333329f * u[1] +
      rho * -0.083333333333333329f * u[2] + rho * 0.027777777777777776f +
      rho * 0.083333333333333329f * (u[1] * u[1]) +
      rho * 0.083333333333333329f * (u[2] * u[2]) + rho * 0.25f * u[1] * u[2];
  pdf_field->getF(&xyz0, uint_t{17u}) =
      rho * -0.083333333333333329f * u[0] +
      rho * -0.083333333333333329f * u[2] + rho * 0.027777777777777776f +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[2] * u[2]) + rho * 0.25f * u[0] * u[2];
  pdf_field->getF(&xyz0, uint_t{18u}) =
      rho * -0.083333333333333329f * u[2] + rho * -0.25f * u[0] * u[2] +
      rho * 0.027777777777777776f + rho * 0.083333333333333329f * u[0] +
      rho * 0.083333333333333329f * (u[0] * u[0]) +
      rho * 0.083333333333333329f * (u[2] * u[2]);
}
} // namespace Equilibrium

namespace Density {
inline float get(GhostLayerField<float, uint_t{19u}> const *pdf_field,
                 Cell const &cell) {
  const float &xyz0 = pdf_field->get(cell, uint_t{0u});
  const float f_0 = pdf_field->getF(&xyz0, uint_t{0u});
  const float f_1 = pdf_field->getF(&xyz0, uint_t{1u});
  const float f_2 = pdf_field->getF(&xyz0, uint_t{2u});
  const float f_3 = pdf_field->getF(&xyz0, uint_t{3u});
  const float f_4 = pdf_field->getF(&xyz0, uint_t{4u});
  const float f_5 = pdf_field->getF(&xyz0, uint_t{5u});
  const float f_6 = pdf_field->getF(&xyz0, uint_t{6u});
  const float f_7 = pdf_field->getF(&xyz0, uint_t{7u});
  const float f_8 = pdf_field->getF(&xyz0, uint_t{8u});
  const float f_9 = pdf_field->getF(&xyz0, uint_t{9u});
  const float f_10 = pdf_field->getF(&xyz0, uint_t{10u});
  const float f_11 = pdf_field->getF(&xyz0, uint_t{11u});
  const float f_12 = pdf_field->getF(&xyz0, uint_t{12u});
  const float f_13 = pdf_field->getF(&xyz0, uint_t{13u});
  const float f_14 = pdf_field->getF(&xyz0, uint_t{14u});
  const float f_15 = pdf_field->getF(&xyz0, uint_t{15u});
  const float f_16 = pdf_field->getF(&xyz0, uint_t{16u});
  const float f_17 = pdf_field->getF(&xyz0, uint_t{17u});
  const float f_18 = pdf_field->getF(&xyz0, uint_t{18u});
  const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const float vel1Term = f_1 + f_11 + f_15 + f_7;
  const float vel2Term = f_12 + f_13 + f_5;
  const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                    vel1Term + vel2Term;
  return rho;
}

inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                float const rho_in, Cell const &cell) {
  const float &xyz0 = pdf_field->get(cell, uint_t{0u});
  const float f_0 = pdf_field->getF(&xyz0, uint_t{0u});
  const float f_1 = pdf_field->getF(&xyz0, uint_t{1u});
  const float f_2 = pdf_field->getF(&xyz0, uint_t{2u});
  const float f_3 = pdf_field->getF(&xyz0, uint_t{3u});
  const float f_4 = pdf_field->getF(&xyz0, uint_t{4u});
  const float f_5 = pdf_field->getF(&xyz0, uint_t{5u});
  const float f_6 = pdf_field->getF(&xyz0, uint_t{6u});
  const float f_7 = pdf_field->getF(&xyz0, uint_t{7u});
  const float f_8 = pdf_field->getF(&xyz0, uint_t{8u});
  const float f_9 = pdf_field->getF(&xyz0, uint_t{9u});
  const float f_10 = pdf_field->getF(&xyz0, uint_t{10u});
  const float f_11 = pdf_field->getF(&xyz0, uint_t{11u});
  const float f_12 = pdf_field->getF(&xyz0, uint_t{12u});
  const float f_13 = pdf_field->getF(&xyz0, uint_t{13u});
  const float f_14 = pdf_field->getF(&xyz0, uint_t{14u});
  const float f_15 = pdf_field->getF(&xyz0, uint_t{15u});
  const float f_16 = pdf_field->getF(&xyz0, uint_t{16u});
  const float f_17 = pdf_field->getF(&xyz0, uint_t{17u});
  const float f_18 = pdf_field->getF(&xyz0, uint_t{18u});
  const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
  const float vel1Term = f_1 + f_11 + f_15 + f_7;
  const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
  const float vel2Term = f_12 + f_13 + f_5;
  const float momdensity_2 =
      f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
  const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                    vel1Term + vel2Term;

  // calculate current velocity (before density change)
  const float conversion = float(1) / rho;
  Vector3<float> velocity;
  velocity[0u] = momdensity_0 * conversion;
  velocity[1u] = momdensity_1 * conversion;
  velocity[2u] = momdensity_2 * conversion;

  Equilibrium::set(pdf_field, velocity, rho_in, cell);
}

inline std::vector<float>
get(GhostLayerField<float, uint_t{19u}> const *pdf_field,
    CellInterval const &ci) {
  std::vector<float> out;
  out.reserve(ci.numCells());
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        const float &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
        const float f_0 = pdf_field->getF(&xyz0, uint_t{0u});
        const float f_1 = pdf_field->getF(&xyz0, uint_t{1u});
        const float f_2 = pdf_field->getF(&xyz0, uint_t{2u});
        const float f_3 = pdf_field->getF(&xyz0, uint_t{3u});
        const float f_4 = pdf_field->getF(&xyz0, uint_t{4u});
        const float f_5 = pdf_field->getF(&xyz0, uint_t{5u});
        const float f_6 = pdf_field->getF(&xyz0, uint_t{6u});
        const float f_7 = pdf_field->getF(&xyz0, uint_t{7u});
        const float f_8 = pdf_field->getF(&xyz0, uint_t{8u});
        const float f_9 = pdf_field->getF(&xyz0, uint_t{9u});
        const float f_10 = pdf_field->getF(&xyz0, uint_t{10u});
        const float f_11 = pdf_field->getF(&xyz0, uint_t{11u});
        const float f_12 = pdf_field->getF(&xyz0, uint_t{12u});
        const float f_13 = pdf_field->getF(&xyz0, uint_t{13u});
        const float f_14 = pdf_field->getF(&xyz0, uint_t{14u});
        const float f_15 = pdf_field->getF(&xyz0, uint_t{15u});
        const float f_16 = pdf_field->getF(&xyz0, uint_t{16u});
        const float f_17 = pdf_field->getF(&xyz0, uint_t{17u});
        const float f_18 = pdf_field->getF(&xyz0, uint_t{18u});
        const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const float vel1Term = f_1 + f_11 + f_15 + f_7;
        const float vel2Term = f_12 + f_13 + f_5;
        const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                          vel1Term + vel2Term;
        out.emplace_back(rho);
      }
    }
  }
  return out;
}

inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                std::vector<float> const &values, CellInterval const &ci) {
  assert(uint_c(values.size()) == ci.numCells());
  auto values_it = values.begin();
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        const float &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
        const float f_0 = pdf_field->getF(&xyz0, uint_t{0u});
        const float f_1 = pdf_field->getF(&xyz0, uint_t{1u});
        const float f_2 = pdf_field->getF(&xyz0, uint_t{2u});
        const float f_3 = pdf_field->getF(&xyz0, uint_t{3u});
        const float f_4 = pdf_field->getF(&xyz0, uint_t{4u});
        const float f_5 = pdf_field->getF(&xyz0, uint_t{5u});
        const float f_6 = pdf_field->getF(&xyz0, uint_t{6u});
        const float f_7 = pdf_field->getF(&xyz0, uint_t{7u});
        const float f_8 = pdf_field->getF(&xyz0, uint_t{8u});
        const float f_9 = pdf_field->getF(&xyz0, uint_t{9u});
        const float f_10 = pdf_field->getF(&xyz0, uint_t{10u});
        const float f_11 = pdf_field->getF(&xyz0, uint_t{11u});
        const float f_12 = pdf_field->getF(&xyz0, uint_t{12u});
        const float f_13 = pdf_field->getF(&xyz0, uint_t{13u});
        const float f_14 = pdf_field->getF(&xyz0, uint_t{14u});
        const float f_15 = pdf_field->getF(&xyz0, uint_t{15u});
        const float f_16 = pdf_field->getF(&xyz0, uint_t{16u});
        const float f_17 = pdf_field->getF(&xyz0, uint_t{17u});
        const float f_18 = pdf_field->getF(&xyz0, uint_t{18u});
        const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
        const float vel1Term = f_1 + f_11 + f_15 + f_7;
        const float momdensity_1 =
            -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
        const float vel2Term = f_12 + f_13 + f_5;
        const float momdensity_2 =
            f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
        const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                          vel1Term + vel2Term;

        // calculate current velocity (before density change)
        const float conversion = float(1) / rho;
        Vector3<float> velocity;
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
inline void set(GhostLayerField<float, uint_t{19u}> *pdf_field,
                GhostLayerField<float, uint_t{3u}> const *force_field,
                Vector3<float> const &u, Cell const &cell) {
  const float &xyz0 = pdf_field->get(cell, uint_t{0u});
  const float f_0 = pdf_field->getF(&xyz0, uint_t{0u});
  const float f_1 = pdf_field->getF(&xyz0, uint_t{1u});
  const float f_2 = pdf_field->getF(&xyz0, uint_t{2u});
  const float f_3 = pdf_field->getF(&xyz0, uint_t{3u});
  const float f_4 = pdf_field->getF(&xyz0, uint_t{4u});
  const float f_5 = pdf_field->getF(&xyz0, uint_t{5u});
  const float f_6 = pdf_field->getF(&xyz0, uint_t{6u});
  const float f_7 = pdf_field->getF(&xyz0, uint_t{7u});
  const float f_8 = pdf_field->getF(&xyz0, uint_t{8u});
  const float f_9 = pdf_field->getF(&xyz0, uint_t{9u});
  const float f_10 = pdf_field->getF(&xyz0, uint_t{10u});
  const float f_11 = pdf_field->getF(&xyz0, uint_t{11u});
  const float f_12 = pdf_field->getF(&xyz0, uint_t{12u});
  const float f_13 = pdf_field->getF(&xyz0, uint_t{13u});
  const float f_14 = pdf_field->getF(&xyz0, uint_t{14u});
  const float f_15 = pdf_field->getF(&xyz0, uint_t{15u});
  const float f_16 = pdf_field->getF(&xyz0, uint_t{16u});
  const float f_17 = pdf_field->getF(&xyz0, uint_t{17u});
  const float f_18 = pdf_field->getF(&xyz0, uint_t{18u});
  const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
  const float vel1Term = f_1 + f_11 + f_15 + f_7;
  const float vel2Term = f_12 + f_13 + f_5;
  const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                    vel1Term + vel2Term;

  const auto x = cell.x();
  const auto y = cell.y();
  const auto z = cell.z();
  const float u_0 =
      -force_field->get(x, y, z, 0) * 0.50000000000000000f / rho + u[0];
  const float u_1 =
      -force_field->get(x, y, z, 1) * 0.50000000000000000f / rho + u[1];
  const float u_2 =
      -force_field->get(x, y, z, 2) * 0.50000000000000000f / rho + u[2];

  Equilibrium::set(pdf_field, Vector3<float>(u_0, u_1, u_2), rho, cell);
}
} // namespace Velocity

namespace MomentumDensity {
inline Vector3<float>
reduce(GhostLayerField<float, uint_t{19u}> const *pdf_field,
       GhostLayerField<float, uint_t{3u}> const *force_field) {
  Vector3<float> momentumDensity(float{0});
  WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
    const float &xyz0 = pdf_field->get(x, y, z, uint_t{0u});
    const float f_0 = pdf_field->getF(&xyz0, uint_t{0u});
    const float f_1 = pdf_field->getF(&xyz0, uint_t{1u});
    const float f_2 = pdf_field->getF(&xyz0, uint_t{2u});
    const float f_3 = pdf_field->getF(&xyz0, uint_t{3u});
    const float f_4 = pdf_field->getF(&xyz0, uint_t{4u});
    const float f_5 = pdf_field->getF(&xyz0, uint_t{5u});
    const float f_6 = pdf_field->getF(&xyz0, uint_t{6u});
    const float f_7 = pdf_field->getF(&xyz0, uint_t{7u});
    const float f_8 = pdf_field->getF(&xyz0, uint_t{8u});
    const float f_9 = pdf_field->getF(&xyz0, uint_t{9u});
    const float f_10 = pdf_field->getF(&xyz0, uint_t{10u});
    const float f_11 = pdf_field->getF(&xyz0, uint_t{11u});
    const float f_12 = pdf_field->getF(&xyz0, uint_t{12u});
    const float f_13 = pdf_field->getF(&xyz0, uint_t{13u});
    const float f_14 = pdf_field->getF(&xyz0, uint_t{14u});
    const float f_15 = pdf_field->getF(&xyz0, uint_t{15u});
    const float f_16 = pdf_field->getF(&xyz0, uint_t{16u});
    const float f_17 = pdf_field->getF(&xyz0, uint_t{17u});
    const float f_18 = pdf_field->getF(&xyz0, uint_t{18u});
    const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
    const float vel1Term = f_1 + f_11 + f_15 + f_7;
    const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
    const float vel2Term = f_12 + f_13 + f_5;
    const float momdensity_2 =
        f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
    const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term +
                      vel1Term + vel2Term;
    const float md_0 =
        force_field->get(x, y, z, 0) * 0.50000000000000000f + momdensity_0;
    const float md_1 =
        force_field->get(x, y, z, 1) * 0.50000000000000000f + momdensity_1;
    const float md_2 =
        force_field->get(x, y, z, 2) * 0.50000000000000000f + momdensity_2;

    momentumDensity[0u] += md_0;
    momentumDensity[1u] += md_1;
    momentumDensity[2u] += md_2;
  });
  return momentumDensity;
}
} // namespace MomentumDensity

namespace PressureTensor {
inline Matrix3<float> get(GhostLayerField<float, uint_t{19u}> const *pdf_field,
                          Cell const &cell) {
  const float &xyz0 = pdf_field->get(cell, uint_t{0u});
  const float f_0 = pdf_field->getF(&xyz0, uint_t{0u});
  const float f_1 = pdf_field->getF(&xyz0, uint_t{1u});
  const float f_2 = pdf_field->getF(&xyz0, uint_t{2u});
  const float f_3 = pdf_field->getF(&xyz0, uint_t{3u});
  const float f_4 = pdf_field->getF(&xyz0, uint_t{4u});
  const float f_5 = pdf_field->getF(&xyz0, uint_t{5u});
  const float f_6 = pdf_field->getF(&xyz0, uint_t{6u});
  const float f_7 = pdf_field->getF(&xyz0, uint_t{7u});
  const float f_8 = pdf_field->getF(&xyz0, uint_t{8u});
  const float f_9 = pdf_field->getF(&xyz0, uint_t{9u});
  const float f_10 = pdf_field->getF(&xyz0, uint_t{10u});
  const float f_11 = pdf_field->getF(&xyz0, uint_t{11u});
  const float f_12 = pdf_field->getF(&xyz0, uint_t{12u});
  const float f_13 = pdf_field->getF(&xyz0, uint_t{13u});
  const float f_14 = pdf_field->getF(&xyz0, uint_t{14u});
  const float f_15 = pdf_field->getF(&xyz0, uint_t{15u});
  const float f_16 = pdf_field->getF(&xyz0, uint_t{16u});
  const float f_17 = pdf_field->getF(&xyz0, uint_t{17u});
  const float f_18 = pdf_field->getF(&xyz0, uint_t{18u});
  const float p_0 =
      f_10 + f_13 + f_14 + f_17 + f_18 + f_3 + f_4 + f_7 + f_8 + f_9;
  const float p_1 = -f_10 - f_7 + f_8 + f_9;
  const float p_2 = -f_13 + f_14 + f_17 - f_18;
  const float p_3 = -f_10 - f_7 + f_8 + f_9;
  const float p_4 =
      f_1 + f_10 + f_11 + f_12 + f_15 + f_16 + f_2 + f_7 + f_8 + f_9;
  const float p_5 = f_11 - f_12 - f_15 + f_16;
  const float p_6 = -f_13 + f_14 + f_17 - f_18;
  const float p_7 = f_11 - f_12 - f_15 + f_16;
  const float p_8 =
      f_11 + f_12 + f_13 + f_14 + f_15 + f_16 + f_17 + f_18 + f_5 + f_6;

  Matrix3<float> pressureTensor;
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