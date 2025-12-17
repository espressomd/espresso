/*
 * Copyright (C) 2021-2025 The ESPResSo project
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

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34,
// sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit
// 007e77e077ad9d22b5eed6f3d3118240993e553c

/*
 * Lattice field accessors.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/a16141524c58ab88386e2a0f8fdd7c63c5edd704/python/lbmpy_walberla/templates/LatticeModel.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Vector3.h>

#include <field/GhostLayerField.h>

#include <array>
#include <cassert>
#include <iterator>
#include <tuple>
#include <vector>

namespace walberla {
namespace ek {
namespace accessor {

namespace Scalar {
inline auto get(GhostLayerField<float, 1u> const *scalar_field,
                Cell const &cell) {
  return scalar_field->get(cell);
}

inline void set(GhostLayerField<float, 1u> *scalar_field, float const &value,
                Cell const &cell) {
  scalar_field->get(cell) = value;
}

inline void add(GhostLayerField<float, 1u> *scalar_field, float const &value,
                Cell const &cell) {
  scalar_field->get(cell) += value;
}

inline void initialize(GhostLayerField<float, 1u> *scalar_field,
                       float const &value) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
      scalar_field, { scalar_field->get(x, y, z) = value; });
}

inline void add_to_all(GhostLayerField<float, 1u> *scalar_field,
                       float const &value) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(
      scalar_field, { scalar_field->get(x, y, z) += value; });
}

inline auto get(GhostLayerField<float, 1u> const *scalar_field,
                CellInterval const &ci) {
  std::vector<float> out;
  out.reserve(ci.numCells());
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        out.emplace_back(scalar_field->get(x, y, z));
      }
    }
  }
  return out;
}

inline void set(GhostLayerField<float, 1u> *scalar_field,
                std::vector<float> const &values, CellInterval const &ci) {
  assert(uint_c(values.size()) == ci.numCells());
  auto values_ptr = values.data();
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        scalar_field->get(x, y, z) = values_ptr[0u];
        std::advance(values_ptr, 1);
      }
    }
  }
}
} // namespace Scalar

namespace Vector {
inline auto get(GhostLayerField<float, uint_t{3u}> const *vec_field,
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

inline void initialize(GhostLayerField<float, uint_t{3u}> *vec_field,
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

inline auto get(GhostLayerField<float, uint_t{3u}> const *vec_field,
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
        std::advance(values_ptr, 3);
      }
    }
  }
}
} // namespace Vector

namespace Flux {
inline auto get(GhostLayerField<float, uint_t{13u}> const *flux_field,
                Cell const &cell) {
  const float &xyz0 = flux_field->get(cell, uint_t{0u});
  std::array<float, 13u> value;
  value[0] = flux_field->getF(&xyz0, uint_t{0u});
  value[1] = flux_field->getF(&xyz0, uint_t{1u});
  value[2] = flux_field->getF(&xyz0, uint_t{2u});
  value[3] = flux_field->getF(&xyz0, uint_t{3u});
  value[4] = flux_field->getF(&xyz0, uint_t{4u});
  value[5] = flux_field->getF(&xyz0, uint_t{5u});
  value[6] = flux_field->getF(&xyz0, uint_t{6u});
  value[7] = flux_field->getF(&xyz0, uint_t{7u});
  value[8] = flux_field->getF(&xyz0, uint_t{8u});
  value[9] = flux_field->getF(&xyz0, uint_t{9u});
  value[10] = flux_field->getF(&xyz0, uint_t{10u});
  value[11] = flux_field->getF(&xyz0, uint_t{11u});
  value[12] = flux_field->getF(&xyz0, uint_t{12u});
  return value;
}

inline void initialize(GhostLayerField<float, uint_t{13u}> *flux_field,
                       std::array<float, 13> const &values) {
  WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flux_field, {
    float &xyz0 = flux_field->get(x, y, z, uint_t{0u});
    flux_field->getF(&xyz0, uint_t{0u}) = values[0u];
    flux_field->getF(&xyz0, uint_t{1u}) = values[1u];
    flux_field->getF(&xyz0, uint_t{2u}) = values[2u];
    flux_field->getF(&xyz0, uint_t{3u}) = values[3u];
    flux_field->getF(&xyz0, uint_t{4u}) = values[4u];
    flux_field->getF(&xyz0, uint_t{5u}) = values[5u];
    flux_field->getF(&xyz0, uint_t{6u}) = values[6u];
    flux_field->getF(&xyz0, uint_t{7u}) = values[7u];
    flux_field->getF(&xyz0, uint_t{8u}) = values[8u];
    flux_field->getF(&xyz0, uint_t{9u}) = values[9u];
    flux_field->getF(&xyz0, uint_t{10u}) = values[10u];
    flux_field->getF(&xyz0, uint_t{11u}) = values[11u];
    flux_field->getF(&xyz0, uint_t{12u}) = values[12u];
  });
}

inline auto get(GhostLayerField<float, uint_t{13u}> const *flux_field,
                CellInterval const &ci) {
  std::vector<float> out;
  out.reserve(ci.numCells() * uint_t(13u));
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        const float &xyz0 = flux_field->get(x, y, z, uint_t{0u});
        out.emplace_back(flux_field->getF(&xyz0, uint_t{0u}));
        out.emplace_back(flux_field->getF(&xyz0, uint_t{1u}));
        out.emplace_back(flux_field->getF(&xyz0, uint_t{2u}));
      }
    }
  }
  return out;
}

inline auto get_vector(GhostLayerField<float, uint_t{13u}> const *flux_field,
                       Cell const &cell) {
  Vector3<float> result = Vector3<float>(0, 0, 0);
  const float &xyz0 = flux_field->get(cell, uint_t{0u});
  std::array<float, 27u> local_value;
  // get fluxes in all directions
  local_value[0] = float(0.0);
  local_value[1] = float(-0.5) * flux_field->getNeighbor(
                                     cell.x(), cell.y(), cell.z(), uint_t{1u},
                                     stencil::Direction(uint_t{1u}));
  local_value[2] = float(0.5) * flux_field->getF(&xyz0, uint_t{1u});
  local_value[3] = float(0.5) * flux_field->getF(&xyz0, uint_t{0u});
  local_value[4] = float(-0.5) * flux_field->getNeighbor(
                                     cell.x(), cell.y(), cell.z(), uint_t{0u},
                                     stencil::Direction(uint_t{4u}));
  local_value[5] = float(-0.5) * flux_field->getNeighbor(
                                     cell.x(), cell.y(), cell.z(), uint_t{2u},
                                     stencil::Direction(uint_t{5u}));
  local_value[6] = float(0.5) * flux_field->getF(&xyz0, uint_t{2u});
  local_value[7] = float(0.5) * flux_field->getF(&xyz0, uint_t{4u});
  local_value[8] = float(-0.5) * flux_field->getNeighbor(
                                     cell.x(), cell.y(), cell.z(), uint_t{3u},
                                     stencil::Direction(uint_t{8u}));
  local_value[9] = float(0.5) * flux_field->getF(&xyz0, uint_t{3u});
  local_value[10] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{4u},
                                      stencil::Direction(uint_t{10u}));
  local_value[11] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{7u},
                                      stencil::Direction(uint_t{11u}));
  local_value[12] = float(0.5) * flux_field->getF(&xyz0, uint_t{8u});
  local_value[13] = float(0.5) * flux_field->getF(&xyz0, uint_t{6u});
  local_value[14] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{5u},
                                      stencil::Direction(uint_t{14u}));
  local_value[15] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{8u},
                                      stencil::Direction(uint_t{15u}));
  local_value[16] = float(0.5) * flux_field->getF(&xyz0, uint_t{7u});
  local_value[17] = float(0.5) * flux_field->getF(&xyz0, uint_t{5u});
  local_value[18] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{6u},
                                      stencil::Direction(uint_t{18u}));
  local_value[19] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{9u},
                                      stencil::Direction(uint_t{19u}));
  local_value[20] = float(0.5) * flux_field->getF(&xyz0, uint_t{12u});
  local_value[21] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{11u},
                                      stencil::Direction(uint_t{21u}));
  local_value[22] = float(0.5) * flux_field->getF(&xyz0, uint_t{10u});
  local_value[23] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{10u},
                                      stencil::Direction(uint_t{23u}));
  local_value[24] = float(0.5) * flux_field->getF(&xyz0, uint_t{11u});
  local_value[25] = float(-0.5) * flux_field->getNeighbor(
                                      cell.x(), cell.y(), cell.z(), uint_t{12u},
                                      stencil::Direction(uint_t{25u}));
  local_value[26] = float(0.5) * flux_field->getF(&xyz0, uint_t{9u});

  result[1] += local_value[1] * 1;

  result[1] += local_value[2] * -1;

  result[0] += local_value[3] * -1;

  result[0] += local_value[4] * 1;

  result[2] += local_value[5] * 1;

  result[2] += local_value[6] * -1;

  result[0] += local_value[7] * -1;
  result[1] += local_value[7] * 1;

  result[0] += local_value[8] * 1;
  result[1] += local_value[8] * 1;

  result[0] += local_value[9] * -1;
  result[1] += local_value[9] * -1;

  result[0] += local_value[10] * 1;
  result[1] += local_value[10] * -1;

  result[1] += local_value[11] * 1;
  result[2] += local_value[11] * 1;

  result[1] += local_value[12] * -1;
  result[2] += local_value[12] * 1;

  result[0] += local_value[13] * -1;
  result[2] += local_value[13] * 1;

  result[0] += local_value[14] * 1;
  result[2] += local_value[14] * 1;

  result[1] += local_value[15] * 1;
  result[2] += local_value[15] * -1;

  result[1] += local_value[16] * -1;
  result[2] += local_value[16] * -1;

  result[0] += local_value[17] * -1;
  result[2] += local_value[17] * -1;

  result[0] += local_value[18] * 1;
  result[2] += local_value[18] * -1;

  result[0] += local_value[19] * 1;
  result[1] += local_value[19] * 1;
  result[2] += local_value[19] * 1;

  result[0] += local_value[20] * -1;
  result[1] += local_value[20] * 1;
  result[2] += local_value[20] * 1;

  result[0] += local_value[21] * 1;
  result[1] += local_value[21] * -1;
  result[2] += local_value[21] * 1;

  result[0] += local_value[22] * -1;
  result[1] += local_value[22] * -1;
  result[2] += local_value[22] * 1;

  result[0] += local_value[23] * 1;
  result[1] += local_value[23] * 1;
  result[2] += local_value[23] * -1;

  result[0] += local_value[24] * -1;
  result[1] += local_value[24] * 1;
  result[2] += local_value[24] * -1;

  result[0] += local_value[25] * 1;
  result[1] += local_value[25] * -1;
  result[2] += local_value[25] * -1;

  result[0] += local_value[26] * -1;
  result[1] += local_value[26] * -1;
  result[2] += local_value[26] * -1;

  return result;
}

inline auto get_vector(GhostLayerField<float, uint_t{13u}> const *flux_field,
                       CellInterval const &ci) {
  std::vector<float> out;
  out.reserve(ci.numCells() * uint_t(3u));
  for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
    for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
      for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
        Vector3<float> result = Vector3<float>(0, 0, 0);
        const float &xyz0 = flux_field->get(x, y, z, uint_t{0u});
        std::array<float, 27u> local_value;
        // get fluxes in all directions
        local_value[0] = float(0.0);
        local_value[1] = float(-0.5) * flux_field->getNeighbor(
                                           x, y, z, uint_t{1u},
                                           stencil::Direction(uint_t{1u}));
        local_value[2] = float(0.5) * flux_field->getF(&xyz0, uint_t{1u});
        local_value[3] = float(0.5) * flux_field->getF(&xyz0, uint_t{0u});
        local_value[4] = float(-0.5) * flux_field->getNeighbor(
                                           x, y, z, uint_t{0u},
                                           stencil::Direction(uint_t{4u}));
        local_value[5] = float(-0.5) * flux_field->getNeighbor(
                                           x, y, z, uint_t{2u},
                                           stencil::Direction(uint_t{5u}));
        local_value[6] = float(0.5) * flux_field->getF(&xyz0, uint_t{2u});
        local_value[7] = float(0.5) * flux_field->getF(&xyz0, uint_t{4u});
        local_value[8] = float(-0.5) * flux_field->getNeighbor(
                                           x, y, z, uint_t{3u},
                                           stencil::Direction(uint_t{8u}));
        local_value[9] = float(0.5) * flux_field->getF(&xyz0, uint_t{3u});
        local_value[10] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{4u},
                                            stencil::Direction(uint_t{10u}));
        local_value[11] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{7u},
                                            stencil::Direction(uint_t{11u}));
        local_value[12] = float(0.5) * flux_field->getF(&xyz0, uint_t{8u});
        local_value[13] = float(0.5) * flux_field->getF(&xyz0, uint_t{6u});
        local_value[14] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{5u},
                                            stencil::Direction(uint_t{14u}));
        local_value[15] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{8u},
                                            stencil::Direction(uint_t{15u}));
        local_value[16] = float(0.5) * flux_field->getF(&xyz0, uint_t{7u});
        local_value[17] = float(0.5) * flux_field->getF(&xyz0, uint_t{5u});
        local_value[18] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{6u},
                                            stencil::Direction(uint_t{18u}));
        local_value[19] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{9u},
                                            stencil::Direction(uint_t{19u}));
        local_value[20] = float(0.5) * flux_field->getF(&xyz0, uint_t{12u});
        local_value[21] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{11u},
                                            stencil::Direction(uint_t{21u}));
        local_value[22] = float(0.5) * flux_field->getF(&xyz0, uint_t{10u});
        local_value[23] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{10u},
                                            stencil::Direction(uint_t{23u}));
        local_value[24] = float(0.5) * flux_field->getF(&xyz0, uint_t{11u});
        local_value[25] = float(-0.5) * flux_field->getNeighbor(
                                            x, y, z, uint_t{12u},
                                            stencil::Direction(uint_t{25u}));
        local_value[26] = float(0.5) * flux_field->getF(&xyz0, uint_t{9u});

        result[1] += local_value[1] * 1;

        result[1] += local_value[2] * -1;

        result[0] += local_value[3] * -1;

        result[0] += local_value[4] * 1;

        result[2] += local_value[5] * 1;

        result[2] += local_value[6] * -1;

        result[0] += local_value[7] * -1;
        result[1] += local_value[7] * 1;

        result[0] += local_value[8] * 1;
        result[1] += local_value[8] * 1;

        result[0] += local_value[9] * -1;
        result[1] += local_value[9] * -1;

        result[0] += local_value[10] * 1;
        result[1] += local_value[10] * -1;

        result[1] += local_value[11] * 1;
        result[2] += local_value[11] * 1;

        result[1] += local_value[12] * -1;
        result[2] += local_value[12] * 1;

        result[0] += local_value[13] * -1;
        result[2] += local_value[13] * 1;

        result[0] += local_value[14] * 1;
        result[2] += local_value[14] * 1;

        result[1] += local_value[15] * 1;
        result[2] += local_value[15] * -1;

        result[1] += local_value[16] * -1;
        result[2] += local_value[16] * -1;

        result[0] += local_value[17] * -1;
        result[2] += local_value[17] * -1;

        result[0] += local_value[18] * 1;
        result[2] += local_value[18] * -1;

        result[0] += local_value[19] * 1;
        result[1] += local_value[19] * 1;
        result[2] += local_value[19] * 1;

        result[0] += local_value[20] * -1;
        result[1] += local_value[20] * 1;
        result[2] += local_value[20] * 1;

        result[0] += local_value[21] * 1;
        result[1] += local_value[21] * -1;
        result[2] += local_value[21] * 1;

        result[0] += local_value[22] * -1;
        result[1] += local_value[22] * -1;
        result[2] += local_value[22] * 1;

        result[0] += local_value[23] * 1;
        result[1] += local_value[23] * 1;
        result[2] += local_value[23] * -1;

        result[0] += local_value[24] * -1;
        result[1] += local_value[24] * 1;
        result[2] += local_value[24] * -1;

        result[0] += local_value[25] * 1;
        result[1] += local_value[25] * -1;
        result[2] += local_value[25] * -1;

        result[0] += local_value[26] * -1;
        result[1] += local_value[26] * -1;
        result[2] += local_value[26] * -1;

        out.emplace_back(result[0u]);
        out.emplace_back(result[1u]);
        out.emplace_back(result[2u]);
      }
    }
  }
  return out;
}
} // namespace Flux

} // namespace accessor
} // namespace ek
} // namespace walberla
