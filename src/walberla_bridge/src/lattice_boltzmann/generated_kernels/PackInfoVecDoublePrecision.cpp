//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file PackInfoVecDoublePrecision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3, lbmpy_walberla/pystencils_walberla from waLBerla commit 04f4adbdfc0af983e2d9b72e244d775f37d77034

#include "PackInfoVecDoublePrecision.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "stencil/Directions.h"

#include <cstddef>

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace walberla {
namespace pystencils {

using walberla::cell::CellInterval;
using walberla::stencil::Direction;

namespace internal_05a1eb9a7382e5e7047cdb22e28b6556 {
static FUNC_PREFIX void pack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE(double *RESTRICT _data_buffer, double *RESTRICT const _data_field, int64_t const _size_field_0, int64_t const _size_field_1, int64_t const _size_field_2, int64_t const _stride_field_0, int64_t const _stride_field_1, int64_t const _stride_field_2, int64_t const _stride_field_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_field_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_field_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_field_0; ctr_0 += 1) {
        _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0] = _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2];
        _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0 + 1] = _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2 + _stride_field_3];
        _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0 + 2] = _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2 + 2 * _stride_field_3];
      }
    }
  }
}
} // namespace internal_05a1eb9a7382e5e7047cdb22e28b6556

namespace internal_1ccccad4ca561e07a0934cadb07d0fc1 {
static FUNC_PREFIX void unpack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE(double *RESTRICT const _data_buffer, double *RESTRICT _data_field, int64_t const _size_field_0, int64_t const _size_field_1, int64_t const _size_field_2, int64_t const _stride_field_0, int64_t const _stride_field_1, int64_t const _stride_field_2, int64_t const _stride_field_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_field_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_field_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_field_0; ctr_0 += 1) {
        _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2] = _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0];
        _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2 + _stride_field_3] = _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0 + 1];
        _data_field[_stride_field_0 * ctr_0 + _stride_field_1 * ctr_1 + _stride_field_2 * ctr_2 + 2 * _stride_field_3] = _data_buffer[3 * _size_field_0 * _size_field_1 * ctr_2 + 3 * _size_field_0 * ctr_1 + 3 * ctr_0 + 2];
      }
    }
  }
}
} // namespace internal_1ccccad4ca561e07a0934cadb07d0fc1

void PackInfoVecDoublePrecision::pack(Direction dir, unsigned char *byte_buffer, IBlock *block) const {
  byte_buffer += sizeof(double) - (reinterpret_cast<std::size_t>(byte_buffer) - (reinterpret_cast<std::size_t>(byte_buffer) / sizeof(double)) * sizeof(double));
  double *buffer = reinterpret_cast<double *>(byte_buffer);

  auto field = block->getData<field::GhostLayerField<double, 3>>(fieldID);

  CellInterval ci;
  field->getSliceBeforeGhostLayer(dir, ci, 1, false);

  switch (dir) {
  case stencil::SW:
  case stencil::BW:
  case stencil::W:
  case stencil::TW:
  case stencil::NW:
  case stencil::BS:
  case stencil::S:
  case stencil::TS:
  case stencil::B:
  case stencil::C:
  case stencil::T:
  case stencil::BN:
  case stencil::N:
  case stencil::TN:
  case stencil::SE:
  case stencil::BE:
  case stencil::E:
  case stencil::TE:
  case stencil::NE: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(field->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(field->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(field->nrOfGhostLayers()))
    double *RESTRICT const _data_field = field->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_field_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_field_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_field_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_field_0 = int64_t(field->xStride());
    const int64_t _stride_field_1 = int64_t(field->yStride());
    const int64_t _stride_field_2 = int64_t(field->zStride());
    const int64_t _stride_field_3 = int64_t(1 * int64_t(field->fStride()));
    internal_05a1eb9a7382e5e7047cdb22e28b6556::pack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE(_data_buffer, _data_field, _size_field_0, _size_field_1, _size_field_2, _stride_field_0, _stride_field_1, _stride_field_2, _stride_field_3);
    break;
  }

  default:
    WALBERLA_ASSERT(false);
  }
}

void PackInfoVecDoublePrecision::unpack(Direction dir, unsigned char *byte_buffer, IBlock *block) const {
  byte_buffer += sizeof(double) - (reinterpret_cast<std::size_t>(byte_buffer) - (reinterpret_cast<std::size_t>(byte_buffer) / sizeof(double)) * sizeof(double));
  double *buffer = reinterpret_cast<double *>(byte_buffer);

  auto field = block->getData<field::GhostLayerField<double, 3>>(fieldID);

  CellInterval ci;
  field->getGhostRegion(dir, ci, 1, false);
  auto communciationDirection = stencil::inverseDir[dir];

  switch (communciationDirection) {
  case stencil::SW:
  case stencil::BW:
  case stencil::W:
  case stencil::TW:
  case stencil::NW:
  case stencil::BS:
  case stencil::S:
  case stencil::TS:
  case stencil::B:
  case stencil::C:
  case stencil::T:
  case stencil::BN:
  case stencil::N:
  case stencil::TN:
  case stencil::SE:
  case stencil::BE:
  case stencil::E:
  case stencil::TE:
  case stencil::NE: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(field->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(field->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(field->nrOfGhostLayers()))
    double *RESTRICT _data_field = field->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_field_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_field_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(field->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_field_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_field_0 = int64_t(field->xStride());
    const int64_t _stride_field_1 = int64_t(field->yStride());
    const int64_t _stride_field_2 = int64_t(field->zStride());
    const int64_t _stride_field_3 = int64_t(1 * int64_t(field->fStride()));
    internal_1ccccad4ca561e07a0934cadb07d0fc1::unpack_SW_BW_W_TW_NW_BS_S_TS_B_C_T_BN_N_TN_SE_BE_E_TE_NE(_data_buffer, _data_field, _size_field_0, _size_field_1, _size_field_2, _stride_field_0, _stride_field_1, _stride_field_2, _stride_field_3);
    break;
  }

  default:
    WALBERLA_ASSERT(false);
  }
}

uint_t PackInfoVecDoublePrecision::size(stencil::Direction dir, const IBlock *block) const {
  auto field = block->getData<field::GhostLayerField<double, 3>>(fieldID);

  CellInterval ci;
  field->getGhostRegion(dir, ci, 1, false);

  uint_t elementsPerCell = 0;

  switch (dir) {
  case stencil::SW:
  case stencil::BW:
  case stencil::W:
  case stencil::TW:
  case stencil::NW:
  case stencil::BS:
  case stencil::S:
  case stencil::TS:
  case stencil::B:
  case stencil::C:
  case stencil::T:
  case stencil::BN:
  case stencil::N:
  case stencil::TN:
  case stencil::SE:
  case stencil::BE:
  case stencil::E:
  case stencil::TE:
  case stencil::NE:
    elementsPerCell = 3;
    break;

  default:
    elementsPerCell = 0;
  }
  return ci.numCells() * elementsPerCell * sizeof(double);
}

} // namespace pystencils
} // namespace walberla
