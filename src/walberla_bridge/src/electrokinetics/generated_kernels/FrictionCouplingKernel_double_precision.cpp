// kernel generated with pystencils v1.0+25.gfe5cece, lbmpy v1.0+16.g030bd5a,
// lbmpy_walberla/pystencils_walberla from commit
// e1fe2ad1dcbe8f31ea79d95e8a5a5cc0ee3691f3

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
//! \\file FrictionCouplingKernel_double_precision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "FrictionCouplingKernel_double_precision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning push
#pragma warning(disable : 1599)
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_828cb3fbc90c26a23ae54639862a1401 {
static FUNC_PREFIX void
frictioncouplingkernel_double_precision_frictioncouplingkernel_double_precision(
    double D, double *RESTRICT _data_f, double *RESTRICT const _data_j,
    int64_t const _size_f_0, int64_t const _size_f_1, int64_t const _size_f_2,
    int64_t const _stride_f_0, int64_t const _stride_f_1,
    int64_t const _stride_f_2, int64_t const _stride_f_3,
    int64_t const _stride_j_0, int64_t const _stride_j_1,
    int64_t const _stride_j_2, int64_t const _stride_j_3, double kT) {
  for (int64_t ctr_2 = 1; ctr_2 < _size_f_2 - 1; ctr_2 += 1) {
    double *RESTRICT _data_f_20_30 = _data_f + _stride_f_2 * ctr_2;
    double *RESTRICT _data_j_2m1_36 =
        _data_j + _stride_j_2 * ctr_2 - _stride_j_2 + 6 * _stride_j_3;
    double *RESTRICT _data_j_2m1_310 =
        _data_j + _stride_j_2 * ctr_2 - _stride_j_2 + 10 * _stride_j_3;
    double *RESTRICT _data_j_2m1_312 =
        _data_j + _stride_j_2 * ctr_2 - _stride_j_2 + 12 * _stride_j_3;
    double *RESTRICT _data_j_20_30 = _data_j + _stride_j_2 * ctr_2;
    double *RESTRICT _data_j_20_310 =
        _data_j + _stride_j_2 * ctr_2 + 10 * _stride_j_3;
    double *RESTRICT _data_j_20_311 =
        _data_j + _stride_j_2 * ctr_2 + 11 * _stride_j_3;
    double *RESTRICT _data_j_20_312 =
        _data_j + _stride_j_2 * ctr_2 + 12 * _stride_j_3;
    double *RESTRICT _data_j_20_33 =
        _data_j + _stride_j_2 * ctr_2 + 3 * _stride_j_3;
    double *RESTRICT _data_j_20_34 =
        _data_j + _stride_j_2 * ctr_2 + 4 * _stride_j_3;
    double *RESTRICT _data_j_20_35 =
        _data_j + _stride_j_2 * ctr_2 + 5 * _stride_j_3;
    double *RESTRICT _data_j_20_36 =
        _data_j + _stride_j_2 * ctr_2 + 6 * _stride_j_3;
    double *RESTRICT _data_j_20_39 =
        _data_j + _stride_j_2 * ctr_2 + 9 * _stride_j_3;
    double *RESTRICT _data_j_21_35 =
        _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 5 * _stride_j_3;
    double *RESTRICT _data_j_21_39 =
        _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 9 * _stride_j_3;
    double *RESTRICT _data_j_21_311 =
        _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 11 * _stride_j_3;
    double *RESTRICT _data_f_20_31 =
        _data_f + _stride_f_2 * ctr_2 + _stride_f_3;
    double *RESTRICT _data_j_2m1_38 =
        _data_j + _stride_j_2 * ctr_2 - _stride_j_2 + 8 * _stride_j_3;
    double *RESTRICT _data_j_20_31 =
        _data_j + _stride_j_2 * ctr_2 + _stride_j_3;
    double *RESTRICT _data_j_20_37 =
        _data_j + _stride_j_2 * ctr_2 + 7 * _stride_j_3;
    double *RESTRICT _data_j_20_38 =
        _data_j + _stride_j_2 * ctr_2 + 8 * _stride_j_3;
    double *RESTRICT _data_j_21_37 =
        _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 7 * _stride_j_3;
    double *RESTRICT _data_f_20_32 =
        _data_f + _stride_f_2 * ctr_2 + 2 * _stride_f_3;
    double *RESTRICT _data_j_20_32 =
        _data_j + _stride_j_2 * ctr_2 + 2 * _stride_j_3;
    double *RESTRICT _data_j_21_32 =
        _data_j + _stride_j_2 * ctr_2 + _stride_j_2 + 2 * _stride_j_3;
    for (int64_t ctr_1 = 1; ctr_1 < _size_f_1 - 1; ctr_1 += 1) {
      double *RESTRICT _data_f_20_30_10 = _stride_f_1 * ctr_1 + _data_f_20_30;
      double *RESTRICT _data_j_2m1_36_10 = _stride_j_1 * ctr_1 + _data_j_2m1_36;
      double *RESTRICT _data_j_2m1_310_11 =
          _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_2m1_310;
      double *RESTRICT _data_j_2m1_312_1m1 =
          _stride_j_1 * ctr_1 - _stride_j_1 + _data_j_2m1_312;
      double *RESTRICT _data_j_20_30_10 = _stride_j_1 * ctr_1 + _data_j_20_30;
      double *RESTRICT _data_j_20_310_10 = _stride_j_1 * ctr_1 + _data_j_20_310;
      double *RESTRICT _data_j_20_311_10 = _stride_j_1 * ctr_1 + _data_j_20_311;
      double *RESTRICT _data_j_20_312_10 = _stride_j_1 * ctr_1 + _data_j_20_312;
      double *RESTRICT _data_j_20_33_10 = _stride_j_1 * ctr_1 + _data_j_20_33;
      double *RESTRICT _data_j_20_34_10 = _stride_j_1 * ctr_1 + _data_j_20_34;
      double *RESTRICT _data_j_20_35_10 = _stride_j_1 * ctr_1 + _data_j_20_35;
      double *RESTRICT _data_j_20_36_10 = _stride_j_1 * ctr_1 + _data_j_20_36;
      double *RESTRICT _data_j_20_39_10 = _stride_j_1 * ctr_1 + _data_j_20_39;
      double *RESTRICT _data_j_20_33_11 =
          _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_20_33;
      double *RESTRICT _data_j_20_34_1m1 =
          _stride_j_1 * ctr_1 - _stride_j_1 + _data_j_20_34;
      double *RESTRICT _data_j_21_35_10 = _stride_j_1 * ctr_1 + _data_j_21_35;
      double *RESTRICT _data_j_21_39_11 =
          _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_21_39;
      double *RESTRICT _data_j_21_311_1m1 =
          _stride_j_1 * ctr_1 - _stride_j_1 + _data_j_21_311;
      double *RESTRICT _data_f_20_31_10 = _stride_f_1 * ctr_1 + _data_f_20_31;
      double *RESTRICT _data_j_2m1_38_11 =
          _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_2m1_38;
      double *RESTRICT _data_j_20_31_10 = _stride_j_1 * ctr_1 + _data_j_20_31;
      double *RESTRICT _data_j_20_37_10 = _stride_j_1 * ctr_1 + _data_j_20_37;
      double *RESTRICT _data_j_20_38_10 = _stride_j_1 * ctr_1 + _data_j_20_38;
      double *RESTRICT _data_j_20_31_11 =
          _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_20_31;
      double *RESTRICT _data_j_21_37_11 =
          _stride_j_1 * ctr_1 + _stride_j_1 + _data_j_21_37;
      double *RESTRICT _data_f_20_32_10 = _stride_f_1 * ctr_1 + _data_f_20_32;
      double *RESTRICT _data_j_20_32_10 = _stride_j_1 * ctr_1 + _data_j_20_32;
      double *RESTRICT _data_j_21_32_10 = _stride_j_1 * ctr_1 + _data_j_21_32;
      for (int64_t ctr_0 = 1; ctr_0 < _size_f_0 - 1; ctr_0 += 1) {
        _data_f_20_30_10[_stride_f_0 * ctr_0] =
            kT *
            (-1.0 * _data_j_20_30_10[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_20_30_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_310_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_311_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_312_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_33_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_33_11[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_20_34_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_34_1m1[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_20_35_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_36_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_39_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_21_311_1m1[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_21_35_10[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_21_39_11[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_2m1_310_11[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_2m1_312_1m1[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_2m1_36_10[_stride_j_0 * ctr_0 + _stride_j_0]) *
            0.5 * ((1.0) / (D));
        _data_f_20_31_10[_stride_f_0 * ctr_0] =
            kT *
            (-1.0 * _data_j_20_310_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_31_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_31_11[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_33_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_33_11[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_20_37_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_38_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_39_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_21_37_11[_stride_j_0 * ctr_0] -
             1.0 * _data_j_21_39_11[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_2m1_310_11[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_2m1_38_11[_stride_j_0 * ctr_0] +
             _data_j_20_311_10[_stride_j_0 * ctr_0] +
             _data_j_20_312_10[_stride_j_0 * ctr_0] +
             _data_j_20_34_10[_stride_j_0 * ctr_0] +
             _data_j_20_34_1m1[_stride_j_0 * ctr_0 + _stride_j_0] +
             _data_j_21_311_1m1[_stride_j_0 * ctr_0 + _stride_j_0] +
             _data_j_2m1_312_1m1[_stride_j_0 * ctr_0 + _stride_j_0]) *
            0.5 * ((1.0) / (D));
        _data_f_20_32_10[_stride_f_0 * ctr_0] =
            kT *
            (-1.0 * _data_j_20_311_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_32_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_35_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_37_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_20_39_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_21_311_1m1[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_21_32_10[_stride_j_0 * ctr_0] -
             1.0 * _data_j_21_35_10[_stride_j_0 * ctr_0 + _stride_j_0] -
             1.0 * _data_j_21_37_11[_stride_j_0 * ctr_0] -
             1.0 * _data_j_21_39_11[_stride_j_0 * ctr_0 + _stride_j_0] +
             _data_j_20_310_10[_stride_j_0 * ctr_0] +
             _data_j_20_312_10[_stride_j_0 * ctr_0] +
             _data_j_20_36_10[_stride_j_0 * ctr_0] +
             _data_j_20_38_10[_stride_j_0 * ctr_0] +
             _data_j_2m1_310_11[_stride_j_0 * ctr_0 + _stride_j_0] +
             _data_j_2m1_312_1m1[_stride_j_0 * ctr_0 + _stride_j_0] +
             _data_j_2m1_36_10[_stride_j_0 * ctr_0 + _stride_j_0] +
             _data_j_2m1_38_11[_stride_j_0 * ctr_0]) *
            0.5 * ((1.0) / (D));
      }
    }
  }
}
} // namespace internal_828cb3fbc90c26a23ae54639862a1401

void FrictionCouplingKernel_double_precision::run(IBlock *block) {
  auto f = block->getData<field::GhostLayerField<double, 3>>(fID);
  auto j = block->getData<field::GhostLayerField<double, 13>>(jID);

  auto &kT = this->kT_;
  auto &D = this->D_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(f->nrOfGhostLayers()));
  double *RESTRICT _data_f = f->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()));
  double *RESTRICT const _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(f->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(f->xSize()) + 2));
  const int64_t _size_f_0 = int64_t(cell_idx_c(f->xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(f->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(f->ySize()) + 2));
  const int64_t _size_f_1 = int64_t(cell_idx_c(f->ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(f->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(f->zSize()) + 2));
  const int64_t _size_f_2 = int64_t(cell_idx_c(f->zSize()) + 2);
  const int64_t _stride_f_0 = int64_t(f->xStride());
  const int64_t _stride_f_1 = int64_t(f->yStride());
  const int64_t _stride_f_2 = int64_t(f->zStride());
  const int64_t _stride_f_3 = int64_t(1 * int64_t(f->fStride()));
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  internal_828cb3fbc90c26a23ae54639862a1401::
      frictioncouplingkernel_double_precision_frictioncouplingkernel_double_precision(
          D, _data_f, _data_j, _size_f_0, _size_f_1, _size_f_2, _stride_f_0,
          _stride_f_1, _stride_f_2, _stride_f_3, _stride_j_0, _stride_j_1,
          _stride_j_2, _stride_j_3, kT);
}

void FrictionCouplingKernel_double_precision::runOnCellInterval(
    const shared_ptr<StructuredBlockStorage> &blocks,
    const CellInterval &globalCellInterval, cell_idx_t ghostLayers,
    IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto f = block->getData<field::GhostLayerField<double, 3>>(fID);
  auto j = block->getData<field::GhostLayerField<double, 13>>(jID);

  auto &kT = this->kT_;
  auto &D = this->D_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(f->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(f->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(f->nrOfGhostLayers()));
  double *RESTRICT _data_f =
      f->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()));
  double *RESTRICT const _data_j =
      j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(f->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 2));
  const int64_t _size_f_0 = int64_t(cell_idx_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(f->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 2));
  const int64_t _size_f_1 = int64_t(cell_idx_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(f->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 2));
  const int64_t _size_f_2 = int64_t(cell_idx_c(ci.zSize()) + 2);
  const int64_t _stride_f_0 = int64_t(f->xStride());
  const int64_t _stride_f_1 = int64_t(f->yStride());
  const int64_t _stride_f_2 = int64_t(f->zStride());
  const int64_t _stride_f_3 = int64_t(1 * int64_t(f->fStride()));
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  internal_828cb3fbc90c26a23ae54639862a1401::
      frictioncouplingkernel_double_precision_frictioncouplingkernel_double_precision(
          D, _data_f, _data_j, _size_f_0, _size_f_1, _size_f_2, _stride_f_0,
          _stride_f_1, _stride_f_2, _stride_f_3, _stride_j_0, _stride_j_1,
          _stride_j_2, _stride_j_3, kT);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif