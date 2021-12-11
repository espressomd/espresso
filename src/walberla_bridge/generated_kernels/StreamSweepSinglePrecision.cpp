// kernel generated with pystencils v0.4.3+4.g30da657, lbmpy v0.4.3+2.g0e17e61,
// lbmpy_walberla/pystencils_walberla from commit
// 88f85eb7a979f81d68e76009811aeed53ec3014e

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
//! \\file StreamSweepSinglePrecision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "StreamSweepSinglePrecision.h"
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

namespace internal_streamsweepsingleprecision_streamsweepsingleprecision {
static FUNC_PREFIX void streamsweepsingleprecision_streamsweepsingleprecision(
    float *RESTRICT const _data_force, float *RESTRICT const _data_pdfs,
    float *RESTRICT _data_pdfs_tmp, float *RESTRICT _data_velocity,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_0,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0,
    int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2,
    int64_t const _stride_pdfs_tmp_3, int64_t const _stride_velocity_0,
    int64_t const _stride_velocity_1, int64_t const _stride_velocity_2,
    int64_t const _stride_velocity_3) {
  for (int64_t ctr_2 = 1; ctr_2 < _size_force_2 - 1; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                        _stride_pdfs_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                       _stride_pdfs_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                         _stride_pdfs_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                         _stride_pdfs_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                         _stride_pdfs_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                         _stride_pdfs_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                        _stride_pdfs_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                        _stride_pdfs_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                        _stride_pdfs_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                        _stride_pdfs_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    float *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    float *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    float *RESTRICT _data_pdfs_tmp_20_30 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2;
    float *RESTRICT _data_pdfs_tmp_20_31 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_32 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_33 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_34 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_35 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_36 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_37 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_38 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_39 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_310 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_311 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_312 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_313 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_314 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_315 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_316 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_317 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3;
    float *RESTRICT _data_pdfs_tmp_20_318 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3;
    for (int64_t ctr_1 = 1; ctr_1 < _size_force_1 - 1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_31_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_32_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_2m1_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_35;
      float *RESTRICT _data_pdfs_21_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_36;
      float *RESTRICT _data_pdfs_20_37_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
      float *RESTRICT _data_pdfs_20_38_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_39_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_310_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_2m1_311_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
      float *RESTRICT _data_pdfs_2m1_312_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
      float *RESTRICT _data_pdfs_2m1_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_313;
      float *RESTRICT _data_pdfs_2m1_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_314;
      float *RESTRICT _data_pdfs_21_315_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
      float *RESTRICT _data_pdfs_21_316_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
      float *RESTRICT _data_pdfs_21_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_317;
      float *RESTRICT _data_pdfs_21_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_318;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      float *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      float *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      float *RESTRICT _data_pdfs_tmp_20_30_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_30;
      float *RESTRICT _data_pdfs_tmp_20_31_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_31;
      float *RESTRICT _data_pdfs_tmp_20_32_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_32;
      float *RESTRICT _data_pdfs_tmp_20_33_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_33;
      float *RESTRICT _data_pdfs_tmp_20_34_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_34;
      float *RESTRICT _data_pdfs_tmp_20_35_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_35;
      float *RESTRICT _data_pdfs_tmp_20_36_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_36;
      float *RESTRICT _data_pdfs_tmp_20_37_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_37;
      float *RESTRICT _data_pdfs_tmp_20_38_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_38;
      float *RESTRICT _data_pdfs_tmp_20_39_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_39;
      float *RESTRICT _data_pdfs_tmp_20_310_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_310;
      float *RESTRICT _data_pdfs_tmp_20_311_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_311;
      float *RESTRICT _data_pdfs_tmp_20_312_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_312;
      float *RESTRICT _data_pdfs_tmp_20_313_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_313;
      float *RESTRICT _data_pdfs_tmp_20_314_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_314;
      float *RESTRICT _data_pdfs_tmp_20_315_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_315;
      float *RESTRICT _data_pdfs_tmp_20_316_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_316;
      float *RESTRICT _data_pdfs_tmp_20_317_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_317;
      float *RESTRICT _data_pdfs_tmp_20_318_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_318;
      for (int64_t ctr_0 = 1; ctr_0 < _size_force_0 - 1; ctr_0 += 1) {
        const float streamed_0 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const float streamed_1 = _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0];
        const float streamed_2 = _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0];
        const float streamed_3 =
            _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const float streamed_4 =
            _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const float streamed_5 = _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const float streamed_6 = _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0];
        const float streamed_7 =
            _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const float streamed_8 =
            _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const float streamed_9 =
            _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const float streamed_10 =
            _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const float streamed_11 =
            _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        const float streamed_12 = _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        const float streamed_13 =
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const float streamed_14 =
            _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const float streamed_15 = _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0];
        const float streamed_16 = _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        const float streamed_17 =
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const float streamed_18 =
            _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const float vel0Term =
            streamed_10 + streamed_14 + streamed_18 + streamed_4 + streamed_8;
        const float vel1Term =
            streamed_1 + streamed_11 + streamed_15 + streamed_7;
        const float vel2Term = streamed_12 + streamed_13 + streamed_5;
        const float rho = streamed_0 + streamed_16 + streamed_17 + streamed_2 +
                          streamed_3 + streamed_6 + streamed_9 + vel0Term +
                          vel1Term + vel2Term;
        const float u_0 =
            (-streamed_13 - streamed_17 - streamed_3 - streamed_7 - streamed_9 +
             vel0Term) /
                rho +
            0.5f * _data_force_20_30_10[_stride_force_0 * ctr_0] / rho;
        const float u_1 =
            (-streamed_10 - streamed_12 - streamed_16 - streamed_2 +
             streamed_8 - streamed_9 + vel1Term) /
                rho +
            0.5f * _data_force_20_31_10[_stride_force_0 * ctr_0] / rho;
        const float u_2 =
            (streamed_11 + streamed_14 - streamed_15 - streamed_16 -
             streamed_17 - streamed_18 - streamed_6 + vel2Term) /
                rho +
            0.5f * _data_force_20_32_10[_stride_force_0 * ctr_0] / rho;
        _data_velocity_20_30_10[_stride_velocity_0 * ctr_0] = u_0;
        _data_velocity_20_31_10[_stride_velocity_0 * ctr_0] = u_1;
        _data_velocity_20_32_10[_stride_velocity_0 * ctr_0] = u_2;
        _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_0;
        _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_1;
        _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_2;
        _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_3;
        _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_4;
        _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_5;
        _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_6;
        _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_7;
        _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_8;
        _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_9;
        _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_10;
        _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_11;
        _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_12;
        _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_13;
        _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_14;
        _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_15;
        _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_16;
        _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_17;
        _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0 * ctr_0] = streamed_18;
      }
    }
  }
}
} // namespace internal_streamsweepsingleprecision_streamsweepsingleprecision

void StreamSweepSinglePrecision::run(IBlock *block) {
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  field::GhostLayerField<float, 19> *pdfs_tmp;
  {
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find(pdfs);
    if (it != cache_pdfs_.end()) {
      pdfs_tmp = *it;
    } else {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_.insert(pdfs_tmp);
    }
  }

  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT _data_velocity = velocity->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->xSize()) + 2));
  const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->ySize()) + 2));
  const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->zSize()) + 2));
  const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_streamsweepsingleprecision_streamsweepsingleprecision::
      streamsweepsingleprecision_streamsweepsingleprecision(
          _data_force, _data_pdfs, _data_pdfs_tmp, _data_velocity,
          _size_force_0, _size_force_1, _size_force_2, _stride_force_0,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0,
          _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0,
          _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3,
          _stride_velocity_0, _stride_velocity_1, _stride_velocity_2,
          _stride_velocity_3);
  pdfs->swapDataPointers(pdfs_tmp);
}

void StreamSweepSinglePrecision::runOnCellInterval(
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

  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  field::GhostLayerField<float, 19> *pdfs_tmp;
  {
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find(pdfs);
    if (it != cache_pdfs_.end()) {
      pdfs_tmp = *it;
    } else {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_.insert(pdfs_tmp);
    }
  }

  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1,
                                -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1,
                                -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1,
                                -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force =
      force->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT const _data_pdfs =
      pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1,
                                -int_c(pdfs_tmp->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1,
                                -int_c(pdfs_tmp->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1,
                                -int_c(pdfs_tmp->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs_tmp =
      pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1,
                                -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1,
                                -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1,
                                -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT _data_velocity =
      velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 2));
  const int64_t _size_force_0 = int64_t(cell_idx_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 2));
  const int64_t _size_force_1 = int64_t(cell_idx_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 2));
  const int64_t _size_force_2 = int64_t(cell_idx_c(ci.zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_streamsweepsingleprecision_streamsweepsingleprecision::
      streamsweepsingleprecision_streamsweepsingleprecision(
          _data_force, _data_pdfs, _data_pdfs_tmp, _data_velocity,
          _size_force_0, _size_force_1, _size_force_2, _stride_force_0,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0,
          _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0,
          _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3,
          _stride_velocity_0, _stride_velocity_1, _stride_velocity_2,
          _stride_velocity_3);
  pdfs->swapDataPointers(pdfs_tmp);
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