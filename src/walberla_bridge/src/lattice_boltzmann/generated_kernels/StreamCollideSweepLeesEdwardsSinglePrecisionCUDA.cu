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
//! \\file StreamCollideSweepLeesEdwardsSinglePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 272d4a09ec35da50685afc9586645e1b9984b423

#include <cmath>

#include "StreamCollideSweepLeesEdwardsSinglePrecisionCUDA.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX __global__

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
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

namespace internal_streamcollidesweepleesedwardssingleprecisioncuda_streamcollidesweepleesedwardssingleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void streamcollidesweepleesedwardssingleprecisioncuda_streamcollidesweepleesedwardssingleprecisioncuda(float *RESTRICT const _data_force, float *RESTRICT const _data_pdfs, float *RESTRICT _data_pdfs_tmp, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3, int64_t lebc_bot_index, int64_t lebc_top_index, float omega_bulk, float omega_even, float omega_odd, float omega_shear, float v_s) {
  if (blockDim.x * blockIdx.x + threadIdx.x + 1 < _size_force_0 - 1 && blockDim.y * blockIdx.y + threadIdx.y + 1 < _size_force_1 - 1 && blockDim.z * blockIdx.z + threadIdx.z + 1 < _size_force_2 - 1) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x + 1;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y + 1;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z + 1;
    const float xi_2 = _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    const float xi_3 = xi_2 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const float xi_4 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const float xi_5 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    const float xi_6 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3];
    const float xi_7 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3];
    const float xi_8 = xi_7 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    const float xi_9 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    const float xi_11 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3];
    const float xi_12 = _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    const float xi_13 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3];
    const float xi_14 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    const float xi_15 = -_data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const float xi_16 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    const float xi_20 = omega_bulk * 0.5f;
    const float xi_21 = 0.16666666666666666f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float xi_22 = 0.083333333333333329f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float xi_33 = 0.16666666666666666f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_34 = 0.083333333333333329f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_39 = 0.16666666666666666f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float xi_40 = 0.083333333333333329f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float xi_47 = omega_shear * 0.041666666666666664f;
    const float xi_51 = omega_bulk * 0.041666666666666664f;
    const float xi_58 = 0.25f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float xi_62 = omega_shear * 0.125f;
    const float xi_63 = xi_62 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float xi_97 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2];
    const float xi_102 = xi_11 + xi_3;
    const float xi_104 = 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3] + 5.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] + 5.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
    const float xi_107 = 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3];
    const float xi_108 = 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const float xi_109 = 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3];
    const float xi_113 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const float xi_114 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3];
    const float xi_115 = -_data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3];
    const float xi_116 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3];
    const float xi_117 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3];
    const float xi_122 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3];
    const float xi_123 = xi_14 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3];
    const float xi_124 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3];
    const float xi_125 = xi_124 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const float xi_126 = xi_122 + xi_123 + xi_125 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    const float xi_127 = omega_odd * 0.25f;
    const float xi_128 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    const float xi_129 = 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    const float xi_130 = -2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3];
    const float xi_131 = xi_125 + xi_128 - xi_129 + xi_130 + xi_5 + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3];
    const float xi_132 = omega_odd * 0.083333333333333329f;
    const float xi_133 = xi_131 * xi_132;
    const float xi_134 = xi_126 * xi_127 + xi_133;
    const float xi_144 = xi_15 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3];
    const float xi_145 = xi_115 + xi_144;
    const float xi_146 = xi_145 - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
    const float xi_147 = -_data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
    const float xi_148 = -xi_129 - xi_130 - xi_145 - xi_147 - xi_7 + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    const float xi_149 = xi_132 * xi_148;
    const float xi_150 = xi_127 * xi_146 + xi_149;
    const float xi_152 = xi_122 + xi_128 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    const float xi_153 = -xi_116 - xi_152 - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3];
    const float xi_154 = -xi_107 - xi_108 + xi_109 + xi_117 + xi_152 + xi_6;
    const float xi_155 = xi_132 * xi_154;
    const float xi_156 = xi_127 * xi_153 + xi_155;
    const float xi_158 = omega_shear * 0.25f;
    const float xi_161 = -xi_149;
    const float xi_173 = omega_odd * 0.041666666666666664f;
    const float xi_174 = xi_154 * xi_173;
    const float xi_175 = omega_odd * 0.125f;
    const float xi_176 = xi_153 * xi_175;
    const float xi_177 = -xi_174 + xi_176;
    const float xi_178 = xi_131 * xi_173;
    const float xi_179 = xi_126 * xi_175;
    const float xi_180 = -xi_178 + xi_179;
    const float xi_181 = xi_178 - xi_179;
    const float xi_188 = xi_146 * xi_175;
    const float xi_189 = xi_148 * xi_173;
    const float xi_190 = -xi_188 + xi_189;
    const float xi_191 = xi_188 - xi_189;
    const float xi_192 = xi_174 - xi_176;
    const float rr_0 = 0.0f;
    const float xi_23 = rr_0 * xi_22;
    const float xi_35 = rr_0 * xi_34;
    const float xi_41 = rr_0 * xi_40;
    const float xi_45 = rr_0 * 0.041666666666666664f;
    const float xi_46 = xi_45 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_50 = xi_45 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float xi_72 = xi_45 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float vel0Term = xi_3 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
    const float vel1Term = xi_4 + xi_5;
    const float vel2Term = xi_6 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3];
    const float delta_rho = vel0Term + vel1Term + vel2Term + xi_8 + xi_9 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2];
    const float rho = delta_rho + 1.0f;
    const float xi_0 = ((1.0f) / (rho));
    const float xi_10 = xi_0 * 0.5f;
    const float u_0 = xi_0 * (vel0Term - xi_11 - xi_8) + xi_10 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_17 = u_0 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_28 = xi_17 * 0.16666666666666666f;
    const float xi_29 = -xi_28;
    const float xi_30 = xi_17 * 0.083333333333333329f;
    const float xi_31 = omega_shear * xi_30 + xi_29;
    const float xi_48 = xi_17 * xi_47 + xi_29;
    const float xi_49 = xi_34 - xi_46 + xi_48;
    const float xi_52 = xi_17 * xi_51;
    const float xi_59 = u_0 * xi_58;
    const float xi_64 = u_0 * xi_63;
    const float xi_68 = -xi_34 + xi_46 + xi_48;
    const float xi_75 = omega_shear * u_0 * -0.083333333333333329f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_85 = u_0 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float xi_86 = xi_85 * 0.25f;
    const float xi_89 = xi_62 * xi_85;
    const float xi_96 = u_0 * u_0;
    const float u_1 = xi_0 * (vel1Term - xi_12 - xi_13 - xi_9) + xi_10 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float xi_18 = u_1 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float xi_26 = xi_18 * 0.16666666666666666f;
    const float xi_36 = omega_shear * u_1 * -0.083333333333333329f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const float xi_42 = -xi_26;
    const float xi_43 = xi_18 * 0.083333333333333329f;
    const float xi_53 = xi_18 * xi_51;
    const float xi_60 = u_1 * 0.25f;
    const float xi_61 = xi_60 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_65 = u_1 * xi_62;
    const float xi_66 = xi_65 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_67 = xi_59 + xi_61 - xi_64 - xi_66;
    const float xi_69 = -xi_59 - xi_61 + xi_64 + xi_66;
    const float xi_77 = xi_60 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float xi_79 = xi_65 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float xi_95 = rho * (u_1 * u_1);
    const float xi_101 = -xi_95;
    const float xi_157 = rho * u_1;
    const float xi_159 = xi_158 * (u_0 * xi_157 + xi_12 - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3]);
    const float xi_160 = -xi_159;
    const float u_2 = xi_0 * (vel2Term - xi_14 - xi_15 - xi_16 - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3]) + xi_10 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float xi_19 = u_2 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float xi_24 = xi_19 * 0.16666666666666666f;
    const float xi_25 = -xi_24;
    const float xi_27 = xi_19 * 0.083333333333333329f;
    const float xi_32 = -omega_shear * xi_26 + omega_shear * xi_27 + xi_18 * 0.33333333333333331f + xi_25 + xi_31;
    const float xi_37 = omega_shear * u_2 * -0.083333333333333329f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const float xi_38 = omega_shear * xi_28 + u_0 * -0.33333333333333331f * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2] + xi_24 + xi_26 + xi_36 + xi_37;
    const float xi_44 = -omega_shear * xi_24 + omega_shear * xi_43 + xi_19 * 0.33333333333333331f + xi_31 + xi_42;
    const float xi_54 = xi_19 * xi_51;
    const float xi_55 = xi_18 * xi_47 + xi_42 + xi_52 + xi_53 + xi_54;
    const float xi_56 = -xi_22 + xi_50 + xi_55;
    const float xi_57 = xi_27 + xi_37 + xi_56;
    const float xi_70 = xi_22 - xi_50 + xi_55;
    const float xi_71 = xi_27 + xi_37 + xi_70;
    const float xi_73 = xi_19 * xi_47 + xi_25;
    const float xi_74 = -xi_40 + xi_72 + xi_73;
    const float xi_76 = xi_30 + xi_56 + xi_75;
    const float xi_78 = u_2 * xi_58;
    const float xi_80 = u_2 * xi_63;
    const float xi_81 = -xi_77 - xi_78 + xi_79 + xi_80;
    const float xi_82 = xi_30 + xi_70 + xi_75;
    const float xi_83 = xi_77 + xi_78 - xi_79 - xi_80;
    const float xi_84 = xi_36 + xi_43 + xi_52 + xi_53 + xi_54 + xi_74;
    const float xi_87 = u_2 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const float xi_88 = xi_87 * 0.25f;
    const float xi_90 = xi_62 * xi_87;
    const float xi_91 = xi_86 + xi_88 - xi_89 - xi_90;
    const float xi_92 = -xi_86 - xi_88 + xi_89 + xi_90;
    const float xi_93 = xi_40 - xi_72 + xi_73;
    const float xi_94 = xi_36 + xi_43 + xi_52 + xi_53 + xi_54 + xi_93;
    const float xi_98 = rho * (u_2 * u_2);
    const float xi_99 = xi_97 + xi_98 * 0.66666666666666663f + 3.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] + 3.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3];
    const float xi_100 = omega_even * (rho * xi_96 * 1.6666666666666667f + xi_95 * 0.66666666666666663f + xi_99 - 3.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] - 3.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] - 3.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] - 3.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3] + 3.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] + 3.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3]);
    const float xi_103 = omega_bulk * (rho * xi_96 - xi_101 - xi_102 - xi_13 - xi_16 - xi_5 - xi_97 + xi_98);
    const float xi_105 = xi_104 + xi_95 * 2.3333333333333335f + xi_99 - 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] - 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3] - 5.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] - 5.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3] - 5.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3] - 5.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const float xi_106 = omega_even * xi_105;
    const float xi_110 = xi_104 + xi_107 + xi_108 + xi_109 + xi_97 + xi_98 * 3.0f - 4.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] - 4.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3] - 7.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] - 7.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] - 7.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] - 7.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3] + 5.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] + 5.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const float xi_111 = omega_even * xi_110;
    const float xi_112 = xi_111 * 0.01984126984126984f;
    const float xi_118 = xi_116 + xi_117 + xi_98 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    const float xi_119 = omega_shear * (-xi_101 - xi_114 - xi_115 - xi_118 - xi_15 - xi_2 - xi_4 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3]);
    const float xi_120 = xi_119 * 0.125f;
    const float xi_121 = -xi_120;
    const float xi_135 = xi_100 * 0.050000000000000003f;
    const float xi_136 = rho * xi_96 * 2.0f - xi_102 - xi_113 - xi_118 - xi_124 - xi_95 - 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] - 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] + 2.0f * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3];
    const float xi_137 = omega_shear * xi_136;
    const float xi_138 = xi_137 * 0.041666666666666664f;
    const float xi_139 = xi_135 + xi_138;
    const float xi_140 = -xi_112;
    const float xi_141 = -xi_138;
    const float xi_142 = -xi_135 + xi_141;
    const float xi_143 = xi_106 * 0.035714285714285712f;
    const float xi_151 = xi_106 * 0.021428571428571429f;
    const float xi_162 = xi_119 * 0.0625f;
    const float xi_163 = xi_111 * 0.013888888888888888f;
    const float xi_164 = xi_103 * 0.041666666666666664f;
    const float xi_165 = xi_137 * 0.020833333333333332f + xi_164;
    const float xi_166 = xi_133 + xi_162 + xi_163 + xi_165;
    const float xi_167 = -xi_133 + xi_162 + xi_163 + xi_165;
    const float xi_168 = xi_111 * -0.003968253968253968f;
    const float xi_169 = xi_106 * -0.0071428571428571426f;
    const float xi_170 = xi_158 * (u_2 * xi_157 + xi_123 + xi_128 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3]);
    const float xi_171 = xi_100 * 0.025000000000000001f;
    const float xi_172 = xi_141 + xi_164 + xi_168 + xi_169 + xi_170 + xi_171;
    const float xi_182 = xi_141 + xi_164 + xi_168 + xi_169 - xi_170 + xi_171;
    const float xi_183 = xi_158 * (rho * u_0 * u_2 + xi_114 + xi_144 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3]);
    const float xi_184 = -xi_183;
    const float xi_185 = -xi_162;
    const float xi_186 = xi_106 * 0.017857142857142856f;
    const float xi_187 = xi_155 + xi_165 + xi_168 + xi_185 + xi_186;
    const float xi_193 = -xi_155 + xi_165 + xi_168 + xi_185 + xi_186;
    const float forceTerm_0 = xi_17 * xi_20 - xi_17 + xi_18 * xi_20 - xi_18 + xi_19 * xi_20 - xi_19;
    const float forceTerm_1 = xi_21 - xi_23 + xi_32;
    const float forceTerm_2 = -xi_21 + xi_23 + xi_32;
    const float forceTerm_3 = -xi_33 + xi_35 - xi_38;
    const float forceTerm_4 = xi_33 - xi_35 - xi_38;
    const float forceTerm_5 = xi_39 - xi_41 + xi_44;
    const float forceTerm_6 = -xi_39 + xi_41 + xi_44;
    const float forceTerm_7 = -xi_49 - xi_57 - xi_67;
    const float forceTerm_8 = -xi_57 - xi_68 - xi_69;
    const float forceTerm_9 = -xi_49 - xi_69 - xi_71;
    const float forceTerm_10 = -xi_67 - xi_68 - xi_71;
    const float forceTerm_11 = -xi_74 - xi_76 - xi_81;
    const float forceTerm_12 = -xi_74 - xi_82 - xi_83;
    const float forceTerm_13 = -xi_49 - xi_84 - xi_91;
    const float forceTerm_14 = -xi_68 - xi_84 - xi_92;
    const float forceTerm_15 = -xi_76 - xi_83 - xi_93;
    const float forceTerm_16 = -xi_81 - xi_82 - xi_93;
    const float forceTerm_17 = -xi_49 - xi_92 - xi_94;
    const float forceTerm_18 = -xi_68 - xi_91 - xi_94;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2] = forceTerm_0 + xi_100 * 0.10000000000000001f + xi_103 * -0.5f + xi_106 * 0.042857142857142858f + xi_111 * 0.023809523809523808f + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3] = forceTerm_1 + omega_even * xi_105 * 0.014285714285714285f - xi_112 - xi_113 - xi_121 - xi_134 - xi_139 + ((ctr_1 >= lebc_top_index) ? (rho * v_s * (u_0 * 2.0f + v_s) * 0.16666666666666666f) : (0.0f));
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3] = forceTerm_2 + xi_106 * 0.014285714285714285f + xi_120 + xi_134 + xi_140 + xi_142 + ((1 >= ctr_1 - lebc_bot_index) ? (rho * v_s * (u_0 * -2.0f + v_s) * 0.16666666666666666f) : (0.0f)) + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3] = forceTerm_3 + xi_137 * 0.083333333333333329f + xi_140 - xi_143 + xi_150 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3] = forceTerm_4 + omega_shear * xi_136 * 0.083333333333333329f - xi_112 - xi_143 - xi_147 - xi_150;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3] = forceTerm_5 + omega_even * xi_110 * 0.015873015873015872f - xi_116 - xi_120 - xi_139 - xi_151 - xi_156;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3] = forceTerm_6 + xi_111 * 0.015873015873015872f + xi_121 + xi_142 - xi_151 + xi_156 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3] = forceTerm_7 + xi_160 + xi_161 + xi_166 + ((ctr_1 >= lebc_top_index) ? (rho * v_s * (u_0 * -2.0f + u_1 * 3.0f - v_s + 1.0f) * 0.083333333333333329f) : (0.0f)) + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3] = forceTerm_8 + xi_149 + xi_159 + xi_166 + ((ctr_1 >= lebc_top_index) ? (rho * v_s * (u_0 * -2.0f + u_1 * -3.0f - v_s - 1.0f) * 0.083333333333333329f) : (0.0f)) + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3] = forceTerm_9 + xi_159 + xi_161 + xi_167 + ((1 >= ctr_1 - lebc_bot_index) ? (rho * v_s * (u_0 * 2.0f + u_1 * 3.0f - v_s - 1.0f) * 0.083333333333333329f) : (0.0f)) + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3] = forceTerm_10 + xi_149 + xi_160 + xi_167 + ((1 >= ctr_1 - lebc_bot_index) ? (rho * v_s * (u_0 * 2.0f + u_1 * -3.0f - v_s + 1.0f) * 0.083333333333333329f) : (0.0f)) + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3] = forceTerm_11 + xi_172 + xi_177 + xi_180 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3] = forceTerm_12 + xi_177 + xi_181 + xi_182 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3] = forceTerm_13 + xi_184 + xi_187 + xi_190 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3] = forceTerm_14 + xi_183 + xi_187 + xi_191 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3] = forceTerm_15 + xi_180 + xi_182 + xi_192 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3] = forceTerm_16 + xi_172 + xi_181 + xi_192 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3] = forceTerm_17 + xi_183 + xi_190 + xi_193 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3] = forceTerm_18 + xi_184 + xi_191 + xi_193 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3];
  }
}
} // namespace internal_streamcollidesweepleesedwardssingleprecisioncuda_streamcollidesweepleesedwardssingleprecisioncuda

void StreamCollideSweepLeesEdwardsSinglePrecisionCUDA::run(IBlock *block, gpuStream_t stream) {

  auto force = block->getData<gpu::GPUField<float>>(forceID);
  auto pdfs = block->getData<gpu::GPUField<float>>(pdfsID);
  gpu::GPUField<float> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  auto &lebc_top_index = this->lebc_top_index_;
  auto &lebc_bot_index = this->lebc_bot_index_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
  auto &v_s = this->v_s_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_bulk = this->omega_bulk_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force->nrOfGhostLayers()))
  float *RESTRICT const _data_force = force->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()))
  float *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  float *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(force->xSize()) + 2))
  const int64_t _size_force_0 = int64_t(int64_c(force->xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(force->ySize()) + 2))
  const int64_t _size_force_1 = int64_t(int64_c(force->ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(force->zSize()) + 2))
  const int64_t _size_force_2 = int64_t(int64_c(force->zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
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
  dim3 _block(uint32_c(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)), uint32_c(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))), uint32_c(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))));
  dim3 _grid(uint32_c(((_size_force_0 - 2) % (((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)) == 0 ? (int64_t)(_size_force_0 - 2) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)) : ((int64_t)(_size_force_0 - 2) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))) + 1)), uint32_c(((_size_force_1 - 2) % (((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))) == 0 ? (int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))) : ((int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) + 1)), uint32_c(((_size_force_2 - 2) % (((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))) == 0 ? (int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))) : ((int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))))) + 1)));
  internal_streamcollidesweepleesedwardssingleprecisioncuda_streamcollidesweepleesedwardssingleprecisioncuda::streamcollidesweepleesedwardssingleprecisioncuda_streamcollidesweepleesedwardssingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, lebc_bot_index, lebc_top_index, omega_bulk, omega_even, omega_odd, omega_shear, v_s);
  pdfs->swapDataPointers(pdfs_tmp);
}

void StreamCollideSweepLeesEdwardsSinglePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block, gpuStream_t stream) {

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto force = block->getData<gpu::GPUField<float>>(forceID);
  auto pdfs = block->getData<gpu::GPUField<float>>(pdfsID);
  gpu::GPUField<float> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  auto &lebc_top_index = this->lebc_top_index_;
  auto &lebc_bot_index = this->lebc_bot_index_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
  auto &v_s = this->v_s_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_bulk = this->omega_bulk_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force->nrOfGhostLayers()))
  float *RESTRICT const _data_force = force->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  float *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
  const int64_t _size_force_0 = int64_t(int64_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
  const int64_t _size_force_1 = int64_t(int64_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
  const int64_t _size_force_2 = int64_t(int64_c(ci.zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
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
  dim3 _block(uint32_c(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)), uint32_c(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))), uint32_c(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))));
  dim3 _grid(uint32_c(((_size_force_0 - 2) % (((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)) == 0 ? (int64_t)(_size_force_0 - 2) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)) : ((int64_t)(_size_force_0 - 2) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))) + 1)), uint32_c(((_size_force_1 - 2) % (((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))) == 0 ? (int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))) : ((int64_t)(_size_force_1 - 2) / (int64_t)(((1024 < ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))) ? 1024 : ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) + 1)), uint32_c(((_size_force_2 - 2) % (((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))) == 0 ? (int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))))) : ((int64_t)(_size_force_2 - 2) / (int64_t)(((64 < ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))) ? 64 : ((_size_force_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2))))))) ? _size_force_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2) * ((_size_force_1 - 2 < 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))) ? _size_force_1 - 2 : 2 * ((int64_t)(128) / (int64_t)(((128 < _size_force_0 - 2) ? 128 : _size_force_0 - 2)))))))))) + 1)));
  internal_streamcollidesweepleesedwardssingleprecisioncuda_streamcollidesweepleesedwardssingleprecisioncuda::streamcollidesweepleesedwardssingleprecisioncuda_streamcollidesweepleesedwardssingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, lebc_bot_index, lebc_top_index, omega_bulk, omega_even, omega_odd, omega_shear, v_s);
  pdfs->swapDataPointers(pdfs_tmp);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
