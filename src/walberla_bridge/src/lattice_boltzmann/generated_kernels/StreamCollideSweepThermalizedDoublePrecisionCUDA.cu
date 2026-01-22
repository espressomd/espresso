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
//! \\file StreamCollideSweepThermalizedDoublePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 272d4a09ec35da50685afc9586645e1b9984b423

#include <cmath>

#include "StreamCollideSweepThermalizedDoublePrecisionCUDA.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include "philox_rand.h"

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

namespace internal_streamcollidesweepthermalizeddoubleprecisioncuda_streamcollidesweepthermalizeddoubleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void streamcollidesweepthermalizeddoubleprecisioncuda_streamcollidesweepthermalizeddoubleprecisioncuda(double *RESTRICT const _data_force, double *RESTRICT const _data_pdfs, double *RESTRICT _data_pdfs_tmp, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double kT, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t seed, uint32_t time_step) {
  if (blockDim.x * blockIdx.x + threadIdx.x + 1 < _size_force_0 - 1 && blockDim.y * blockIdx.y + threadIdx.y + 1 < _size_force_1 - 1 && blockDim.z * blockIdx.z + threadIdx.z + 1 < _size_force_2 - 1) {
    const int64_t ctr_0 = blockDim.x * blockIdx.x + threadIdx.x + 1;
    const int64_t ctr_1 = blockDim.y * blockIdx.y + threadIdx.y + 1;
    const int64_t ctr_2 = blockDim.z * blockIdx.z + threadIdx.z + 1;

    double random_7_0{};
    double random_7_1{};
    if (kT > 0.) {
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 7, seed, random_7_0, random_7_1);
    }

    double random_6_0{};
    double random_6_1{};
    if (kT > 0.) {
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 6, seed, random_6_0, random_6_1);
    }

    double random_5_0{};
    double random_5_1{};
    if (kT > 0.) {
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 5, seed, random_5_0, random_5_1);
    }

    double random_4_0{};
    double random_4_1{};
    if (kT > 0.) {
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 4, seed, random_4_0, random_4_1);
    }

    double random_3_0{};
    double random_3_1{};
    if (kT > 0.) {
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, random_3_0, random_3_1);
    }

    double random_2_0{};
    double random_2_1{};
    if (kT > 0.) {
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, random_2_0, random_2_1);
    }

    double random_1_0{};
    double random_1_1{};
    if (kT > 0.) {
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, random_1_0, random_1_1);
    }

    double random_0_0{};
    double random_0_1{};
    if (kT > 0.) {
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, random_0_0, random_0_1);
    }
    const double xi_2 = _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    const double xi_3 = xi_2 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const double xi_4 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const double xi_5 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    const double xi_6 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3];
    const double xi_7 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3];
    const double xi_8 = xi_7 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    const double xi_9 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    const double xi_11 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3];
    const double xi_12 = _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    const double xi_13 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3];
    const double xi_14 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    const double xi_15 = -_data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const double xi_16 = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    const double xi_20 = omega_bulk * 0.5;
    const double xi_21 = 0.16666666666666666 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const double xi_22 = 0.083333333333333329 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const double xi_33 = 0.16666666666666666 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_34 = 0.083333333333333329 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_39 = 0.16666666666666666 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double xi_40 = 0.083333333333333329 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double xi_47 = omega_shear * 0.041666666666666664;
    const double xi_51 = omega_bulk * 0.041666666666666664;
    const double xi_58 = 0.25 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const double xi_62 = omega_shear * 0.125;
    const double xi_63 = xi_62 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const double xi_97 = 3.7416573867739413;
    const double xi_98 = random_6_0 - 0.5;
    const double xi_100 = 5.4772255750516612;
    const double xi_101 = random_7_0 - 0.5;
    const double xi_103 = random_2_1 - 0.5;
    const double xi_104 = 2.4494897427831779;
    const double xi_107 = 8.3666002653407556;
    const double xi_108 = random_6_1 - 0.5;
    const double xi_112 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2];
    const double xi_118 = xi_11 + xi_3;
    const double xi_121 = 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3] + 5.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] + 5.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
    const double xi_124 = 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3];
    const double xi_125 = 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const double xi_126 = 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3];
    const double xi_129 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const double xi_133 = random_0_1 - 0.5;
    const double xi_137 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3];
    const double xi_138 = -_data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3];
    const double xi_139 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3];
    const double xi_140 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3];
    const double xi_145 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3];
    const double xi_146 = xi_14 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3];
    const double xi_147 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3];
    const double xi_148 = xi_147 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const double xi_149 = xi_145 + xi_146 + xi_148 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    const double xi_150 = omega_odd * 0.25;
    const double xi_151 = random_5_1 - 0.5;
    const double xi_155 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    const double xi_156 = -_data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3];
    const double xi_157 = 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    const double xi_158 = -2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3];
    const double xi_159 = xi_148 + xi_155 + xi_156 - xi_157 + xi_158 + xi_5 + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    const double xi_160 = omega_odd * 0.083333333333333329;
    const double xi_161 = xi_159 * xi_160;
    const double xi_162 = random_3_0 - 0.5;
    const double xi_173 = 1.7320508075688772;
    const double xi_174 = random_0_0 - 0.5;
    const double xi_185 = xi_15 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3];
    const double xi_186 = xi_138 + xi_185;
    const double xi_187 = xi_186 - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
    const double xi_188 = random_4_1 - 0.5;
    const double xi_189 = -_data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
    const double xi_190 = -xi_157 - xi_158 - xi_186 - xi_189 - xi_7 + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    const double xi_191 = xi_160 * xi_190;
    const double xi_192 = random_4_0 - 0.5;
    const double xi_198 = xi_145 + xi_155 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    const double xi_199 = -xi_139 - xi_198 - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3];
    const double xi_200 = random_5_0 - 0.5;
    const double xi_201 = -xi_124 - xi_125 + xi_126 + xi_140 + xi_198 + xi_6;
    const double xi_202 = xi_160 * xi_201;
    const double xi_203 = random_3_1 - 0.5;
    const double xi_217 = omega_shear * 0.25;
    const double xi_223 = omega_odd * 0.041666666666666664;
    const double xi_224 = xi_159 * xi_223;
    const double xi_226 = omega_odd * 0.125;
    const double xi_227 = xi_149 * xi_226;
    const double xi_233 = xi_201 * xi_223;
    const double xi_234 = xi_199 * xi_226;
    const double xi_254 = xi_187 * xi_226;
    const double xi_255 = xi_190 * xi_223;
    const double rr_0 = 0.0;
    const double xi_23 = rr_0 * xi_22;
    const double xi_35 = rr_0 * xi_34;
    const double xi_41 = rr_0 * xi_40;
    const double xi_45 = rr_0 * 0.041666666666666664;
    const double xi_46 = xi_45 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_50 = xi_45 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const double xi_72 = xi_45 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double vel0Term = xi_3 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
    const double vel1Term = xi_4 + xi_5;
    const double vel2Term = xi_6 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3];
    const double delta_rho = vel0Term + vel1Term + vel2Term + xi_8 + xi_9 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2];
    const double rho = delta_rho + 1.0;
    const double xi_95 = kT * rho;
    const double xi_96 = pow(xi_95 * (1.0 - ((-omega_even + 1.0) * (-omega_even + 1.0))), 0.5);
    const double xi_99 = xi_96 * xi_97 * xi_98;
    const double xi_102 = xi_100 * xi_101 * xi_96;
    const double xi_105 = pow(xi_95 * (1.0 - ((-omega_bulk + 1.0) * (-omega_bulk + 1.0))), 0.5);
    const double xi_106 = xi_103 * xi_104 * xi_105;
    const double xi_109 = xi_107 * xi_108 * xi_96;
    const double xi_131 = xi_99 * 0.11904761904761904;
    const double xi_134 = pow(xi_95 * (1.0 - ((-omega_shear + 1.0) * (-omega_shear + 1.0))), 0.5);
    const double xi_135 = xi_134 * 0.5;
    const double xi_136 = xi_133 * xi_135;
    const double xi_152 = pow(xi_95 * (1.0 - ((-omega_odd + 1.0) * (-omega_odd + 1.0))), 0.5);
    const double xi_153 = xi_152 * 1.4142135623730951;
    const double xi_154 = xi_153 * 0.5;
    const double xi_163 = xi_104 * xi_152;
    const double xi_164 = xi_163 * 0.16666666666666666;
    const double xi_165 = xi_162 * xi_164;
    const double xi_166 = xi_161 + xi_165;
    const double xi_167 = xi_149 * xi_150 + xi_151 * xi_154 + xi_166;
    const double xi_169 = xi_102 * 0.10000000000000001;
    const double xi_175 = xi_134 * xi_173 * xi_174;
    const double xi_176 = xi_175 * 0.16666666666666666;
    const double xi_184 = xi_109 * 0.071428571428571425;
    const double xi_193 = xi_164 * xi_192;
    const double xi_194 = xi_191 + xi_193;
    const double xi_195 = xi_150 * xi_187 + xi_154 * xi_188 + xi_194;
    const double xi_197 = xi_109 * 0.042857142857142858;
    const double xi_204 = xi_164 * xi_203;
    const double xi_205 = xi_202 + xi_204;
    const double xi_206 = xi_150 * xi_199 + xi_154 * xi_200 + xi_205;
    const double xi_207 = xi_133 * xi_134 * 0.25;
    const double xi_210 = xi_99 * 0.083333333333333329;
    const double xi_214 = -xi_191 - xi_193;
    const double xi_215 = xi_135 * (random_1_0 - 0.5);
    const double xi_222 = xi_135 * (random_2_0 - 0.5);
    const double xi_228 = xi_163 * 0.083333333333333329;
    const double xi_229 = xi_162 * xi_228;
    const double xi_230 = xi_153 * 0.25;
    const double xi_231 = xi_151 * xi_230;
    const double xi_235 = xi_203 * xi_228;
    const double xi_236 = xi_200 * xi_230;
    const double xi_237 = -xi_233 + xi_234 - xi_235 + xi_236;
    const double xi_239 = xi_109 * 0.014285714285714285;
    const double xi_241 = xi_99 * 0.023809523809523808;
    const double xi_244 = xi_233 - xi_234 + xi_235 - xi_236;
    const double xi_246 = -xi_207;
    const double xi_249 = xi_109 * 0.035714285714285712;
    const double xi_251 = xi_135 * (random_1_1 - 0.5);
    const double xi_256 = xi_188 * xi_230;
    const double xi_257 = xi_192 * xi_228;
    const double xi_258 = -xi_254 + xi_255 - xi_256 + xi_257;
    const double xi_260 = xi_254 - xi_255 + xi_256 - xi_257;
    const double xi_0 = ((1.0) / (rho));
    const double xi_10 = xi_0 * 0.5;
    const double u_0 = xi_0 * (vel0Term - xi_11 - xi_8) + xi_10 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_17 = u_0 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_28 = xi_17 * 0.16666666666666666;
    const double xi_29 = -xi_28;
    const double xi_30 = xi_17 * 0.083333333333333329;
    const double xi_31 = omega_shear * xi_30 + xi_29;
    const double xi_48 = xi_17 * xi_47 + xi_29;
    const double xi_49 = xi_34 - xi_46 + xi_48;
    const double xi_52 = xi_17 * xi_51;
    const double xi_59 = u_0 * xi_58;
    const double xi_64 = u_0 * xi_63;
    const double xi_68 = -xi_34 + xi_46 + xi_48;
    const double xi_75 = omega_shear * u_0 * -0.083333333333333329 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_85 = u_0 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double xi_86 = xi_85 * 0.25;
    const double xi_89 = xi_62 * xi_85;
    const double xi_111 = u_0 * u_0;
    const double u_1 = xi_0 * (vel1Term - xi_12 - xi_13 - xi_9) + xi_10 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const double xi_18 = u_1 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const double xi_26 = xi_18 * 0.16666666666666666;
    const double xi_36 = omega_shear * u_1 * -0.083333333333333329 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3];
    const double xi_42 = -xi_26;
    const double xi_43 = xi_18 * 0.083333333333333329;
    const double xi_53 = xi_18 * xi_51;
    const double xi_60 = u_1 * 0.25;
    const double xi_61 = xi_60 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_65 = u_1 * xi_62;
    const double xi_66 = xi_65 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_67 = xi_59 + xi_61 - xi_64 - xi_66;
    const double xi_69 = -xi_59 - xi_61 + xi_64 + xi_66;
    const double xi_77 = xi_60 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double xi_79 = xi_65 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double xi_110 = rho * (u_1 * u_1);
    const double xi_117 = -xi_110;
    const double xi_216 = rho * u_1;
    const double xi_218 = xi_217 * (u_0 * xi_216 + xi_12 - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3]);
    const double xi_219 = -xi_215 - xi_218;
    const double xi_220 = xi_215 + xi_218;
    const double u_2 = xi_0 * (vel2Term - xi_14 - xi_15 - xi_16 - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3]) + xi_10 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double xi_19 = u_2 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double xi_24 = xi_19 * 0.16666666666666666;
    const double xi_25 = -xi_24;
    const double xi_27 = xi_19 * 0.083333333333333329;
    const double xi_32 = -omega_shear * xi_26 + omega_shear * xi_27 + xi_18 * 0.33333333333333331 + xi_25 + xi_31;
    const double xi_37 = omega_shear * u_2 * -0.083333333333333329 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3];
    const double xi_38 = omega_shear * xi_28 + u_0 * -0.33333333333333331 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2] + xi_24 + xi_26 + xi_36 + xi_37;
    const double xi_44 = -omega_shear * xi_24 + omega_shear * xi_43 + xi_19 * 0.33333333333333331 + xi_31 + xi_42;
    const double xi_54 = xi_19 * xi_51;
    const double xi_55 = xi_18 * xi_47 + xi_42 + xi_52 + xi_53 + xi_54;
    const double xi_56 = -xi_22 + xi_50 + xi_55;
    const double xi_57 = xi_27 + xi_37 + xi_56;
    const double xi_70 = xi_22 - xi_50 + xi_55;
    const double xi_71 = xi_27 + xi_37 + xi_70;
    const double xi_73 = xi_19 * xi_47 + xi_25;
    const double xi_74 = -xi_40 + xi_72 + xi_73;
    const double xi_76 = xi_30 + xi_56 + xi_75;
    const double xi_78 = u_2 * xi_58;
    const double xi_80 = u_2 * xi_63;
    const double xi_81 = -xi_77 - xi_78 + xi_79 + xi_80;
    const double xi_82 = xi_30 + xi_70 + xi_75;
    const double xi_83 = xi_77 + xi_78 - xi_79 - xi_80;
    const double xi_84 = xi_36 + xi_43 + xi_52 + xi_53 + xi_54 + xi_74;
    const double xi_87 = u_2 * _data_force[_stride_force_0 * ctr_0 + _stride_force_1 * ctr_1 + _stride_force_2 * ctr_2];
    const double xi_88 = xi_87 * 0.25;
    const double xi_90 = xi_62 * xi_87;
    const double xi_91 = xi_86 + xi_88 - xi_89 - xi_90;
    const double xi_92 = -xi_86 - xi_88 + xi_89 + xi_90;
    const double xi_93 = xi_40 - xi_72 + xi_73;
    const double xi_94 = xi_36 + xi_43 + xi_52 + xi_53 + xi_54 + xi_93;
    const double xi_113 = rho * (u_2 * u_2);
    const double xi_114 = xi_112 + xi_113 * 0.66666666666666663 + 3.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] + 3.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3];
    const double xi_115 = rho * xi_111 * 1.6666666666666667 + xi_110 * 0.66666666666666663 + xi_114 - 3.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] - 3.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] - 3.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] - 3.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3] + 3.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] + 3.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const double xi_116 = omega_even * xi_115;
    const double xi_119 = rho * xi_111 - xi_112 + xi_113 - xi_117 - xi_118 - xi_13 - xi_16 - xi_5;
    const double xi_120 = omega_bulk * xi_119;
    const double xi_122 = xi_110 * 2.3333333333333335 + xi_114 + xi_121 - 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] - 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3] - 5.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3] - 5.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3] - 5.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3] - 5.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    const double xi_123 = omega_even * xi_122;
    const double xi_127 = xi_112 + xi_113 * 3.0 + xi_121 + xi_124 + xi_125 + xi_126 - 4.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3] - 4.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3] - 7.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] - 7.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] - 7.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] - 7.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3] + 5.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] + 5.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
    const double xi_128 = omega_even * xi_127;
    const double xi_130 = xi_128 * 0.01984126984126984;
    const double xi_132 = xi_130 + xi_131;
    const double xi_141 = xi_113 + xi_139 + xi_140 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    const double xi_142 = omega_shear * (-xi_117 - xi_137 - xi_138 - xi_141 - xi_15 - xi_2 - xi_4 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3]);
    const double xi_143 = xi_142 * 0.125;
    const double xi_144 = -xi_136 - xi_143;
    const double xi_168 = xi_116 * 0.050000000000000003;
    const double xi_170 = rho * xi_111 * 2.0 - xi_110 - xi_118 - xi_129 - xi_141 - xi_147 - 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] - 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3] + 2.0 * _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3] - _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3];
    const double xi_171 = omega_shear * xi_170;
    const double xi_172 = xi_171 * 0.041666666666666664;
    const double xi_177 = xi_172 + xi_176;
    const double xi_178 = xi_168 + xi_169 + xi_177;
    const double xi_179 = -xi_130 - xi_131;
    const double xi_180 = xi_136 + xi_143;
    const double xi_181 = -xi_172 - xi_176;
    const double xi_182 = -xi_168 - xi_169 + xi_181;
    const double xi_183 = xi_123 * 0.035714285714285712;
    const double xi_196 = xi_123 * 0.021428571428571429;
    const double xi_208 = xi_142 * 0.0625;
    const double xi_209 = xi_128 * 0.013888888888888888;
    const double xi_211 = xi_106 * 0.083333333333333329 + xi_120 * 0.041666666666666664;
    const double xi_212 = xi_171 * 0.020833333333333332 + xi_175 * 0.083333333333333329 + xi_211;
    const double xi_213 = xi_166 + xi_207 + xi_208 + xi_209 + xi_210 + xi_212;
    const double xi_221 = -xi_161 - xi_165 + xi_207 + xi_208 + xi_209 + xi_210 + xi_212;
    const double xi_225 = xi_217 * (u_2 * xi_216 + xi_146 + xi_155 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3]);
    const double xi_232 = xi_222 - xi_224 + xi_225 + xi_227 - xi_229 + xi_231;
    const double xi_238 = xi_123 * 0.0071428571428571426;
    const double xi_240 = xi_128 * 0.003968253968253968;
    const double xi_242 = -xi_240 - xi_241;
    const double xi_243 = xi_102 * 0.050000000000000003 + xi_116 * 0.025000000000000001 + xi_181 + xi_211 - xi_238 - xi_239 + xi_242;
    const double xi_245 = omega_bulk * xi_119 * -0.041666666666666664 + omega_even * xi_115 * -0.025000000000000001 + xi_100 * xi_101 * xi_96 * -0.050000000000000003 + xi_103 * xi_104 * xi_105 * -0.083333333333333329 + xi_177 + xi_238 + xi_239 + xi_240 + xi_241;
    const double xi_247 = -xi_208;
    const double xi_248 = xi_123 * 0.017857142857142856;
    const double xi_250 = xi_205 + xi_212 + xi_242 + xi_246 + xi_247 + xi_248 + xi_249;
    const double xi_252 = xi_217 * (rho * u_0 * u_2 + xi_137 + xi_185 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3]);
    const double xi_253 = -xi_251 - xi_252;
    const double xi_259 = xi_251 + xi_252;
    const double xi_261 = xi_222 + xi_224 + xi_225 - xi_227 + xi_229 - xi_231;
    const double xi_262 = -xi_202 - xi_204 + xi_212 + xi_242 + xi_246 + xi_247 + xi_248 + xi_249;
    const double forceTerm_0 = xi_17 * xi_20 - xi_17 + xi_18 * xi_20 - xi_18 + xi_19 * xi_20 - xi_19;
    const double forceTerm_1 = xi_21 - xi_23 + xi_32;
    const double forceTerm_2 = -xi_21 + xi_23 + xi_32;
    const double forceTerm_3 = -xi_33 + xi_35 - xi_38;
    const double forceTerm_4 = xi_33 - xi_35 - xi_38;
    const double forceTerm_5 = xi_39 - xi_41 + xi_44;
    const double forceTerm_6 = -xi_39 + xi_41 + xi_44;
    const double forceTerm_7 = -xi_49 - xi_57 - xi_67;
    const double forceTerm_8 = -xi_57 - xi_68 - xi_69;
    const double forceTerm_9 = -xi_49 - xi_69 - xi_71;
    const double forceTerm_10 = -xi_67 - xi_68 - xi_71;
    const double forceTerm_11 = -xi_74 - xi_76 - xi_81;
    const double forceTerm_12 = -xi_74 - xi_82 - xi_83;
    const double forceTerm_13 = -xi_49 - xi_84 - xi_91;
    const double forceTerm_14 = -xi_68 - xi_84 - xi_92;
    const double forceTerm_15 = -xi_76 - xi_83 - xi_93;
    const double forceTerm_16 = -xi_81 - xi_82 - xi_93;
    const double forceTerm_17 = -xi_49 - xi_92 - xi_94;
    const double forceTerm_18 = -xi_68 - xi_91 - xi_94;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2] = forceTerm_0 + xi_102 * 0.20000000000000001 - xi_106 + xi_109 * 0.085714285714285715 + xi_116 * 0.10000000000000001 + xi_120 * -0.5 + xi_123 * 0.042857142857142858 + xi_128 * 0.023809523809523808 + xi_99 * 0.14285714285714285 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3] = forceTerm_1 + omega_even * xi_122 * 0.014285714285714285 + xi_107 * xi_108 * xi_96 * 0.028571428571428571 - xi_129 - xi_132 - xi_144 - xi_167 - xi_178;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3] = forceTerm_2 + xi_109 * 0.028571428571428571 + xi_123 * 0.014285714285714285 + xi_167 + xi_179 + xi_180 + xi_182 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3] = forceTerm_3 + xi_171 * 0.083333333333333329 + xi_175 * 0.33333333333333331 + xi_179 - xi_183 - xi_184 + xi_195 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3] = forceTerm_4 + omega_shear * xi_170 * 0.083333333333333329 - xi_132 + xi_134 * xi_173 * xi_174 * 0.33333333333333331 - xi_183 - xi_184 - xi_189 - xi_195;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3] = forceTerm_5 + omega_even * xi_127 * 0.015873015873015872 - xi_139 - xi_178 - xi_180 - xi_196 - xi_197 - xi_206 + xi_96 * xi_97 * xi_98 * 0.095238095238095233;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3] = forceTerm_6 + xi_128 * 0.015873015873015872 + xi_144 + xi_182 - xi_196 - xi_197 + xi_206 + xi_99 * 0.095238095238095233 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3] = forceTerm_7 + xi_213 + xi_214 + xi_219 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3] = forceTerm_8 + xi_194 + xi_213 + xi_220 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3] = forceTerm_9 + xi_214 + xi_220 + xi_221 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3] = forceTerm_10 + xi_194 + xi_219 + xi_221 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3] = forceTerm_11 + xi_232 + xi_237 + xi_243 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3] = forceTerm_12 - xi_156 - xi_232 - xi_244 - xi_245;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3] = forceTerm_13 + xi_250 + xi_253 + xi_258 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3] = forceTerm_14 + xi_250 + xi_259 + xi_260 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3] = forceTerm_15 - xi_145 - xi_237 - xi_245 - xi_261;
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3] = forceTerm_16 + xi_243 + xi_244 + xi_261 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3] = forceTerm_17 + xi_258 + xi_259 + xi_262 + _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3];
    _data_pdfs_tmp[_stride_pdfs_tmp_0 * ctr_0 + _stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3] = forceTerm_18 + xi_253 + xi_260 + xi_262 + _data_pdfs[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3];
  }
}
} // namespace internal_streamcollidesweepthermalizeddoubleprecisioncuda_streamcollidesweepthermalizeddoubleprecisioncuda

void StreamCollideSweepThermalizedDoublePrecisionCUDA::run(IBlock *block, gpuStream_t stream) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  auto pdfs = block->getData<gpu::GPUField<double>>(pdfsID);
  auto force = block->getData<gpu::GPUField<double>>(forceID);
  gpu::GPUField<double> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  auto &omega_bulk = this->omega_bulk_;
  auto &seed = this->seed_;
  auto &omega_even = this->omega_even_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_shear = this->omega_shear_;
  auto &time_step = this->time_step_;
  auto &block_offset_1 = this->block_offset_1_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &kT = this->kT_;
  auto &block_offset_2 = this->block_offset_2_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force->nrOfGhostLayers()))
  double *RESTRICT const _data_force = force->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, -1, 0);
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
  internal_streamcollidesweepthermalizeddoubleprecisioncuda_streamcollidesweepthermalizeddoubleprecisioncuda::streamcollidesweepthermalizeddoubleprecisioncuda_streamcollidesweepthermalizeddoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
  pdfs->swapDataPointers(pdfs_tmp);
}

void StreamCollideSweepThermalizedDoublePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block, gpuStream_t stream) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto pdfs = block->getData<gpu::GPUField<double>>(pdfsID);
  auto force = block->getData<gpu::GPUField<double>>(forceID);
  gpu::GPUField<double> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  auto &omega_bulk = this->omega_bulk_;
  auto &seed = this->seed_;
  auto &omega_even = this->omega_even_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_shear = this->omega_shear_;
  auto &time_step = this->time_step_;
  auto &block_offset_1 = this->block_offset_1_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &kT = this->kT_;
  auto &block_offset_2 = this->block_offset_2_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force->nrOfGhostLayers()))
  double *RESTRICT const _data_force = force->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
  internal_streamcollidesweepthermalizeddoubleprecisioncuda_streamcollidesweepthermalizeddoubleprecisioncuda::streamcollidesweepthermalizeddoubleprecisioncuda_streamcollidesweepthermalizeddoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
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
