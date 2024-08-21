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
//! \\file DiffusiveFluxKernelThermalized_single_precision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3, lbmpy_walberla/pystencils_walberla from waLBerla commit b0842e1a493ce19ef1bbb8d2cf382fc343970a7f

#include <cmath>

#include "DiffusiveFluxKernelThermalized_single_precision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include "philox_rand.h"

#define FUNC_PREFIX

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

namespace internal_13067439141d91943f586adb1c937d5b {
static FUNC_PREFIX void diffusivefluxkernelthermalized_single_precision_diffusivefluxkernelthermalized_single_precision(float D, float *RESTRICT const _data_j, float *RESTRICT const _data_rho, int64_t const _size_j_0, int64_t const _size_j_1, int64_t const _size_j_2, int64_t const _stride_j_0, int64_t const _stride_j_1, int64_t const _stride_j_2, int64_t const _stride_j_3, int64_t const _stride_rho_0, int64_t const _stride_rho_1, int64_t const _stride_rho_2, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, uint32_t field_size_0, uint32_t field_size_1, uint32_t field_size_2, uint32_t seed, uint32_t time_step) {
  {
    {
      {
        if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {

          float random_3_0;
          float random_3_1;
          float random_3_2;
          float random_3_3;
          philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

          float random_2_0;
          float random_2_1;
          float random_2_2;
          float random_2_3;
          philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

          float random_1_0;
          float random_1_1;
          float random_1_2;
          float random_1_3;
          philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

          float random_0_0;
          float random_0_1;
          float random_0_2;
          float random_0_3;
          philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

          _data_j[_stride_j_0 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0] - _data_rho[_stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0] + 0.5f * _data_rho[_stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 * ctr_0 + 12 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0]), 0.5f) * 1.5025119784898082f;
          }
        }
        if (0 < _size_j_1 - 1 && 0 < _size_j_2 - 1) {

          float random_3_0;
          float random_3_1;
          float random_3_2;
          float random_3_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

          float random_2_0;
          float random_2_1;
          float random_2_2;
          float random_2_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

          float random_1_0;
          float random_1_1;
          float random_1_2;
          float random_1_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

          float random_0_0;
          float random_0_1;
          float random_0_2;
          float random_0_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

          _data_j[_stride_j_0 * (_size_j_0 - 1) + 12 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * (_size_j_0 - 1)]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1)]), 0.5f) * 1.5025119784898082f;
        }
      }
      for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
        {
          {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + 6 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2]) * 0.11520472029718914f + (random_1_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && 1 < _size_j_0 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + 8 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1]) * 0.11520472029718914f + (random_2_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + 6 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * 0.11520472029718914f + (random_1_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_0 < _size_j_0 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + 8 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1]) * 0.11520472029718914f + (random_2_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
          {
            if (ctr_1 > 0 && 0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + 6 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]) * 0.11520472029718914f + (random_1_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && 0 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (0 < _size_j_2 - 1 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
        }
      }
      {
        {
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1 && 1 < _size_j_0 - 1) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + 8 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)]) * 0.11520472029718914f + (random_2_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)]), 0.5f) * 1.6628028407278295f;
          }
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] - _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1)] + 0.5f * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
          }
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1 && ctr_0 < _size_j_0 - 1) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + 8 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)]) * 0.11520472029718914f + (random_2_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)]), 0.5f) * 1.6628028407278295f;
          }
          if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1)] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
          }
        }
        if (_size_j_1 - 1 > 0 && 0 < _size_j_2 - 1) {

          float random_3_0;
          float random_3_1;
          float random_3_2;
          float random_3_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

          float random_2_0;
          float random_2_1;
          float random_2_2;
          float random_2_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

          float random_1_0;
          float random_1_1;
          float random_1_2;
          float random_1_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

          float random_0_0;
          float random_0_1;
          float random_0_2;
          float random_0_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

          _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1)] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
        }
      }
    }
    for (int64_t ctr_2 = 1; ctr_2 < _size_j_2 - 1; ctr_2 += 1) {
      {
        {
          {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_1_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_1_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
          {
            if (ctr_2 > 0 && 0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_1_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_2 > 0 && 0 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (0 < _size_j_1 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
        }
        for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
          {
            {
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
              }
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
              }
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_0_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_1_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.11520472029718914f + (random_1_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.11520472029718914f + (random_1_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              if (ctr_1 > 0 && ctr_2 > 0 && 1 < _size_j_0 - 1) {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.11520472029718914f + (random_1_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              if (ctr_1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.11520472029718914f + (random_2_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
            }
            for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 2 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_0_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_1_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.11520472029718914f + (random_1_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.11520472029718914f + (random_1_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.11520472029718914f + (random_1_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.11520472029718914f + (random_2_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
            }
            {
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_0_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 4 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_1_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 5 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.11520472029718914f + (random_1_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 6 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.11520472029718914f + (random_1_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
              {

                float random_3_0;
                float random_3_1;
                float random_3_2;
                float random_3_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

                float random_2_0;
                float random_2_1;
                float random_2_2;
                float random_2_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

                float random_1_0;
                float random_1_1;
                float random_1_2;
                float random_1_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

                float random_0_0;
                float random_0_1;
                float random_0_2;
                float random_0_3;
                philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

                _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * ctr_2 + 12 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_3_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
              }
            }
          }
        }
        {
          {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_0_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && 1 < _size_j_0 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.11520472029718914f + (random_1_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (_size_j_1 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.11520472029718914f + (random_2_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.16292407789368385f + (random_0_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.977416969040271f;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_0_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 7 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.11520472029718914f + (random_1_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (_size_j_1 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 8 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.11520472029718914f + (random_2_0 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
          {
            if (ctr_2 > 0 && _size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 3 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]) * 0.11520472029718914f + (random_0_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_2 > 0 && _size_j_1 - 1 > 0) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (_size_j_1 - 1 > 0 && ctr_2 < _size_j_2 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((block_offset_2 + ctr_2) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * ctr_2 + 10 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]) * 0.09406426022938992f + (random_2_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * ctr_2] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * ctr_2 + _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
        }
      }
    }
    {
      {
        if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {

          float random_3_0;
          float random_3_1;
          float random_3_2;
          float random_3_3;
          philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

          float random_2_0;
          float random_2_1;
          float random_2_2;
          float random_2_3;
          philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

          float random_1_0;
          float random_1_1;
          float random_1_2;
          float random_1_3;
          philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

          float random_0_0;
          float random_0_1;
          float random_0_2;
          float random_0_3;
          philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

          _data_j[_stride_j_0 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
          }
        }
        if (_size_j_2 - 1 > 0 && 0 < _size_j_1 - 1) {

          float random_3_0;
          float random_3_1;
          float random_3_2;
          float random_3_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

          float random_2_0;
          float random_2_1;
          float random_2_2;
          float random_2_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

          float random_1_0;
          float random_1_1;
          float random_1_2;
          float random_1_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

          float random_0_0;
          float random_0_1;
          float random_0_2;
          float random_0_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

          _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
        }
      }
      for (int64_t ctr_1 = 1; ctr_1 < _size_j_1 - 1; ctr_1 += 1) {
        {
          {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 2 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]) * 0.16292407789368385f + (random_0_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]), 0.5f) * 1.977416969040271f;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.11520472029718914f + (random_1_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.11520472029718914f + (random_1_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
          for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 2 * _stride_j_3] = D * (-_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]) * 0.16292407789368385f + (random_0_2 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)]), 0.5f) * 1.977416969040271f;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.11520472029718914f + (random_1_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.11520472029718914f + (random_1_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
          {
            if (ctr_1 > 0 && _size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 5 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.11520472029718914f + (random_1_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
            }
            if (ctr_1 > 0 && _size_j_2 - 1 > 0) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
            if (_size_j_2 - 1 > 0 && ctr_1 < _size_j_1 - 1) {

              float random_3_0;
              float random_3_1;
              float random_3_2;
              float random_3_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

              float random_2_0;
              float random_2_1;
              float random_2_2;
              float random_2_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

              float random_1_0;
              float random_1_1;
              float random_1_2;
              float random_1_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

              float random_0_0;
              float random_0_1;
              float random_0_2;
              float random_0_3;
              philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((block_offset_1 + ctr_1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

              _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * ctr_1 + _stride_j_2 * (_size_j_2 - 1) + 11 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * ctr_1 + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * ctr_1 + _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
            }
          }
        }
      }
      {
        {
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0 && 1 < _size_j_0 - 1) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.11520472029718914f + (random_1_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
          }
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
          }
        }
        for (int64_t ctr_0 = 2; ctr_0 < _size_j_0 - 1; ctr_0 += 1) {
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0 && ctr_0 < _size_j_0 - 1) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 7 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.11520472029718914f + (random_1_3 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.6628028407278295f;
          }
          if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {

            float random_3_0;
            float random_3_1;
            float random_3_2;
            float random_3_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

            float random_2_0;
            float random_2_1;
            float random_2_2;
            float random_2_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

            float random_1_0;
            float random_1_1;
            float random_1_2;
            float random_1_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

            float random_0_0;
            float random_0_1;
            float random_0_2;
            float random_0_3;
            philox_float4(time_step, ((block_offset_0 + ctr_0) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

            _data_j[_stride_j_0 * ctr_0 + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * ctr_0 + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * ctr_0 - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
          }
        }
        if (_size_j_1 - 1 > 0 && _size_j_2 - 1 > 0) {

          float random_3_0;
          float random_3_1;
          float random_3_2;
          float random_3_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);

          float random_2_0;
          float random_2_1;
          float random_2_2;
          float random_2_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);

          float random_1_0;
          float random_1_1;
          float random_1_2;
          float random_1_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);

          float random_0_0;
          float random_0_1;
          float random_0_2;
          float random_0_3;
          philox_float4(time_step, ((_size_j_0 + block_offset_0 - 1) % (field_size_0)), ((_size_j_1 + block_offset_1 - 1) % (field_size_1)), ((_size_j_2 + block_offset_2 - 1) % (field_size_2)), 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);

          _data_j[_stride_j_0 * (_size_j_0 - 1) + _stride_j_1 * (_size_j_1 - 1) + _stride_j_2 * (_size_j_2 - 1) + 9 * _stride_j_3] = D * (_data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] - _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]) * 0.09406426022938992f + (random_2_1 - 0.5f) * powf(D * (0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) + _stride_rho_1 * (_size_j_1 - 1) + _stride_rho_2 * (_size_j_2 - 1)] + 0.5f * _data_rho[_stride_rho_0 * (_size_j_0 - 1) - _stride_rho_0 + _stride_rho_1 * (_size_j_1 - 1) - _stride_rho_1 + _stride_rho_2 * (_size_j_2 - 1) - _stride_rho_2]), 0.5f) * 1.5025119784898082f;
        }
      }
    }
  }
}
} // namespace internal_13067439141d91943f586adb1c937d5b

void DiffusiveFluxKernelThermalized_single_precision::run(IBlock *block) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  auto rho = block->getData<field::GhostLayerField<float, 1>>(rhoID);
  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);

  auto &field_size_0 = this->field_size_0_;
  auto &block_offset_1 = this->block_offset_1_;
  auto &seed = this->seed_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &time_step = this->time_step_;
  auto &field_size_2 = this->field_size_2_;
  auto &field_size_1 = this->field_size_1_;
  auto &D = this->D_;
  auto &block_offset_2 = this->block_offset_2_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(j->nrOfGhostLayers()))
  float *RESTRICT const _data_j = j->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho->nrOfGhostLayers()))
  float *RESTRICT const _data_rho = rho->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(int64_c(j->xSize()) + 2))
  const int64_t _size_j_0 = int64_t(int64_c(j->xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(int64_c(j->ySize()) + 2))
  const int64_t _size_j_1 = int64_t(int64_c(j->ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(int64_c(j->zSize()) + 2))
  const int64_t _size_j_2 = int64_t(int64_c(j->zSize()) + 2);
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  internal_13067439141d91943f586adb1c937d5b::diffusivefluxkernelthermalized_single_precision_diffusivefluxkernelthermalized_single_precision(D, _data_j, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2, block_offset_0, block_offset_1, block_offset_2, field_size_0, field_size_1, field_size_2, seed, time_step);
}

void DiffusiveFluxKernelThermalized_single_precision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto rho = block->getData<field::GhostLayerField<float, 1>>(rhoID);
  auto j = block->getData<field::GhostLayerField<float, 13>>(jID);

  auto &field_size_0 = this->field_size_0_;
  auto &block_offset_1 = this->block_offset_1_;
  auto &seed = this->seed_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &time_step = this->time_step_;
  auto &field_size_2 = this->field_size_2_;
  auto &field_size_1 = this->field_size_1_;
  auto &D = this->D_;
  auto &block_offset_2 = this->block_offset_2_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(j->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(j->nrOfGhostLayers()))
  float *RESTRICT const _data_j = j->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho->nrOfGhostLayers()))
  float *RESTRICT const _data_rho = rho->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(j->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
  const int64_t _size_j_0 = int64_t(int64_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
  const int64_t _size_j_1 = int64_t(int64_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(j->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
  const int64_t _size_j_2 = int64_t(int64_c(ci.zSize()) + 2);
  const int64_t _stride_j_0 = int64_t(j->xStride());
  const int64_t _stride_j_1 = int64_t(j->yStride());
  const int64_t _stride_j_2 = int64_t(j->zStride());
  const int64_t _stride_j_3 = int64_t(1 * int64_t(j->fStride()));
  const int64_t _stride_rho_0 = int64_t(rho->xStride());
  const int64_t _stride_rho_1 = int64_t(rho->yStride());
  const int64_t _stride_rho_2 = int64_t(rho->zStride());
  internal_13067439141d91943f586adb1c937d5b::diffusivefluxkernelthermalized_single_precision_diffusivefluxkernelthermalized_single_precision(D, _data_j, _data_rho, _size_j_0, _size_j_1, _size_j_2, _stride_j_0, _stride_j_1, _stride_j_2, _stride_j_3, _stride_rho_0, _stride_rho_1, _stride_rho_2, block_offset_0, block_offset_1, block_offset_2, field_size_0, field_size_1, field_size_2, seed, time_step);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
